__version__ = '0.2.8'
import pysam
import re
import os
import json
import queue
import pandas as pd
import numpy as np
from scipy.stats import entropy
from itertools import chain, product
from collections import Counter,defaultdict
from svcaller import SvTyper
from utils import loadConfig, getLogger, BamRegionViewer, BamViewer_Error

        
class StarTyper:
    '''
    TODO annot
    '''
    
    def __init__( self, config_yaml, logFile='star_typer.log', callMode='wgs',
                  logLevel='INFO', verbose=False, randSeed=42, grayscale=False ):
        '''
        TODO: fill this out
        '''
        self.config    = loadConfig( config_yaml )
        self.reference = self.config[ 'reference' ]
        self.minCov    = self.config[ 'minCov' ][ callMode ]
        self.minFreq   = self.config[ 'minFreq' ][ callMode ]
        self.warnings  = queue.Queue()
        self.log       = getLogger( self.__repr__(), logFile, self.warnings, stdout=verbose, logLevel=logLevel )
        self.svtyper   = SvTyper( self.config, self.log )  
        self.bamRegion = BamRegionViewer( self.config, minCov=self.minCov, minFreq=self.minFreq, logger=self.log, callMode=callMode )
        self.coreVar   = pd.read_csv( self.config[ 'coreVariants' ], index_col=0 )
        self.coreMeta  = self._loadCoreMeta( self.config[ 'coreMetaData' ] )
        self.diplotype = Diplotype( self.config, self.log, self.coreMeta, grayscale )
        self.randGen   = np.random.default_rng( randSeed )
        self.grayscale = grayscale
        self.callMode = callMode
        self.report   = []

    def run( self, sampleBam, vcf=None ):
        '''
        TODO: annot
        '''
        self.log.info( f'Processing {sampleBam}' )
        # clear any previous result
        self.diplotype.reset()
        # load sample data
        try:
            self.bamRegion.loadSample( sampleBam, vcf=vcf )
            self.sampleVars = self.bamRegion.sampleVars
            self.readMeta   = self.bamRegion.readMeta
        except BamViewer_Error as e:
            self.log.error( e.args[0] )
            return None
        # print ("xxxx")
        # check for common parameter mistake -- if no reads match ccs naming convention, might be consensus
        if not self.readMeta.index.get_level_values('hifi_read').str.match(r'^.*/\d+/ccs$').any() \
           and self.callMode != "consensus":
            self.log.warn( 'Query names do not look like HiFi reads -- use \"consensus\" mode for consensus' )
        # call sv
        # print ("self.sampleVars", self.sampleVars)
        self.log.info( 'Identifying SV signatures' )
        self.svtyper.call_sv( self.sampleVars, self.readMeta )
        
        # print ("self.diplotype.calls", self.diplotype.calls)
        self.log.info( 'Phasing and Calling Alleles' )
        # call the star alleles
        for lbl in ['noSV'] + list( self.readMeta.SV.dropna().unique() ): #do noSV first
            self.log.debug( f'Calling {lbl}' )
            self.log.info( f'Calling {lbl}' )
            # print (">>>1", self.diplotype.calls)
            calls = self.call_group( lbl, self.callers[ lbl ] )
            # print (">>>2", self.diplotype.calls)
            if calls not in [ None, {} ]:
                self.diplotype.calls[ lbl ] = calls
            # print ("<<<",lbl, calls, self.diplotype.calls)
        if len( self.diplotype.calls ) == 0:
            self.log.error( 'No Calls found!' )
            return None
        # print ("self.diplotype.calls", self.diplotype.calls)
        self.log.debug( f'Alleles Identified: {dict( self.diplotype.calls )}' )
        # make diplotypeG
        self.log.info( 'Calling Diplotype' )
        self.diplotype._make_diplotype()
        self.log.info( f'Diplotype called: {self.diplotype}' )
        self.report.append( self.makeReport() )
        return self.diplotype.diplotype
    
    def return_read_meta(self):
        # return self.svtyper.readMeta
        return self.readMeta

    def _loadCoreMeta( self, metaFile ):
        evalStr = ( lambda v: pd.Series( eval( v ) ) )
        return pd.read_csv( metaFile, 
                            index_col=0, 
                            converters={ 'varPos' : lambda v: pd.Series( eval( v ) ),
                                         'rsID'   : lambda v: eval(v) if v else [] } )

    def call_group( self, label, cfunc, updateHaplotype=True ):
        mask = self.readMeta.SV.fillna( 'noSV' ) == label 
        if label == 'noSV':
            mask &= self.readMeta[ 'cov_CYP2D6' ] > 0
        grpIdxs = self.readMeta[ mask ].index
        nrecs   = len(grpIdxs)
        # check that there is at least minCov reads to remove spurious stuff
        if nrecs < self.minCov:
            self.log.warning( f'Low Coverage call: {label} ({nrecs}<{self.minCov})' )
            return {}
        if updateHaplotype:
            clusters = None
            if nrecs >= ( 2 * self.minCov ):
                self.log.debug( f'Attempting to phase {label} reads' )
                clusters = self.bamRegion.updateHaplotypes( grpIdxs, label=label, 
                                                            window=self.config[ 'phaseWindows' ][ label ] )
            else:
                self.log.info( f'Not enough reads to phase {label} reads' )
            # print ("clusters", clusters)
            # if no phasing, set HP to {label}_0
            if clusters is None:
                self.readMeta.loc[ grpIdxs, 'HP' ] = f'{label}_0'                

        return cfunc( self, self.readMeta.loc[ grpIdxs ] )
        #TODO warn of low cov issues
        #if len( lowcov ):
        #    total = len( self.coreMatchesSpanning )
        #    lc    = ', '.join( lowcov.apply( lambda r: f'{r.call}({r.nreads}/{total} reads)', axis=1 ) )
        #    self.log.warning( f'{label}: low-coverage call {lc}' )
        #return list( counts.query( 'nreads >= @self.minCov' ).call.values )

    def _get_del_calls( self, reads ):
        if len( reads ) < self.minCov:
            self.log.warning( f'{len( reads )} deletion (*5) reads found!' )
            return { self.readMeta.HP.loc[ reads.index[0] ] : ('(*5)',) }
        else:
            return { hp : ('*5',) for hp in self.readMeta.loc[ reads.index ].HP.unique() }

    def _get_dup_calls( self, reads ):
        round1_dups = pd.Series( self.matchStarAlleles( reads.index ) )
        #for amp/consensus we don't need to resort -- just return original calls
        if self.callMode in [ 'amplicon','consensus'] :
            return round1_dups.to_dict()
        # after round one, we merge same allele dup/nosv and re-sort reads into first/last cyp2d6
        dupLabels = round1_dups.groupby( round1_dups ).apply( lambda d: list( d.index ) ).to_dict()
        # print ("dupLabels", dupLabels)
        noSV = self.diplotype.calls[ 'noSV' ]     
        for lbl,alleles in noSV.items():
            # print ("lbl, alleles", lbl, alleles)
            if alleles in dupLabels:
                dupLabels[ alleles ].append( lbl )
        # print ("dupLabels2", dupLabels)
        uclipLoc = self.config['duplicates']['dup_upstream']['rstart'][1]
        dclipLoc = self.config['duplicates']['dup_downstream']['rstop'][1]
        window   = self.config[ 'phaseWindows' ][ 'duplicate' ]
        # self.readMeta.to_csv("/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/NA19207/readMeta4.csv")
        dupcalls = {}        
        for call,labels in dupLabels.items():
            # print ("call, labels", call, labels)
            creads = self.readMeta.query( 'HP in @labels' )
            #supplemetary reads labeled as "dup_upstream" relabeld to downstream (other alignment)
            sidx = self.readMeta.loc[ creads.index ].query( 'SVlabel == "dup_upstream" and is_primary == False' ).index
            self.readMeta.loc[ sidx, 'SVlabel' ] = 'dup_downstream'
            #unlabeled reads starting before first clip are downstream copy (downstream is left in hg38)
            didx = self.readMeta.loc[ creads.index ].query( 'SVlabel.isnull() & rstart < @uclipLoc' ).index
            self.readMeta.loc[ didx, 'SVlabel' ] = 'dup_downstream'
            #unlabeled reads that map past downstream clips and have no deletion are upstream    
            uidx = self.readMeta.loc[ creads.index ].query( 'SVlabel.isnull() & rstop > @dclipLoc & maxDeletion < 1000' ).index
            self.readMeta.loc[ uidx, 'SVlabel' ] = 'dup_upstream'
            #leftovers are randomly assigned at this point
            #TODO properly sort leftovers by phasing them independently and assigning them to up-/downstream
            ridx = self.readMeta.loc[ creads.index ].query( 'SVlabel.isnull()' ).index
            self.readMeta.loc[ ridx, 'SVlabel' ] = self.randGen.choice( ['dup_rand_upstream', 'dup_rand_downstream'], size=len( ridx ) )
        
            #now re-phase reads, if possible and reassign calls
            #for capture data, we restict possible splits to variants on reads with dup SV sigs
            for position, atype in [ ( 'upstream', 'noSV' ),
                                     ( 'downstream', 'duplicate') ]:
                # print ("position, atype", position, atype)
                idxs     = self.readMeta[ self.readMeta.SVlabel.fillna('').str.endswith( position ) 
                                         & self.readMeta.HP.isin( labels ) ].index
                candIdx = self.readMeta.loc[ idxs ].query( 'SV.notnull()' ).index \
                          if ( self.callMode == 'capture' ) \
                             and ( position == 'upstream' ) \
                             and ( len( set( noSV.values() ) ) >= 2 ) \
                          else None #none = all
                # print (candIdx, idxs)
                prevCall = self.diplotype.calls.get( atype, {} ) if atype == 'noSV' else { **dupcalls, **round1_dups.to_dict() }
                # print("prevCall", prevCall)
                offset   = max( [int( hp.split('_')[-1] ) for hp in prevCall.keys() ], default=-1 )  + 1
                clusters = self.bamRegion.updateHaplotypes( idxs, label=atype, window=window, 
                                                            offset=offset, round=2, candidateSubset=candIdx )
                # print ("new clusters", clusters)
                #remove old noSV labels
                if atype == 'noSV':
                    for hp in labels:
                        if hp.startswith( atype ):
                            del self.diplotype.calls[ atype ][ hp ]
                #add in new ones
                vals = [offset] if clusters is None else clusters.unique()
                for n in vals:
                    lbl = f'{atype}_{n}'
                    # print ("lbl", lbl, call)
                    if atype == 'noSV':
                        self.diplotype.calls[ atype ][ lbl ] = call
                    else:
                        dupcalls[ lbl ] = call

                    # print (lbl, clusters, idxs)        
                if clusters is None:
                    self.readMeta.loc[ idxs, 'HP' ] = lbl
        # self.readMeta.to_csv("/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/NA19207/readMeta3.csv")
        # print ("dupcalls", dupcalls)
        return dupcalls

    def _get_hyb_calls( self, reads ):
        byRead =  reads.SVlabel.map( self.svtyper.config[ 'starAlleles' ][ 'hybrid' ] )\
                       .fillna( '' ).map( tuple ).rename( 'call' )
        res = {}
        for (hp,candidates),rds in self.readMeta.loc[ byRead.index ].join( byRead ).groupby( ['HP','call'] ):
            if len( rds ) < self.minCov:
                self.log.warning( f'Low-coverage hybrid called: {candidates} ({len(rds)}<{self.minCov})' )
                continue
            if len( candidates ) > 1: #identify defining variants
                #will this fail if some hybrid has no variants in the db?
                candDf = self.coreMeta.loc[ map( self.diplotype.hapSortKey, candidates ) ]
                candDf.index = candidates
                candVars = pd.concat( candDf.varPos.values, keys=candDf.index )
                ind6 = candVars[ ~candVars.duplicated( keep=False ) ].reset_index(level=1, drop=True)
                indicators = pd.concat( [ ind6, ind6.map( self.config['d6_d7_homologous_pos'] ) ], axis=1 )
                indicators.index.name = 'allele'
                #consensus at indicators
                matchFound = False
                for col in indicators.columns:
                    cons = self.bamRegion.getConsensusVariants( rds.index, positions=indicators[col].dropna() )\
                                         .assign( allele=indicators.index )\
                                         .set_index( 'allele' )
                    #matching all indicators
                    matching = cons.groupby( 'allele' ).filter( lambda v: (v.VAR != '.').all() )
                    if not matching.empty:
                        newCand = (matching.index[0],)
                        matchFound = True
                        break
                if not matchFound:
                    #select any allele(s) with no indicators as the new candidate list
                    newCand = tuple( candDf.index.difference( indicators.index ) )
                    if len( newCand ) == 0:
                        self.log.warning( f'Hybrid does not match any known allele!' )
                        newCand = candidates
            else:
                newCand = candidates
            res[ hp ] = newCand
        #if only a single "noSV" allele exists, try to split by looking at 
        #reads spanning both d6 and d7 -- is the d7 region "normal" or hybrid
        noSValleles = self.diplotype.calls.get( 'noSV', {} ) 
        if len( noSValleles ) == 1:
            #reads that span both loci
            zmws = self.readMeta.query( 'cov_CYP2D6 > 0 & cov_CYP2D7 > 0' ).index.get_level_values( 'hifi_read' ).unique()
            # dict of locus -> reads, indicating which reads have hybrids mapped where
            hybcalls = defaultdict(set)
            for locus,read in self.svtyper.hybrids.index.droplevel( 'is_primary' ):
                hybcalls[ locus ].add( read )
            # reads with hyb in d7 position
            upstreamHyb = zmws[ ~zmws.isin( hybcalls['CYP2D6'] ) & zmws.isin( hybcalls['CYP2D7'] ) ] 
            # reads with no hyb
            upstreamD7 = zmws[ ~zmws.isin( hybcalls['CYP2D6'] ) & ~zmws.isin( hybcalls['CYP2D7'] ) ]
            if ( len( upstreamHyb ) >= self.minCov ) and ( len( upstreamD7 ) >= self.minCov ):
                self.log.info( 'Separating alleles by upstream hybrid presence/absence' )
                old_noSV, allele = list( noSValleles.items() )[0]
                offset = int( old_noSV.split('_')[-1] )
                #reassign reads
                new_noSV = [ f'noSV_{offset+i}' for i in [1,2] ]
                self.readMeta.loc[ upstreamD7, 'HP' ] = new_noSV[0]
                self.readMeta.loc[ upstreamHyb, 'HP' ] = new_noSV[1] 
                #randomly assign rest of noSV reads that do not have d7 coverage
                undetIdx = self.readMeta.query( 'HP==@old_noSV' ).index
                self.readMeta.loc[ undetIdx, 'HP' ] = self.randGen.choice( new_noSV, len( undetIdx ) )
                self.diplotype.calls[ 'noSV' ] = { lbl: allele for lbl in new_noSV }
        return res

    def _get_noSV_calls( self, reads ):
        # TODO include long reads with normal/hyb or normal/dup together
        return self.matchStarAlleles( reads.index )

    def matchStarAlleles( self, idxs ):
        # catch empty index
        if idxs.empty:
            return {}
        start, stop = self.config[ 'genes' ][ 'CYP2D6' ]
        cov = self._getConsensus( idxs, start, stop, suppressWarnings=True )
        if ( cov.coverage < self.minCov ).any():
            self.log.warning( f'LOW COVERAGE VARIANTS CALLED' ) 
        res         = cov.VAR.unstack( level=0 ).fillna( '.' )
        variants    = res[ ( res != '.' ).any( axis=1 ) ]
        coreMatches = self.coreVar.drop( columns='contig' )\
                                  .set_index( 'pos' )\
                                  .join( variants ).dropna()
        getStar = self._starMatches( ( start, stop ) )
        return {
              hp : getStar( coreMatches[ [ 'coreAllele', 'VAR', hp ] ].query( f'VAR == {hp}' ) )
              for hp in coreMatches.columns[ 2: ]
        }

    def makeReport( self ):
        idxs = self.readMeta.query( 'HP != "-1"' ).index
        start, stop = self.config[ 'genes' ][ 'CYP2D6' ]
        cons = self._getConsensus( idxs, start, stop ) 
        stats = cons.groupby( 'HP' ).coverage.describe()
        count = sum( 1 if not len( calls ) 
                       else (  2 if 'x2' in calls[0] 
                                 else ( 0 if calls[0] == '*5' else 1 ) )
                     for _,alleles in self.diplotype.calls.items() 
                     for _,calls in alleles.items() )
        warnings = []
        while not self.warnings.empty():
            warnings.append( self.warnings.get() )
        res = dict( input=self.bamRegion.sampleBam, 
                    diplotype=str( self.diplotype ),
                    copynumber=count,
                    haplotypes=[],
                    warnings=warnings )
        ## output the self.readMeta to a file
        # self.readMeta.to_csv("/mnt/d/HLAPro_backup/Nanopore_optimize/cyp_results/NA19207/readMeta.csv")
        # print (self.diplotype.haplotypes)
        for i, (hap, hpTags) in enumerate( self.diplotype.haplotypes ):
            # print (hap, hpTags)
            hapData = dict( call=hap, alleles=[] )
            for j, tag in enumerate( hpTags ):
                allele = self.diplotype.callMap[ tag ].split()[-1]
                # check if tag is in stats
                if tag not in stats.index:
                    self.log.warning("No reads assigned to the label %s for %s" % (tag,allele))
                    alleleData = dict(
                        call=      allele,
                        num_reads= 0,
                        meanCover= 0,
                        maxCover=  0,
                        minCover=  0,
                        rsIDs   =  []
                    )
                else:
                    alleleData = dict(
                        call=      allele,
                        num_reads= self.readMeta.query('HP==@tag').index.get_level_values('hifi_read').nunique(),
                        meanCover= round( stats.loc[ tag ]["mean"], 2 ),
                        maxCover=  int( stats.loc[ tag ]["max"] ),
                        minCover=  int( stats.loc[ tag ]["min"] ),
                        rsIDs   =  self.diplotype.rsIDs[ hap ].get( allele, [] )
                    )
                hapData[ 'alleles' ].append( alleleData )
            res[ 'haplotypes' ].append( hapData )
        return res

    def _getConsensus( self, idxs, start, stop, suppressWarnings=False ):
        consFunc = lambda hp: self.bamRegion.getConsensusVariants( hp.index, start=start, stop=stop, suppressWarnings=suppressWarnings )
        return self.readMeta.loc[ idxs ]\
                            .groupby( 'HP' )\
                            .apply( consFunc )

    def _hasAllCoreVars( self, region=None ):
        ' region is a tuple of coordinates to restrict counting (defaults to chr22) '
        def checkCount( varset ):
            allele   = int( varset.coreAllele.iloc[0][1:] )
            expCount = self.coreMeta.loc[ allele ].varPos.between( *region ).sum()
            return len( varset ) == expCount
        return checkCount

    def _alleleSort( self, allele ):
        'sort by impact decr, nvars decr, corenumber incr'
        # change core alleles to integers
        a    = int( allele[1:] )
        meta = self.coreMeta.loc[ a ]
        return -meta.impact, -meta.nVars, a  

    def _starMatches( self, region ):
        def _getstar( matches ):
            # TODO: *1 heuristic check
            if matches.empty: # and add a coverage check?
                return ( '*1', )
            else:
                res = sorted( matches.groupby( 'coreAllele' )\
                                             .apply( self._hasAllCoreVars( region ) )\
                                             .where( lambda x:x ).dropna().index,
                              key=self._alleleSort )
                if len(res) == 0: # only partial hits, but some core variants here
                    part = ','.join( matches.coreAllele.unique() )
                    self.log.warn( f"Partial hits on core alleles {part}, defaulting to *1")
                    return ('*1',)
                else:                    
                    return tuple( res )
        return _getstar

    # caller list 
    callers = dict( deletion  = _get_del_calls, 
                    duplicate = _get_dup_calls, 
                    hybrid    = _get_hyb_calls,
                    noSV      = _get_noSV_calls )

    
class CYP2D6_CoreFormater:

    def __init__( self, alignments, reference, impactXls, 
                  haplotypeTsv, config_yaml, logFile, verbose=True ):
        self.alignments = alignments
        self.reference  = reference
        self.impactXls  = impactXls
        self.hapTsv     = haplotypeTsv
        self.config     = loadConfig( config_yaml )
        self.log        = getLogger( self.__repr__(), logFile, stdout=verbose )
        self.region     = self.config[ 'region' ]
        self.loader     = BamRegionViewer( self.config, logger=self.log )
        self.coreVar    = self._loadCoreVar()
        self.coreMeta   = self._loadCoreMeta()

    def saveCoreFiles( self, version, outdir, backup=False ):
        '''
        Save compressed core data files to outdir.
        if backup=True, save a backup of existing files first (<file>.bak)
        '''
        for corefile, df in zip( [ 'coreVariants', 'coreMetaData' ], [ 'coreVar', 'coreMeta' ] ):
            if os.path.isfile( self.config[ corefile ] ):
                os.rename( self.config[ corefile ], self.config[ corefile ] + '.bak' )
            fname = f'{outdir}/{df.lower()}.{version}.csv.gz'
            getattr( self, df ).to_csv( fname )
            self.log.info( f'Saving core data {fname} -- UPDATE your config.yaml!' )

    def _loadCoreVar( self ):
        # load sample variants over region
        core = self.loader._makeDf( self.alignments, self.reference, self.region, qryCol='coreAllele' )
        # remove non-core alleles (alleles with dots in the name, *DD.DDD)
        core = core[~core.coreAllele.str.contains('\.')]
        # rename allele by removing the gene name
        core.coreAllele = core.coreAllele.str.extract( 'CYP2D6(\*\d+)' )
        return core
    
    def _loadCoreMeta( self ):
        import warnings
        # suppress format warning
        with warnings.catch_warnings( record=True ):
            warnings.simplefilter( "always" )
            coreDef = pd.read_excel( self.impactXls,
                                     sheet_name='Allele Function',
                                     header=1)

        alleleField   = 'Allele/cDNA/rsID'
        functionField = 'Allele Clinical Functional Status (Required)'
        
        # drop tandems
        coreDef = coreDef[ ~coreDef[ alleleField ].str.contains( 'x' ) ]
        # change core allele to integer
        coreDef[ alleleField ] = coreDef[ alleleField ].str[1:].astype( int )
        # change column name to something easier
        coreDef.rename( columns={ alleleField : 'coreAllele' }, inplace=True )
        
        coreImpacts = coreDef.set_index( 'coreAllele' )\
                             [ functionField ].map( self.config[ 'impactScores' ] )\
                             .rename( 'impact' )
        
        # make any adjustments to impact vals
        #for allele, impact in self.config[ 'impactAdjustments' ].items():
        #    coreImpacts[ allele ] = impact            
        # get count and list of var positions
        coreVars = pd.concat( [ self.coreVar.groupby( 'coreAllele' ).size().rename( 'nVars' ),
                                self.coreVar.groupby( 'coreAllele' ).pos.apply( sorted ).rename( 'varPos' ) ],
                              axis=1 )
        #change index to integers
        coreVars.index = coreVars.index.str[1:].astype(int)

        #add in rs numbers where known
        haplotypes = pd.read_csv( self.hapTsv, 
                                  sep='\t',
                                  skiprows=[0], 
                                  index_col=0 )
        core_rsID = haplotypes[ ~haplotypes.index.str.contains('\.') & haplotypes.rsID.str.contains('rs')]\
                              .dropna()\
                              .groupby( level=0 )\
                              .rsID.apply(lambda v: [*{*v}])
        core_rsID.index = core_rsID.index.str.split('\*').str[-1].astype(int)
        
        return pd.concat( [ coreImpacts, coreVars, core_rsID ], axis=1 )\
                 .fillna( dict( impact=0, nVars=0, varPos=-1 ) )\
                 .astype( dict( impact=int, nVars=int ) )


class Diplotype:
    '''
    Hold info for star diplotype
    '''
    def __init__( self, config, log, coreMeta, grayscale=False ):
        self.config     = config
        self.log        = log
        self.coreMeta   = coreMeta
        self.grayscale  = grayscale
        self.calls      = defaultdict( dict )
        self.haplotypes = []
        self.diplotype  = None

    def __repr__( self ):
        gene = 'CYP2D6'
        return f'{gene} {self.diplotype}'

    def pharmCat( self ):
        gene = 'CYP2D6'
        def getPharmCat( hap ):
            if '+' in hap:
                return f'[{" + ".join( hap.split("+") )}]'
            else:
                return hap
        # when pangu returns a single call, we need to pad it out for pharmcat
        haps = self.haplotypes + [("[]","")] if len(self.haplotypes) == 1 else self.haplotypes
        return f'{gene}\t{"/".join( getPharmCat( hap ) for hap,_ in haps )}'
    
    def reset( self ):
        self.__init__( self.config, self.log, self.coreMeta, self.grayscale )

    def _make_diplotype( self ):
        def getDupHap( labels, ndups=1 ):
            '''no variant comparison, just return pairs of parent/child.
               returns  tuple( list(pair), list(rest) )   
            '''
            lbls = sorted( labels )
            # print ("getDupHap, lbls", lbls, lbls[:ndups] + lbls[-1:], lbls[ndups:-1])
            return lbls[:ndups] + lbls[-1:], lbls[ndups:-1]
        # print ("self.calls",self.calls)
        # haplotype counters
        hybs             = Counter( self.calls.get( 'hybrid', {} ).values() )
        nosv, dups, dels = [ Counter( [ calls[0] for calls in self.calls.get( atype, {} ).values() ] )
                             for atype in ['noSV', 'duplicate', 'deletion' ] ] 
        #collect labels for each call
        calls = defaultdict(list)
        for _,subset in self.calls.items():
            # print ("subset", subset)
            for hp,call in subset.items():
                a = call[0] if len( call ) else 'unknown'
                calls[ a ].append( hp )
        self.haplotypes = []
        # merge dup haps
        for allele, dcount in dups.items():
            pcount = nosv[ allele ]
            if dcount > pcount: #(multi-dups are rare[?] -> warning)
                self.log.warning( 'Multi-tandem DUP call -- check reads' )
            found = False 
            while pcount and dcount:
                if ( pcount == 1 ): # if one parent, all remaining dups go here
                    mult = dcount
                else:
                    if ( dcount > 2 ):  # if pcount != 1 and hcount > 2, split across parents, starting with 2
                        mult = 2
                    elif ( dcount <= 2 ): # one or less on each parent
                        mult = 1
                    else: # catch rest
                        mult = dcount
                suffix = mult + 1
                haplbls, rest = getDupHap( calls.pop(allele), mult )
                self.haplotypes.append( ( f'{allele}x{suffix}', haplbls ) )
                calls[ allele ] = rest
                nosv[ allele ] -= 1
                pcount -= 1
                dups[ allele ] -= mult
                dcount -= mult
                found = True
            if not found: #no parent -> issue/problem call
                self.log.error( f'No DUP parent found for tandem {allele}' )
                self.haplotypes.append( ( f'{allele}x{dcount}', calls.pop(allele) ) )
                dups[ allele ] -= dcount

        # merge hybrid haps
        self._hybridSelections = {}
        for alleles,hcount in hybs.items():
            found = False
            # potential parents, when not already assigned to some other SV match
            candidateParents = [ allele for allele,count in nosv.items() if count > 0 ]
            for allele in alleles:
                expParents = self.config[ 'starAlleles' ][ 'hybrid_parent' ].get( allele.split('x')[0], [] )
                pMask      = np.in1d( candidateParents, expParents )
                if pMask.any():
                    parent = candidateParents[ np.where( pMask )[ 0 ][ 0 ] ] #first matching candidate
                    pcount = nosv[ parent ]
                    if hcount > pcount: #(multi-hybs are rare[?] -> warning)
                        self.log.warning( 'Multi-tandem HYBRID call -- check reads' )
                    while pcount and hcount:
                        if ( pcount == 1 ): 
                            mult = hcount
                        else:
                            if ( hcount > 2 ):
                                mult = 2
                            elif ( hcount <= 2 ):
                                mult = 1
                            else:
                                mult = hcount
                        upstream = allele if ( mult == 1 ) else f'{allele}x{mult}'
                        lbls = []
                        for a in [parent,allele]:
                            a_lbls = calls.pop(a)
                            lbls.extend( a_lbls[:mult] )
                            calls[a] = a_lbls[mult:]
                        self.haplotypes.append( ( f'{upstream}+{parent}', lbls ) )
                        nosv[ parent ] -= 1
                        pcount -= 1
                        hybs[ alleles ] -= mult
                        hcount -= mult
                        self._hybridSelections[ alleles ] = allele
                    found  = True
                    break
            if not found:
                if len( alleles ) == 0:  # unknown breakpoint, removing here with warning
                            self.log.warning( 'UNKNOWN hybrid breakpoint' )
                else:
                            self.log.warning( f'NO HYBRID parent found for {alleles} (singleton)' )

        # add in remaining alleles to haps
        self.haplotypes.extend( ( call, [lbl] ) for call,lbls in calls.items() for lbl in lbls )
        self.haplotypes = sorted( self.haplotypes, key=lambda h: self.hapSortKey( h[0]) ) 

        if len( self.haplotypes ) > 2:
            self.log.error( f'Extra haplotypes called: {self.haplotypes}' )
    
        self.diplotype = self.joinHaplotypes( self.haplotypes )
        self.rsIDs     = self.get_rsIDs( self.haplotypes )
        return self.diplotype

    def hapSortKey( self, hap ):
        if not type( hap ) == str: # singleton hyb
            hap = hap[0]
        match = re.search( '\*(\d+)[x+.)]?', hap )
        return 999 if match is None else int( match.groups()[0] )

    def get_rsIDs( self, haps ):
        spatt = re.compile( '\*(\d+)' )
        return { hap: { a.group() : self.coreMeta.rsID.get( int( a.groups()[0] ), [] ) 
                      for a in spatt.finditer( hap if type( hap ) == str else hap[0] ) }
                for hap,labels in haps }
    
    def joinHaplotypes( self, haps ):
        return '/'.join( h if type( h ) == str else h[0]
                         for h,lbls in haps )
    
    @property
    def callMap( self ):
        def getCall( allele, kind ):
            if type( allele ) == str:
                return allele
            if allele in self._hybridSelections:
                return self._hybridSelections[ allele ]
            if not len( allele ):
                return f'unknown {kind}'
            return allele[0]

        tag2hap = { lbl : f'haplotype {i+1}{chr(97+j) if len(labels)>1 else ""}'
                    for i, ( call,labels ) in enumerate( self.haplotypes ) 
                    for j,lbl in enumerate( labels ) }

        cm = defaultdict( lambda: 'No Label' )
        cm.update( { hp : f'{tag2hap[hp]}: {getCall(allele,kind)}'
                     for kind,alleles in self.calls.items() 
                     for hp,allele in alleles.items() } )
        return cm

    @property
    def colorMap( self ):
        colors = self.config[ 'grayscale' ] if self.grayscale else self.config[ 'color' ]
        def getColor( lbls ):
            if len( lbls ) == 1:
                kind = lbls[0].split('_')[0]
            else:
                try:
                    kind = ( { lbl.split('_')[0] for lbl in lbls } - {'noSV'} ).pop()
                except KeyError: 
                    # catch case where all alleles are the same, but some sv (eg *10+*10/*10)
                    # this is just a patch to prevent the program from failing to color the outputs
                    kind = '-1'
            return kind, colors[ kind ]
        default = ','.join( map( str, colors[ '-1' ][ 'base' ] ) )
        cmap = defaultdict( lambda: default )
        usedColors = Counter()
        for hap,lbls in self.haplotypes:
            kind,color = getColor( lbls )
            usedColors[ kind ] += 1
            for lbl in lbls:
                base = np.array( color['base'] )
                offset = ( usedColors[kind] - 1 ) * np.array(color['step'])
                cmap[ lbl ] = ','.join( map( str, base + offset ) )
        return cmap

class StarTyper_Error( Exception ):
    pass

class Diplotype_Error( Exception ):
    pass
