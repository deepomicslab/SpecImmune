__version__ = '0.2.8'
import re
import pandas as pd
import numpy as np
import mappy as mp
from collections import Counter
            
class SvTyper:
    '''
    port of Xiao's sv typer, using sampleVar and read metrics gathered by star typer
    '''
    minClip = 100

    def __init__( self, config, log ):
        self.config          = config
        self.log             = log
        self._aligner        = Aligner( self.config[ 'reference' ] )
        self.diff_sites      = self._loadDiffs( config[ 'diff_sites' ] )
        self._hybridPatterns = list( self._makeHybPatterns() )
        self._call_functions = [ self.call_del_and_dup, self.call_hybrid ]

    def call_sv( self, sampleVars, readMeta ):
        # clear values
        self.sampleVars = sampleVars
        self.readMeta   = readMeta
        self.readMeta[ 'SV' ] = None
        self.readMeta[ 'SVlabel' ] = None
        self.hybrids    = None
        # align clips
        self.clips      = self._alignClips()
        # print ("self.clips", self.clips)
        # map dsnps and call source sequence
        self.dsnpVectors = self._mapDiffs( self.sampleVars )
        # call sv categories
        for f in self._call_functions: f()
        # check for un-labeled large alignment features
        self._checkUnlabeled()
        # output the self.readMeta into a table
        # csv_file = "//mnt/d/HLAPro_backup/Nanopore_optimize/data/CYP2D6/pangu-main/demo_data/HG01190.readMeta.csv"
        # self.readMeta.to_csv(csv_file)
        # for lbl, patt in self._hybridPatterns:
        #     print ("lbl", lbl, "patt", patt)

    def _loadDiffs( self, dFile ):
        diffs = pd.read_csv( dFile )
        for loc in self.config['loci'].keys():
            diffs[ f'base_{loc}' ] = diffs[ f'base_{loc}' ].str.upper()
            # shift pos for zero indexing
            diffs[ f'pos_{loc}' ] -= 1
        return diffs

    def _mapDiffs( self, varset ):
        '''
        Extract variants at key differentiating positions and mark gene source/no-data
        '''
        result = {}
        # unstack variants to dataframe
        pileup = varset.set_index( [ 'pos', 'hifi_read', 'is_primary' ] ).VAR.unstack( 'pos' )
        for gene,s in self.config['loci'].items():
            
            replace = { ',|\.' : str( s ),
                        '\*'   : 'x' }
            othergene = ( set( self.config['loci'].keys() ) - {gene} ).pop()
            # get vector of key sites in gene
            diffsites = self.diff_sites.set_index( f'pos_{gene}' )[ f'base_{othergene}' ]
            # null values = no info (x) 
            # refcalls (,.) assigned gene number
            keyvars = pileup.reindex( columns=diffsites.index )\
                            .fillna( 'x' )\
                            .replace( replace, regex=True )
            # realign arrays for comparison broadcast
            # assign "othergene" number where bases match
            a,b = keyvars.align( diffsites, axis=1, copy=False, fill_value='x' ) 
            # print ("gene", gene, keyvars)
            keyvars.mask( a==b, str( self.config[ 'loci' ][ othergene ]  ), inplace=True  ) 
            # fill in x for anything else left over
            keyvars.mask( ~keyvars.isin( [ '1', '2', 'x' ] ), 'x', inplace=True )
            result[ gene ] = keyvars.apply( ''.join, axis=1 )
        return pd.concat( result, names=['locus'] ).rename( 'dsnp' )

    def _alignClips( self ):
        columns = [ 'hifi_read', 'end', 'rstart', 'rstop', 'strand' ]
        recs = []
        for end in [ '5p', '3p' ]:
            for ( read, isprimary ), seq in self.readMeta[ f'clipSeq_{end}' ].dropna().items():
                aln = self._aligner(seq)
                if aln:
                    recs.append( [ read, end, aln.r_st, aln.r_en, aln.strand ] )
        return pd.DataFrame( recs, columns=columns ).set_index( 'hifi_read' )

    def call_del_and_dup( self ):
        ' Updates readMeta SV and SVlabel columns'
        candidates = pd.merge( self.readMeta.xs( True, level='is_primary' ), 
                               self.clips, on='hifi_read', how='left', suffixes=[ '', '_clip' ] )
        for sv in [ 'deletions', 'duplicates' ]:
            for feat, crits in self.config[ sv ].items():
                # print ("feat", feat)
                mask = np.ones(len(candidates), dtype=bool)
                
                for fld, (op, val) in crits.items():
                    # print ("fld", fld, "val", val, candidates[fld], op(candidates[fld], val).values)
                    ## handle each element separately
                    # print (fld, val, op)
                    mask &= op(candidates[fld], val).values
                
                svReads    = candidates[ mask ].index
                # if sv == 'deletions':
                #     print (sv, feat, svReads)
                self.readMeta.loc[ svReads, [ 'SV', 'SVlabel' ] ] = [ sv[:-1], feat ]    
            # break  
                # print (self.readMeta)   
        ## print each row in self.readMeta
        # for index, row in self.readMeta.iterrows():
        #    if row['SV'] == 'deletions':
        #         print (index, row)   

    def call_hybrid( self ):
        # filter for hapstrings with both genes
        minChain = self.config[ 'hybrids' ][ 'minChain' ]
        mask     = ( self.dsnpVectors.str.count( '1' ) >= minChain ) & ( self.dsnpVectors.str.count( '2' ) >= minChain )
        # print ("mask", mask)
        print ("self.dsnpVectors", self.dsnpVectors[mask])
        # find breakpoints and update readmeta

        self.hybrids = self.dsnpVectors[ mask ]\
                           .apply( self._findBp( start=3, end=-3 ) )\
                           .dropna() # drop reads that are in the mask, but failed to match a bp
        # print ("self.hybrids", self.hybrids)
        idxs         = self.hybrids.reset_index( 'locus', drop=True ).index
        self.readMeta.loc[ idxs, 'SV' ] = 'hybrid'
        self.readMeta.loc[ idxs, 'SVlabel' ] = self.hybrids.values 
        self.readMeta[ 'hybrid_maploc' ] = None
        self.readMeta.loc[ idxs, 'hybrid_maploc' ] = self.hybrids.reset_index( 'locus' )\
                                                         .groupby( ['hifi_read','is_primary'], group_keys=False )\
                                                         .locus.apply( lambda loci: '+'.join( loci.sort_values() ) )

        # for index, row in self.readMeta.iterrows():
        #    if row['SV'] == 'hybrid':
        #         print (index, row)  

    def _makeHybPatterns( self ):
        label = list( '12' )
        maxChain = 6
        minChain = self.config[ 'hybrids' ][ 'minChain' ]
        for s in range( maxChain, minChain - 1, -1 ):
            q = [ '1' ] * s + [ '2' ] * s 
            for o in [ 1, -1 ]:
                yield self._makeHybLabel( label[ ::o ] ), re.compile( ''.join( q )[ ::o ] )
            for e in range( 2 * s ): #allow one x
                qe = list( q ) 
                qe[ e ] = 'x'
                for o in [ 1, -1 ]:
                    yield self._makeHybLabel( label[ ::o ] ), re.compile( ''.join( qe )[ ::o ] )

    def _findBp( self, start=0, end=None ):
        ' use breakpoint search pattern with error to locate bp (optionally skipping chars on end) '
        def fbp( hap ):
            h = hap[ start : end ]
            for lbl, patt in self._hybridPatterns:
                # print ("lbl", lbl, "patt", patt, h)
                m = patt.search( h )
                if m: 
                    print ("lbl", lbl, "patt", patt, h)
                    bp = m.start() + len( m.group() ) // 2 - 1 + start
                    # make sure its actually a hyb breakpoint and not an island
                    left, right = Counter( hap[ :bp+1 ].strip( 'x' ) ), Counter( hap[ bp+1: ].strip( 'x' ) )
                    if sorted( [ left.most_common(1)[0][0], right.most_common(1)[0][0] ] ) == [ '1','2' ]:
                        # print ("bp", bp, "lbl", lbl, left, right, self.diff_sites.loc[ bp ].annotation)
                        return '_'.join( [ lbl, self.diff_sites.loc[ bp ].annotation ] )
            return None
        return fbp

    def _makeHybLabel( self, lbl ):
        orient = self.config[ 'hybrids' ][ 'exonOrder' ]
        return '_'.join( [ self.config[ 'abbreviation' ][ int( i ) ] for i in lbl ][ ::orient ] ) 

    def _checkUnlabeled( self, minClip=100, minSameClip=2, tolerance=50 ):
        '''TODO: better set min clip; document'''
        def findUnlabeled( reads ):
            return ( ( reads.clip_5p > self.minClip ).any() \
                    | ( reads.clip_3p > self.minClip ).any() ) \
                   &  reads.SV.isnull().all()

        unlabeledSV = self.readMeta.groupby( level='hifi_read' ).filter( findUnlabeled ) 
        if len( unlabeledSV ):
            # only report if more than one read has the same clip and is unlabeled
            unlabeledGroup = pd.DataFrame( [ ( np.abs( np.subtract.outer( unlabeledSV[ end ].values, 
                                                                          unlabeledSV[ end ].values ) ) <= tolerance ).sum( axis=0 ) >= minSameClip 
                                             for end in [ 'rstart','rstop' ] ] ).any( axis=0 )
            if unlabeledGroup.any():
                reads = "\n".join( unlabeledSV[ unlabeledGroup.values ].index.get_level_values( 'hifi_read' ).unique() )
                self.log.warning( f'Unlabeled Clipped Reads\n{reads}' ) 

class Aligner:
    '''
    Wrapper around mappy/minimap2 aligner, allows for parameter setting and failure handling
    Best to use the same alignment parameters as input bam for consistent calls
    '''
    presets = {
                'splice'    : { 'preset'  : 'splice' },
                'map-hifi'  : { 'preset'  : 'map-hifi' }, #minimap2 hifi preset
                'pb-hifi'   : { 'scoring' : (2,5,5,4,56,1) }, #pbmm2 hifi settings
                'cyp2d6'    : { 'scoring' : ( 1,4,10,2,26,1) } #recommended for cyp2d6 to align hybrid exon9
              } 

    def __init__( self, reference, preset='pb-hifi' ):
        self.kwargs = {'fn_idx_in'   : reference,
                       'best_n'      : 1,
                       'k'           : 19,
                       'w'           : 19,
                       'min_dp_score': 200}
        self.kwargs.update( self.presets[ preset ] )
        self._aligner = mp.Aligner( **self.kwargs )
    
    def __call__( self, seq, skipFailed=True ):
        try:
            return list( filter( lambda a: a.is_primary, self._aligner.map( seq=seq ) ) )[0]
        except IndexError:
            if skipFailed:
                return None
            else:
                raise SVTyper_Error( f'Unable to align {rec.name}' )

class SVTyper_Error( Exception ):
   pass 

