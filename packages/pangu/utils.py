__version__ = '0.2.8'
import logging
import sys
import yaml
import pysam
import re
import numpy as np
import pandas as pd
# from importlib.resources import files
from scipy.stats import entropy
from scipy.linalg import svd, diagsvd
from sklearn.cluster import KMeans

if sys.version_info >= (3, 9):
    import importlib.resources as resources
else:
    import importlib_resources as resources




class QueuingHandler(logging.Handler):
    """A thread safe logging.Handler that writes messages into a queue object.
    """
    def __init__(self, *args, message_queue, **kwargs):
        """Initialize by copying the queue and sending everything else to superclass."""
        logging.Handler.__init__(self, *args, **kwargs)
        self.message_queue = message_queue

    def emit(self, record):
        """Add the formatted log message (sans newlines) to the queue."""
        self.message_queue.put(self.format(record).rstrip('\n'))


def getLogger( name, filename, message_queue, stdout=False, logLevel='DEBUG'):
    logger = logging.getLogger( name )
    level  = getattr( logging, logLevel )
    formatter = logging.Formatter( '%(asctime)s | %(levelname)s | %(message)s' )
    fhandler  = logging.FileHandler( filename=filename, mode='w')
    fhandler.setLevel( level )
    fhandler.setFormatter( formatter )
    logger.addHandler( fhandler )
    # save warnings in queue for exporting later
    warnqhandler = QueuingHandler( message_queue=message_queue, level=logging.WARNING )
    wformat = logging.Formatter( '%(message)s' )
    warnqhandler.setFormatter( wformat )
    logger.addHandler( warnqhandler ) 
    if stdout:
        shandler = logging.StreamHandler( sys.stdout )
        shandler.setLevel( level )
        shandler.setFormatter( formatter )
        logger.addHandler( shandler )
    logger.setLevel( level )
    return logger

def isclose( atol=20, rtol=0 ):
    def _isclose( a, b ):
        return pd.Series( np.isclose( a.fillna( np.nan ), b, atol=atol, rtol=rtol ) )
    return _isclose

OPS = {
       '='  : np.equal,
       '!=' : np.not_equal,
       '<'  : np.less,
       '>'  : np.greater,
       '<=' : np.less_equal,
       '>=' : np.greater_equal,
      }

def loadConfig( configYaml ):
    # import importlib.resources as resources
    # data_files = resources.files('pangu.data.CYP2D6')
    data_files = sys.path[0] + '/data/CYP2D6'
    with open( configYaml, "r" ) as stream:
        try:
            config = yaml.safe_load( stream )
        except yaml.YAMLError as e:
            raise Config_Error( f'Error reading {configYaml}: \n\n{e}' )
    #set data locations
    for key,fname in config[ 'data_files' ].items():
        # config[ key ] = str( data_files.joinpath( fname ) )
        config[ key ] = str( data_files) + "/" + fname
        # print (data_files, config[ key ])
    # replace sv rule operator labels with op functions
    # and include "isclose" operator with config tolerance
    ops = { **OPS, **{ '~' : isclose( atol=config['tolerance'] ) } }
    for sv in [ 'deletions', 'duplicates' ]:
        for feat, crit in config[ sv ].items():
            for field, ( op, val ) in crit.items():
                config[ sv ][ feat ][ field ][0] = ops[ op ]
    return config

class BamRegionViewer:
    cigarregx = dict(
                        clip_5p      = re.compile( r"^\d+S" ),
                        clip_3p      = re.compile( r"\d+S$" ),
                        maxDeletion  = re.compile( r"\d+D" ),
                        maxInsertion = re.compile( r"\d+I" ) 
                    )
    
    def __init__( self, config, resetHP=True, minCov=3, minFreq=0.1, logger=None, logFile='bamloader.log', callMode='wgs', verbose=False ):
        self.config    = config
        self.reference = config[ 'reference' ]
        self.region    = config[ 'region' ]
        self.genes     = config[ 'genes' ]
        self.haptag    = config[ 'phaseTag' ]
        self.minClip   = config[ 'minClip' ]
        self.maxClip   = config[ 'maxClip' ]
        self.minCov    = minCov
        self.minFreq   = minFreq
        self.callMode  = callMode
        self.hpindex   = self._loadHP( config[ 'homopolymers' ] )
        self.resetHP   = resetHP
        self.vcf       = None
        self.log       = logger if logger is not None else getLogger( self.__repr__(), logFile, stdout=verbose )
        self._reset()

    def __repr__( self ):
        return f'BamRegionViewer: {self.region}: ( { " / ".join( [ g for g in self.genes.keys() ] ) } )'

    def _reset( self ):
        # initialize empty properties
        self.sampleBam, self.sampleVars, self.readMeta, self.pileup = None, None, None, None

    def loadSample( self, sampleBam, vcf=None ):
        #load bam pileup  and remove strand info from sample variants
        # load non-primary alignments separately so overlapping alignments don't blend
        #self.log.info( f'Loading {sampleBam}' )
        # print ("loadSample")
        self._reset()
        self.sampleBam = sampleBam
        sampledf = pd.concat( [ self._makeDf( sampleBam, self.reference, self.region, varsOnly=False, qryCol='hifi_read', **kws )\
                                    .assign( is_primary=primary )
                                for primary,kws in zip( [ True, False ], [ { 'exclude':0xD00 }, { 'exclude': 0x400, 'require':0x900 } ] ) ] )
        # check that we found some data
        if sampledf.empty:
            raise BamViewer_Error( 'No Data!' )
        sampledf.VAR = sampledf.VAR.str.replace( '.', ',', regex=False )
        # get some read metadata, incl secondary/supp alignments
        readMeta     = self._getReadMeta( sampleBam ) 
        
        # flag coverage fractions
        for gene, ( gstart, gstop ) in self.genes.items():
            genesize = gstop - gstart
            overlap  = ( readMeta.rstop.clip( upper=gstop ) - readMeta.rstart.clip( lower=gstart ) ).clip( lower=0 )
            readMeta[ f'cov_{gene}' ] = ( overlap / genesize ).replace( np.inf, 0 ).round( 3 )

        self.sampleVars = sampledf
        self.readMeta   = readMeta
        if self.resetHP:
            self.readMeta.HP = str( -1 )
        self.pileup     = self._pileup()
        if vcf is not None:
            self.vcf = Vcf( vcf, region=self.region, geneRegion=self.config[ 'coreRegion' ] )
        # check for low coverage
        coverage = self.pileup.loc[ :, slice( *self.genes[ 'CYP2D6' ] ) ].notnull().sum()
        if coverage.mean() < self.config[ 'lowCoverageWarning' ][ self.callMode ]:
            self.log.warning( f'LOW COVERAGE! MeanCov over gene: {coverage.mean():.3}' )

    def _loadHP( self, hpbed, tol=1 ):
        homopoly = pd.read_csv( hpbed, sep='\t', comment='#', names=[ 'CHR', 'start', 'end', 'label' ] ) 
        return pd.IntervalIndex.from_arrays( homopoly.start - tol, homopoly.end + tol, closed='both' )

    def _checkhp( self, p ):
        try:
            i = self.hpindex.get_loc( p )
            return True
        except KeyError:
            return False

    def fetch( self, indices, start=None, stop=None, positions=None, dropnull=True ):
        ''' Fetch a subset of the pileup using row indices and optional positional bounds.

            indices: vector of (hifi_read,is_primary) indices
            start/stop/positions:  use either start/stop or a vect of positions [optional]
            dropnull: remove reads/positions where all elem are null, if a 2-bool list passed, 
                      positions relate to rxc [ dropnull rows, dropnull cols ]
            Returns: dataframe reads x postion with elem = pileup vals'''
        if positions is None:
            res = self.pileup.loc[ indices, slice( start, stop ) ]
        else:
            res = self.pileup.loc[ indices ].reindex( columns=positions ).dropna( axis=1, how='all' )
        if dropnull == True: # drop both
            r, c = True, True
        elif dropnull == False:
            return res
        else: # assume 2-mer of bool is passed
            r, c = dropnull
        for ax,val in enumerate( [ r, c ] ):
            if val:
                res = res.dropna( axis=ax, how='all' )
        return res

    def _makeDf( self, bamfile, reference, 
                region, truncate=True, varsOnly=True, 
                exclude=0, require=0, minBaseQ=0, 
                qryCol='query_name', dropAllRef=True ):
        bam = pysam.AlignmentFile( bamfile, 'r' )
        ref = pysam.FastaFile( reference )
        #refcall = list(',.') if dropAllRef else ''
        # load pileup
        df  = pd.DataFrame( { 
                              ( column.reference_name,
                                column.reference_pos  )  : dict( zip( column.get_query_names(),
                                                                 column.get_query_sequences( mark_matches=True,
                                                                                             add_indels=True ) ) )
                            for column in bam.pileup( flag_filter=exclude,
                                                      flag_require=require,
                                                      fastafile=ref,
                                                      region=region,
                                                      truncate=truncate,
                                                      min_base_quality=minBaseQ,
                                                      compute_baq=False)
#                            if ~pd.Series( column.get_query_sequences(mark_matches=True,add_indels=True )).isin( refcall ).all()
                            } )\
                .unstack()\
                .rename( 'VAR' )
        if len( df ):
            df = df.str.upper() # remove strand info (case)
            df.index.names = [ 'contig', 'pos', qryCol ]
            if varsOnly:
                df = df[ (df != '.') & (df != ',') & (df != '*') & (df.notna()) ]
            return df.reset_index()
        else:
            return pd.DataFrame()

    def _pileup( self ):
        return self.sampleVars.set_index( ['hifi_read','is_primary','pos'] )\
                              .VAR.replace( ',', 'ref' )\
                              .unstack( 'pos' )

    @staticmethod
    def _parseCigar( cigar ):
        return { 
                name : max( [ int( m[:-1] ) for m in patt.findall( cigar ) ], default=0 )
                 for name,patt in BamRegionViewer.cigarregx.items() 
               }

    def _getReadMeta( self, bam ):
        # define what subset of data we need here
        def makeRow( rec ):
            data = dict(
                        hifi_read  = rec.query_name,
                        is_primary = not ( rec.is_secondary | rec.is_supplementary ), 
                        rstart     = rec.reference_start,
                        rstop      = rec.reference_end,
                        qlen       = rec.query_length
                   )
            data.update( self._getHapTag( rec ) )
            data.update( self._parseCigar( rec.cigarstring ) )
            data.update( self._getClipped( rec, data ) )
            return data

        recgen = filter( lambda rec: not 0x400 & rec.flag, pysam.AlignmentFile( bam, 'r' ).fetch( region=self.region ) )
        idxCol = [ 'hifi_read', 'is_primary' ]
        return  pd.DataFrame( map( makeRow, recgen ) )\
                  .drop_duplicates( subset=idxCol, keep='first' )\
                  .set_index( idxCol )
    
    def _getClipped( self, rec, data ):
        # set a null return to generate columns in cases with no clips
        res = { f'clipSeq_{i}p' : None for i in [ 3, 5 ] }
        if rec.flag & 0x900:
            return res # no need to re-map primary
        if data[ 'clip_5p' ] >= self.minClip:
            # use maxClip to minimize compute and allow for calling on assemblies
            start = rec.qstart - min( data[ 'clip_5p' ], self.maxClip )
            res[ 'clipSeq_5p' ] = rec.query_sequence[ start : rec.qstart ]
        if data[ 'clip_3p' ] >= self.minClip:
            stop  = rec.qend + min( data[ 'clip_3p' ], self.maxClip )
            res[ 'clipSeq_3p' ] = rec.query_sequence[ rec.qend : stop ]
        return res

    def updateHaplotypes( self, idxs, method='svt', label=None, round=1, candidateSubset=None, window=None, offset=0 ):
        methods = { 'svt' : self._svt }
        if idxs.empty:
            return None
        # reset values
        lbl = f'{ len( idxs ) } { label + " " if label else "" }'
        self.log.debug( f'Reseting HP for {lbl} reads' )
        self.readMeta.loc[ idxs, 'HP' ] = label #None
        if label == 'hybrid':
            clusters = self._split_hybrids( idxs )
            if clusters is not None:
                clusters += offset
                self.readMeta.update( label + '_' + clusters.astype( str ) )
            return clusters
        else:
            # use candidate subset for getting variant pos if defined
            vIdx = idxs if candidateSubset is None else candidateSubset
            if self.vcf is not None:
                varPos = self._getVariablePosFromVCF( vIdx, minGroups=2, window=window )
            else:
                varPos = self._getVariablePos( vIdx, minGroups=2, minEntropy=0.2, window=window )
            if varPos.empty:
                self.log.debug( f'No variable positions found for {lbl}reads. Exiting updateHaplotypes' )
                return None
            self.log.debug( f'Positions used for read clustering: {list( map( int, varPos.index ) )}' )
            clusters = methods[ method ]( varPos ) 
        # validate result
        if self.validateClusters( clusters, varPos, label=label, round=round ):
            clusters += offset
            self.readMeta.update( label + '_' + clusters.astype( str ) )
            return clusters
        return None  

    def _split_hybrids( self, idxs ):
        '''First split by hybrid, if different, then Only split same-hybrids if they are in tandem or mode=consensus'''
        reads  = self.readMeta.loc[ idxs ]
        offset = 0 #for unique cluster numbers
        res    = []
        for lbl,rds in reads.groupby( 'SVlabel' ):
            tandemReads = rds.hybrid_maploc.str.contains('\+')
            if tandemReads.sum() >= 2:
                self.log.info( f'Found {tandemReads.sum()} reads with direct support for tandem hybrid' )
                doubles = pd.Series( np.random.randint( 0, 2, tandemReads.sum() ), index=rds.index[ tandemReads ], name='HP' )
                singles = pd.Series( rds[ ~tandemReads ].hybrid_maploc.factorize( sort=True )[0], index=rds.index[ ~tandemReads ], name='HP' )
                res.append( pd.concat( [ singles, doubles ] )  + offset )
                offset += 2
            else:
                if self.callMode == 'consensus' and len(rds) > 1:  # assume each pbaa consensus is unique
                    incr = len( rds )
                    vals = [ i + offset for i in range(incr) ]
                else:
                    incr = 1
                    vals = offset
                res.append( pd.Series( vals, index=rds.index, name='HP' ) )
                offset += incr
        return pd.concat( res )


    def validateClusters( self, clusters, varPos, minFrac=0.9, label=None, round=1 ):
        "Checks on cluster outputs"
        # check coverage
        reads     = self.readMeta.loc[ clusters.index ]
        rstart    = reads.rstart.min()
        rstop     = reads.rstop.max()
        covRegion = rstop - rstart
        coverage  = pd.DataFrame( np.zeros( ( len( clusters ), covRegion ), dtype=bool ),
                                  index=clusters.index, columns=range( rstart, rstop ) )
        for i, row in reads.iterrows():
            coverage.loc[ i, row.rstart : row.rstop ] = True
        covsum = coverage.groupby( clusters ).sum()
        # coverage over gene region 
        gstart, gstop = self.config[ 'genes' ][ 'CYP2D6' ]
        uncovered = ( ( covsum.loc[ :, gstart:gstop ] >= self.minCov ).sum( axis=1 ) / ( gstop - gstart ) ) < minFrac 
        if uncovered.any() and label != 'deletion':
            # reject haplotypes if both subsets don't cover minfrac of the gene
            self.log.debug( f'Rejecting phasing because {uncovered.sum()} group(s) cover < {minFrac} of the gene at >={self.minCov} reads' )
            return False
        #coverage at variant positions used to cluser
        # if either cluster covers no positions used, then reject
        uncovVarpos = ( covsum[ varPos.index ] < self.minCov ).all( axis=1 )
        if uncovVarpos.any():
            self.log.debug( f'Rejecting phasing due to {uncovVarpos.sum()} uncovered variants in one/both groups' )
            return False
        # reject tandem-dup splits unless we have minCov reads with SV signatures in each subset
        #if label == 'duplicate' and round==2:
        if label == 'duplicate' and round==2:
            counts = self.readMeta.reindex( clusters.index )\
                         .query( 'SV=="duplicate"' )\
                         .groupby( clusters ).size()
            if len( counts ) == 1 or ( counts < self.minCov ).any():
                self.log.debug( f'Rejecting dup phasing because at least one subset has <{self.minCov} reads with SV signature' )
                return False
        return True

    def _svt( self, varPos, rank=2, clusterSeed=42 ):
        # shrinkage operator
        def D( M, r=2 ):
            U, s, Vh = svd( M )
            Ur = U[ :, :r ]
            Sr = diagsvd( s[ :r ], r, Vh.shape[ 0 ] )
            return Ur @ Sr @ Vh.T
        # observed matrx R in { ref: 1, alt:-1, nodata:0 }
        cats = { 'ref': 1, 'nc': 0 }
        R = pd.concat([ p.calls.map(cats).where( p.calls.isin( cats.keys() ), -1 ) for p in varPos.sort_index() ], axis=1 ).fillna(0)
        # include dummy constant feature
        R = pd.concat( [ R, pd.Series( 1, index=R.index, name='dummy' ) ], axis=1 )
        # down project to rank
        shrunk = D( R, rank )
        # get prediction from kmeans clustering of filled (shrunken) matrix
        #prediction = KMeans( n_clusters=rank, random_state=clusterSeed ).fit( shrunk ).predict( shrunk ) 
        prediction = np.argmin( KMeans( n_clusters=rank, random_state=clusterSeed ).fit_transform( shrunk ), axis=1 )
        return pd.Series( prediction, name='HP', index=R.index )

    def _filter_pcolumns( self, minGroups=2, maxIndel=50 ):
        '''Filters: 
            variable pos ; no 1-bp indel next to hp ; 
            avoid low-complexity rgn ; indel less than maxIndel ;
            no 1-indel < 30% of reads (all), or < minfrac for amp/capture '''
        def filt( p ):
            return   bool(
                               ( p._ngrps >= minGroups ) \
                            & ~( ( p.vtype == '1indel' ) & ( self._checkhp( p.pos ) ) )\
                            & ~( ( p.pos >= 42132023 ) & ( p.pos <= 42132049 ) ) \
                            & ~( ( p.vtype in ['ins','del'] ) & ( p.vlen > maxIndel ) )\
                            & ~( ( p.vtype == '1indel' ) & ( p.altFrac < 0.30  ) )\
                            & ( p.altFrac >= self.minFreq )
                         )
        return filt
    
    def _getVariablePosFromVCF( self, idxs, minGroups=2, window=None, maxIndel=50 ):
        positions = self.vcf.positions if window is None \
                    else self.vcf.positions[ self.vcf.positions.between( *window, inclusive='both' ) ]
        varPos    =   self.fetch( idxs, positions=positions, dropnull=[ False, True ] )\
                          .apply( lambda c: PileupColumn( c, minCov=self.minCov, minFreq=self.minFreq, dropna=False ) )
        return varPos[ varPos.map( self._filter_pcolumns( minGroups=minGroups, maxIndel=maxIndel ) ) ]

    def _getVariablePos( self, idxs, minGroups=2, window=None, maxIndel=50, minEntropy=0 ):
        ' return vect of Position objects where >= ( minCov - 1 ) reads are non-reference calls '
        # refcall = 'ref'
        start,stop = ( None,None ) if window is None else window
        pileup = self.fetch( idxs, start=start, stop=stop, dropnull=[ False, True ] )
        nonRef = ( ( pileup != 'ref' ) & ( pileup != '*' ) & pileup.notnull() ).sum( axis=0 ) >= self.minCov
        pcolumn = pileup.loc[ :, nonRef ].apply( lambda c: PileupColumn( c, minCov=self.minCov, minFreq=self.minFreq, dropna=False ) )
        return pcolumn[ pcolumn.map( self._filter_pcolumns( minGroups=minGroups, maxIndel=maxIndel ) ) ]

    def _getHapTag( self, rec ):
        try:
            res = rec.get_tag( self.haptag )
        except KeyError: #not set
            res = 'notPhased'
        return { 'HP' : str( res ) }

    def getConsensusVariants( self, idxs, positions=None, start=None, stop=None, suppressWarnings=False ):
        if positions is not None:
            columns=positions
        elif start is not None and stop is not None:
            columns=self.pileup.columns[ pd.Series( self.pileup.columns ).between( start, stop ) ]
        else:
            raise BamViewer_Error( 'Must pass either vector of positions or start/stop' )
        pile      = self.pileup.reindex( index=idxs, columns=columns )
        coverage  = pile.notnull().sum()
        counts    = pile.apply( lambda p: p.value_counts() )
        plurality = counts.apply( lambda c: c.idxmax() ) if len( counts ) \
                    else pd.Series( None, index=counts.columns, dtype=int )
        support   = counts.apply( lambda c: c.max() ).fillna( 0 ).astype( int )
        variants  = plurality[ plurality != 'ref' ]
        res       = pd.concat( [ coverage, variants, support ], axis=1 ).fillna( '.' )
        res.columns = [ 'coverage', 'VAR', 'support' ]
        # some warnings
        dropped = res.coverage == 0
        lowcov  = res.coverage < self.minCov
        if dropped.any() and not suppressWarnings and not self.callMode == 'consensus':
            self.log.warn( f'{dropped.sum()} bases in variant region with no coverage' )
        if lowcov.any() and not suppressWarnings and not self.callMode == 'consensus':
            self.log.warn( f'{lowcov.sum()} bases in variant region with coverage < {self.minCov}' )
        return res 

    def exportSubset( self, idxs, outbam, updateHPtag=True, labelMap=None, colorMap=None ):
        indices = self.readMeta.index if idxs is None else idxs
        lblmap = ( lambda x: 'No Label' ) if labelMap is None else ( lambda x: labelMap[ x ] )
        with pysam.AlignmentFile( self.sampleBam, 'r' ) as ibam:
            with pysam.AlignmentFile( outbam, 'wb', template=ibam ) as obam:
                for rec in ibam.fetch( region=self.region ):
                    ix = ( rec.query_name, not bool( rec.flag & 0x900 ) )
                    if ix in indices:
                        if updateHPtag:
                            hp  = self.readMeta.loc[ ix ].HP
                            rec.set_tag( 'HP', lblmap( hp ), 'Z' )
                            if colorMap is not None:
                                rec.set_tag( 'YC', colorMap[ hp ], 'Z' )
                        obam.write( rec )
        pysam.index( outbam )
        return outbam

    def _jaccardSimilarity( self, positions=None ):
        pos = positions if positions is not None else slice( None )
        def jaccard( idxs ):
            paired = self.pileup.reindex( idxs )[ pos ] \
                           .dropna(axis=1).apply( lambda p: p[0]==p[1] )
            if paired.shape == ( 2, 0 ): #no overlap
                return None
            else:
                return paired.sum() / len( paired )
        return jaccard
    
class PileupColumn:
    def __init__( self, callVector, minCov=3, minFreq=0.1, refCall='ref', chrom=None, dropna=True ):
        self.chr     = chrom
        self.pos     = callVector.name
        self.minCov  = minCov
        self.minFreq = minFreq
        self.refCall = refCall
        self.dropna  = dropna
        self._getCalls( callVector )

    def __repr__( self ):
        return str( self._counts.to_dict() )

    def __lt__( self, other ):
        ' For sorting by entropy '
        return self.entropy < other.entropy

    def __contains__( self, read ):
        return read in self.calls.index

    def _getCalls( self, callVector ):
        nocall = 'nc'
        def _kind( variant ):
            if variant == self.refCall: return 'refCall'
            if variant == nocall: return 'nocall'
            if len( variant ) == 4: return '1indel' # single-base indel
            o = variant.strip(',')[0]
            if o == '-': return 'del'
            if o == '+': return 'ins'
            return 'snp'

        calls = callVector.dropna() if self.dropna else callVector.fillna( nocall )

        self.calls   = calls
        self._counts = calls.value_counts()
        alts         = self._counts[ self._counts >= self.minCov ].drop( [ self.refCall, nocall, '*' ], errors='ignore' )
        self._ngrps  = int( self._counts.get( self.refCall, 0 ) >= self.minCov ) + len( alts )
        self.vtype   = None if alts.empty else _kind( alts.idxmax() )
        self.altFrac = 0 if alts.empty else alts.max() / self._counts.drop( nocall, errors='ignore' ).sum() 
        self.vlen    = sum( 1 for b in re.finditer( '[AGCT]', alts.idxmax(), re.I ) ) if self.vtype in ['ins','del'] else 1
        self.reads   = { vnt : rows.index for vnt,rows in calls.groupby( calls ) }

    @property
    def entropy( self ):
        return entropy( self._counts[ self._counts >= self.minCov ] )

    @property
    def bifrac( self ):
        ' Fraction of total reads represented by two most common groups '
        return -1 if len( self._counts ) < 2 \
                 else self._counts[ :2 ].sum() / self._counts.sum()

    @property
    def coverage( self ):
        return self._counts.drop( 0 ).sum()

class Vcf:
    def __init__( self, vcf_file, region=None, geneRegion=None, toZeroIdx=True ):
        self.vcf        = vcf_file
        self.region     = region
        self.geneRegion = geneRegion
        self.toZeroIdx = toZeroIdx
   
        with pysam.VariantFile( vcf_file ) as inVcf:
            recGen = inVcf if region is None else inVcf.fetch( region=region )
            self.positions = pd.Series( [ rec.pos - int( toZeroIdx )
                                          for rec in recGen ] )
        self.geneRegion = self.positions if geneRegion is None \
                          else self.positions[ self.positions.between( *self.geneRegion, inclusive='both' ) ]

    def __contains__( self, item ):
        return item in self.positions

class Config_Error( Exception ):
    pass

class BamViewer_Error( Exception ):
    pass
