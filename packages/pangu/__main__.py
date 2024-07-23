#! /usr/bin/env python

__version__ = '0.2.8'
import sys
import os
import logging
import json
from pathlib import Path
from itertools import repeat
import warnings
from cyp2d6_typer import StarTyper

def main( parser ):
    args = parser.parse_args()

    # silence pandas v2 future warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    # set up output prefix
    path = Path( args.prefix )
    if path.is_dir():
        prefix = lambda name: f'{path / name}' 
    elif path.parent.is_dir():
        prefix = lambda name: f'{path}_{name}'
    else:
        raise CYP2D6_typer_Error( f'Cannot create output for prefix {path}' )

    logfile = prefix( 'caller.log' ) if args.logFile is None else args.logFile

    typer = StarTyper( args.config_yaml, 
                       logFile=logfile, 
                       verbose=args.verbose, 
                       logLevel=args.logLevel,
                       callMode=args.mode,
                       randSeed=args.seed,
                       grayscale=args.grayscale ) 

    # make input gnerators
    if len( args.inBam ) >= 1:
        bamgen = args.inBam
        if args.vcf is not None:
            if len( args.inBam ) > 1:
                raise CYP2D6_typer_Error( '--vcf option only valid for single-bam input' )
            vcfgen = [ args.vcf ]
        else:
            vcfgen = repeat( None )
 
    elif args.bamFofn is not None:
        bamgen = open( args.bamFofn, 'r' ).read().strip().split('\n')
        if args.vcfFofn is not None:
            vcfgen = open( args.vcfFofn, 'r' ).read().strip().split('\n')
        else:
            vcfgen = repeat( None )
    else:
        raise CYP2D6_typer_Error( 'No BAM input' )

    if len( args.inBam ) > 1 or args.bamFofn is not None:
        #phamrCat requires one sample per file
        #use bam input as name
        pharmcat_file = ( lambda bam: prefix( f'{Path( bam ).stem}.pharmcat.tsv'  ) )
    else:
        pharmcat_file = ( lambda bam: prefix( 'pharmcat.tsv' ) )
    
    for bamfile, vcffile in zip( bamgen, vcfgen ):
        typer.run( bamfile, vcffile )
        # export pharmcat callfile
        open( pharmcat_file( bamfile ), 'w' ).write( f'{typer.diplotype.pharmCat()}\n' )
        # export reads
        if args.exportLabeledReads:
            ofile = prefix( 'labeledReads.bam' )
            typer.log.info( f'Exporting {ofile}' )
            typer.bamRegion.exportSubset( None, ofile, 
                                          labelMap=typer.diplotype.callMap,
                                          colorMap=typer.diplotype.colorMap )
    # write report
    with open( prefix( 'report.json' ), 'w' ) as report:
        report.write( json.dumps( typer.report, indent=4 ) )

    return typer

class CYP2D6_typer_Error( Exception ):
    pass

def main_cli():
    import argparse
    import pkgutil
    from datetime import datetime
    # from importlib.resources import files

    now  = datetime.now().strftime( '%Y-%m-%d_%H%M%S' )

    parser = argparse.ArgumentParser( prog='pangu', description='Call CYP2D6 star alleles from HiFi WGS data' )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument( 'inBam', metavar='inBam', nargs='*', type=str, help='Aligned BAM file(s) of HiFi WGS reads' )

    inputp = parser.add_argument_group( 'Input Options' )
    inputp.add_argument( '--vcf', dest='vcf', default=None,
                         help='Vcf file name for single bam input.  Default None' )
    inputp.add_argument( '--bamFofn', dest='bamFofn', default=None,
                         help='Text file of bam file names.  Default None' )
    inputp.add_argument( '--vcfFofn', dest='vcfFofn', default=None,
                         help='Text file of vcf file names in same order as bamFofn.  Default None' )
    # inputp.add_argument( '--config', dest='config_yaml', default=files('pangu.data.CYP2D6').joinpath('CYP2D6.yaml'),
    #                      help='Override installed configuration yaml file' ) 
    inputp.add_argument( '--config', dest='config_yaml', default=f'{sys.path[0]}/data/CYP2D6/CYP2D6.yaml',
                         help='Override installed configuration yaml file' ) 
    inputp.add_argument( '-p','--prefix', dest='prefix', default='./',
                         help='Prefix for output files.  Default cwd' ) 
    inputp.add_argument( '-m','--mode', dest='mode', choices=['wgs','amplicon','capture','consensus'], default='wgs',
                         help='Calling mode by HiFi input type.  Default wgs' ) 
    inputp.add_argument( '-s','--seed', dest='seed', default=42, type=int,
                         help='Seed for random generator.  Default 42' ) 
 
    outputp = parser.add_argument_group( 'Output Options' )
    outputp.add_argument( '-x','--exportLabeledReads', dest='exportLabeledReads', action='store_true', default=False,
                         help='Write labeled reads to output prefix' )
    outputp.add_argument( '-g','--grayscale', dest='grayscale', action='store_true', default=False,
                         help='Use grayscale to annotate haplotypes in export bam.  Default use color' )
    outputp.add_argument( '-v','--verbose', dest='verbose', action='store_true', default=False,
                         help='Print logging info to stdout' )
    outputp.add_argument( '--logFile', dest='logFile', default=None,
                         help='Log file.  Default {prefix}[_/]caller.log' )
    outputp.add_argument( '--logLevel', dest='logLevel', choices=logging._levelToName.values(), default='INFO',
                         help='Logging level. Default INFO' )
    
    # print ("xxxx")
    main( parser )
    # try:
    #     starTyper = main( parser )
    # except CYP2D6_typer_Error as e:
    #     print( f'\nERROR: {e}\n' )
    #     sys.exit( 1 )

main_cli()

