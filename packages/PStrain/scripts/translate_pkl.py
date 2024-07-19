"""
convert the *.pkl file in the Metaphlan database into human-readable file.


Example: 
python ~/softwares/PStrain/scripts/translate_pkl.py mpa_vOct22_CHOCOPhlAnSGB_202403.pkl mpa_vOct22_CHOCOPhlAnSGB_202403.species_markers.txt

Output is like:

UniRef90_UPI000E65AEFC|1__27|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E646519|1__24|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E647771|1__22|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64864D|1__22|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E646DCF|1__21|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E65281B|1__21|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E65A109|1__21|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E646212|1__20|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E6514C5|1__19|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E653299|1__18|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64897A|1__18|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64FDA9|1__18|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E653B3C|1__18|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E649898|1__18|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E648AF2|5__21|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64B60E|1__17|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64636E|1__17|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64D6B9|1__17|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E646A07|1__17|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E65E13F|11__26|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E659221|1__16|SGB32561 s__Calidifontimicrobium_sediminis
UniRef90_UPI000E64ADB9|1__16|SGB32561 s__Calidifontimicrobium_sediminis
"""

import pickle
import bz2
import sys


def read_metaphlan_pkl(metaphalan_database_pkl, output_marker_gene_species):


    f = open(output_marker_gene_species, 'w')
    with bz2.BZ2File(metaphalan_database_pkl, 'r' ) as a:
        mpa_pkl = pickle.load( a )

    for k, p in mpa_pkl['markers'].items():
        gene = k
        species = p['taxon'].split("|")[6]
        print (gene, species, file = f)
    
    f.close()


if __name__ == "__main__":
    metaphalan_database_pkl = sys.argv[1]
    output_marker_gene_species = sys.argv[2]
    read_metaphlan_pkl(metaphalan_database_pkl, output_marker_gene_species)