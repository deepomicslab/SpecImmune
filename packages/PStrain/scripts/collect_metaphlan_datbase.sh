#!/bin/bash

###
### Collect the Metaphlan database
### 
###
### Usage:
###   bash collect_metaphlan_datbase.sh -x mpa_vOct22_CHOCOPhlAnSGB_202403 -m 4 -d ../
###   or
###   bash collect_metaphlan_datbase.sh -x mpa_v31_CHOCOPhlAn_201901 -m 3 -d ../
###
### Options:
###   -x        metaphlan database index. e.g., mpa_vOct22_CHOCOPhlAnSGB_202403 <required> 
###   -d        folder to save the metaphlan database. <required> 
###   -m        metaphlan version (3 or 4). <required>
###   -h        Show this message.

help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
    help
    exit 1
fi

while getopts "x:d:m:" opt; do
  case $opt in
    x) prefix="$OPTARG"
    ;;
    d) db="$OPTARG"
    ;;
    m) version="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


db_prefix=${prefix:-mpa_vOct22_CHOCOPhlAnSGB_202403}
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
db_dir=${db:-"$script_dir/../"}
metaphlan_version=${version:-4}


echo "Metaphlan database prefix: " $db_prefix
echo "Metaphlan database folder: " $db_dir
echo "Metaphlan version: " $metaphlan_version
echo "Metaphlan database source: " http://cmprod1.cibio.unitn.it/biobakery$metaphlan_version/metaphlan_databases/$db_prefix.tar


mkdir $db_dir
cd $db_dir

wget http://cmprod1.cibio.unitn.it/biobakery$metaphlan_version/metaphlan_databases/$db_prefix.tar
tar -xvf $db_dir/$db_prefix.tar

bzip2 -d $db_dir/$db_prefix.fna.bz2
bowtie2-build $db_dir/$db_prefix.fna $db_dir/$db_prefix


# script_dir=$(dirname "$0")
# mkdir $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403

# cd $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403

# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202403.tar
# tar -xvf $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403.tar

# bzip2 -d $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403.fna.bz2
# bowtie2-build $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403.fna $script_dir/../mpa_vJun23_CHOCOPhlAnSGB_202403



# script_dir=$(dirname "$0")

# mkdir $script_dir/../mpa_v31_CHOCOPhlAn_201901

# cd $script_dir/../mpa_v31_CHOCOPhlAn_201901

# wget http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v31_CHOCOPhlAn_201901.tar
# tar -xvf $script_dir/../mpa_v31_CHOCOPhlAn_201901.tar

# bzip2 -d $script_dir/../mpa_v31_CHOCOPhlAn_201901.fna.bz2
# bowtie2-build $script_dir/../mpa_v31_CHOCOPhlAn_201901.fna $script_dir/../mpa_v31_CHOCOPhlAn_201901





