#!/bin/bash

if [ "$1" = "" ]; then
    echo "Usage: bash generate_rna_filetypes.sh /path/to/RNAStructure/directory /path/to/sequence/file_seq.txt [seed]"
    echo "    Seed defaults to 1234 if not set"
    echo "    This script assumes that RNAprofile is in a parallel folder"
    echo "    For example: ~/testing/ipam-wbio-scripts/generate_rna_filetypes.sh and ~/testing/src/RNAprofile"
    echo "    It creates a test_files folder which contains subfolders per sequence and per seed"
    exit
fi

# Get the directory of this script, rather than the current working directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo Sequence:
echo $2

SEED=$3
if [ "$3" = "" ]; then
    SEED=1234
fi

BASE_NAME=$(basename $2)
FILE_NAME=${BASE_NAME%.*}
FILE_NAME=${FILE_NAME%_seq*}

# make the output folder for the results of sampling and profiling
mkdir test_files
mkdir test_files/$FILE_NAME
mkdir test_files/$FILE_NAME/$SEED

# Run RNAStructure to sample an ensemble of structures using the given seed
RNAStructure=$1
export DATAPATH=$RNAStructure/data_tables
$RNAStructure/exe/stochastic $2 test_files/$FILE_NAME/$SEED/$FILE_NAME.ct --sequence --seed $SEED

# Convert the RNAStructure output from .ct to .gtboltz files
python3 $DIR/RNAStructure_to_gtboltzmann.py test_files/$FILE_NAME/$SEED/$FILE_NAME.ct test_files/$FILE_NAME/$SEED/$FILE_NAME.gtboltz

# Adjust the sequence file location to work properly when calling RNAprofile, if it is a relative path
ADJUSTED_LOC=../../../$2
if [[ "$2" = /* ]]; then
    ADJUSTED_LOC=$2
fi

cd test_files/$FILE_NAME/$SEED/
$DIR/../src/RNAprofile -v -o output -e $FILE_NAME.gtboltz $ADJUSTED_LOC > output.txt
