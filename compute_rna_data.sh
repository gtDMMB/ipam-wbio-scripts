#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo Sequence:
echo $2

if [ "$2" = "" ]; then
    echo "Usage: generate_rna_filetypes.sh /path/to/RNAStructure/directory /path/to/sequence/file_seq.txt [seed]"
    exit
fi

SEED=$3
if [ "$3" = "" ]; then
    SEED=1234
fi

BASE_NAME=$(basename $2)
FILE_NAME=${BASE_NAME%.*}
FILE_NAME=${FILE_NAME%_seq*}

mkdir test_files
mkdir test_files/$FILE_NAME
mkdir test_files/$FILE_NAME/$SEED

RNAStructure=$1
export DATAPATH=$RNAStructure/data_tables
$RNAStructure/exe/stochastic $2 test_files/$FILE_NAME/$SEED/$FILE_NAME.ct --sequence --seed $SEED

python3 $DIR/RNAStructure_to_gtboltzmann.py test_files/$FILE_NAME/$SEED/$FILE_NAME.ct test_files/$FILE_NAME/$SEED/$FILE_NAME.gtboltz

ADJUSTED_LOC=../../../$2
if [[ "$2" = /* ]]; then
    ADJUSTED_LOC=$2
fi

cd test_files/$FILE_NAME/$SEED/
../../../../src/RNAprofile -v -o output -e $FILE_NAME.gtboltz $ADJUSTED_LOC > output.txt
