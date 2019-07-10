#!/bin/bash

cd `dirname $0`

###########################################################################

GOLD_STANDARDS=$(cd gold_standard; ls *.txt | sed -e 's/.txt//g' | sort)
VERSION=0.93

###########################################################################

IMPL=impl/cluster_one-${VERSION}.jar
PYTHON=`which python`
JAVA=`which java`

echo "Using ClusterONE ${VERSION}..."
echo ""

for GOLD_STANDARD in ${GOLD_STANDARDS}; do
    echo "Using gold standard: ${GOLD_STANDARD}"
    echo ""

    for dataset in datasets/*.txt; do
	    echo `basename $dataset .txt`
	    ${JAVA} -jar ${IMPL} $dataset >results.txt 2>/dev/null
	    ${PYTHON} scripts/match_standalone.py -q -n ${dataset} \
		    -m frac -m cws -m ppv -m acc -m mmr \
		    gold_standard/${GOLD_STANDARD}.txt results.txt | awk '{ printf "  "; print $0 }'
	    echo ""
	    rm results.txt
	done
done
