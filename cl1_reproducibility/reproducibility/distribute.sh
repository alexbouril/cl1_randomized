#!/bin/bash

cd `dirname $0`
CWD=`pwd`
cd ..
BASENAME=`basename ${CWD}`

rm -f $BASENAME/cl1_$BASENAME.zip
zip -r $BASENAME/cl1_$BASENAME.zip $BASENAME --exclude $BASENAME/scripts/mwmatching.pyc