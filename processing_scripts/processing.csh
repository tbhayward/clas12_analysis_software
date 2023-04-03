#!/bin/csh

pwd
cd clasqaDB/
source env.csh
cd ..
coatjava/bin/run-groovy $1 $2 $3 $4

