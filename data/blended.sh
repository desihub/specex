#!/bin/bash

echo "looking for indicernable lines in list"
dist=0.5 # angs

ifile=lamplines-specex.par

cat $ifile | grep ^arclineid | awk 'BEGIN{l1=0}{l2=$2;if((l1-l2)<'$dist' && (l1-l2)>'-$dist') print "PB WITH ",l1,l2,"USE",(l1+l2)/2. ;l1=l2}'
