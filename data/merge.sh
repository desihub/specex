#!/bin/bash


elements=`ls nist_*.txt | sed 's/nist_//' | sed 's/\.txt//'`
#elements=`echo HI`

tmp="tmp.txt"
rm -f $tmp

for element in $elements ; do

    echo $element
    file=nist_$element.txt
 
    cat $file | grep -v '\-\-\-' | sed 's/ //g' | sed 's/||/|0|/g' | grep  '^.[0-9]*' | sed 's/|/ /g' | sed 's/(//g' | sed 's/)//g' | sed 's/\*//g' | awk '{if(($1-$2)<10 && ($1-$2)>-10 && $1>3000. && $1<100000.) printf("%09.4f %g '$element'\n",$1,$3)}' >> $tmp
    

   

done

cat $tmp | sort -n -k 1 | grep -v "HgII" | awk '{if($2>100) print $0}'> lamplines.txt

