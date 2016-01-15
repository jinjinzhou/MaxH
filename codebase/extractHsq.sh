#!/bin/sh

#  extractHsq.sh
#  
#
#  Created by Jin Zhou on 1/14/16.
#

for entry in "."/*.hsq
do
hsq=$(grep -r 'V(G)'  $entry)
x='./'
hsq=${hsq//$x/''}
n=$(grep -A 1 'Pval' $entry | grep  'n')
echo $hsq" "$n
done