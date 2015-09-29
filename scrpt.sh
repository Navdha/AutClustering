#!/bin/bash
while read -r line
do 
for((c=1;c<="$1";c++))
do 
eval "$line"
done
done<coms.txt;
