#!/bin/bash
for i in {1..20}
do
   echo "Executing $i program"
   g++ main.cpp -o main && ./main > Experimento3_$i.csv
done