#!/bin/bash
# This program runs ARACNE-AP (https://sourceforge.net/p/aracne-ap/wiki/Home/)

# calculate threshold with a fixed seed and saves in a new directory named "outputFolder"
java -Xmx5G -jar Aracne.jar -e $DATA/expression_lung.txt  -o outputFolder --tfs $DATA/regulators.txt --pvalue 1E-8 --seed 1 --calculateThreshold

# run 100 reproducible bootstraps and saves in the output folder
for i in {1..100}; do
    java -Xmx5G -jar Aracne.jar -e  $DATA/expression_lung.txt  -o outputFolder --tfs $DATA/regulators.txt --pvalue 1E-8 --threads 2 --seed $i
done

# consolidate bootstraps in the output folder
java -Xmx5G -jar Aracne.jar -o outputFolder --consolidate
