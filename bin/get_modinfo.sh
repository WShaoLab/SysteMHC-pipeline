#!/bin/sh


grep -v '#' $1 | awk '{if (length($1)>=8 && length($1)<=14) print $0}' | grep '/' | awk -F '|' '{print $1}'> modify.peps;

grep -v '#' $1 | awk '{if (length($1)>=8 && length($1)<=14) print $0}' | grep '/' | awk -F '|' '{print $2}' | awk -F ',' '{print $2,$3}' | awk -F '/' '{print $1}'  > modify.name;

grep -v '#' $1 | awk '{if (length($1)>=8 && length($1)<=14) print $0}' | grep '/' | awk -F '|' '{print $2}' | awk -F ',' '{print $1}' | awk -F '/' '{print $2}' > modify.location;

paste modify.peps modify.name modify.location | awk '{print $1$2,$0,$1$2$3$4$5}' | sort -k 5 | uniq | awk '{print $1,$2,$3,$4,$5,$6}'> modify.csv;

grep -v '#' $1 | awk '{if (length($1)>=8 && length($1)<=14) print $0}' | grep -v '/' | awk -F '|' '{print $1}' | awk '{print $1$2,$0}' | sort -k 1 | uniq > unmodify.csv;

grep -v '#' $1 | awk '{if (length($1)>=8 && length($1)<=14) print $1}' | sort | uniq > peptides.csv

rm modify.peps;
rm modify.name;
rm modify.location;
