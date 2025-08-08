#!/bin/bash

mkdir fastqc

for i in *fastq.gz; do
	fastqc $i -o ./fastqc
done 



