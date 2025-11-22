#!/bin/bash

for size in 32 64 128 256 512 1024; do
    echo "Times for $size:" >> times.txt
    for i in $(seq 1 3); do
        ./Smith_Waterman_Cuda.exe resources/Homo_sapiens.fasta resources/Pan_troglodytes.fasta $size >> times.txt
    done
    echo "" >> times.txt # salto de linea
done
