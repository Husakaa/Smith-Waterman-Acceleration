#!/bin/bash

rm ./times.txt
for i in $(seq 1 10); do
	echo "" >> times.txt
	echo "Tiempo $i para Secuencial" >> times.txt
	./Smith_Waterman_fasta.exe resources/Homo_sapiens.fasta resources/Pan_troglodytes.fasta >> times.txt
	echo "" >> times.txt
	echo "Tiempo $i para OpenMP" >> times.txt
	./Smith_Waterman_Cpu.exe resources/Homo_sapiens.fasta resources/Pan_troglodytes.fasta >> times.txt
    echo "Tiempo $i para CUDA" >> times.txt
	./Smith_Waterman_Cuda.exe resources/Homo_sapiens.fasta resources/Pan_troglodytes.fasta 128 >> times.txt
	echo "" >> times.txt
done
