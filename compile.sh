#!/bin/bash

# =======================================
# Script de Compilación - Smith-Waterman 
# =======================================

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}[INFO] Iniciando compilación del proyecto...${NC}"

# 1. Compilar Versión Secuencial
# FIX: Se añade -fopenmp porque el código usa omp_get_wtime() para medir tiempos
if [ -f "src/Smith_Waterman.c" ]; then
    echo -e "Compilando versión Secuencial..."
    gcc -O3 -fopenmp src/Smith_Waterman.c -o Smith_Waterman_fasta.exe
    if [ $? -eq 0 ]; then echo -e "${GREEN}✔ Secuencial compilado correctamente.${NC}"; else echo -e "${RED}✘ Error compilando Secuencial.${NC}"; fi
else
    echo -e "${RED}⚠ No se encontró src/Smith_Waterman.c${NC}"
fi

# 2. Compilar Versión OpenMP
if [ -f "src/Smith_Waterman_OpenMP.c" ]; then
    echo -e "Compilando versión OpenMP..."
    gcc -O3 -fopenmp src/Smith_Waterman_OpenMP.c -o Smith_Waterman_Cpu.exe
    if [ $? -eq 0 ]; then echo -e "${GREEN}✔ OpenMP compilado correctamente.${NC}"; else echo -e "${RED}✘ Error compilando OpenMP.${NC}"; fi
else
    echo -e "${RED}⚠ No se encontró src/Smith_Waterman_OpenMP.c${NC}"
fi

# 3. Compilar Versión CUDA
if command -v nvcc &> /dev/null; then
    if [ -f "src/Smith_Waterman_Cuda.cu" ]; then
        echo -e "Compilando versión CUDA..."
        nvcc -O3 src/Smith_Waterman_Cuda.cu -o Smith_Waterman_Cuda.exe
        if [ $? -eq 0 ]; then echo -e "${GREEN}✔ CUDA compilado correctamente.${NC}"; else echo -e "${RED}✘ Error compilando CUDA.${NC}"; fi
    else
        echo -e "${RED}⚠ No se encontró src/Smith_Waterman_Cuda.cu${NC}"
    fi
else
    echo -e "${RED}[SKIP] Compilador nvcc no detectado o no compatible. Saltando versión GPU.${NC}"
fi

echo -e "${GREEN}[INFO] Proceso finalizado.${NC}"
