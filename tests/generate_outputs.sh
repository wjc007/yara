#!/bin/bash
#
# Masai reference output generation. Version v0.6.1 was used.

PATH=~/Code/seqan-builds/Release-Gcc/bin
INDEXER=$PATH/cuda_indexer
MAPPER=$PATH/cuda_mapper

# ============================================================
# Run Indexer
# ============================================================

# Run with different organismnisms.
#for organism in celegans hg18; do
for organism in celegans; do
    ${INDEXER} datasets/$organism/genome.fasta -xp datasets/$organism/genome.fm &> logs/indexer.$organism.stdout
done

# ============================================================
# Run Single-End Mapper
# ============================================================

# Run on CPU with 1 thread, on CPU with 8 threads, on GPU.
MAPPER_ARGS=("--no-cuda --threads 1" "--no-cuda --threads 8" "")
LOG_SUFFIX=("t1" "t8" "cuda")

# Run with different organismnisms.
for organism in celegans hg18; do
    # Run with different arguments.
    for ((a=0; a<${#MAPPER_ARGS[@]}; a++)); do
        ${MAPPER} datasets/$organism/genome.fasta -xp datasets/$organism/genome.fm datasets/$organism/reads.fastq ${MAPPER_ARGS[$a]} &> logs/mapper.$organism.${LOG_SUFFIX[$a]}.stdout
    done
done
