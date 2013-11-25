#!/bin/bash
#
# Masai reference output generation. Version v0.6.1 was used.

PATH=~/Code/seqan-builds/Release-Gcc/bin/
INDEXER=$PATH/cuda_indexer
MAPPER=$PATH/cuda_mapper

# ============================================================
# Run Indexer
# ============================================================

# Run with different organisms.
for i in celegans; do #hg18
    ${INDEXER} datasets/$i/genome.fasta &> logs/indexer.$i.stdout
done

# ============================================================
# Run Single-End Mapper
# ============================================================

MAPPER_ARGS=("--no-cuda --threads 1" "--no-cuda --threads 8" "")

# Run with different organisms.
for i in celegans hg18; do
    # Run with different seed length.
#    for sl in 20 33; do
        # Run on CPU with 1 thread,  on CPU with 8 threads, on GPU.        
        for args in "${MAPPER_ARGS[@]}"; do
            ${MAPPER} datasets/$i/genome.fasta datasets/$i/reads.fastq ${args} &> logs/mapper.$i.stdout
        done
#    done
done
