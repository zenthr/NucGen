#!/bin/bash
#SBATCH --job-name=DRAW\!
#SBATCH --output=stdout.%j
#SBATCH --account=122804261495
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --mail-type=NONE
#SBATCH --mail-user=zenthrelin@gmail.com

./FFTW.exe