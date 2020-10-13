#!/bin/bash
#conda activate Johannes_turing
#for i in {1..8}
#do
cd "./ScriptsToRun/" 
    #mpirun -n 50 python3 mainFile.py
export OPENBLAS_NUM_THREADS=1
python3 mainFile.py 
rm -r __pycache__ 
cd "../" 



