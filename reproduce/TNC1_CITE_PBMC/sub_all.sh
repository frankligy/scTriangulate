#!/bin/bash

#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 100000
#BSUB -J TNC1
#BSUB -R "span[hosts=1]"
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

./run.py











