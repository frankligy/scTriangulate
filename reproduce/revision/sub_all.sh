#BSUB -W 10:00
#BSUB -M 250G
#BSUB -n 1
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J revision
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./simulation.py







