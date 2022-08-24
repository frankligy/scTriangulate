#BSUB -W 20:00
#BSUB -M 196000
#BSUB -n 4
#BSUB -J run_python
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./run.py







