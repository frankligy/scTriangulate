#BSUB -W 15:00
#BSUB -M 196000
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -J scTriangulate
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./run.py










