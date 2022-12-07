#BSUB -W 15:00
#BSUB -M 196000
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -J scTriangulate
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./run1.py










