#BSUB -W 10:00
#BSUB -M 100000
#BSUB -q gpu-v100
#BSUB -gpu "num=1"
#BSUB -J tea
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


./total_vi.py







