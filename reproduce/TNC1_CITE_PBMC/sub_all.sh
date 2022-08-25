#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 100000
#BSUB -J "total_vi"
#BSUB -R "span[hosts=1]"
#BSUB -q "gpu-v100"
#BSUB -gpu "num=1"
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# ./run.py
sleep 5
#./total_vi.py










