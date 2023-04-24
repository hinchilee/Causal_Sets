#PBS -lwalltime=23:59:00
#PBS -lselect=1:ncpus=32:mem=64gb

module load anaconda3/personal
python3 $HOME/Causal_Sets/min_time.py $PBS_O_WORKDIR 