import subprocess

script_path = 'basic_model.py'
cmd = 'bsub -q new-short -R "rusage[mem=100]" -u /dev/null python ' + script_path

#runs_missing = [(18, 53), (54, 37), (36, 58), (38, 50)]

for run in range(100):
  for param in range(72):
#for param, run in runs_missing:
    paramstr = ' '.join(map(str,[param, run]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)