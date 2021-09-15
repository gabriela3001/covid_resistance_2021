import subprocess

script_path = 'script_immunity_waning.py'
cmd = 'bsub -q new-medium -R "rusage[mem=100]" -u /dev/null python ' + script_path
#cmd = 'bsub -q new-long -R "rusage[mem=1000]" python ' + script_path


for param in range(72):
  for run in range(100):
    paramstr = ' '.join(map(str,[param, run]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)