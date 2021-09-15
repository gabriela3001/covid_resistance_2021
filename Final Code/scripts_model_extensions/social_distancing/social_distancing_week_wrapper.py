import subprocess

script_path = 'social_distancing_week_script.py'
cmd = 'bsub -q new-short -R "rusage[mem=100]" -u /dev/null python ' + script_path

#cmd = 'bsub -q new-long python ' + script_path

for param in range(72):
  for run in range(100):
    paramstr = ' '.join(map(str,[param, run]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)