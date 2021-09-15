import subprocess

script_path = 'script_realdata_vax_hesitancy.py'
cmd = 'bsub -q new-short -R "rusage[mem=100]" -u /dev/null python ' + script_path


for param in range(96):
  for run in range(100):
    paramstr = ' '.join(map(str,[param, run]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)