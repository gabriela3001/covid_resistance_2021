import subprocess

script_path = 'basic_model.py'
cmd = 'bsub -q new-medium -R "rusage[mem=100]" -u /dev/null python ' + script_path

for run in range(10):
    for param in range(5000):
        paramstr = ' '.join(map(str,[param, run]))
        _cmd = ' '.join([cmd, paramstr])
        subprocess.run(_cmd, shell=True)