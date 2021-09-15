import subprocess

script_path = 'two_step_mutant.py'
cmd = 'bsub -q new-medium -R "rusage[mem=100]" -u /dev/null python ' + script_path


for run in range(100):
    for param in range(72):
        paramstr = ' '.join(map(str,[param, run]))
        _cmd = ' '.join([cmd, paramstr])
        subprocess.run(_cmd, shell=True)