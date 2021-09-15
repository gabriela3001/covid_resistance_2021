import subprocess

script_path = 'basic_sensitivity_analysis.py'
cmd = 'bsub -q new-short -R "rusage[mem=100]" -u /dev/null python ' + script_path

paramfiles = ['paramgrid_appearance_sensitivity_analysis_sthreshold.txt',
              'paramgrid_appearance_sensitivity_analysis_death.txt',
              'paramgrid_appearance_sensitivity_analysis_recoveryrate.txt',
              'paramgrid_appearance_sensitivity_analysis_reproductionrate.txt'
]


for p in range(3):
    for r in range(100):
        for pf in paramfiles:
            paramstr = ' '.join(map(str,[p, r, pf]))
            _cmd = ' '.join([cmd, paramstr])
            subprocess.run(_cmd, shell=True)
