import subprocess

with open("/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem/commandlines.txt") as f:
    commands = [line.strip() for line in f if line.strip()]

processes = [subprocess.Popen(cmd, shell=True) for cmd in commands]

for p in processes:
    p.wait()