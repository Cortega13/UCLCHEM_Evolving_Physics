import os

working_path = "/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem/"#"/work/09338/carlos9/vista/UCLCHEM/benchmarks/turbulence/"

tracer_script_path = os.path.join(working_path, 'tracer_chemical_evolution.py')
tracer_folder_path = os.path.join(working_path, 'data', 'turbulence_tracers_csv')
commandlines_path = os.path.join(working_path, 'commandlines.txt')
finished_folder = os.path.join(working_path, 'data', 'chemistry_tracers_csv')


if not os.path.exists(finished_folder):
    os.makedirs(finished_folder)
tracer_files = set(os.listdir(tracer_folder_path))
finished_files = set(os.listdir(finished_folder))

next_files = tracer_files - finished_files

with open(commandlines_path, 'w') as f:
    for tracer_file in next_files:
        tracer_filepath = os.path.join(tracer_folder_path, tracer_file)
        f.write(f"python {tracer_script_path} {tracer_filepath}\n")
f.close()
print("Command lines written to commandlines.txt")
print("Total number of tracers: ", len(tracer_files))
print("Total number of completed tracers: ", len(finished_files))