import numpy as np 
import pandas as pd
import os

np.random.seed(42)

fraction_selected_tracers = 0.005 # Randomly select this many tracers
data_working_path = '/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem/data' # '/scratch/09338/carlos9/turbulence_original'
final_save_path = os.path.join(data_working_path, 'turbulence_tracers_csv')

model_counter = 0

benchmarks = {
    # "M150_1": {
    #     "path": os.path.join(data_working_path, "M150_seed1_trace_cells.npy"),
    #     "timestep": 16.559, # in kyrs
    #     "clip": 200, # only take first 200 timesteps
    #     "discretization": 4, # take every timestep
    # },
    # # "M150_42": {
    # #     "path": os.path.join(data_working_path, "M150_seed42_trace_cells.npy"),
    # #     "timestep": 16.559,
    # #     "clip": 200,
    # #     "discretization": 1,
    # # },
    "M600_1": {
        "path": os.path.join(data_working_path, "M600_seed1_trace_cells.npy"),
        "timestep": 8.299,
        "clip": 400,
        "discretization": 1,
    },
    # # "M600_42": {
    # #     "path": os.path.join(data_working_path, "M600_seed42_trace_cells.npy"),
    # #     "timestep": 8.299,
    # #     "clip": 400,
    # #     "discretization": 2,
    # # },
    # "M2400_1": {
    #     "path": os.path.join(data_working_path, "M2400_seed1_trace_cells.npy"),
    #     "timestep": 4.145,
    #     "clip": 800,
    #     "discretization": 4,
    # },
    # "M2400_42": {
    #     "path": os.path.join(data_working_path, "M2400_seed42_trace_cells.npy"),
    #     "timestep": 4.145,
    #     "clip": 800,
    #     "discretization": 4,
    # },
}

def density_to_number_density(density):
    hydrogen_mass = 1.66053906660e-24
    mean_molecular_mass = 1.4168138025 # Calculated separately. Should be ~ constant.
    return density / (mean_molecular_mass * hydrogen_mass)

def kalman_filter(series, process_var=1e-2, meas_var=1e-1, init_error=1e-3):
    series = np.asarray(series)
    n = len(series)
    xhat = np.zeros(n)
    P = np.zeros(n)
    K = np.zeros(n)
    xhat[0] = series[0]
    P[0] = init_error
    for k in range(1, n):
        xhat_minus = xhat[k-1]
        P_minus = P[k-1] + process_var
        K[k] = P_minus / (P_minus + meas_var)
        xhat[k] = xhat_minus + K[k] * (series[k] - xhat_minus)
        P[k] = (1 - K[k]) * P_minus
    return xhat    

for benchmark_name, benchmark_info in benchmarks.items():
    benchmark_path = benchmark_info["path"]
    timestep = benchmark_info["timestep"]
    clip = benchmark_info["clip"]
    discretization = benchmark_info["discretization"]
    
    print(f"Processing benchmark {benchmark_name} with timestep {timestep} kyrs")
    
    data = np.load(benchmark_path)
    print(f"Loaded data shape: {data.shape}")
    
    # random_indices = np.random.choice(data.shape[1], size=int(fraction_selected_tracers * data.shape[1]), replace=False)
    random_indices = np.random.choice(data.shape[1], size=100, replace=False)
    # random_indices = [10]
    selected_data = data[:clip:discretization, random_indices, :]

    os.makedirs(final_save_path, exist_ok=True)

    for tracer_index in range(selected_data.shape[1]):
        current_tracer = selected_data[:, tracer_index, :]
        current_index = random_indices[tracer_index]
        
        df = pd.DataFrame(current_tracer, columns=[
            'density', 'gasTemp', 'av', 'PI_Rad', 'radField', 'NUV_Rad', 'NIR_Rad', 'IR_Rad'
        ])
        
        df['time'] = np.arange(len(current_tracer)) * timestep * discretization  # in kyrs
        df['tracer'] = current_index  # original tracer index
        df['benchmark'] = benchmark_name
        df['model'] = model_counter
        model_counter += 1

        df = df[['model', 'benchmark', 'tracer', 'time', 'gasTemp', 'density', 'av', 'radField']]
        
        df['density'] = density_to_number_density(df['density'].values)
        
        df['radField'] = kalman_filter(np.log10(df['radField'].values))
        df['radField'] = np.power(10, df['radField'])
                
        csv_filename = os.path.join(final_save_path, f'{benchmark_name}_{discretization}_Tracer_{current_index}.csv')
        df.to_csv(csv_filename, index=False)
        
        print(f"Saved tracer {current_index} to {csv_filename}")

print("All benchmarks processed and saved.")