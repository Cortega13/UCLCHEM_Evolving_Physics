import numpy as np

file_paths = [
    "/scratch1/09338/carlos9/turbulence_chemistry/M150_seed1_trace_cells.npy",
    "/scratch1/09338/carlos9/turbulence_chemistry/M150_seed42_trace_cells.npy",
    "/scratch1/09338/carlos9/turbulence_chemistry/M2400_seed1_trace_cells.npy",
    "/scratch1/09338/carlos9/turbulence_chemistry/M2400_seed42_trace_cells.npy",
    "/scratch1/09338/carlos9/turbulence_chemistry/M600_seed1_trace_cells.npy",
    "/scratch1/09338/carlos9/turbulence_chemistry/M600_seed42_trace_cells.npy",
]


for path in file_paths:
    data = np.load(path)
    print(f"{path}: {data.shape}")
    del data