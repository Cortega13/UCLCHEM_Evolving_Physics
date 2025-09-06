import numpy as np
import pandas as pd
import uclchem
import os
import re
import argparse

working_path = '/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem' #"/work/09338/carlos9/vista/UCLCHEM/benchmarks/turbulence"
save_path = '/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem/data/chemistry_tracers_csv' #"/scratch/09338/carlos9/turbulence_chemistry/chemistry_tracers_csv"
initial_abundances_path = os.path.join(working_path, 'initial_abundances.csv')
df_init = pd.read_csv(initial_abundances_path)


def rename_columns(columns):
    """Renames column names containing chemical species using substring replacement."""
    name_mapping = {
        'H2COH+': 'H3CO+',
        'H2COH': 'H3CO',
        'H2CSH+': 'H3CS+',
        'SISH+': 'HSIS+',
        'HOSO+': 'HSO2+',
        'OCSH+': 'HOCS+',
        'HCOO': 'HCO2',
        'HCOOH': 'H2CO2',
        'CH2CO': 'C2H2O',
        'CH2OH': 'CH3O',
        'CH3CCH': 'C3H4',
        'CH3CHO': 'C2H4O',
        'CH3CN': 'C2H3N',
        'CH3CNH': 'C2H4N',
        'CH3OH': 'CH4O',
        'CH3OH2+': 'CH5O+',
        'CH3CNH+': 'C2H4N+',
        'NH2CHO': 'CH3NO',
        'HCO2H': 'H2CO2',
        'HCNH': 'H2CN',
        'NCCN': 'N2C2',
    }

    sorted_mapping = sorted(name_mapping.items(), key=lambda x: -len(x[0]))

    columns = [col.strip() for col in columns]
    new_columns = []
    for col in columns:
        new_col = col
        for old, new in sorted_mapping:
            if old in new_col:
                new_col = new_col.replace(old, new)
        new_columns.append(new_col)

    return new_columns


def log_linear_interpolation(timestep, prev, curr):
    intercept = prev
    slope = (curr - intercept) / (timestep * 3.16e7)
    return (intercept, slope)

def process_row(row, starting_chemistry_init, timestep, templinear, avlinear, radfieldlinear, densitylinear):
    starting_chemistry = starting_chemistry_init.copy()
    
    param_dict = {
        "initialDens": densitylinear[0],  # starting density
        "initialTemp": templinear[0],  # temperature of gas
        "initialBaseAv": avlinear[0],  # visual extinction at cloud edge.
        'initialRadfield': radfieldlinear[0],  # initial radiation field strength
            
        'tempRate': templinear[1],
        'densRate': densitylinear[1],
        'baseAvRate': avlinear[1],
        'radfieldRate': radfieldlinear[1],

        'evolvingphysicalparams': True,
        'useavdirectly': True,
        "currentTime": 0,
        "endatfinaldensity": False,
        "freefall": False,
        "zeta": 1.2213740458,
        "finalTime": timestep,
        "rout": 1.5,
        "abstol_min": 1e-22,
    }
    result = uclchem.model.cloud(
        param_dict=param_dict,
        return_dataframe=True,
        starting_chemistry=starting_chemistry
    )
    df_temp = pd.concat([result[0], result[1]], axis=1).iloc[2:].reset_index(drop=True)
    df_temp["Time"] = row['time'] + df_temp['Time'] - timestep
    df_temp['tracer'] = row['tracer']
    return df_temp, result[3], result[4]

def process_csv_file(csv_path, out_path):
    
    tracer_df = pd.read_csv(csv_path)
    tracer_df["radField"] = tracer_df["radField"] * 1.7
    tracer_df['time'] = tracer_df['time'] * 1000

    prev_av = tracer_df['av'].iloc[0]
    prev_radfield = tracer_df['radField'].iloc[0]
    prev_temp = tracer_df['gasTemp'].iloc[0]
    prev_density = tracer_df['density'].iloc[0]

    results = [df_init]
    starting_chemistry = df_init.to_numpy().flatten()
    timestep = tracer_df['time'].iloc[1]

    for i, row in tracer_df.iterrows():
        if i == 0:
            continue
        print(f"Processing at time {row['time']/1000} kyr")

        templinear = log_linear_interpolation(timestep, prev_temp, row['gasTemp'])
        avlinear = log_linear_interpolation(timestep, prev_av, row['av'])
        radfieldlinear = log_linear_interpolation(timestep, prev_radfield, row['radField'])
        densitylinear = log_linear_interpolation(timestep, prev_density, row['density'])

        df_result, starting_chemistry, success_flag = process_row(row, starting_chemistry, timestep, templinear, avlinear, radfieldlinear, densitylinear)
        if success_flag < 0:
            print("Too many fails error. Exiting tracer.")
            return
        if i == 1:
            print(df_result['Time'])
        results.append(df_result)
        prev_av = row['av']
        prev_radfield = row['radField']
        prev_temp = row['gasTemp']
        prev_density = row['density']

    result_df = pd.concat(results, ignore_index=True).drop(columns=["dstep", "dustTemp", "zeta"])
    result_df.rename(columns={
        'Time': 'time',
        'Av': 'av',
        'radfield': 'radField',
        'Density': 'density',
    }, inplace=True)
    result_df.loc[result_df.index[0], 'time'] = 0.0
    result_df.loc[result_df.index[0], 'av'] = tracer_df['av'].iloc[0]
    result_df.loc[result_df.index[0], 'radField'] = tracer_df['radField'].iloc[0]
    result_df.loc[result_df.index[0], 'density'] = tracer_df['density'].iloc[0]
    result_df.loc[result_df.index[0], 'gasTemp'] = tracer_df['gasTemp'].iloc[0]
    
    os.makedirs(save_path, exist_ok=True)
    result_df.to_csv(out_path, index=False)

if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description="Generate chemical abundances from physical parameter evolution in a csv file.")
    # parser.add_argument("csv_path", type=str, help="The path of the csv file containing the physical parameter evolution.")
    # args = parser.parse_args()
    # csv_path = args.csv_path
    csv_path = "/mnt/c/users/carlo/projects/UCLCHEM_Evolving_Physics/benchmarks/turbulence_uclchem/data/turbulence_tracers_csv/M600_1_1_Tracer_10.csv"
    
    info = csv_path.split("/")[-1].split("_Tracer_")
    benchmark = info[0]
    tracer_num = info[1].split(".")[0]
    out_path = os.path.join(save_path, f"{benchmark}_Tracer_{tracer_num}.csv")
    
    if False: #os.path.exists(out_path):
        print(f"Output file {out_path} already exists. Skipping processing.")
    else:
        process_csv_file(csv_path, out_path)
