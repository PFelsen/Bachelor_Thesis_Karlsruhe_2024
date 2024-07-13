# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:29:55 2024

@author: Philipp Felsen
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Dateipfade der neuen pickle-Dateien
new_file_paths = [
    'df_alpha.pickle',
    'df_length.pickle',
    'df_sigma.pickle'
]

# Laden der neuen Datenrahmen
df_alpha = pd.read_pickle(new_file_paths[0])
df_length = pd.read_pickle(new_file_paths[1])
df_sigma = pd.read_pickle(new_file_paths[2])

# Funktion zum Bereinigen der Datenrahmen
def clear_data_frames(df):
    if 'trunc' in df.columns:
        df = df[df["trunc"] == 4]
    
    columns_to_remove = ['MSE', 'RMSE', 'sVar', 'Samples']
    
    if 'Cost' in df.columns:
        columns_to_remove.append('Cost')
    if 'Estimator' in df.columns:
        columns_to_remove.append('Estimator')
    
    df = df.drop(columns=columns_to_remove, errors='ignore')
    
    return df

# Bereinigen der neuen Datenrahmen
df_alpha = clear_data_frames(df_alpha)
df_length = clear_data_frames(df_length)
df_sigma = clear_data_frames(df_sigma)

print(df_alpha)

# Funktion zum Plotten und Schätzen des Parameters a
def plot_and_estimate(df, stoch_levels, param_name):
    unique_params = df[param_name].unique()
    
    for param in unique_params:
        for stoch_level in stoch_levels:
            filtered_data = df[(df[param_name] == param) & (df['stochLevel'] == stoch_level)]
            filtered_data = filtered_data.reset_index(drop=True)
            
            
            # Hinzufügen der SpaceLevel-Spalte
            filtered_data['SpaceLevel'] = filtered_data.index + 1

            print(filtered_data)
            """
            # Berechnung der Differenzen
            diffs = []
            for i in range(1, len(filtered_data)):    
                diffs.append(abs(filtered_data["Mean"][i] - filtered_data["Mean"][i - 1]))

            # Plotten der Differenzen und der Mean-Werte auf derselben Y-Achse mit SpaceLevel auf der x-Achse
            plt.figure(figsize=(12, 8))

            # Plot der Differenzen mit logarithmischer Skala
            plt.plot(filtered_data["SpaceLevel"][1:], diffs, marker='o', linestyle='--', color='blue', label=r'log(E[|Q_{h_l} - Q_{h_{l-1}}|])')
            plt.yscale('log')

            # Plot der Mean-Werte
            plt.plot(filtered_data["SpaceLevel"], filtered_data['Mean'], marker='x', linestyle='-', color='red', label=r'log(E[|Q_{h_l}|])')

            # Achsenbeschriftungen und Titel
            plt.xlabel('SpaceLevel')
            plt.ylabel('log(Wert)')
            
            # LaTeX-formatierter Titel mit param_name und stoch_level
            plt.title(r'$\log\left(\left|\mathbb{E}[Q_{h_l} - Q_{h_{l-1}}]\right|\right)$ and $\log\left(\left|\mathbb{E}[Q_{h_l}]\right|\right)$ for ' 
                      + param_name + f'={param}, stochLevel={stoch_level}', fontsize=14)
            
            plt.xticks(ticks=filtered_data["SpaceLevel"])

            # Legende und Gitter
            plt.legend()
            plt.grid(True, which="both", ls="--")

            # Speichern des Plots
            plt.savefig(f'{param_name}_convergence_level_{stoch_level}_param_{param}.png')
            plt.show()

            # Schätzung des Parameters a
            hl = np.array([2**-l for l in filtered_data["SpaceLevel"][2:]])
            log_diffs = np.log(diffs[2:])
            log_hl = np.log(hl)
            
            # Sicherstellen, dass die Arrays die gleiche Länge haben
            min_length = min(len(log_diffs), len(log_hl))
            log_diffs = log_diffs[:min_length]
            log_hl = log_hl[:min_length]
            
            # Nutzung von numpy.polyfit zur Schätzung des Parameters a
            a_estimated, _ = np.polyfit(log_hl, log_diffs, 1)
            
            print(f"Estimated parameter a for {param_name}={param}, stochLevel={stoch_level}: {a_estimated}")
            """
            
            
# Beispielaufruf der Funktion für alle neuen Datenrahmen und stochLevels 3 und 6
stoch_levels = [3, 6]
param_names = ['alpha', 'length', 'sigma']
dataframes = [df_alpha, df_length, df_sigma]

for df, param_name in zip(dataframes, param_names):
    plot_and_estimate(df, stoch_levels, param_name)


