# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:06:13 2024

@author: Philipp Felsen
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths for the pickle files
file_paths = [
    'DataFrames/df_sc_1.pickle',
    'DataFrames/df_sc_2.pickle',
    'DataFrames/df_sc_3.pickle',
    'DataFrames/df_sc_4.pickle',
    'DataFrames/df_sc_5.pickle'
]

# Load the dataframes from the pickle files
data_frames = [pd.read_pickle(file) for file in file_paths]

# Function to clean the dataframes
def clear_data_frames(df):
    # Check if the 'trunc' column exists and filter rows where 'trunc' is 4
    if 'trunc' in df.columns:
        df = df[df["trunc"] == 4]
    
    # Define columns to be removed
    columns_to_remove = ['MSE', 'RMSE']
    
    # Conditionally add 'Cost' and 'Estimator' to the columns to remove if they exist
    if 'Cost' in df.columns:
        columns_to_remove.append('Cost')
    if 'Estimator' in df.columns:
        columns_to_remove.append('Estimator')
    
    # Drop the specified columns, ignoring errors if columns don't exist
    df = df.drop(columns=columns_to_remove, errors='ignore')
    
    return df

# Apply the cleaning function to each dataframe
cleared_data_frames = [clear_data_frames(df) for df in data_frames]

# Function to plot mean convergence for a given stochastic level
def plot_mean_convergence(stoch_level):
    label_fontsize = 20
    tick_fontsize = 16

    # Filter dataframes for the given stochastic level
    filtered_data_frames = [df[df['stochLevel'] == stoch_level].copy() for df in cleared_data_frames]

    # Add a 'SpaceLevel' column to each dataframe and merge them
    for i, df in enumerate(filtered_data_frames):
        df.loc[:, 'SpaceLevel'] = i + 1

    combined_df_with_spacelevel = pd.concat(filtered_data_frames, ignore_index=True)

    # Calculate the differences between consecutive Mean values
    diffs = []
    for i in range(1, len(combined_df_with_spacelevel)):    
        diffs.append(abs(combined_df_with_spacelevel["Mean"][i] - combined_df_with_spacelevel["Mean"][i - 1]))

    # Plot the differences and Mean values on the same Y-axis with SpaceLevel on the X-axis
    plt.figure(figsize=(12, 8))

    # Plot the differences with a logarithmic scale
    plt.plot(combined_df_with_spacelevel["SpaceLevel"][1:], diffs, marker='o', linestyle='--', color='blue', label=r'$\log\left(\left|\mathbb{E}[Q_{h_l} - Q_{h_{l-1}}]\right|\right)$')
    plt.yscale('log')

    # Plot the Mean values
    plt.plot(combined_df_with_spacelevel["SpaceLevel"], combined_df_with_spacelevel['Mean'], marker='x', linestyle='-', color='red', label=r'$\log\left(\left|\mathbb{E}[Q_{h_l}]\right|\right)$.')

    # Set axis labels and title
    plt.xlabel('SpaceLevel', fontsize=label_fontsize)
    plt.ylabel(r'$\log\left(\left|\mathbb{E}[\cdot]\right|\right)$', fontsize=label_fontsize)
    plt.title(r'$\log\left(\left|\mathbb{E}[Q_{h_l} - Q_{h_{l-1}}]\right|\right)$ and $\log\left(\left|\mathbb{E}[Q_{h_l}]\right|\right)$ for StochLevel=' + str(stoch_level), fontsize=label_fontsize)

    plt.xticks(ticks=combined_df_with_spacelevel["SpaceLevel"], fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(tick_fontsize)

    # Add legend and grid
    plt.legend(fontsize=label_fontsize)
    plt.grid(True)

    # Save the plot as a PNG file
    plt.savefig(f'mean_convergence_level_{stoch_level}.png', dpi=500)
    plt.show()
    
    # Estimate the parameter a
    hl = np.array([2**-l for l in combined_df_with_spacelevel["SpaceLevel"][2:]])
    log_diffs = np.log(diffs[1:])
    log_hl = np.log(hl)
    
    # Use numpy.polyfit to estimate the parameter a
    a_estimated, _ = np.polyfit(log_hl, log_diffs, 1)
    
    print(f"StochLevel {stoch_level}, Estimated parameter a: {a_estimated}")

# Example function calls to plot mean convergence for different stochastic levels
plot_mean_convergence(1)
plot_mean_convergence(2)
plot_mean_convergence(3)
plot_mean_convergence(4)
plot_mean_convergence(5)
plot_mean_convergence(6)
plot_mean_convergence(7)


