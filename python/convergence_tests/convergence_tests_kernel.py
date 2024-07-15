# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:29:55 2024

@author: Philipp Felsen
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths for the pickle files
sigma_paths = [
    'DataFrames/df_sigma_1.pickle',
    'DataFrames/df_sigma_2.pickle',
    'DataFrames/df_sigma_3.pickle',
    'DataFrames/df_sigma_4.pickle',
    'DataFrames/df_sigma_5.pickle'
]

alpha_paths = [
    'DataFrames/df_alpha_1.pickle',
    'DataFrames/df_alpha_2.pickle',
    'DataFrames/df_alpha_3.pickle',
    'DataFrames/df_alpha_4.pickle',
    'DataFrames/df_alpha_5.pickle'
]

length_paths = [
    'DataFrames/df_length_1.pickle',
    'DataFrames/df_length_2.pickle',
    'DataFrames/df_length_3.pickle',
    'DataFrames/df_length_4.pickle',
    'DataFrames/df_length_5.pickle'
]

# Load the new dataframes
df_sigma = [pd.read_pickle(path) for path in sigma_paths]
df_alpha = [pd.read_pickle(path) for path in alpha_paths]
df_length = [pd.read_pickle(path) for path in length_paths]

def clear_data_frames(df):
    """
    Cleans the provided dataframe by removing specific columns and filtering rows based on the 'trunc' column.

    Parameters:
    df (pd.DataFrame): The dataframe to clean.

    Returns:
    pd.DataFrame: The cleaned dataframe.
    """
    # Remove rows where the 'trunc' column is not equal to 4
    if 'trunc' in df.columns:
        df = df[df["trunc"] == 4]
    
    # Specify columns to remove
    columns_to_remove = ['MSE', 'RMSE']
    
    # Additional columns to remove if they exist
    if 'Cost' in df.columns:
        columns_to_remove.append('Cost')
    if 'Estimator' in df.columns:
        columns_to_remove.append('Estimator')
    
    # Drop the specified columns
    df = df.drop(columns=columns_to_remove, errors='ignore')
    
    return df

# Organize the dataframes into dictionaries for easier processing
dfs_sigma = {i+1: clear_data_frames(df) for i, df in enumerate(df_sigma)}
dfs_alpha = {i+1: clear_data_frames(df) for i, df in enumerate(df_alpha)}
dfs_length = {i+1: clear_data_frames(df) for i, df in enumerate(df_length)}

def clear_and_concatenate(dfs):
    """
    Cleans and concatenates a dictionary of dataframes, adding a 'SpaceLevel' column to each.

    Parameters:
    dfs (dict): Dictionary of dataframes to be concatenated.

    Returns:
    pd.DataFrame: The concatenated dataframe.
    """
    cleared_dfs = []
    for level, df in dfs.items():
        # Add 'SpaceLevel' column
        df['SpaceLevel'] = level
        cleared_dfs.append(df)
    # Concatenate all cleaned dataframes
    return pd.concat(cleared_dfs, ignore_index=True)

# Clear and concatenate the dataframes for each parameter
combined_df_sigma = clear_and_concatenate(dfs_sigma)
combined_df_alpha = clear_and_concatenate(dfs_alpha)
combined_df_length = clear_and_concatenate(dfs_length)

def plot_and_estimate(df, stoch_levels, param_name, other_params):
    """
    Plots the convergence of the parameter estimates and estimates the parameter 'a' for different stochastic levels.

    Parameters:
    df (pd.DataFrame): The dataframe containing the data to be plotted and analyzed.
    stoch_levels (list): List of stochastic levels to be analyzed.
    param_name (str): Name of the parameter being analyzed.
    other_params (list): List of other parameters for labeling purposes.
    """
    label_fontsize = 20
    tick_fontsize = 16
    
    colors = ['red', 'blue', 'orange']
    
    latex_names = {'alpha': r'$\alpha$', 'length': r'$\lambda$', 'sigma': r'$\sigma$'}
    
    for stoch_level in stoch_levels:
        plt.figure(figsize=(12, 8))
        for idx, param in enumerate(df[param_name].unique()):
            # Filter data for the current parameter and stochastic level
            filtered_data = df[(df[param_name] == param) & (df['stochLevel'] == stoch_level)]
            filtered_data = filtered_data.reset_index(drop=True)
            
            # Add 'SpaceLevel' column
            filtered_data['SpaceLevel'] = filtered_data.index + 1

            # Calculate the differences
            diffs = []
            for i in range(1, len(filtered_data)):    
                diffs.append(abs(filtered_data["Mean"][i] - filtered_data["Mean"][i - 1]))

            # Create labels for the plot
            label_diffs = rf'{latex_names[param_name]}={param}'
            label_means = rf'{latex_names[param_name]}={param}'
            for other_param in other_params:
                if other_param in df.columns:
                    label_diffs += rf', {latex_names[other_param]}={filtered_data[other_param].iloc[0]}'
                    label_means += rf', {latex_names[other_param]}={filtered_data[other_param].iloc[0]}'
            
            color = colors[idx % len(colors)]
            # Plot the differences and mean values on the same y-axis with SpaceLevel on the x-axis
            plt.plot(filtered_data["SpaceLevel"][1:], diffs, marker='o', linestyle='--', label=label_diffs, color=color)
            plt.plot(filtered_data["SpaceLevel"], filtered_data['Mean'], marker='x', linestyle='-', color=color)

            plt.yscale('log')
            
            # Estimate the parameter 'a'
            hl = np.array([2**-l for l in filtered_data["SpaceLevel"][1:]])
            log_diffs = np.log(diffs)
            log_hl = np.log(hl)
            
            # Use numpy.polyfit to estimate the parameter 'a'
            a_estimated, _ = np.polyfit(log_hl, log_diffs, 1)
            
            print(f"Estimated parameter a for {param_name}={param}, stochLevel={stoch_level}: {a_estimated}")

        # Set axis labels and title
        plt.xlabel('SpaceLevel', fontsize=label_fontsize)
        plt.ylabel(r'$\log\left(\left|\mathbb{E}[\cdot]\right|\right)$', fontsize=label_fontsize)
        
        # LaTeX-formatted title with param_name and stoch_level
        plt.title(r'--$\log\left(\left|\mathbb{E}[Q_{h_l} - Q_{h_{l-1}}]\right|\right)$ and -$\log\left(\left|\mathbb{E}[Q_{h_l}]\right|\right)$ for StochLevel=' + str(stoch_level), fontsize=label_fontsize)
        
        plt.xticks(ticks=filtered_data["SpaceLevel"], fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        
        ax = plt.gca()
        ax.yaxis.get_offset_text().set_fontsize(tick_fontsize)

        # Add legend and grid
        plt.legend(fontsize=label_fontsize)
        plt.grid(True)

        # Save the plot
        plt.savefig(f'{param_name}_convergence_level_{stoch_level}.png')
        plt.show()

# Example call of the function for all new dataframes and stochLevels 2 and 5
stoch_levels = [2, 4]
param_names = ['sigma', 'alpha', 'length']
combined_dfs = [combined_df_sigma, combined_df_alpha, combined_df_length]

for df, param_name in zip(combined_dfs, param_names):
    other_params = [name for name in param_names if name != param_name]
    plot_and_estimate(df, stoch_levels, param_name, other_params)
