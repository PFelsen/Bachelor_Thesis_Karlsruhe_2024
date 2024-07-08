# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:47:28 2024

@author: phili
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

def estimate_coefficients(stoch_levels, values):
    """
    Estimate the coefficients for the inequality |value| <= c * h_l^a using linear regression.
    
    Parameters:
    - stoch_levels: The levels of stochastic (l).
    - values: The absolute differences of the values between consecutive levels.
    
    Returns:
    - a: The estimated coefficient a.
    """
    # Filter out zero values
    non_zero_indices = values != 0
    h_l = 2 ** -stoch_levels[non_zero_indices]
    filtered_values = values[non_zero_indices]

    # Perform linear regression on log-log scale using numpy.polyfit
    X = np.log(h_l)
    y = np.log(filtered_values)

    a, _ = np.polyfit(X, y, 1)

    return a


def plot_mean_and_variance(file, parameter, parameter_values, plotting_bool=True):
    """
    Function to plot mean and variance values and their absolute differences for given parameter values.
    
    Parameters:
    - file: Path to the input pickle file containing the DataFrame.
    - parameter: The parameter name to filter the DataFrame on.
    - parameter_values: List of parameter values to plot.
    """
    
    label_fontsize = 20
    tick_fontsize = 16
    
    # Load the new DataFrame with alpha values
    df = pd.read_pickle(file)

    # Define the fixed list of colors excluding red and green
    colors = ['blue', 'yellow', 'purple', 'orange', 'red', 'green', 'cyan', 'magenta', 'pink', 'brown'][:len(parameter_values)] # Ensure color list matches number of values

    # Initialize the plot
    plt.figure(figsize=(24, 8))

    # Subplot for Mean
    plt.subplot(1, 2, 1)
    for idx, value in enumerate(parameter_values):
        # Filter the DataFrame for the current value 
        df_filtered = df[df[parameter] == value].sort_values(by='stochLevel')
        df_filtered = df_filtered[df_filtered['stochLevel'] != 6]

        # Calculate the absolute differences of Mean values between consecutive stochLevels
        mean_diff = abs(df_filtered['Mean'].diff().iloc[1:])

        # Plot the absolute differences of mean
        plt.plot(df_filtered['stochLevel'].iloc[1:], mean_diff, marker='o', linestyle='--', label=f'{parameter} {value}', color=colors[idx])

        # Plot the mean values
        plt.plot(df_filtered['stochLevel'], df_filtered['Mean'], marker='o', linestyle='-', label=f'{parameter} {value}', color=colors[idx])
        
        stoch_levels = df_filtered['stochLevel'].iloc[1:]
        a_mean = estimate_coefficients(stoch_levels, mean_diff)
        print(f"For {parameter}= {value}: Estimated a (Mean) = {a_mean:.2f}")

    # Set the labels and title for Mean subplot
    plt.xlabel('stochLevel l', fontsize = label_fontsize)
    plt.ylabel('Mean', fontsize = label_fontsize)
    plt.yscale('log')
    plt.title(r"-$\log(|\mathbb{E}[Q_{h_l}]|)$ and --$\log(|\mathbb{E}[Q_{h_l}-Q_{h_{l-1}}]|)$", fontsize=label_fontsize)
    plt.legend(fontsize = tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.yticks(fontsize=tick_fontsize)
    plt.grid(True)

    # Subplot for Variance
    plt.subplot(1, 2, 2)
    for idx, value in enumerate(parameter_values):
        # Filter the DataFrame for the current value
        df_filtered = df[df[parameter] == value].sort_values(by='stochLevel')
        df_filtered = df_filtered[df_filtered['stochLevel'] != 6]

        # Calculate the absolute differences of variance values between consecutive stochLevels
        var_diff = abs(df_filtered['sVar'].diff().iloc[1:])

        # Plot the absolute differences of variance
        plt.plot(df_filtered['stochLevel'].iloc[1:], var_diff, marker='o', linestyle='--', label=f'{parameter} {value}', color=colors[idx])

        # Plot the variance values
        plt.plot(df_filtered['stochLevel'], df_filtered['sVar'], marker='o', linestyle='-', label=f'{parameter} {value}', color=colors[idx])
        
        stoch_levels = df_filtered['stochLevel'].iloc[1:]
        b_var = estimate_coefficients(stoch_levels, var_diff)
        print(f"For {parameter}= {value}: Estimated b (Variance) = {b_var:.2f}")
        

    # Set the labels and title for Variance subplot
    plt.xlabel('stochLevel l', fontsize = label_fontsize)
    plt.ylabel('Variance', fontsize = label_fontsize)
    plt.yscale('log')
    plt.title(r"-$\log(|\operatorname{Var}[Q_{h_l}]|)$ and --$\log(|\operatorname{Var}[Q_{h_l}-Q_{h_{l-1}}]|)$", fontsize=label_fontsize)
    plt.legend(fontsize = tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(tick_fontsize)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.grid(True)
    plt.tight_layout()
    
    if plotting_bool==True:
        name = parameter + "_convergence_plot.png"
        plt.savefig(name, dpi=500)
    plt.show()
    
  


sigma_values = [1.0, 1.2, 1.4, 1.6]
plot_mean_and_variance('df_sigma.pickle', "sigma", sigma_values)

alpha_values = [1.0, 1.2, 1.4, 1.6]
plot_mean_and_variance('df_alpha.pickle', "alpha", alpha_values)


length_values = [0.4, 0.6, 0.8, 1.0]
plot_mean_and_variance('df_length.pickle', "length", length_values)
