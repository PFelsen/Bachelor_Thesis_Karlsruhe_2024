# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 19:30:38 2024

@author: phili
"""


import numpy as np
#from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import time

label_fontsize = 20
tick_fontsize = 16


# Definition of phi function
def phi(x):
    """
    Computes the value of the piecewise linear function phi(x).
    
    Parameters:
        x (float): The input value.
        
    Returns:
        float: The value of phi(x).
    """
    return np.where(np.logical_and(-1 <= x, x <= 1), 1 - np.abs(x), 0)


# Definition of phi_j_k function
def phi_j_k(x, j, k, a, b):
    """
    Computes the value of the scaled and translated piecewise linear function phi(x).
    
    Parameters:
        x (float): The input value.
        j (int): The scaling factor.
        k (int): The translation factor.
        a (float): Left bound of the intervall.
        b (float): Right bound of the intervall.
        
    Returns:
        float: The value of phi(2^j * x - k).
    """
    return phi(((x-a)/(b-a))*2 ** j  - k)


# Integral Function
def int_phi(l,a,b):
    """
    Computes the integral values for the function phi_{l,i}(x).
    
    Parameters:
        l (float): Level.
        a (float): Left bound of the interval.
        b (float): Right bound of the interval.
        
    Returns:
        float: The value of the integral.
    """
   
    if (l==0):
        return (b-a)/2
        
    if (l>=1):
        return (b-a)*2**(-l)


# Computes hierarchical surplus
def delta(f,l_1, l_2, i_1, i_2, a, b, c, d):
    """
    Computes the hierarchical surplus for a given function at specific levels and indices.
    
    Parameters:
        f (function): The function for which to compute the surplus.
        l_1 (int): The level in the first dimension.
        l_2 (int): The level in the second dimension.
        i_1 (int): The index in the first dimension.
        i_2 (int): The index in the second dimension.
        a (float): Left bound of the interval in x-direction.
        b (float): Right bound of the interval in x-direction.
        c (float): Left bound of the interval in y-direction.
        d (float): Right bound of the interval in y-direction.
        
    Returns:
        float: The hierarchical surplus value.
    """
    erg = 0.0
        
    
    h_l_1 = (b - a) * 2**(-l_1)
    h_l_2 = (d - c) * 2**(-l_2)
    
    x = a + i_1 * h_l_1
    y = c + i_2 * h_l_2
    
    
    if l_1==0 and l_2==0:
        return f(x,y)
    
    if l_1==0 and l_2 >=1:
        return f(x, y) -0.5 * (f(x, y - h_l_2) + f(x, y + h_l_2))
    
    if l_2==0 and l_1 >=1:
        return f(x, y) -0.5 * (f(x - h_l_1, y) + f(x + h_l_1, y))
    
        
    erg += f(x, y) -0.5 * (f(x, y - h_l_2) + f(x, y + h_l_2))
    erg += -0.5*(f(x - h_l_1 , y) -0.5 * (f(x - h_l_1, y - h_l_2) + f(x - h_l_1, y + h_l_2)))
    erg += -0.5*(f(x + h_l_1 , y) -0.5 * (f(x + h_l_1, y - h_l_2) + f(x + h_l_1, y + h_l_2)))

    return erg


# Returns an np.array with index sets
def J(k):
    """
    Generates the index set for a given level k.
    
    Parameters:
        k (int): The level for which to generate the index set.
        
    Returns:
        list: The index set for level k.
    """
    erg = []
    if k==0:
        return np.array([0,1])
    
    for i in range(1,2**k):
        if i%2 !=0:
            erg.append(i)
    return erg
   

def sparse_grid_integration(f, L, X, Y):
    """
    Computes the integral for a given function in two dimensions over [a,b]x[c,d].
    
    Parameters:
        f (function): The function to be integrated. Must accept two arguments (x, y).
        L (int): Maximum level of refinement in each dimension.
        X (array): An array containing the interval [a,b] in x-direction.
        Y (array): An array containing the interval [c,d] in y-direction.
        
    Returns:
        float: The value of the integral.
        
    """
    erg = 0
    a = X[0]
    b = X[1]
    c = Y[0]
    d = Y[1]
    
    l = L + 2 -1 # L + d - 1 sparse grid construction
      
    
    for l_1 in range(0,l + 1):
        for l_2 in range(0, l - l_1 + 1):

            for i_1 in J(l_1):
                for i_2 in J(l_2):

                    erg += delta(f, l_1, l_2, i_1, i_2, a = a, b = b, c = c, d = d) * int_phi(l_1, a, b) * int_phi(l_2, c, d) 
    return erg


def full_grid_integration(f, L, X, Y):
    """
    Computes the integral for a given function in two dimensions over [a,b]x[c,d].
    
    Parameters:
        f (function): The function to be integrated. Must accept two arguments (x, y).
        L (int): Maximum level of refinement in each dimension.
        X (array): An array containing the interval [a,b] in x-direction.
        Y (array): An array containing the interval [c,d] in y-direction.
        
    Returns:
        float: The value of the integral.
        
    """
    erg = 0
    a = X[0]
    b = X[1]
    c = Y[0]
    d = Y[1]
    
    
    l = L + 2 -1 # L + d - 1 sparse grid construction
      
    
    for l_1 in range(0, l + 1):
        for l_2 in range(0, l + 1):

            for i_1 in J(l_1):
                for i_2 in J(l_2):

                    erg += delta(f, l_1, l_2, i_1, i_2, a = a, b = b, c = c, d = d) * int_phi(l_1, a, b) * int_phi(l_2, c, d) 
    return erg










def random_polynomial(degree):
    """
    Generates random polynomial coefficients.
    
    Parameters:
        degree (int): The degree of the polynomial.
        
    Returns:
        array: An array of random coefficients.
    """
    return np.random.uniform(-10, 10, degree + 1)

def polynomial_function(coefficients):
    """
    Creates a polynomial function based on given coefficients.
    
    Parameters:
        coefficients (array): Array of polynomial coefficients.
        
    Returns:
        function: A function representing the polynomial.
    """
    def f(x, y):
        return sum(coefficients[i] * (x**(len(coefficients) - i - 1)) * (y**j) 
                   for i in range(len(coefficients)) for j in range(len(coefficients) - i))
    return f

def exact_integral(coefficients, a, b, c, d):
    """
    Computes the exact integral of a polynomial over a given interval.
    
    Parameters:
        coefficients (array): Array of polynomial coefficients.
        a (float): Lower bound of interval in x-direction.
        b (float): Upper bound of interval in x-direction.
        c (float): Lower bound of interval in y-direction.
        d (float): Upper bound of interval in y-direction.
        
    Returns:
        float: The exact value of the integral.
    """
    integral = 0
    for i in range(len(coefficients)):
        for j in range(len(coefficients) - i):
            integral += coefficients[i] * (b**(len(coefficients) - i) - a**(len(coefficients) - i)) * (d**(j+1) - c**(j+1)) / ((len(coefficients) - i) * (j+1))
    return integral

def polynom_test(n, L, saving = True):
    """
    Tests polynomial integration using sparse and full grid methods for various polynomial degrees and levels of refinement.
    
    Parameters:
        n (int): Maximum degree of polynomial.
        L (list): List of refinement levels.
        saving (boolean): If True the plots will be saved as .png.
        
    Returns:
        DataFrame: A DataFrame containing the results of the integration tests.
    """
    X = [0, 1]  # Example interval [0, 1]
    Y = [0, 1]  # Example interval [0, 1]
    results = []

    for degree in range(n + 1):
        coefficients = random_polynomial(degree)
        f = polynomial_function(coefficients)
        exact_value = exact_integral(coefficients, X[0], X[1], Y[0], Y[1])
        
        for level in L:
            start_time = time.time()
            sparse_value = sparse_grid_integration(f, level, X, Y)
            sparse_time = time.time() - start_time

            start_time = time.time()
            full_value = full_grid_integration(f, level, X, Y)
            full_time = time.time() - start_time

            sparse_error = abs(sparse_value - exact_value)
            full_error = abs(full_value - exact_value)

            results.append({
                'Polynomial Degree': degree,
                'Level': level,
                'Exact Value': exact_value,
                'Sparse Value': sparse_value,
                'Full Value': full_value,
                'Sparse Error': sparse_error,
                'Full Error': full_error,
                'Sparse Time': sparse_time,
                'Full Time': full_time
            })
        
        print("Degree =", degree, "done.")

    df = pd.DataFrame(results)

    # Plot average absolute error for each level L over all polynomials
    avg_sparse_errors = df.groupby('Level')['Sparse Error'].mean()
    avg_full_errors = df.groupby('Level')['Full Error'].mean()
    avg_sparse_times = df.groupby('Level')['Sparse Time'].mean()
    avg_full_times = df.groupby('Level')['Full Time'].mean()

    plt.figure(figsize=(12, 8))
    plt.plot(avg_sparse_errors.index, avg_sparse_errors.values, label='Sparse Grid Error')
    plt.plot(avg_full_errors.index, avg_full_errors.values, label='Full Grid Error')
    plt.xlabel('Level', fontsize=20)
    plt.ylabel('Average Absolute Error', fontsize=20)
    plt.xticks(L, fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Average Absolute Error by Level', fontsize=20)
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.grid(True)
    if saving == True:
        plt.savefig("Average_Absolute_Error_by_Level_polynomial_test.png", dpi=500)
    plt.show()

    plt.figure(figsize=(12, 8))
    plt.plot(avg_sparse_times.index, avg_sparse_times.values, label='Sparse Grid Time')
    plt.plot(avg_full_times.index, avg_full_times.values, label='Full Grid Time')
    plt.xlabel('Level', fontsize=20)
    plt.ylabel('Average Time (s)', fontsize=20)
    plt.xticks(L, fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Average Time by Level', fontsize=20)
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.grid(True)
    if saving == True:
        plt.savefig("Average_Time_by_Level_polynomial_test.png", dpi=500)
    plt.show()

    max_level = max(L)
    max_level_data = df[df['Level'] == max_level]
    plt.figure(figsize=(24, 8))
    plt.plot(max_level_data['Polynomial Degree'], max_level_data['Sparse Error'], label='Sparse Grid Error')
    plt.plot(max_level_data['Polynomial Degree'], max_level_data['Full Error'], label='Full Grid Error')
    plt.xlabel('Polynomial Degree', fontsize=20)
    plt.ylabel('Absolute Error', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16)
    plt.title(f'Error by Polynomial Degree at Level {max_level}', fontsize=20)
    ax = plt.gca()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.grid(True)
    if saving == True:
        plt.savefig("Error_by_Polynomial_Degree_at_max_level_polynomial_test.png", dpi=500)
    plt.show()

    return df


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)


# Setting the seed value to 42
np.random.seed(42)


df = polynom_test(10, [2, 3, 4, 5, 6, 7, 8])
print(df)