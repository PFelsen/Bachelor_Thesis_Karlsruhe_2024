# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:43:51 2024

@author: Philipp Felsen
"""

"""
Performs quadrature based on hierarchical interpolation on a sparse grid for functions f:[a,b]x[c,d] \to \R.

"""

import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

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
        a (float): Left bound of the interval.
        b (float): Right bound of the interval.
        
    Returns:
        float: The value of phi(2^j * x - k).
    """
    return phi(((x - a) / (b - a)) * 2 ** j - k)

# Integral function
def int_phi(l, a, b):
    """
    Computes the integral values for the function phi_{l,i}(x).
    
    Parameters:
        l (int): Level.
        a (float): Left bound of the interval.
        b (float): Right bound of the interval.
        
    Returns:
        float: The value of the integral.
    """
    if l == 0:
        return (b - a) / 2
    if l >= 1:
        return (b - a) * 2 ** (-l)

# Computes hierarchical surplus
def delta(f, l_1, l_2, i_1, i_2, a, b, c, d):
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
    
    h_l_1 = (b - a) * 2 ** (-l_1)
    h_l_2 = (d - c) * 2 ** (-l_2)
    
    x = a + i_1 * h_l_1
    y = c + i_2 * h_l_2
    
    if l_1 == 0 and l_2 == 0:
        return f(x, y)
    
    if l_1 == 0 and l_2 >= 1:
        return f(x, y) - 0.5 * (f(x, y - h_l_2) + f(x, y + h_l_2))
    
    if l_2 == 0 and l_1 >= 1:
        return f(x, y) - 0.5 * (f(x - h_l_1, y) + f(x + h_l_1, y))
    
    erg += f(x, y) - 0.5 * (f(x, y - h_l_2) + f(x, y + h_l_2))
    erg += -0.5 * (f(x - h_l_1, y) - 0.5 * (f(x - h_l_1, y - h_l_2) + f(x - h_l_1, y + h_l_2)))
    erg += -0.5 * (f(x + h_l_1, y) - 0.5 * (f(x + h_l_1, y - h_l_2) + f(x + h_l_1, y + h_l_2)))
    
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
    if k == 0:
        return np.array([0, 1])
    
    for i in range(1, 2 ** k):
        if i % 2 != 0:
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
    
    l = L + 2 - 1 # L + d - 1 sparse grid construction
    
    for l_1 in range(0, l + 1):
        for l_2 in range(0, l - l_1 + 1):
            for i_1 in J(l_1):
                for i_2 in J(l_2):
                    erg += delta(f, l_1, l_2, i_1, i_2, a=a, b=b, c=c, d=d) * int_phi(l_1, a, b) * int_phi(l_2, c, d)
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
    
    l = L + 2 - 1 # L + d - 1 sparse grid construction
    
    for l_1 in range(0, l + 1):
        for l_2 in range(0, l + 1):
            for i_1 in J(l_1):
                for i_2 in J(l_2):
                    erg += delta(f, l_1, l_2, i_1, i_2, a=a, b=b, c=c, d=d) * int_phi(l_1, a, b) * int_phi(l_2, c, d)
    return erg

def data_output(f, L_values, X, Y, integral_value):
    """
    Generates the output data for the sparse and full grid integration.
    
    Parameters:
        f (function): The function to be integrated.
        L_values (array): Array of levels for the integration.
        X (array): Array containing the interval [a,b] in x-direction.
        Y (array): Array containing the interval [c,d] in y-direction.
        integral_value (float): The true value of the integral.
        
    Returns:
        tuple: Arrays containing the absolute errors for the sparse and full grid.
    """
    data = []
    
    abs_sparse_array = np.array([])
    abs_full_array = np.array([])
    
    for L in tqdm(L_values):
        sparse_value = sparse_grid_integration(f, L, X, Y)
        full_value = full_grid_integration(f, L, X, Y)
        
        error_abs_sparse = abs(sparse_value - integral_value)
        abs_sparse_array = np.append(abs_sparse_array, error_abs_sparse)
        
        error_abs_full = abs(full_value - integral_value)
        abs_full_array = np.append(abs_full_array, error_abs_full)
        
        error_rel_sparse = error_abs_sparse / integral_value
        error_rel_full = error_abs_full / integral_value
        
        # Add results to list
        data.append({
            "Level": L,
            "Sparse Grid Value": sparse_value,
            "Sparse Abs. Error": error_abs_sparse,
            "Sparse Rel. Error": error_rel_sparse,
            "Full Grid Value": full_value,
            "Full Abs. Error": error_abs_full,
            "Full Rel. Error": error_rel_full
        })
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Adjust DataFrame display settings
    pd.set_option('display.max_columns', None)  # Show all columns
    pd.set_option('display.width', 1000)        # Adjust the width of the output
    pd.set_option('display.colheader_justify', 'center')  # Center column headers
    
    # Print placeholder
    print(" ")
    # Print DataFrame
    print(df.to_string(index=False))
    print(" ")
    print("Full Grid abs. Error:", abs_full_array)
    print("Sparse Grid abs. Error:", abs_sparse_array)
    
    return abs_sparse_array, abs_full_array

def convergence_rate(L_values, abs_sparse_array, abs_full_array):
    """
    Computes the convergence rate for sparse and full grid integration.
    
    Parameters:
        L_values (array): Array of levels for the integration.
        abs_sparse_array (array): Array of absolute errors for the sparse grid.
        abs_full_array (array): Array of absolute errors for the full grid.
        
    Returns:
        tuple: Convergence rates for full and sparse grid.
    """
    # Log of the errors
    log_error_full = np.log(abs_full_array)
    log_error_sparse = np.log(abs_sparse_array)
    
    # Linear regression log(E(L)) = a*L + b
    coefficients_full = np.polyfit(L_values, log_error_full, 1)
    coefficients_sparse = np.polyfit(L_values, log_error_sparse, 1)
    
    # Slope of the line a (corresponds to -p * log(2))
    p_log_2_full = coefficients_full[0]
    p_log_2_sparse = coefficients_sparse[0]
    # Convergence order p
    p_full = -p_log_2_full / np.log(2)
    p_sparse = -p_log_2_sparse / np.log(2)
    # Output
    print("Full Grid Convergence Order p = :", p_full)
    print("Sparse Grid Convergence Order p = :", p_sparse)
    return p_full, p_sparse

# Functions to be tested ###
def sin_sin(x, y):
    return np.sin(x) * np.sin(y)

def rosenbrock(x, y):
    return (1 - x) ** 2 + 100 * (y - x ** 2) ** 2

def g(x, y):
    if 0.5 <= x <= 1:
        return 1.0
    else:
        return 0

def polynomial(x, y):
    return x ** 500 * y ** 17


"""
Begin of test cases.

"""


L_values = np.array([2, 3])
X = [0, 1]
Y = [0, 1]

print("sin(x)*sin(y)")
abs_sparse_array_sin_sin, abs_full_array_sin_sin = data_output(sin_sin, L_values, X, Y, (1 - np.cos(1)) ** 2)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_sin_sin, abs_full_array_sin_sin)

print("Rosenbrock")
abs_sparse_array_rose, abs_full_array_rose = data_output(rosenbrock, L_values, X, Y, 61 / 3)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_rose, abs_full_array_rose)

print("g")
abs_sparse_array_g, abs_full_array_g = data_output(g, L_values, X, Y, 1 / 2)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_g, abs_full_array_g)

print("x**500 * y**17")
abs_sparse_array_poly, abs_full_array_poly = data_output(polynomial, L_values, X, Y, 1 / 501 * 1 / 18)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_poly, abs_full_array_poly)

# Plot for sin(x)*sin(y)
plt.figure(figsize=(10, 8))
plt.plot(L_values, abs_sparse_array_sin_sin, label='Sparse Grid')
plt.plot(L_values, abs_full_array_sin_sin, label='Full Grid')
plt.yscale('log')
plt.title(r'$h(x,y) = \sin(x) \cdot \sin(y)$', fontsize=20)
plt.xlabel('L values', fontsize=20)
plt.ylabel('Absolute Error)', fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.show()

# Plot for Rosenbrock
plt.figure(figsize=(10, 8))
plt.plot(L_values, abs_sparse_array_rose, label='Sparse Grid')
plt.plot(L_values, abs_full_array_rose, label='Full Grid')
plt.yscale('log')
plt.title(r'Rosenbrock $f(x, y)$', fontsize=20)
plt.xlabel('L values', fontsize=20)
plt.ylabel('Absolute Error', fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.show()

# Plot for g
plt.figure(figsize=(10, 8))
plt.plot(L_values, abs_sparse_array_g, label='Sparse Grid')
plt.plot(L_values, abs_full_array_g, label='Full Grid')
plt.yscale('log')
plt.title(r'$g(x,y)$', fontsize=20)
plt.xlabel('L values', fontsize=20)
plt.ylabel('Absolute Error', fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.show()

# Plot for x**500 * y**17
plt.figure(figsize=(10, 8))
plt.plot(L_values, abs_sparse_array_poly, label='Sparse Grid')
plt.plot(L_values, abs_full_array_poly, label='Full Grid')
plt.yscale('log')
plt.title(r'$p(x,y) = x^{500} \cdot y^{17}$', fontsize=20)
plt.xlabel('L values', fontsize=20)
plt.ylabel('Absolute Error', fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.show()

# Additional: Chessboard function
L_values = np.array([2, 3, 4, 5])

def create_chessboard_pattern(n):
    """
    Creates a chessboard pattern of size 2^n x 2^n.
    
    Parameters:
        n (int): The exponent determining the size of the chessboard (2^n x 2^n).
        
    Returns:
        numpy.ndarray: A 2^n x 2^n array with a chessboard pattern where 1 represents a white square and 0 represents a black square.
    """
    size = 2 ** n
    chessboard = np.zeros((size, size), dtype=int)
    
    for i in range(size):
        for j in range(size):
            if (i + j) % 2 == 0:
                chessboard[i, j] = 1
            else:
                chessboard[i, j] = 0
                
    return chessboard

def chessboard_2(x, y):
    """
    Determines the value at coordinates (x, y) on a 4x4 chessboard.
    
    Parameters:
        x (float): The x-coordinate in the interval [0, 1].
        y (float): The y-coordinate in the interval [0, 1].
        
    Returns:
        int: The value at the coordinates (x, y) on the chessboard pattern (1 for white square, 0 for black square).
    """
    n = 2
    chessboard = create_chessboard_pattern(n)
    size = 2 ** n
    
    i = np.floor(x * size).astype(int) - 1
    j = np.floor(y * size).astype(int) - 1
    
    return chessboard[i, j]


# Function for n=4
def chessboard_4(x, y):
    """
    Determines the value at coordinates (x, y) on a 16x16 chessboard.
    
    Parameters:
        x (float): The x-coordinate in the interval [0, 1].
        y (float): The y-coordinate in the interval [0, 1].
        
    Returns:
        int: The value at the coordinates (x, y) on the chessboard pattern (1 for white square, 0 for black square).
    """
    n = 4
    chessboard = create_chessboard_pattern(n)
    size = 2 ** n
    
    i = np.floor(x * size).astype(int) - 1
    j = np.floor(y * size).astype(int) - 1
    
    return chessboard[i, j]



print("Chessboard: n=2")
abs_sparse_array_chess_2, abs_full_array_chess_2 = data_output(chessboard_2, L_values, X, Y, 1 / 2)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_chess_2, abs_full_array_chess_2)

print("Chessboard: n=4")
abs_sparse_array_chess_4, abs_full_array_chess_4 = data_output(chessboard_4, L_values, X, Y, 1 / 2)
p_full, p_sparse = convergence_rate(L_values, abs_sparse_array_chess_4, abs_full_array_chess_4)


