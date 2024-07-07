# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 00:08:28 2024

@author: Philipp Felsen
"""

"""
Performs quadrature based on hierarchical interpolation on a sparse grid for functions f:[a,b]x[c,d] \to \R.

"""

import numpy as np
import random

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

# Integral Function
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

def function_test():
    """
    Tests the sparse grid integration with random polynomial, exponential, and trigonometric functions in two dimensions.
    """
    L_values = [8, 10, 12]
    
    # Random coefficients
    a = random.randint(-4, 4)
    while a == 0:  # Ensure a is not zero
        a = random.randint(-4, 4)
    b = random.randint(-4, 4)
    while b == 0:  # Ensure b is not zero
        b = random.randint(-4, 4)
    
    # Random exponents
    k = random.randint(1, 3)
    n = random.randint(1, 3)
    
    # Random integration domain
    step_x = random.randint(1, 2)
    step_y = random.randint(1, 2)
    
    left_x = random.randint(0, 5)
    left_y = random.randint(0, 5)
    
    right_x = left_x + step_x
    right_y = left_y + step_y
    
    X = [left_x, right_x]
    Y = [left_y, right_y]
    
    print("Integration domain:")
    print("X-axis:", X)
    print("Y-axis:", Y)
    
    # Random polynomial function
    def f(x, y):
        x = float(x)
        y = float(y)
        return a * x ** k * b * y ** n
    
    print("f(x,y) =", a, "x**", k, "*", b, "y**", n)
    
    # Compute the exact integral value for the polynomial function
    real_value_f = (a / (k + 1) * (right_x ** (k + 1) - left_x ** (k + 1)) *
                    b / (n + 1) * (right_y ** (n + 1) - left_y ** (n + 1)))
    
    print("Exact integral value:", real_value_f)
    
    for L in L_values:
        quad = sparse_grid_integration(f, L, X, Y)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_f - quad))
    
    # Test with exponential function
    def g(x, y):
        return a * np.exp(2 * x)
    
    # Compute the exact integral value for the exponential function
    real_value_g = (a / 2 * (np.exp(2 * right_x) - np.exp(2 * left_x)) *
                    (right_y - left_y))
    
    print("g(x,y) =", a, "*", "exp(2x)")
    print("Exact integral value:", real_value_g)
    for L in L_values:
        quad = sparse_grid_integration(g, L, X, Y)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_g - quad))
    
    # Test with trigonometric function
    def h(x, y):
        return np.cos(a * x) * np.cos(b * y)
    
    # Compute the exact integral value for the trigonometric function
    real_value_h = (1 / a * (np.sin(a * right_x) - np.sin(a * left_x)) *
                    1 / b * (np.sin(b * right_y) - np.sin(b * left_y)))
    
    print("h(x,y) = cos(", a, "x)*cos(", b, "y)")
    print("Exact integral value:", real_value_h)
    for L in L_values:
        quad = sparse_grid_integration(h, L, X, Y)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_h - quad))

function_test()
