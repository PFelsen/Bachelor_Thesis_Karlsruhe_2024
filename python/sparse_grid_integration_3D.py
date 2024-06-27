# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 16:45:02 2024

@author: Philipp Felsen
"""

"""
Performs quadrature based on hierarchical interpolation on a sparse grid for functions f:[a,b]x[c,d]x[e,g] \to \R.

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
def delta(f, l_1, l_2, l_3, i_1, i_2, i_3, a, b, c, d, e, g):
    """
    Computes the hierarchical surplus for a given function at specific levels and indices.
    
    Parameters:
        f (function): The function for which to compute the surplus.
        l_1 (int): The level in the first dimension.
        l_2 (int): The level in the second dimension.
        l_3 (int): The level in the third dimension.
        i_1 (int): The index in the first dimension.
        i_2 (int): The index in the second dimension.
        i_3 (int): The index in the third dimension.
        a (float): Left bound of the interval in x-direction.
        b (float): Right bound of the interval in x-direction.
        c (float): Left bound of the interval in y-direction.
        d (float): Right bound of the interval in y-direction.
        e (float): Left bound of the interval in z-direction.
        g (float): Right bound of the interval in z-direction.
        
    Returns:
        float: The hierarchical surplus value.
    """
    result = 0.0
    
    h_l_1 = (b - a) * 2 ** (-l_1)
    h_l_2 = (d - c) * 2 ** (-l_2)
    h_l_3 = (g - e) * 2 ** (-l_3)
    
    x = a + i_1 * h_l_1
    y = c + i_2 * h_l_2
    z = e + i_3 * h_l_3
    
    if l_1 == 0 and l_2 == 0 and l_3 == 0:
        return f(x, y, z)
    
    if l_1 >= 1 and l_2 == 0 and l_3 == 0:
        return f(x, y, z) - 0.5 * (f(x - h_l_1, y, z) + f(x + h_l_1, y, z))
    
    if l_1 == 0 and l_2 >= 1 and l_3 == 0:
        return f(x, y, z) - 0.5 * (f(x, y - h_l_2, z) + f(x, y + h_l_2, z))
    
    if l_1 == 0 and l_2 == 0 and l_3 >= 1:
        return f(x, y, z) - 0.5 * (f(x, y, z - h_l_3) + f(x, y, z + h_l_3))
    
    if l_1 >= 1 and l_2 >= 1 and l_3 == 0:
        result += f(x, y, z) - 0.5 * (f(x, y - h_l_2, z) + f(x, y + h_l_2, z))
        result += -0.5 * (f(x - h_l_1, y, z) - 0.5 * (f(x - h_l_1, y - h_l_2, z) + f(x - h_l_1, y + h_l_2, z)))
        result += -0.5 * (f(x + h_l_1, y, z) - 0.5 * (f(x + h_l_1, y - h_l_2, z) + f(x + h_l_1, y + h_l_2, z)))
        return result
    
    if l_1 == 0 and l_2 >= 1 and l_3 >= 1:
        result += f(x, y, z) - 0.5 * (f(x, y - h_l_2, z) + f(x, y + h_l_2, z))
        result += -0.5 * (f(x, y, z - h_l_3) - 0.5 * (f(x, y - h_l_2, z - h_l_3) + f(x, y + h_l_2, z - h_l_3)))
        result += -0.5 * (f(x, y, z + h_l_3) - 0.5 * (f(x, y - h_l_2, z + h_l_3) + f(x, y + h_l_2, z + h_l_3)))
        return result
    
    if l_1 >= 1 and l_2 == 0 and l_3 >= 1:
        result += f(x, y, z) - 0.5 * (f(x, y, z - h_l_3) + f(x, y, z + h_l_3))
        result += -0.5 * (f(x - h_l_1, y, z) - 0.5 * (f(x - h_l_1, y, z - h_l_3) + f(x - h_l_1, y, z + h_l_3)))
        result += -0.5 * (f(x + h_l_1, y, z) - 0.5 * (f(x + h_l_1, y, z - h_l_3) + f(x + h_l_1, y, z + h_l_3)))
        return result
    
    if l_1 >= 1 and l_2 >= 1 and l_3 >= 1:
        result += f(x, y, z) - 0.5 * (f(x, y - h_l_2, z) + f(x, y + h_l_2, z))
        result += -0.5 * (f(x, y, z - h_l_3) - 0.5 * (f(x, y - h_l_2, z - h_l_3) + f(x, y + h_l_2, z - h_l_3)))
        result += -0.5 * (f(x, y, z + h_l_3) - 0.5 * (f(x, y - h_l_2, z + h_l_3) + f(x, y + h_l_2, z + h_l_3)))
        
        result += -0.5 * (f(x - h_l_1, y, z) - 0.5 * (f(x - h_l_1, y - h_l_2, z) + f(x - h_l_1, y + h_l_2, z)))
        result += 0.25 * (f(x - h_l_1, y, z - h_l_3) - 0.5 * (f(x - h_l_1, y - h_l_2, z - h_l_3) + f(x - h_l_1, y + h_l_2, z - h_l_3)))
        result += 0.25 * (f(x - h_l_1, y, z + h_l_3) - 0.5 * (f(x - h_l_1, y - h_l_2, z + h_l_3) + f(x - h_l_1, y + h_l_2, z + h_l_3)))
        
        result += -0.5 * (f(x + h_l_1, y, z) - 0.5 * (f(x + h_l_1, y - h_l_2, z) + f(x + h_l_1, y + h_l_2, z)))
        result += 0.25 * (f(x + h_l_1, y, z - h_l_3) - 0.5 * (f(x + h_l_1, y - h_l_2, z - h_l_3) + f(x + h_l_1, y + h_l_2, z - h_l_3)))
        result += 0.25 * (f(x + h_l_1, y, z + h_l_3) - 0.5 * (f(x + h_l_1, y - h_l_2, z + h_l_3) + f(x + h_l_1, y + h_l_2, z + h_l_3)))
        
        return result

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

def sparse_grid_integration(f, L, X, Y, Z):
    """
    Computes the integral for a given function in three dimensions over [a,b]x[c,d]x[e,g].
    
    Parameters:
        f (function): The function to be integrated. Must accept three arguments (x, y, z).
        L (int): Global level of refinement.
        X (array): An array containing the interval [a,b] in x-direction.
        Y (array): An array containing the interval [c,d] in y-direction.
        Z (array): An array containing the interval [e,g] in z-direction.
      
    Returns:
        float: The value of the integral.
    """
    erg = 0
    
    a = X[0]
    b = X[1]
    c = Y[0]
    d = Y[1]
    e = Z[0]
    g = Z[1]
    
    l = L + 3 - 1 # L + d - 1 sparse grid construction
    
    for l_1 in range(0, l + 1):
        for l_2 in range(0, l - l_1 + 1):
            for l_3 in range(0, l - l_1 - l_2 + 1):
                for i_1 in J(l_1):
                    for i_2 in J(l_2):
                        for i_3 in J(l_3):
                            erg += delta(f, l_1, l_2, l_3, i_1, i_2, i_3, a=a, b=b, c=c, d=d, e=e, g=g) * int_phi(l_1, a, b) * int_phi(l_2, c, d) * int_phi(l_3, e, g)
    return erg


def function_test():
    """
    Tests the sparse grid integration with random polynomial, exponential, and trigonometric functions.
    """
    L_values = [6, 7, 8, 9]
    
    # Random coefficients
    a = random.randint(-4, 4)
    b = random.randint(-4, 4)
    c = random.randint(-4, 4)
    
    # Random exponents
    k = random.randint(1, 3)
    n = random.randint(1, 3)
    j = random.randint(1, 3)
    
    # Random integration domain
    step_x = random.randint(1, 2)
    step_y = random.randint(1, 2)
    step_z = random.randint(1, 2)
    
    left_x = random.randint(0, 5)
    left_y = random.randint(0, 5)
    left_z = random.randint(0, 5)
    
    right_x = left_x + step_x
    right_y = left_y + step_y
    right_z = left_z + step_z
    
    X = [left_x, right_x]
    Y = [left_y, right_y]
    Z = [left_z, right_z]
    
    print("Integration domain:")
    print("X-axis:", X)
    print("Y-axis:", Y)
    print("Z-axis:", Z)
    
    # Random polynomial function
    def f(x, y, z):
        x = float(x)
        y = float(y)
        z = float(z)
        return a * x ** k * b * y ** n * c * z ** j
    
    print("f(x,y,z) =", a, "x**", k, "*", b, "y**", n, "*", c, "z**", j)
    
    # Compute the exact integral value for the polynomial function
    real_value_f = (a / (k + 1) * (right_x ** (k + 1) - left_x ** (k + 1)) *
                    b / (n + 1) * (right_y ** (n + 1) - left_y ** (n + 1)) *
                    c / (j + 1) * (right_z ** (j + 1) - left_z ** (j + 1)))
    
    print("Exact integral value:", real_value_f)
    
    for L in L_values:
        quad = sparse_grid_integration(f, L, X, Y, Z)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_f - quad))
    
    # Test with exponential function
    def g(x, y, z):
        return a * np.exp(2 * x)
    
    # Compute the exact integral value for the exponential function
    real_value_g = (a / 2 * (np.exp(2 * right_x) - np.exp(2 * left_x)) *
                    (right_y - left_y) * (right_z - left_z))
    
    print("g(x,y,z) =", a, "*", "exp(2x)")
    print("Exact integral value:", real_value_g)
    for L in L_values:
        quad = sparse_grid_integration(g, L, X, Y, Z)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_g - quad))
    
    # Test with trigonometric function
    def h(x, y, z):
        return np.cos(a * x) * np.cos(b * y) * np.cos(c * z)
    
    # Compute the exact integral value for the trigonometric function
    real_value_h = (1 / a * (np.sin(a * right_x) - np.sin(a * left_x)) *
                    1 / b * (np.sin(b * right_y) - np.sin(b * left_y)) *
                    1 / c * (np.sin(c * right_z) - np.sin(c * left_z)))
    
    print("h(x,y,z) = cos(", a, "x)*cos(", b, "y)*cos(", c, "z)")
    print("Exact integral value:", real_value_h)
    for L in L_values:
        quad = sparse_grid_integration(h, L, X, Y, Z)
        print("L=", L, " Quadrature:", quad, " Abs. error:", abs(real_value_h - quad))

function_test()
