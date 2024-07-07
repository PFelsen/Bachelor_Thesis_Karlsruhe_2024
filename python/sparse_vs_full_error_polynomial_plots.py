# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 21:11:56 2024

@author: phili
"""


import numpy as np
import matplotlib.pyplot as plt

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
        L (tuple): A tuple containing the maximum levels of refinement in each dimension (L[0] for x, L[1] for y).
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



def integral_polynom(a, b):
    if a <= -1 or b <= -1:
        raise ValueError("a and b must be greater than -1")
    integral_value = 1 / ((a + 1) * (b + 1))
    return integral_value

# Set up parameters
L = 7
polynomial_degrees = range(2,101)
X = [0, 1]
Y = [0, 1]

# Initialize lists to store errors
sparse_errors = []
full_errors = []

# Compute errors
for degree in polynomial_degrees:
    f = lambda x, y: (x * y) ** degree
    exact_value = integral_polynom(degree, degree)
    sparse_value = sparse_grid_integration(f, L, X, Y)
    full_value = full_grid_integration(f, L, X, Y)
    sparse_errors.append(abs(sparse_value - exact_value))
    full_errors.append(abs(full_value - exact_value))

# Plot errors
plt.figure(figsize=(24, 8))

plt.plot(polynomial_degrees, sparse_errors, marker='o', label='Sparse Grid')
plt.plot(polynomial_degrees, full_errors, marker='x', label='Full Grid')

plt.xlabel('Polynomial Degree', fontsize=label_fontsize)
plt.ylabel('Absolute Quadrature Error (log)', fontsize=label_fontsize)
plt.yscale('log')  # Set y-axis to logarithmic scale
plt.title('Quadrature Error for (x*y)^degree with L=7', fontsize=label_fontsize)
plt.legend(fontsize=label_fontsize)
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(tick_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.grid(True)
plt.savefig("polynomial_error_analysis.png", dpi=500)
plt.show()
