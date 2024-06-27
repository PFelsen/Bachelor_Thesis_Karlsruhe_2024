
# Sparse Grid Integration

This repository contains Python programs to test and compare the integration of various mathematical functions using sparse grid and full grid quadrature methods. The scripts provided are designed to handle integration in both two-dimensional and three-dimensional spaces.

## File Descriptions

### 1. `2D_sparse_grid_integration.py`

This script performs sparse grid integration for random polynomial, exponential, and trigonometric functions in two dimensions. The main function in this script is `function_test()`.

- **Random Coefficients and Exponents**: The script generates random coefficients and exponents for the functions.
- **Integration Domain**: The integration domain is randomly chosen within specified bounds.
- **Functions Tested**:
  - Polynomial: \( f(x, y) = a \cdot x^k \cdot b \cdot y^n \)
  - Exponential: \( g(x, y) = a \cdot \exp(2x) \)
  - Trigonometric: \( h(x, y) = \cos(a \cdot x) \cdot \cos(b \cdot y) \)
- **Integration Method**: Sparse grid integration is used to compute the integral for different levels of discretization \( L \).
- **Output**: The script prints the integration domain, the exact integral values, and the computed quadrature values along with their absolute errors for each function.

### 2. `3D_sparse_grid_integration.py`

This script extends the sparse grid integration to three dimensions. Similar to the 2D script, it tests random polynomial, exponential, and trigonometric functions but in three dimensions.

- **Random Coefficients and Exponents**: Generates random coefficients and exponents for the functions in three dimensions.
- **Integration Domain**: The integration domain is chosen randomly within specified bounds.
- **Functions Tested**:
  - Polynomial: \( f(x, y, z) = a \cdot x^k \cdot b \cdot y^n \cdot c \cdot z^m \)
  - Exponential: \( g(x, y, z) = a \cdot \exp(2x + 3y + 4z) \)
  - Trigonometric: \( h(x, y, z) = \cos(a \cdot x) \cdot \cos(b \cdot y) \cdot \cos(c \cdot z) \)
- **Integration Method**: Uses sparse grid integration to compute the integral for different levels of discretization \( L \).
- **Output**: Prints the integration domain, exact integral values, computed quadrature values, and their absolute errors for each function.

### 3. `sparse_grid_vs_full_grid_quadrature.py`

This script compares the efficiency and accuracy of sparse grid quadrature against full grid quadrature.

- **Integration Methods**:
  - Sparse Grid Quadrature: Uses sparse grids for multidimensional integration.
  - Full Grid Quadrature: Uses full tensor grids for multidimensional integration.
- **Functions Tested**: Can be configured to test various polynomial, exponential, and trigonometric functions.
- **Performance Metrics**:
  - Accuracy: Compares the exact integral value with the computed quadrature values.
  - Efficiency: Evaluates the computational time and resource usage for both methods.
- **Output**: Provides a detailed comparison of quadrature values, errors, and performance metrics for both sparse grid and full grid methods.

