#!/usr/bin/env/ python3
"""
Read the binary temperature field and produce a contour plot. 
File format:
    - 2 x int32 : nx, ny
    - nx*ny * float64 : temperature data (column-major order)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def read_temperature(filename):
    with open(filename, 'rb') as f:
        # Read Header
        nx = np.fromfile(f, dtype=np.int32, count=1)[0]
        ny = np.fromfile(f, dtype=np.int32, count=1)[0]
        print(f"Grid size: {nx} x {ny}")

        # Read data
        data = np.fromfile(f, dtype=np.float64, count=nx*ny)

        # Reshape to 2D (column-major: first index is x, second is y)
        # Fortran order: (nx, ny) --> reshape with order = 'F'
        T = data.reshape((nx,ny), order='F')

        # Transpose so that y is rows , x is columns for plotting
        T = T.T     # now shpae (ny, nx)
        return T, nx, ny

def plot_temperature(T, nx, ny, output_file=None):
    fig, ax = plt.subplots(figsize=(8,6))
    # Create coordinate arrays
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    # Filled contour plot
    c = ax.contourf(x,y, T, levels=20, cmap='hot')
    fig.colorbar(c, label="Temperature")

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Steady-State Temperature Field (Heat Equation)')

    if output_file:
        fig.savefig(output_file, dpi=150)
        print(f"Plot saved to {output_file}")
    else: 
        plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 2: 
        filename = 'temperature.bin'
    else:
        filename = sys.argv[1]


    T, nx, ny = read_temperature(filename)
    plot_temperature(T, nx, ny, output_file='temperature.png')
