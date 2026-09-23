"""
2D Short-Channel MOSFET TCAD Simulation & DIBL Analysis
Author: Ebad Naeem
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def solve_mosfet_2d(V_ds: float, V_gs: float = 1.0) -> np.ndarray:
    Nx, Ny = 60, 40
    dx, dy = 1.0e-9, 1.0e-9
    
    # Permittivity matrix setup
    eps_0 = 8.854e-12
    eps_r = np.full((Nx, Ny), 11.7)  # Silicon
    eps_r[:, 30:] = 3.9               # Oxide layer
    
    N = Nx * Ny
    A = lil_matrix((N, N))
    b = np.zeros(N)

    def get_idx(i, j):
        return i * Ny + j

    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            idx = get_idx(i, j)
            
            # Boundary conditions (Contacts)
            if j == Ny - 1:  # Gate contact
                A[idx, idx] = 1.0
                b[idx] = V_gs
            elif i < 10 and j < 20:  # Source contact
                A[idx, idx] = 1.0
                b[idx] = 0.0
            elif i > Nx - 10 and j < 20:  # Drain contact
                A[idx, idx] = 1.0
                b[idx] = V_ds
            else:
                # 5-Point Laplacian Stencil Operator
                A[idx, idx] = -4.0
                A[idx, get_idx(i + 1, j)] = 1.0
                A[idx, get_idx(i - 1, j)] = 1.0
                A[idx, get_idx(i, j + 1)] = 1.0
                A[idx, get_idx(i, j - 1)] = 1.0
                b[idx] = 0.0

    phi_sparse = spsolve(A.tocsr(), b)
    return phi_sparse.reshape((Nx, Ny))


def run_simulation():
    phi_low = solve_mosfet_2d(V_ds=0.05, V_gs=1.0)
    phi_high = solve_mosfet_2d(V_ds=1.20, V_gs=1.0)

    os.makedirs('images', exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    c1 = ax1.contourf(phi_low.T, cmap='magma', levels=20)
    fig.colorbar(c1, ax=ax1)
    ax1.set_title('Electrostatic Potential @ V_DS = 0.05V', fontweight='bold')
    ax1.set_xlabel('X Domain (nm)')
    ax1.set_ylabel('Y Domain (nm)')

    # Conduction channel profile
    y_channel = 15
    ax2.plot(phi_low[:, y_channel], label='V_DS = 0.05 V (Low Bias)', color='blue', linewidth=2)
    ax2.plot(phi_high[:, y_channel], label='V_DS = 1.20 V (High Bias - Barrier Lowered)', color='red', linestyle='--', linewidth=2)
    ax2.set_title('Channel Barrier Lowering (DIBL)', fontweight='bold')
    ax2.set_xlabel('Channel Distance (nm)', fontweight='bold')
    ax2.set_ylabel('Potential Phi (V)', fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()

    plt.tight_layout()
    plt.savefig('images/viz3_mosfet_2d_dibl_landscape.png', dpi=300)
    print("MOSFET TCAD Simulation Completed. Figure saved to images/viz3_mosfet_2d_dibl_landscape.png")


if __name__ == '__main__':
    run_simulation()