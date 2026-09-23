"""
Dynamic Electrostatic MEMS Comb-Drive Actuator & Pull-In Analysis
Author: Ebad Naeem
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def mems_ode(t: float, state: list, p: dict, V_app: float) -> list:
    """State-Space representation [x, v] of parallel-plate / comb-drive actuator."""
    x, v = state
    g0 = p['g0']
    
    # Avoid mathematical singularity at contact
    gap = max(g0 - x, 1.0e-9)
    
    F_spring = p['k'] * x
    F_damping = p['b'] * v
    F_elec = (p['eps0'] * p['Area'] * V_app**2) / (2.0 * gap**2)
    
    dxdt = v
    dvdt = (F_elec - F_spring - F_damping) / p['m']
    return [dxdt, dvdt]


def run_simulation():
    p = {
        'eps0': 8.854e-12, 'Area': 100.0e-8, 'g0': 3.0e-6,
        'k': 0.5, 'm': 1.0e-6, 'b': 2.0e-4
    }

    # Theoretical Pull-In Voltage
    V_pi = np.sqrt((8.0 / 27.0) * (p['k'] * p['g0']**3) / (p['eps0'] * p['Area']))
    print(f"Calculated Theoretical Pull-In Voltage: {V_pi:.2f} V")

    t_span = (0.0, 2.0e-3)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)

    sol_stable = solve_ivp(
        fun=lambda t, y: mems_ode(t, y, p, V_app=0.8 * V_pi),
        t_span=t_span, y0=[0.0, 0.0], t_eval=t_eval, method='RK45'
    )
    
    sol_pullin = solve_ivp(
        fun=lambda t, y: mems_ode(t, y, p, V_app=1.15 * V_pi),
        t_span=t_span, y0=[0.0, 0.0], t_eval=t_eval, method='RK45'
    )

    os.makedirs('images', exist_ok=True)
    plt.figure(figsize=(8, 5))
    plt.plot(sol_stable.t * 1e3, sol_stable.y[0] * 1e6, label=f'Sub-Critical (0.8 V_pi)', color='blue', linewidth=2)
    plt.plot(sol_pullin.t * 1e3, sol_pullin.y[0] * 1e6, label=f'Super-Critical (1.15 V_pi - Pull-In)', color='red', linestyle='--', linewidth=2)
    plt.axhline(y=(p['g0'] / 3.0) * 1e6, color='black', linestyle=':', label='Pull-In Limit (g0/3)')

    plt.title('MEMS Comb-Drive Dynamic Step Response & Snap-Through', fontsize=12, fontweight='bold')
    plt.xlabel('Time (ms)', fontweight='bold')
    plt.ylabel('Displacement x (µm)', fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig('images/viz2_mems_comb_drive_step.png', dpi=300)
    print("MEMS Simulation Completed. Figure saved to images/viz2_mems_comb_drive_step.png")


if __name__ == '__main__':
    run_simulation()