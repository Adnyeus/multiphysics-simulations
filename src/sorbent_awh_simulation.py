"""
Solar-Thermal Sorbent Bed Dynamics for Atmospheric Water Harvesting (AWH)
Author: Ebad Naeem
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def tetens_p_sat(T_k: float) -> float:
    """Calculates saturation vapor pressure (Pa) using Tetens Equation."""
    T_c = T_k - 273.15
    return 610.78 * np.exp((17.27 * T_c) / (T_c + 237.3))


def awh_system_ode(t: float, y: list, p: dict) -> list:
    """Coupled non-linear ODEs for Sorbent Mass Loading q(t) and Temperature T(t)."""
    q, T = y
    
    # Environment profile (Overnight Adsorption vs. Daytime Desorption)
    if t < p['t_desorb_start']:
        T_amb = p['T_amb_night']
        RH = p['RH_night']
        I_solar = 0.0
    else:
        T_amb = p['T_amb_day']
        RH = p['RH_day']
        I_solar = p['I_solar_max'] * np.sin(np.pi * (t - p['t_desorb_start']) / 8.0)

    P_sat = tetens_p_sat(T)
    P_v = RH * tetens_p_sat(T_amb)
    
    # Equilibrium uptake (Langmuir approximation)
    q_eq = p['q_max'] * (p['b_0'] * P_v) / (1.0 + p['b_0'] * P_v)
    
    # Mass balance (LDF model)
    dqdt = p['k_m'] * (q_eq - q)
    
    # Energy balance
    m_bed = p['m_sorbent']
    c_p_bed = p['cp_sorbent']
    Q_solar = p['alpha_abs'] * p['A_bed'] * I_solar
    Q_conv = p['h_conv'] * p['A_bed'] * (T - T_amb)
    Q_ads = m_bed * p['delta_H_ads'] * dqdt
    
    dTdt = (Q_solar - Q_conv + Q_ads) / (m_bed * c_p_bed)
    return [dqdt, dTdt]


def run_simulation():
    params = {
        'q_max': 0.45, 'b_0': 1.0e-3, 'k_m': 1.2e-4, 'm_sorbent': 1.5,
        'cp_sorbent': 920.0, 'delta_H_ads': 2.8e6, 'A_bed': 1.0,
        'alpha_abs': 0.92, 'h_conv': 12.0, 'T_amb_night': 293.15,
        'RH_night': 0.75, 'T_amb_day': 308.15, 'RH_day': 0.30,
        'I_solar_max': 800.0, 't_desorb_start': 8.0 * 3600.0
    }

    t_span = (0.0, 16.0 * 3600.0)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    y0 = [0.05, params['T_amb_night']]

    sol = solve_ivp(
        fun=lambda t, y: awh_system_ode(t, y, params),
        t_span=t_span, y0=y0, t_eval=t_eval, method='RK45'
    )

    t_hours = sol.t / 3600.0
    q_water = sol.y[0]
    T_bed = sol.y[1] - 273.15

    # Plot results
    os.makedirs('images', exist_ok=True)
    fig, ax1 = plt.subplots(figsize=(8, 5))

    color = 'tab:blue'
    ax1.set_xlabel('Time (Hours)', fontweight='bold')
    ax1.set_ylabel('Sorbent Moisture Loading q (kg/kg)', color=color, fontweight='bold')
    ax1.plot(t_hours, q_water, color=color, linewidth=2.5, label='Moisture Uptake q(t)')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--', alpha=0.5)

    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Bed Temperature (°C)', color=color, fontweight='bold')
    ax2.plot(t_hours, T_bed, color=color, linewidth=2.0, linestyle='--', label='Temperature T(t)')
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('AWH Diurnal Adsorption-Desorption Kinetics', fontsize=12, fontweight='bold')
    fig.tight_layout()
    plt.savefig('images/viz1_solar_awh_kinetics.png', dpi=300)
    print("AWH Simulation Completed. Figure saved to images/viz1_solar_awh_kinetics.png")


if __name__ == '__main__':
    run_simulation()