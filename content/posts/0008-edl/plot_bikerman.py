import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const


# constants
e = const.e # (C)
k_B = const.k # (J/K)
T = 298.15  # (K)
eps_0 = const.epsilon_0 # (F/m)
N_A = const.N_A # (mol^-1)

z = 1  # 1:1 electrolyte
c0_M = 0.01 # Bulk concentration (M)
v = 4/3*np.pi*(0.8e-9)**3 # Volume of the ion (m^3)

# Permittivity
eps_bulk = eps_0 * 78.4 # Bulk water

# Derived thermodynamic quantities
n0 = c0_M * 1000 * N_A  # Number density (m^-3)
V_T = (k_B * T) / (z * e)
kappa = np.sqrt((2 * z**2 * e**2 * n0) / (eps_bulk * k_B * T)) # (m^-1)
C_GC_0 = eps_bulk * kappa  # Baseline GC capacitance (F/m^2)
Phi = 2 * v * n0 # Packing parameter Phi_0 (total volume fraction of ions in the bulk)

psi_0_vec = np.linspace(-0.6, 0.6, 1000)
u_vec = psi_0_vec / V_T

print(k_B * T / (z * e) * np.log(1/(v*n0)))

# Gouy-Chapman Capacitance
C_gouy_chapman = (C_GC_0 * np.cosh(u_vec / 2)) * 1e2  # Convert to uF/cm^2

#  Bikerman Capacitance
C_bikerman = np.zeros_like(psi_0_vec)

safe_mask = np.abs(u_vec) > 1e-5 # We must avoid division by zero at the PZC (u = 0)
u_safe = u_vec[safe_mask]

denom_inner = 1.0 + Phi * (np.cosh(u_safe) - 1.0)
numerator = np.sqrt(Phi / 2.0) * np.abs(np.sinh(u_safe))
denominator = denom_inner * np.sqrt(np.log(denom_inner))

C_bik_F_m2 = C_GC_0 * (numerator / denominator)
C_bikerman[safe_mask] = C_bik_F_m2 * 1e2  # Convert to uF/cm^2
C_bikerman[~safe_mask] = C_GC_0 * 1e2

# plot
plt.figure(figsize=(8, 6))

plt.plot(psi_0_vec, C_gouy_chapman, 'g-.', lw=2.5, label='Gouy-Chapman')
plt.plot(psi_0_vec, C_bikerman, 'b-', lw=2.5, label=f'Bikerman (Steric Volume $v={v*1e27:.3f}$ nm$^3$)')

plt.axvline(0, color='gray', linestyle=':', label='PZC ($\\psi_0 = 0$)')

plt.xlabel('Electrode Potential $\\psi_0 - \\psi_{\\mathrm{PZC}}$ (V)')
plt.ylabel('Differential Capacitance $C$ ($\\mu\\mathrm{F/cm}^2$)')

plt.xlim(-0.6, 0.6)
plt.ylim(0, max(np.max(C_bikerman) * 1.5, 30))
plt.legend()

plt.tight_layout()
plt.savefig('bikerman.svg')
