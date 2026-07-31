import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import scipy.constants as const

# constants
e = const.e # (C)
k_B = const.k # (J/K)
T = 298.15  # (K)
eps_0 = const.epsilon_0 # (F/m)
N_A = const.N_A # (mol^-1)

z = 1  # 1:1 electrolyte
c0_M = 0.01 # Bulk concentration (M)

# Permittivities and Distances
eps_bulk = eps_0 * 78.4 # Bulk water
eps_stern = eps_0 * 15.0  # Stern layer effective permittivity
d_helm = 0.5e-9  # Stern layer thickness (m)

# Derived thermodynamic quantities
n0 = c0_M * 1000 * N_A  # Number density (m^-3)
V_T = (k_B * T) / (z * e)
kappa = np.sqrt((2 * z**2 * e**2 * n0) / (eps_bulk * k_B * T)) # (m^-1)
C_GC_0 = eps_bulk * kappa  # Baseline GC capacitance (F/m^2)

psi_0_vec = np.linspace(-1, 1, 800)

def sigma_GC(psi):
    """Diffuse layer surface charge density (C/m^2)"""
    return np.sqrt(8 * eps_bulk * n0 * k_B * T) * np.sinh(psi / (2 * V_T))

def C_GC(psi):
    """Diffuse layer differential capacitance (F/m^2)"""
    return C_GC_0 * np.cosh(psi / (2 * V_T))

# A. Helmholtz (Constant)
C_H_val = (eps_stern / d_helm) * 1e2  # Convert F/m^2 to uF/cm^2
print(C_H_val)
C_helmholtz = np.full_like(psi_0_vec, C_H_val)

# B. Gouy-Chapman (Pure Diffuse)
C_gouy_chapman = C_GC(psi_0_vec) * 1e2

# C. Stern Model
C_stern = np.zeros_like(psi_0_vec)
for i, psi_0 in enumerate(psi_0_vec):
    # Charge balance: sigma_compact = sigma_diffuse
    def stern_balance(psi_OHP):
        sigma_comp = (eps_stern / d_helm) * (psi_0 - psi_OHP)
        return sigma_comp - sigma_GC(psi_OHP)
    
    # Root finding bracketed widely to ensure convergence
    res = root_scalar(stern_balance, bracket=[-1.0, 1.0])
    psi_OHP = res.root
    
    # Series formula
    C_s = 1.0 / ( (1.0 / (eps_stern / d_helm)) + (1.0 / C_GC(psi_OHP)) )
    C_stern[i] = C_s * 1e2

# D. Grahame model
C_grahame = np.zeros_like(psi_0_vec)

def K_comp(dpsi):
    base = 25.0 * 1e2  
    asym = 4.0 * 1e2 * np.tanh(4.0 * dpsi)
    hump = 12.0 * 1e2 * np.exp(-((dpsi - 0.15) / 0.10)**2)
    return base + asym + hump

def C_comp_diff(dpsi):
    K = K_comp(dpsi)
    # dK/d(dpsi)
    dK = (4.0 * 1e2 * 4.0 / (np.cosh(4.0 * dpsi)**2) - 12.0 * 1e2 * 2 * (dpsi - 0.15) / (0.10**2) * np.exp(-((dpsi - 0.15) / 0.10)**2))
    return K + dpsi * dK

for i, psi_0 in enumerate(psi_0_vec):
    def grahame_balance(psi_OHP):
        dpsi = psi_0 - psi_OHP
        sigma_comp = K_comp(dpsi) * dpsi
        return sigma_comp - sigma_GC(psi_OHP)

    res_g = root_scalar(grahame_balance, bracket=[-1.0, 1.0])
    psi_OHP = res_g.root
    dpsi = psi_0 - psi_OHP
    
    # Calculate exact differential compact capacitance at this potential drop
    C_compact_val = C_comp_diff(dpsi)
    
    # True Grahame Series Formula: 1/C_G = 1/C_compact + 1/C_GC(psi_OHP)
    C_g = 1.0 / ( (1.0 / C_compact_val) + (1.0 / C_GC(psi_OHP)) )
    C_grahame[i] = C_g * 1e-2

# Plot
plt.figure(figsize=(8, 6))

plt.plot(psi_0_vec, C_helmholtz, 'r--', lw=2, label='Helmholtz')
plt.plot(psi_0_vec, C_gouy_chapman, 'g-.', lw=2, label='Gouy-Chapman')
plt.plot(psi_0_vec, C_stern, 'b-', lw=2.5, label='Stern Model')
plt.plot(psi_0_vec, C_grahame, 'm-', lw=2.5, label=r'Grahame')

plt.axvline(0, color='gray', linestyle=':', label='PZC ($\\psi_0 = 0$)')

plt.xlabel('Electrode Potential $\\psi_0 - \\psi_{\\mathrm{PZC}}$ (V)')
plt.ylabel('Differential Capacitance $C$ ($\\mu\\mathrm{F/cm}^2$)')

plt.xlim(-1, 1)
plt.ylim(0, 50) # let GC shoot off-screen naturally
plt.legend()

plt.tight_layout()
plt.savefig('capa.svg')
