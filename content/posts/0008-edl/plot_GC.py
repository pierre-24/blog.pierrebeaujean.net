import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
import scipy.constants as const

# constants
e = const.e # (C)
k_B = const.k # (J/K)
T = 298.15 # (K)
eps_0 = const.epsilon_0  # (F/m)
eps_r = 78.4
eps = eps_0 * eps_r
N_A = const.N_A # (mol^-1)

z = 1  # Valency
c0_M = 0.01 # (mol/L)
psi_0 = 0.050 # (V)

# Derived quantities
n0 = c0_M * 1000 * N_A # Bulk number density (ions/m^3)
V_T = (k_B * T) / (z * e)
kappa = np.sqrt((2 * z**2 * e**2 * n0) / (eps * k_B * T)) # (m^-1)
lambda_D = 1 / kappa # (m)
psi_tilde_0 = psi_0 / V_T

print('lambda_D = {} nm'.format(lambda_D * 1e9))

# numerical solver
def pb_equation(x, y):
    return np.vstack((y[1], np.sinh(y[0])))

def boundary_conditions(ya, yb):
    return np.array([ya[0] - psi_tilde_0, yb[0]])

x_tilde = np.linspace(0, 10, 500)
y_guess = np.zeros((2, x_tilde.size))
y_guess[0] = psi_tilde_0 * np.exp(-x_tilde)
y_guess[1] = -psi_tilde_0 * np.exp(-x_tilde)

sol = solve_bvp(pb_equation, boundary_conditions, x_tilde, y_guess)
if not sol.success:
    raise RuntimeError("BVP solver failed: " + sol.message)

# Data processing
x_nm = (x_tilde * lambda_D) * 1e9          
psi_exact = sol.sol(x_tilde)[0] * V_T      
psi_dh = psi_0 * np.exp(-x_tilde)          

c_plus_exact = c0_M * np.exp(-psi_exact / V_T)   
c_minus_exact = c0_M * np.exp(psi_exact / V_T)   
c_plus_dh = c0_M * (1 - psi_dh / V_T)
c_minus_dh = c0_M * (1 + psi_dh / V_T)

# Plotting Results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), sharex=True)

# Plot 1: Electric Potential vs. Distance
ax1.plot(x_nm, psi_exact * 1000, 'b-', lw=2.5, label='Exact (Poisson-Boltzmann)')
ax1.plot(x_nm, psi_dh * 1000, 'r--', lw=2, label='Approx (Debye-Hückel)')
ax1.set_ylabel('Potential $\\psi(x)$ (mV)')
ax1.set_xlim(0, 10)
ax1.legend()

# Plot 2: Ion Concentration vs. Distance
ax2.plot(x_nm, c_minus_exact, 'b-', lw=2.5, label='Anions (-) Exact')
ax2.plot(x_nm, c_minus_dh, 'b--', lw=2, label='Anions (-) DH Approx')

ax2.plot(x_nm, c_plus_exact, 'r-', lw=2.5, label='Cations (+) Exact')
ax2.plot(x_nm, c_plus_dh, 'r--', lw=2, label='Cations (+) DH Approx')

ax2.axhline(c0_M, color='k', linestyle=':', label='Bulk Concentration ($n_0$)')
ax2.axhline(0, color='gray', linestyle='-') 

ax2.set_xlabel('Distance from electrode $x$ (nm)')
ax2.set_ylabel('Concentration (mol/L)')
ax2.set_xlim(0, 5 * lambda_D * 1e9)
ax2.set_ylim(-0.03, np.max(c_minus_exact) * 1.1) 
ax2.legend()

plt.tight_layout()
plt.savefig('GC.svg')
