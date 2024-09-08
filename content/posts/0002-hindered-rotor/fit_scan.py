import pandas
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy

data = pandas.read_csv('ethane_scan.csv')

figure = plt.figure(figsize=(10, 5))
ax1 = figure.subplots()

def V(theta: float, barrier: float, sigma: int = 3) -> float:
    return barrier/2 * (1-numpy.cos(sigma*theta))

popt, _ = curve_fit(lambda theta, barrier: V(theta,barrier), numpy.radians(data['angle'] + 60), data['rel_energy'])

ax1.plot(data['angle'] + 60, data['rel_energy'], 'ob')
T = numpy.arange(0, 360)
ax1.plot(T, V(numpy.radians(T), popt[0]), color='blue', label='V($\\theta)$, with $V_0$={:.2f} kJ mol$^{{-1}}$'.format(popt[0]))

ax1.set_xlabel('Angle (°)')
ax1.set_ylabel('Energy (kJ mol$^{-1}$)')
ax1.set_ylim(0,13)

ax1.legend()

plt.tight_layout()
plt.show()
figure.savefig('ethane_scan.svg')
