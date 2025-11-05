import numpy
import matplotlib.pyplot as plt

Tt_A = 1.2 # in T/T_0
Tt_B = 1 # in T/T_0
S = 7 # R

figure = plt.figure(figsize=(5, 10))
ax1, ax2 = figure.subplots(2)

def dGt_A(T: float) -> float:
    return Tt_A*S - T * S  # in RT/T_0
    
def dGt_B(T: float) -> float:
    return Tt_B*S - T * S  # in RT/T_0
    
def G_AB_alpha(x: float, T:float) ->float:
    return dGt_B(T) + (dGt_A(T)-dGt_B(T))*x + T*(x*numpy.log(x)+(1-x)*numpy.log(1-x)) # in RT/T_0
    
def G_AB_beta(x: float, T:float) ->float:
    return T*(x*numpy.log(x)+(1-x)*numpy.log(1-x)) # in RT/T_0

# free energy at Tx
Tx = 1.05 # T/T_0

x = numpy.linspace(0, 1, 201)

ax1.plot(x, G_AB_alpha(x, Tx), 'k-')
ax1.plot(x, G_AB_beta(x, Tx), 'r-')

# phase diagram
T = numpy.linspace(1.0, 1.2, 21)
x_A_beta = (numpy.exp(dGt_B(T))-1)/(numpy.exp(dGt_B(T))-numpy.exp(dGt_A(T)))
x_A_alpha = x_A_beta * numpy.exp(dGt_A(T))


ax2.plot(x_A_alpha, T)
ax2.plot(x_A_beta, T)

ax2.set_xlim(0, 1)
ax2.set_xlabel("$x_A$")
ax2.set_ylabel("T/T$_0$")

plt.tight_layout()
plt.show()
