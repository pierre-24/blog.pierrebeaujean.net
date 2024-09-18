import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

Tt_A = 1.2 # T/T_0
Tt_B = 1.0 # T/T_0

S = 7 # R

def Gt_A(T: float) -> float:
    return Tt_A * S - T * S

def Gt_B(T: float) -> float:
    return Tt_B * S - T * S

def G_AB_alpha(x: float, T: float) -> float:
    return T*(x*numpy.log(x)+(1-x)*numpy.log(1-x))

def dG_AB_alpha(x: float, T: float) -> float:
    return T * numpy.log(x / (1-x))

def tg_G_AB_alpha(x: float, T: float, xint: float):
    return dG_AB_alpha(xint, T) * x + (G_AB_alpha(xint,T) - dG_AB_alpha(xint,T) * xint)

def G_AB_beta(x, T):
    return  Gt_B(T) + (Gt_A(T)-Gt_B(T)) * x + T*(x*numpy.log(x)+(1-x)*numpy.log(1-x))

def dG_AB_beta(x: float, T: float) -> float:
    return Gt_A(T)-Gt_B(T) + T * numpy.log(x / (1-x))


figure = plt.figure(figsize=(5, 10))
axes = figure.subplots(3, sharex=True)

x = numpy.linspace(0, 1, 401)

T = numpy.array([1.3, 1.12, 1.06])
x_A_beta = (numpy.exp(Gt_B(T))-1)/(numpy.exp(Gt_B(T)) - numpy.exp(Gt_A(T)))
x_A_alpha = x_A_beta * numpy.exp(Gt_A(T))

for i in range(3): 
    axes[i].plot(x, G_AB_alpha(x, T[i]), 'r-')
    axes[i].plot(x, G_AB_beta(x, T[i]), 'b-')
    axes[i].set(yticklabels=[])
    axes[i].tick_params(left=False)  # remove the ticks
    axes[i].text(.5, .9, 'T$_{}$'.format(i+1), transform=axes[i].transAxes, fontsize=16)
    
    axes[i].text(0, 0, '$\mu_B^{0,\\alpha}$', fontsize=16, ha='right', color='red')
    axes[i].text(1, 0, '$\mu_A^{0,\\alpha}$', fontsize=16, ha='left', color='red')
    
    axes[i].text(0, Gt_B(T[i]), '$\mu_B^{0,\\beta}$', fontsize=16, ha='right', color='blue')
    axes[i].text(1, Gt_A(T[i]), '$\mu_A^{0,\\beta}$', fontsize=16, ha='left', color='blue')
    
    if i == 0:
        axes[i].text(.2, G_AB_alpha(.2, T[i]), ' $G_{AB}^{\\alpha}$', fontsize=16, color='r')
        axes[i].text(.8, G_AB_beta(.8, T[i]), '   $G_{AB}^{\\beta}$', fontsize=16, color='b')
    
    if i > 0:  # tangent
        axes[i].plot(x, tg_G_AB_alpha(x, T[i], x_A_alpha[i]), 'k--', linewidth=.75)
        axes[i].plot([x_A_alpha[i], x_A_beta[i]], [G_AB_alpha(x_A_alpha[i], T[i]), G_AB_beta(x_A_beta[i], T[i])], 'ko')
        axes[i].text(x_A_alpha[i],  G_AB_alpha(x_A_alpha[i], T[i]), "\n$Q_{}$".format(i), fontsize=16, va='center')
        axes[i].text(x_A_beta[i],  G_AB_beta(x_A_beta[i], T[i]), "\n$P_{}$".format(i), fontsize=16, va='center')
    
        axes[i].text(0, tg_G_AB_alpha(0, T[i], x_A_alpha[i]), '$\mu_B^{\\beta}=\mu_B^{\\alpha}$', fontsize=12, ha='right')
        axes[i].text(1, tg_G_AB_alpha(1, T[i], x_A_alpha[i]), '$\mu_A^{\\beta}=\mu_A^{\\alpha}$', fontsize=12, ha='left')

axes[2].set_xlim(0, 1)
axes[2].set_xlabel('molar fraction, $x_A$')

plt.tight_layout()

figure.savefig('twophase_T.svg')

# ---

figure = plt.figure(figsize=(5, 7))
ax1, ax2 = figure.subplots(2, sharex=True)

Tx=T[2]
ax1.plot(x, G_AB_alpha(x, Tx), 'r-')
ax1.plot(x, G_AB_beta(x, Tx), 'b-')
ax1.plot(x, tg_G_AB_alpha(x, Tx, x_A_alpha[2]), 'k--', linewidth=.75)
ax1.text(.2, G_AB_alpha(.2, Tx), ' $G_{AB}^{\\alpha}$', fontsize=16, color='r')
ax1.text(.8, G_AB_beta(.8, Tx), '$G_{AB}^{\\beta}$  ', fontsize=16, color='b', ha='right')

xy0=x_A_alpha[2], G_AB_alpha(x_A_alpha[2], Tx)
xy1=x_A_beta[2], G_AB_beta(x_A_beta[2], Tx)
ax1.plot(*xy0, 'ko')
ax1.plot(*xy1, 'ko')
ax1.text(*xy0, '$Q_2$', fontsize=16, va='bottom')
ax1.text(*xy1, '$P_2$', fontsize=16, va='bottom')

ax1.set(yticklabels=[])
ax1.tick_params(left=False)  # remove the ticks
ax1.set_ylabel('Free energy')

ax2.set_xlim(0, 1)
ax2.set_xlabel('molar fraction, $x_A$')
ax2.set_ylabel('Potential, $E^0$')

ax2.set(yticklabels=[])
ax2.tick_params(left=False)  # remove the ticks

x0 = numpy.linspace(0, xy1[0])
ax2.plot(x0, -dG_AB_beta(x0, Tx), 'k-')

x1 = numpy.linspace(xy1[0], xy0[0])
ax2.plot(x1, [-dG_AB_alpha(xy0[0], Tx)]*x1.shape[0], 'k-')

ax2.plot(xy0[0], -dG_AB_alpha(xy0[0], Tx), 'ko')
ax2.plot(xy1[0], -dG_AB_alpha(xy0[0], Tx), 'ko')

ax2.add_artist(ConnectionPatch(xyA=xy0, coordsA=ax1.transData, xyB=(xy0[0], -dG_AB_alpha(xy0[0], Tx)), coordsB=ax2.transData))
ax2.add_artist(ConnectionPatch(xyA=xy1, coordsA=ax1.transData, xyB=(xy1[0], -dG_AB_alpha(xy0[0], Tx)), coordsB=ax2.transData))

x2 = numpy.linspace(xy0[0], 1)
ax2.plot(x2, -dG_AB_alpha(x2, Tx), 'k-')

plt.tight_layout()
figure.savefig('twophase_E0.svg')

#------
figure = plt.figure(figsize=(5, 7))
ax1, ax2 = figure.subplots(2, sharex=True)

ax1.plot(x, G_AB_alpha(x, Tx), 'r-')
ax1.plot(x, G_AB_beta(x, Tx), 'b-')
ax1.plot(x, tg_G_AB_alpha(x, Tx, x_A_alpha[2]), 'k--', linewidth=.75)
ax1.text(.2, G_AB_alpha(.2, Tx), ' $G_{AB}^{\\alpha}$', fontsize=16, color='r')
ax1.text(.8, G_AB_beta(.8, Tx), '$G_{AB}^{\\beta}$  ', fontsize=16, color='b', ha='right')

ax1.plot(*xy0, 'ko')
ax1.plot(*xy1, 'ko')
ax1.text(*xy0, '$Q_2$', fontsize=16, va='bottom')
ax1.text(*xy1, '$P_2$', fontsize=16, va='bottom')

ax1.set(yticklabels=[])
ax1.tick_params(left=False)  # remove the ticks
ax1.set_ylabel('Free energy')

T = numpy.linspace(Tt_A, Tt_B, int((Tt_A - Tt_B) * 100) + 2)
x_A_beta = (numpy.exp(Gt_B(T))-1)/(numpy.exp(Gt_B(T)) - numpy.exp(Gt_A(T)))
x_A_alpha = x_A_beta * numpy.exp(Gt_A(T))

ax2.plot(x_A_alpha, T, 'r-')
ax2.plot(x_A_beta, T, 'b-')
ax2.plot(x, [Tx] * x.shape[0], 'k--', linewidth=.75)

ax2.add_artist(ConnectionPatch(xyA=xy0, coordsA=ax1.transData, xyB=(xy0[0], Tx), coordsB=ax2.transData))
ax2.add_artist(ConnectionPatch(xyA=xy1, coordsA=ax1.transData, xyB=(xy1[0], Tx), coordsB=ax2.transData))
ax2.plot(xy0[0], Tx, 'ko')
ax2.plot(xy1[0], Tx, 'ko')
ax2.text(xy0[0], Tx, '$Q_2$', fontsize=16, va='top')
ax2.text(xy1[0], Tx, '$P_2$', fontsize=16, va='top')

ax2.set_xlim(0, 1)
ax2.set_xlabel('molar fraction, $x_A$')
ax2.set_ylabel('Temperature')

ax2.set(yticklabels=[])
ax2.tick_params(left=False)  # remove the ticks

ax2.text(0,Tt_B, '$T_B^{\\alpha\\rightarrow\\beta}$', ha='right')
ax2.text(1,Tt_A, '$T_A^{\\alpha\\rightarrow\\beta}$', ha='left')

ax2.text(.1, .9, '$\\beta$', transform=ax2.transAxes, fontsize=16, color='blue')
ax2.text(.5, .5, '$\\alpha+\\beta$', transform=ax2.transAxes, fontsize=16)
ax2.text(.9, .1, '$\\alpha$', transform=ax2.transAxes, fontsize=16, color='red')


plt.tight_layout()
figure.savefig('twophase_Tx.svg')
