---
title: "Hindered rotor correction in thermochemistry"
date: 2024-09-05T19:00:00+02:00
draft: true
tldr: while small, a correct treatment of hindered rotation might be required in some situations.
tags: 
- lang:en
- quantum chemistry
---

When computing thermochemical values for molecules in gas (or, with some approximation, condensed) phase, we generally split the $3N$ degree of freedom into 3 rotational modes, 3 translational modes (or 2 if the molecule is linear) and the remaining $3N-6$ [vibrational](https://en.wikipedia.org/wiki/Molecular_vibration) modes (or $3N-5$ if the molecule is linear).
The laters are generally treated using [(quantum) harmonic oscillators](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator) (QHOs).
This induces different issues, such as the lack of hanarmonicity, or the one I'm interested here, the fact that it describes a single minimum for this mode.

Indeed, we learn pretty quickly in organic chemistry lectures that given a certain amount of energy, a single bond can rotate from one minima to another.
For example, in ethane, there are 3 stable minima, or staggered conformation.

![](https://switkes.chemistry.ucsc.edu/teaching/CHEM1B/RotationalConformation/images/1.2_ethane.gif)

**Figure:** the different conformations of ethane ([source](https://switkes.chemistry.ucsc.edu/teaching/CHEM1B/RotationalConformation/Conformations_of_ethane.html)).

When the barrier is lower or comparable to the [characteristic energy](https://en.wikipedia.org/wiki/KT_(energy)) ($k_BT$), the system can (more or less freely) move from one minima to another, and describing this "vibration" with a QHO (in which the system is confined in a single minima) is thus incorrect.

This is the case for most C-C bonds as presented above, but it might also be the case for other molecules, such as the one I was interested in when I bump in this topic, [ferrocene](https://en.wikipedia.org/wiki/Ferrocene).
Indeed, it possess an low-lying internal rotation mode located around 40 cm⁻¹.

![](vib_Fc.gif)

**Figure:** visualization of the he rotational mode of Ferrocene (located at 31.87 cm⁻¹) as computed at the ωB97X-D/6-311G* level in DMSO (SMD).

While there are different article on the subject, few of them goes into the nitty gritty of treating this correction.
In this blog post, I explore the approach proposed in [10.1007/s00214-007-0376-5](https://dx.doi.org/10.1007/s00214-007-0376-5) to correctly account for [hindered rotation](https://doi.org/10.1351/goldbook.F02520) (when the barrier is comparable to $k_BT$).

Note that this has been treated by Pitzer and Gwinn in 1942, which provided reference tables (see [10.1063/1.1723744](https://dx.doi.org/10.1063/1.1723744)), used since then in many articles and programs.

## Theory

### The 1D-hindered rotor (1-DHR) model

Like most problems in quantum physics, the one of treating internal rotations is *hard*, which means that we can only provide approximate solutions.
One of such is to assume uncoupled rotations, and to model the problem by assuming the rotation of two counter-rotating [tops](https://en.wikipedia.org/wiki/Spinning_top) around an axis (generally, a single bond), treated in practice with an effective one-dimensional Hamiltonian:

$$\tag{1}\left[-\frac{\hbar^2}{2I_r}\frac{d^2}{d\theta^2} + V(\theta)\right]\Theta(\theta)  = E\ \Theta(\theta)$$

where $I_r$ is the reduced moment of inertia of the rotating tops (see below), and $V(\theta)$ is the rotational hindrance potential, *e.g.*, the potential energy surface for the rotation, given by $\theta\in[0,2\pi]$.

Note that if $V(\theta)=0$, this correspond to the case of a [particle in a ring](https://en.wikipedia.org/wiki/Particle_in_a_ring). 
Provided that we use the boundary conditions $\Theta(\theta) = \Theta(\theta + \frac{2\pi}{\sigma})$, the solutions are given by:

$$\tag{2}\Theta_{m}(\theta)=\frac{1}{\sqrt{2\pi}} e^{im\theta} \text{ and } E_m = \frac{m^2\sigma^2\hbar^2}{2I_r},$$

where $m\in\mathbb{Z}$ is a quantum number, and $\sigma\in\mathbb{N}_0$ is the rotational symmetry number, which accounts for orientation that interchanges identical atoms.
For example, in ethane, a rotation of 120° around the C-C left the molecule unchanged, so $\sigma=3$.

### The reduced moment of inertia, $I_r$

The [moment of inertia](https://en.wikipedia.org/wiki/Moment_of_inertia) is the relevant quantity when dealing with rotation, as it is a measure on how "difficult" it is to rotate a body around a given axis.
It is thus the equivalent of the mass in the linear movements (*i.e.*, [inertia](https://en.wikipedia.org/wiki/Inertia)). 

Following section 6.6 of [10.1007/978-3-030-52006-9](https://dx.doi.org/10.1007/978-3-030-52006-9), if we consider the rotational motion of two symmetric tops (one at the left, $L$, and the other at the right, $R$) around a common rotational axis (*e.g.*, the methyl groups in ethane), we can write an expression for the kinetic energy,

$$T_{rot} = \frac{1}{2}(I_L\omega_L^2+I_R\omega_R^2),$$

where the $I$'s are the moment of inertia of each tops, and $\omega$'s are their angular velocities. 
This can be rewritten in term of the rotational kinetic energy of the center of mass and relative rotational motion, as:

$$T_{rot} = \frac{1}{2}(I_{tot}\Omega^2+I_r\omega^2),$$

where $\Omega = \frac{1}{I_{tot}}(I_L\omega_L+I_R\omega_R)$ is the center of mass [angular velocity](https://en.wikipedia.org/wiki/Angular_velocity) and $\omega = \omega_R-\omega_L$ is the relative angular velocity.
Furthermore, we have the total moment of inertia, $I_{tot} = I_L + I_R$, and

$$I_r = \frac{I_LI_R}{I_L+I_R},$$

the reduced moment of inertia for the relative motion.
Since we want to treat internal rotation, it is the later quantity which is relevant. 
This formula corresponds to the $I^{(2,n)}$ approximation of [10.1063/1.473958](https://dx.doi.org/10.1063/1.473958), since the proper coupling with the molecular rotations is not accounted for.

The inertia moment for each group $G=L, R$ is computed as:

$$I_{G} = \sum_{i\in G} m_id_i^2,$$

where $m_i$ is the mass of the atom $i$ in group $G$, and $d_i$ is the distance between atom $i$ and the rotation axis (see [there](https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line)).
The $n$ in $I^{(2,n)}$ refers to the definition of the rotation axis:

+ if $n=1$, the rotation axis corresponds to the twisting bond, 
+ if $n=2$, the rotation axis is parallel to the bond, but pass through the center of mass of the rotating group in question, and
+ if $n=3$, the rotation axis pass through the center of mass of both group (and is, thus, not necessarily parallel to the bond).

More accurate $I^{(m,n)}$ with $m>2$ exists, see [10.1063/1.473958](https://dx.doi.org/10.1063/1.473958) and the seminal work**s** of Pitzer, *e.g.* [10.1063/1.1932193](https://dx.doi.org/10.1063/1.1932193), which involve Coriolis coupling (see the "Proof" section).
Another approximation comes from the fact that $I_L$ and $I_R$ (generally evaluated with the geometry of a minimum) are considered constant in the rotational motion, which may not be the case (think about, *e.g.*, rotations around C-C bond in butane).
The impact of this approximation, among others, is discussed in [10.1021/acs.jced.6b00757](https://dx.doi.org/10.1021/acs.jced.6b00757).

### The potential, $V(\theta)$

This can be obtained either:

1. by performing a rigid scan around the rotation axis,
2. by performing a relaxed scan around the rotation axis, or
3. by performing as scan, followed by saddle point optimization for the maxima.

Note that in Gaussian 16, such relaxed scan can be performed thanks to [Generalized Internal Coordinates](https://gaussian.com/gic/) (GIC) with:

```txt
#P wB97XD/6-311G* opt geom=(ModRedundant,GIC)

Scan over 120° around the C-C bond, by step of 10°.

0  1
C          -2.01593         0.66575         0.00003
C          -0.48083         0.73954        -0.00003
H          -0.07239         0.17360        -0.84025
H          -0.14934         1.77672        -0.08652
H          -0.07942         0.32314         0.92656
H          -2.41733         1.08215        -0.92659
H          -2.34739        -0.37145         0.08652
H          -2.42440         1.23167         0.84021

Dih1(StepSize=10.0,NSteps=12)=Dihedral(6,1,2,3)
```

In general, the obtained profile is fitted to a [Fourrier series](https://en.wikipedia.org/wiki/Fourier_series).

### Thermochemistry

For a given system that can exists in different microstates defined by energy levels $\\{\varepsilon_n\\}$, we can compute the [partition function](https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)#Canonical_partition_function) as:

$$Q(T) = \sum_n e^{-\beta\ \varepsilon_n},$$

where $\beta = (k_BT)^{-1}$, the so called [thermodynamic beta](https://en.wikipedia.org/wiki/Thermodynamic_beta). 
The sum extend over the set of all energy levels (a degeneracy factor can be used to consider the set of unique energy levels).
Then, we just use the expressions coming from statistical mechanics.
The internal (thermal) energy, $U$, at a given temperature is given by:

$$ U = \braket{E} = \frac{1}{Q(T)}\ \sum_n \varepsilon_n\ e^{-\varepsilon_n\beta}  = k_BT^2\ \frac{\partial}{\partial T}\ln[Q(t)] ,$$

while the entropy, $S$, is given by:

$$S = k_B\ln[Q(T)] + \frac{U - U_0}{T},$$

where $U_0$ is obtained by setting $T=0$ in the expression for $U$.

Combining those two quantities gives the [Helmholtz free energy](https://en.wikipedia.org/wiki/Helmholtz_free_energy) given by $A = U -TS = -k_BT\ln[Q(T)]$.
If one assume the pressure change to be negligible, this is equivalent to the [Gibbs free energy](https://en.wikipedia.org/wiki/Helmholtz_free_energy).

For example, let's consider a 1D **free rotor** (FR), using the solution given in Eq. (2). 
The partition function is then given by:

$$Q_{FR}(T) = \sum_{m\in\mathbb{Z}} e^{-B\sigma^2\beta m^2}, \text{ with } B = \frac{\hbar^2}{2I_r}.$$

where $B$ is the [rotational constant](https://en.wikipedia.org/wiki/Rigid_rotor#Quantum_mechanical_linear_rigid_rotor) of the rotor.
According to [10.1016/0009-2614(94)87058-6](https://dx.doi.org/10.1016/0009-2614(94)87058-6), if $B\sigma^2\beta$ is sufficiently small (large temperature and/or large $I_r$), one can approximate the partition function by computing the integral over $m$:

$$Q_{FR}(T) =  \int_{-\infty}^\infty e^{-B\sigma^2\beta m^2}\ dm = \sqrt{\frac{\pi}{B\sigma^2\beta}} = \sqrt{\frac{2\pi k_B T\ I_r}{\sigma^2\hbar^2}}.$$

Therefore,

$$\tag{3} U_{FR} = \frac{k_BT}{2} \text{  and } S_{FR} = k_B\ \left[\ln[Q(T)] + \frac{1}{2}\right].$$

For a single QHO characterized by a vibrational energy $\varepsilon_{vib} = h\nu_{vib}$ (where $\nu$ is the vibrational frequency), one has:

$$Q_{QHO}(T) = \frac{e^{-\frac{1}{2}\varepsilon_{vib}\beta}}{1-e^{-\varepsilon_{vib}\beta}},$$

so that:

$$\tag{4} U_{QHO} = \varepsilon_{vib}\ \left[\frac{1}{2}+\frac{1}{1-e^{-\varepsilon_{vib}\beta}}\right], \text{ and }
S_{QHO} = k_B\left[\frac{\varepsilon_{vib}\beta}{1-e^{-\varepsilon_{vib}\beta}} - \ln\left(1-e^{-\varepsilon_{vib}\beta}\right)\right].$$

Note that here, the [zero-point energy](https://en.wikipedia.org/wiki/Zero-point_energy) have been directly included in both $Q(T)$ and $U$ (referred to as the "bottom of the well" convention in some text and programs, such as [Gaussian](https://gaussian.com/wp-content/uploads/dl/thermo.pdf)).
This correction is sometimes added separately.

## Example: a simple $\sigma$-fold cosine potential

In the **very** simple cases, such as ethane (or ferrocene, for that matter), instead of FR or QHO, one can approximate $V(\theta)$ with a $\sigma$-fold cosine potential:

$$V(\theta)=\frac{V_0}{2}\ [1 - \cos(\sigma\theta)],$$

where $\sigma$ thus corresponds to the number of minima/maxima in the potential.
Note that what is following could be easily extended to a more general Fourrier series, as done, *e.g.*, in [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py).

In that case, according to [10.1063/1.475616](https://dx.doi.org/10.1063/1.475616), an approximate barrier can be estimated from the vibrational frequency of the mode, $\nu$ (in s⁻¹).
Using $\theta = \theta_0\cos(2\pi\nu t)$ and assuming a [torsional harmonic oscillator](https://en.wikipedia.org/wiki/Torsion_spring#Torsional_harmonic_oscillators) (*i.e.*, $I_r\ddot{\theta}+\kappa\theta =0$), one gets:

$$I_r\ \frac{d^2\theta}{dt^2}+ \frac{V_0\sigma^2}{2} \theta = 0 \Leftrightarrow V_0 = \frac{2\nu^2I_r}{\sigma^2\hbar^2}.$$

### Solving the 1-DHR

Note that here, we are only interested in the energies from Eq. (1) in order to compute the thermodynamical quantities.
While there actually are analytical solutions for this case, given by [Mathieu functions](https://en.wikipedia.org/wiki/Mathieu_function), I would like to stress another method (or, rather, "remind", as it is generally addressed in quantum physics lectures) which, while being approximate, can virtually apply to any $V(\theta)$: the **secular determinant**.
One can find its application for this case in [10.1021/ed1000137](https//dx.doi.org/10.1021/ed1000137) (which covers different cases in its Supporting Information).

Reminder: it is based on the good old **variational principle** (more precisely the [Rayleigh-Ritz method](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Ritz_method)).
Provided a set of basis functions $\\{\phi_i\\}$, we know from the variational principle that the energy of a trial function $\Psi = \sum_i c_i\phi_i$ is always larger or equal to the true ground state energy, $\varepsilon_0$,

$$\varepsilon_0 \leq \varepsilon, \text{ with } \varepsilon =  \frac{\braket{\Psi|H|\Psi}}{\braket{\Psi|\Psi}},$$

so that we can find the set of coefficients $\\{c_i\\}$ that minimize the energy of our trial function by setting $\frac{d\varepsilon}{dc_i} = 0$.
If we do that, we eventually obtain a set of secular equations of the form:

$$\forall k: \sum_i c_i(H_{ki}-\varepsilon_k S_{ki}) = 0,$$

where $S_{ki}=\braket{\phi_k|\phi_i}$ (overlap matrix) and $H_{ki} = \braket{\phi_k|\hat H|\phi_i}$ (Hamiltonian matrix).
This is a [generalized eigenvalue problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem).
However, if the set of functions is orthogonal, $S_{ki} = \delta_{k,i}$ and this becomes an (normal) eigenvalue problem, for which we "just" need to diagonalize $H$ to get the energies.

Any set of [orthogonal function](https://en.wikipedia.org/wiki/Orthogonal_functions) would work, but depending on the system, some are easier than others.
In this case, the solution for the FR, $\\{\Theta_i\\}$, given in Eq. (1), will do fine.
The elements of $H$ are given by:

$$\begin{aligned}
H_{mn} &=\braket{\Theta_m|\hat H|\Theta_n}\\\\ 
&= -B\ \int_0^{2\pi} \Theta_m^\star\frac{d^2}{d\theta^2}\Theta_n\ d\theta + \frac{V_0}{2}\int_0^{2\pi} \Theta_m^\star[1 - \cos(\sigma\theta)]\Theta_n\ d\theta  \\\\
&= \frac{Bn^2}{4\pi}\int_0^{2\pi} e^{i(n-m)\theta}\ d\theta + \frac{V_0}{4\pi}\int_0^{2\pi}  e^{i(n-m)\theta}\ [1 - \cos(\sigma\theta)]\ d\theta \\\\
&= \left(Bn^2 + \frac{V_0}{2}\right)\delta_{m,n} - \frac{V_0}{4}\delta_{m,n\pm\sigma}.
\end{aligned}$$

Here, $m, n \in\mathbb{Z}$ are not the position in the matrix, but the quantum numbers of the basis functions.
In practice, we will chose $-M \leq m,n \leq M$, so that $H$ contains $2M+1$ basis functions, and thus columns and rows.

## Some literature

+ *Statistical Thermodynamics for Pure and Applied Sciences* (book), [10.1007/978-3-030-52006-9](https://dx.doi.org/10.1007/978-3-030-52006-9)
+ *The hindered rotor theory: A review*, [10.1002/wcms.1583](https://doi.org/10.1002/wcms.1583).
+ *Solution of the Schrödinger Equation for One-Dimensional Anharmonic Potentials: An Undergraduate Computational Experiment*, [10.1021/ed1000137](https://dx.doi.org/10.1021/ed1000137).
+ *Hindered Translator and Hindered Rotor Models for Adsorbates: Partition Functions and Entropies*, [10.1021/acs.jpcc.5b11616](https://dx.doi.org/10.1021/acs.jpcc.5b11616).
