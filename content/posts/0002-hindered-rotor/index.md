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
The laters are generally treated using [(quantum) harmonic oscillators](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator).
This induces different issues, such as the lack of hanarmonicity, or the one I'm interested here, the fact that it describes a single minimum for this mode.

Indeed, we learn pretty quickly in organic chemistry lectures that given a certain amount of energy, a single bond can rotate from one minima to another.
For example, in ethane, there are 3 stable minima, or staggered conformation.

![](https://switkes.chemistry.ucsc.edu/teaching/CHEM1B/RotationalConformation/images/1.2_ethane.gif)

**Figure:** the different conformations of ethane ([source](https://switkes.chemistry.ucsc.edu/teaching/CHEM1B/RotationalConformation/Conformations_of_ethane.html)).

When the barrier is lower or comparable to the [characteristic energy](https://en.wikipedia.org/wiki/KT_(energy)) ($k_BT$), the system can (more or less freely) move from one minima to another, and using an harmonic oscillator(in which the system is confined in a single minima) is thus incorrect.

This is the case in the ethane molecule (and its derivatives) presented above, but it might also be the case for other molecules, such as the one I was interested in when I bump in this topic, [Ferrocene](https://en.wikipedia.org/wiki/Ferrocene).
Indeed, it possess an low-lying internal rotation mode located around 40 cm⁻¹.

![](vib_Fc.gif)

**Figure:** visualization of the he rotational mode of Ferrocene (located at 31.87 cm⁻¹) as computed at the ωB97X-D/6-311G* level in DMSO (SMD).

While there are different article on the subject, few of them goes into the nitty gritty of treating this correction.
In this blog post, I explore the approach proposed in [10.1007/s00214-007-0376-5](https://dx.doi.org/10.1007/s00214-007-0376-5) to correctly account for [hindered rotation](https://doi.org/10.1351/goldbook.F02520) (when the barrier is comparable to $k_BT$).

## Theory

### The 1D-hindered rotor (1-DHR) model

Like most problems in quantum physics, the one of treating internal rotations is *hard*, which means that we can only provide approximate solutions.
One of such is to assume uncoupled rotations, and to model the problem by assuming the rotation of two counter-rotating [tops](https://en.wikipedia.org/wiki/Spinning_top) around an axis (generally, a single bond), treated in practice with an effective one-dimensional Hamiltonian:

$$\tag{1}\left[-\frac{\hbar^2}{2I_r}\frac{d^2}{d\theta^2} + V(\theta)\right]\Psi(\theta)  = E\ \Psi(\theta)$$

where $I_r$ is the reduced moment of inertia of the rotating top (see below), and $V(\theta)$ is the rotational hindrance potential, *e.g.*, the potential energy surface for the rotation, given by $\theta\in[0,2\pi]$.

Note that if $V(\theta)=0$, this correspond to the case of a [particle in a ring](https://en.wikipedia.org/wiki/Particle_in_a_ring), for which the solutions are given by:

$$\tag{2}\Psi_{m}(\theta)=\frac{1}{\sqrt{2\pi}} e^{im\theta} \text{ and } E_m = \frac{m^2\hbar^2}{2I_r},$$

where $m\in\mathbb{Z}$.

### The reduced moment of inertia, $I_r$

xxx

### The potential, $V(\theta)$

This can be obtained either:

1. by performing a rigid scan around the rotation axis,
2. by performing a relaxed scan around the rotation axis, or
3. by performing as scan, followed by saddle point optimization for the maxima.

Note that in Gaussian 16, such relaxed scan can be performed thanks to [Generalized Internal Coordinates](https://gaussian.com/gic/) (GIC) with:

```txt
#P wB97XD/6-311G* opt geom=(ModRedundant,GIC)

Ethane

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

If we are able to solve (1) to get a set of level $\{\varepsilon_n\}$, we can compute the [partition function](https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)#Canonical_partition_function) as:

$$Q(T) = \sum_n g_n\ e^{-\beta\ \varepsilon_n},$$

where $\beta = (k_BT)^{-1}$, the so called [thermodynamic beta](https://en.wikipedia.org/wiki/Thermodynamic_beta). 
The sum extend over the set of **unique** energy levels, and thus $g_m$ is the degeneracy (this would be the case for the particle in a ring).
Then, we just use the expressions coming from statistical mechanics.
The internal (thermal) energy, $U$, at a given temperature is given by:

$$ U = \braket{E} = \frac{1}{Q(T)}\ \sum_n g_n\ \varepsilon_n\ e^{-\beta(T)\ \varepsilon_n}  = -\frac{\partial}{\partial\beta}\ln[Q(t)],$$

while the entropy, $S$, is given by:

$$S = k_B\ln[Q(T)] + \frac{U - U_0}{T},$$

where $U_0$ is obtained for $T=0$.

Combining those two quantities gives the [Helmholtz free energy](https://en.wikipedia.org/wiki/Helmholtz_free_energy) given by $A = U -TS = -k_BT\ln[Q(T)]$.
If one assume the pressure change to be negligible, this is equivalent to the [Gibbs free energy](https://en.wikipedia.org/wiki/Helmholtz_free_energy).

## Appplication: a simple $k$-fold potential

In the **very** simple cases, such as ethane or Ferrocene, this is a $k$-fold potential described by:

$$V(\theta)=\frac{V_0}{2}\ [1 + \cos(k\theta)],$$

where $k$ is the number of minima/maxima in the potential, also referred to as the **symmetry number**.
Note that what is following could be easily extended to a more general Fourrier series, as done, *e.g.*, in [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py).

In that case, according to [10.1063/1.475616](https://dx.doi.org/10.1063/1.475616), an approximate barrier can be estimated from the vibrational frequency of the mode, $\nu$ (in s⁻¹).
Using $\theta = \theta_0\cos(2\pi\nu t)$ and assuming a [torsional harmonic oscillator](https://en.wikipedia.org/wiki/Torsion_spring#Torsional_harmonic_oscillators) (*i.e.*, $I_r\ddot{\theta}+\kappa\theta =0$), one gets:

$$I_r\ \frac{d^2\theta}{dt^2}+ \frac{V_0k^2}{2} \theta = 0 \Leftrightarrow V_0 = \frac{2\nu^2I_r}{\hbar^2k^2}.$$

### Solving the 1-DHR

xxx

## Some literature

+ *The hindered rotor theory: A review*, [10.1002/wcms.1583](https://doi.org/10.1002/wcms.1583).
+ *Solution of the Schrödinger Equation for One-Dimensional Anharmonic Potentials: An Undergraduate Computational Experiment*, [10.1021/ed1000137](https://dx.doi.org/10.1021/ed1000137).
+ *Hindered Translator and Hindered Rotor Models for Adsorbates: Partition Functions and Entropies*, [10.1021/acs.jpcc.5b11616](https://dx.doi.org/10.1021/acs.jpcc.5b11616).
