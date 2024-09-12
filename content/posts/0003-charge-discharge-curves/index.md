---
title: "Understanding the charge/discharge curve of batteries"
date: 2024-09-12T09:03:20+02:00
draft: true
tldr: These curves give information about the cathode/anode materials.
tags: 
- lang:en
---

## Introduction

Following [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5), in first approximation, a battery can actually be modeled as a simple reaction between two metallic materials, $A$ and $B$, reacting as $A + B \rightarrow AB$.

However, in a battery, $A$ and $B$ are electrode materials separated by an **electrolyte**. 
This electrolyte is specifically chosen to allow the transport of ionic species while acting as an electrical insulator. 

Therefore, for $A$ and $B$ to react, one of the species must either gain or lose electrons to be transported. 
These electrons, however, cannot travel directly through the electrolyte.
Instead, they move through an external electrical circuit to reach the other electrode.
This movement of electrons through the circuit is what powers devices and represents the working principle of a battery.
Note that ionic species and electrons must move at the same rate to maintain charge balance.

Let's assume that $A$ is the species being transported, and thus $A$ is the **anode**. 
At the interface of electrode $A$, the following reaction occurs: $A \rightarrow A^+ + e^-$. 
The ion $A^+$ then moves through the electrolyte toward electrode $B$. 
At electrode $B$, the electron meet $A^+$ and a second reaction takes place: $A^+ + e^- \rightarrow A$, followed by $A + B \rightarrow AB$.
Schematically, one has:

```goat
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+----------+
|  A  | electrolyte |  B  | -->  | A | electrolyte | AB | B | -->  | electrolyte |    AB    |
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+----------+
```

Note that the diffusion of $A$, $B$, and $e^-$ through $AB$ is necessary for the reaction to occur. 
Therefore, $AB$ must function as both an electrical and ionic conductor.

When the circuit is open, no reaction takes place, meaning the ionic species stop moving between the electrodes. 
This indicates that the chemical driving force acting on the ions is counterbalanced by another force, specifically, an electrostatic force.
The chemical driving force arises from the difference in chemical potential between the two electrodes, which corresponds to the $\Delta G_r$ of the reaction $A + B \rightarrow AB$ (the fact that $A$ gets transported from one side to the other is actually irrelevant here). 
The electrostatic force derived from the energy of a set of charge, by $-z\mathcal{F} E^0$, where $z$ is the charge number of the ionic species, $\mathcal{F}$ is the [Faraday constant](https://en.wikipedia.org/wiki/Faraday_constant), and $E^0$ is the (absolute) potential difference between the two electrodes.
The energy balance is expressed as:

$$\tag{1} \Delta G_r = -z\mathcal{F} E^0.$$

This is (one form of) the [Nernst equation](https://en.wikipedia.org/wiki/Nernst_equation).

Following $E^0$ alongside the charge/discharge process gives some very interesting information about the nature of the electrodes.
In this post, I will explore the thermodynamic reasons for that.

## A bit of thermodynamics

### The chemical potential

The [chemical potential](https://en.wikipedia.org/wiki/Chemical_potential) of a given species, $\mu_i$, represents the energy absorbed or released due to a change in the number of particles of that species:

$$\mu_i = \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n\_{j\neq i}}.$$

Particles tend to move from regions of high chemical potential to regions of low chemical potential, as this reduces the system's (free) energy. 

Consider a system with two compartments, 1 and 2, separated by a membrane, containing species X at different concentrations. 
The movement of particles from compartment 1 to compartment 2 occurs such that $-dn_1 = dn_2$. 
The corresponding change in Gibbs free energy is:

$$dG = d(G_1+ G_2) = \mu_1dn_1 + \mu_2dn_2 = -(\mu_1-\mu_2)dn_2,$$

where $\mu_i$ is the chemical potential in compartment $i$.
If $\mu_1 > \mu_2$, then $dG < 0$, indicating that the transfer from compartment 1 (high chemical potential) to compartment 2 (low chemical potential) decreases the total Gibbs free energy. 
Additionally, as particles move, the chemical potentials of the compartments adjust: $\mu_1$ decreases in compartment 1, and $\mu_2$ increases in compartment 2. At equilibrium, $\mu_1 = \mu_2$.

This thought experiment leads to a fundamental principle: **when two (or more) phases are at equilibrium, the chemical potential of any components is equal in all phases**.

Furthermore, in a solution, the chemical potential of $i$ only depends on temperature and pressure. 
This leads to:

$$\mu_i = \mu_i^0+RT\ln a_i,$$

where $a_i = \gamma_i x_i$ is the activity (an effective concentration), with $\gamma_i$ the [activity coefficient](https://en.wikipedia.org/wiki/Activity_coefficient), $x_i =n_i / \sum_j n_j$ the molar fraction of $i$ in solution and $\mu_i^0$ is the standard state chemical potential.
Therefore, $a_i\in[0,1]$ represent how "pure" compound $i$ is in solution.

From the definition of the chemical potential, it follows that **when two (or more) phases are at equilibrium, the activity of any components is the same in all phases**, provided that the activity is expressed with respect to the same standard state.

Finally, one can re-write the Nernst equation as:

$$\tag{2}  E^0 = -\frac{\mu_{cathode}-\mu_{anode}}{z\mathcal F}.$$

## Sources

+ *Advanced Batteries (Materials Science Aspects)* (book), [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5).
