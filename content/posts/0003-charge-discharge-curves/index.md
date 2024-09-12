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

Let's assume that $A$ is the species being transported. 
At the interface of electrode $A$, the following reaction occurs: $A \rightarrow A^+ + e^-$. 
The ion $A^+$ then moves through the electrolyte toward electrode $B$. 
At electrode $B$, the electron meet $A^+$ and a second reaction takes place: $A^+ + e^- \rightarrow A$, followed by $A + B \rightarrow AB$.
Schematically, one has:

```goat
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+--------+
|  A  | electrolyte |  B  | -->  | A | electrolyte | AB | B | -->  | electrolyte |   AB   |
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+--------+
```

Note that the diffusion of $A$, $B$, and $e^-$ through $AB$ is necessary for the reaction to occur. 
Therefore, $AB$ must function as both an electrical and ionic conductor.

When the circuit is open, no reaction takes place, meaning the ionic species stop moving between the electrodes. 
This indicates that the chemical driving force acting on the ions is counterbalanced by another force, specifically, an electrostatic force.
The chemical driving force arises from the difference in chemical potential between the two electrodes, which corresponds to the $\Delta G_r$ of the reaction $A + B \rightarrow AB$ (the fact that $A$ gets transported from one side to the other is actually irrelevant here). 
The electrostatic force derived from the energy of a set of charge, by $-z\mathcal{F} E^0$, where $z$ is the charge number of the ionic species, $\mathcal{F}$ is the [Faraday constant](https://en.wikipedia.org/wiki/Faraday_constant), and $E^0$ is the (absolute) potential difference between the two electrodes.
The energy balance is expressed as:

$$\Delta G_r = -z\mathcal{F} E^0.$$

This is (one form of) the [Nernst equation](https://en.wikipedia.org/wiki/Nernst_equation).

## Sources

+ *Advanced Batteries (Materials Science Aspects)* (book), [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5).
