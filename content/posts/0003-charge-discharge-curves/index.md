---
title: "Understanding the charge/discharge curve of batteries"
date: 2024-09-18T11:40:00+02:00
draft: false
tldr: These curves give information about the cathode/anode materials in the charge/discharge process.
tags: 
- lang:en
- quantum chemistry
- batteries
---

The charge/discharge curves are invaluable tools for investigating batteries and are part of the "standard set" of measurements when comparing a new, promising battery to existing ones. 
However, as a theoretical chemist not deeply immersed in electrochemistry, I initially found this concept a bit confusing. 
So, I decided to dig a little deeper.


## Introduction

Following [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5), in first approximation, a battery can actually be modeled as a simple reaction between two metallic materials, $A$ and $B$, reacting for example as $A + B \rightarrow AB$.

However, in a battery, $A$ and $B$ are electrode materials separated by an **electrolyte**. 
This electrolyte is specifically chosen to allow the transport of ionic species while acting as an electrical insulator. 

Therefore, for $A$ and $B$ to react, one of the species must either gain or lose electrons to be transported. 
These electrons, however, cannot travel directly through the electrolyte.
Instead, they move through an external electrical circuit to reach the other electrode.
This movement of electrons through the circuit is what powers devices and represents the working principle of a battery.
Note that ionic species and electrons must move at the same rate to maintain charge balance.

Let's assume that $A$ is the species being transported. The pathway is the following:

1. At the interface of electrode $A$, say that the following reaction occurs: $A \rightarrow A^+ + e^-$. $A$ is thus the **anode**.
2. The ion $A^+$ then moves through the electrolyte toward electrode $B$. 
3. At electrode $B$, the electron meet $A^+$ and a second reaction takes place: $A^+ + e^- \rightarrow A$, followed by $A + B \rightarrow AB$. $B$ is thus the **cathode**.

Schematically, one has:

```goat
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+----------+
|  A  | electrolyte |  B  | -->  | A | electrolyte | AB | B | -->  | electrolyte |    AB    |
+-----+-------------+-----+      +---+-------------+----+---+      +-------------+----------+
```

During the charge/discharge process, at the cathode side, there are two phases: one containing $AB$ and another with $B$.

Note that the diffusion of $A$, $B$, and $e^-$ through $AB$ is necessary for the reaction to occur. 
Therefore, $AB$ must function as both an electrical and ionic conductor.

When the circuit is open, no reaction takes place, meaning the ionic species stop moving between the electrodes. 
This indicates that the chemical driving force acting on the ions is counterbalanced by another force, specifically, an electrostatic force.
The chemical driving force arises from the difference in chemical potential between the two electrodes, which corresponds to the $\Delta G_r$ of the reaction $A + B \rightarrow AB$ (the fact that $A$ gets transported from one side to the other is actually irrelevant here). 
The electrostatic force derived from the energy of a set of charge, by $-z\mathcal{F} E^0$, where $z$ is the charge number of the ionic species, $\mathcal{F}$ is the [Faraday constant](https://en.wikipedia.org/wiki/Faraday_constant), and $E^0$ is the (absolute) potential difference between the two electrodes.
The energy balance is expressed as:

$$\tag{1} \Delta G_r = -z\mathcal{F} E^0.$$

This is (one form of) the [Nernst equation](https://en.wikipedia.org/wiki/Nernst_equation).

Note that the formation reaction, $A + B \rightarrow AB$, where a new phase $AB$ forms alongside $B$, is not the only mechanism in battery operation.
Other mechanisms include:

+ **Displacement Reactions**: $A + BX \rightarrow AX + B$. 
  This involves the creation of a new phase containing $B$, for which the driving force is that the phase $AX$ is more stable than $BX$. 
  An example of this mechanism can be found in [CuS batteries](https://dx.doi.org/10.1002/aenm.202002394).
+ **Insertion Reactions**: $xA + BX \rightarrow A_xBX$. 
  In this case, $A$ is inserted into unoccupied sites within the structure of $BX$, forming a **solid solution** phase that exhibits a continuous range of compositions, *i.e.*, a range of $x$ values. 
  This is the working principle of many modern batteries, such as [lithium-ion batteries](https://en.wikipedia.org/wiki/Lithium-ion_battery) (originally utilizing cobalt oxide as the cathode, thought other oxides have since been explored).

```goat
+-----+-------------+-----+        +---+-------------+--+--+---+
|  A  | electrolyte | BX  | -+-->  | A | electrolyte |B |AX|BX |  (diplacement)
+-----+-------------+-----+  |     +---+-------------+--+--+---+
                             |
                             |     +---+-------------+-------+
                             +-->  | A | electrolyte | AₓBX  |    (insertion)
                                   +---+-------------+-------+
```

The potential along the charge/discharge process depends on the mechanism in place, and thus following $E^0$ alongside the charge/discharge process gives some information.
In this post, I will explore the (thermodynamic) reasons for that.

## A bit of thermodynamics

### The chemical potential

The [chemical potential](https://en.wikipedia.org/wiki/Chemical_potential) of a given species, $\mu_i$, represents the energy absorbed or released due to a change in the number of particles of that species:

$$\mu_i = \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n\_{j\neq i}}.$$

Particles tend to move from regions of high chemical potential to regions of low chemical potential, as this reduces the system's (free) energy. 

Consider a system with two compartments, 1 and 2, separated by a membrane, containing species X at different concentrations. 
The movement of particles from compartment 1 to compartment 2 occurs such that $-dn_1 = dn_2$. 
The corresponding change in Gibbs free energy is:

$$dG = d(G_1+ G_2) = \mu_1dn_1 + \mu_2dn_2 = -(\mu_1-\mu_2)dn_2,$$

where $\mu_i$ is the chemical potential in compartment $i$ and we have used the relationship:

$$dG = \sum_i \mu_idn_i.$$

If $\mu_1 > \mu_2$, then $dG < 0$, indicating that the transfer from compartment 1 (high chemical potential) to compartment 2 (low chemical potential) decreases the total Gibbs free energy. 
Additionally, as particles move, the chemical potentials of the compartments adjust: $\mu_1$ decreases in compartment 1, and $\mu_2$ increases in compartment 2. At equilibrium, $\mu_1 = \mu_2$.

This thought experiment leads to a fundamental principle: **when two (or more) phases are at equilibrium, the chemical potential of any components is equal in all phases**.

Finally, in a solution, the chemical potential of $i$ only depends on temperature and pressure. 
This leads to:

$$\mu_i = \mu_i^0+RT\ln a_i,$$

where $a_i = \gamma_i x_i$ is the activity (an effective concentration), with $\gamma_i$ the [activity coefficient](https://en.wikipedia.org/wiki/Activity_coefficient), $x_i =n_i / \sum_j n_j$ the molar fraction of $i$ in solution and $\mu_i^0$ is the standard state chemical potential.
Therefore, $a_i\in[0,1]$ represent how "pure" compound $i$ is in solution.

From the definition of the chemical potential, it follows that **when two (or more) phases are at equilibrium, the activity of any components is the same in all phases**, provided that the activity is expressed with respect to the same standard state.

### The electrochemical potential

As noted in [10.1021/acsenergylett.0c02443](https://dx.doi.org/10.1021/acsenergylett.0c02443), the chemical potential does not account for the electrostatic contribution. 
Thus, the **electrochemical potential**, $\bar\mu_i$, is defined as:

$$\bar\mu_i = \mu_i + z_i \mathcal{F} \phi,$$

where $z_i$ is the charge of species $i$, and $\phi$ is the electrical potential in a given phase.

We can extend the principle of chemical equilibrium by stating that **when two (or more) phases are at equilibrium, the electrochemical potential of each component is equal across all phases**.
This principle is particularly useful in electrochemistry, where the electrochemical potential of electrons is often the focus. 
For example, the (relative) electrode or cell potential is defined by the difference between the electrochemical potential of the working electrode, $\bar\mu_e^{w}$, and that of a reference electrode, $\bar\mu_e^{ref}$:

$$E^0_{w} = -\frac{\bar\mu_e^{w} - \bar\mu_e^{ref}}{\mathcal{F}}.$$

Moreover, if the working electrode is in equilibrium with a redox couple $O/R$, where the reaction is $O + ne^- \rightarrow R$, then the electrochemical potentials of all species involved must be equal:

$$\bar\mu_e^{w} = \bar\mu_{O/R} \Rightarrow \bar\mu_e^s = \frac{\bar\mu_R - \bar\mu_O}{n},$$

where $\bar\mu_e^s$ is the electrochemical potential of the solution. 
This gives rise to the *solution* potential:

$$E^0 = -\frac{\bar\mu_e^s}{\mathcal{F}} = -\frac{\bar\mu_R - \bar\mu_O}{n\mathcal{F}},$$

which is another form of the Nernst equation (Eq. 1). 
In practice, since $\bar\mu_e^{w}$ is measured relative to a reference, $E^0$ is also defined with respect to a reference electrode.

With those tools in hand, let's address our different charge/discharge mechanisms.

### Thermodynamics of mixing in a single phase

Let's first consider the case of a single phase containing two components, $A$ and $B$. 
This could represent either a liquid phase or a [solid solution](https://en.wikipedia.org/wiki/Solid_solution) (often referred to as an [alloy](https://en.wikipedia.org/wiki/Alloy) in metallurgy).

The molar fraction of $A$ is denoted as $x_A$, while that of $B$ is $x_B = 1 - x_A$. 
The free energy of this system in a given phase, $\alpha$, is expressed as:

$$\tag{2} G^\alpha_{AB} = x_A\mu_A^\alpha + x_B\mu_B^\alpha.$$

Note that using the electrochemical potential, $\bar\mu_i$, instead of the chemical potential, $\mu_i$, would not change this result.

The free energy of mixing is the difference between the free energy of the individual components when separated and when mixed together as phase $\alpha$. This can be written as:

$$\Delta G_{mix} = G^\alpha_{AB} - G_A - G_B = x_A \Delta\mu_A^\alpha + x_B \Delta\mu_B^\alpha = RT(x_A \ln a^\alpha_A + x_B \ln a^\alpha_B),$$

where we have used the fact that $G_A = x_A \mu^0_A$, since the activity of pure $A$ is equal to 1.

![](mix.png)

**Figure:** Evolution of the free energy of phase $\alpha$ with $x_A$ (inspired by [10.1016/B978-0-444-53770-6.00003-4](https://doi.org/10.1016/B978-0-444-53770-6.00003-4)). 
The tangent (dashed line) illustrates that at a given $x$ (here $x=0.4$), $G^\alpha_{AB} = x_A\mu_A^\alpha + x_B\mu_B^\alpha$, and the intercept at $x=0$ and $x=1$ provides $\mu_A^\alpha$ and $\mu_B^\alpha$, respectively. 

In an ideal solution (or Raoultian solution, where [Raoult's law](https://en.wikipedia.org/wiki/Raoult's_law) holds), we assume that $a_i = x_i$, so:

$$\Delta G_{mix}^{ideal} = RT[x_A \ln x_A + (1 - x_A) \ln (1 - x_A)]$$

In this case, the mixing process is purely entropic, since $\Delta H_{mix}^{ideal} = \Delta G_{mix}^{ideal} + T \Delta S_{mix}^{ideal} = 0$, with:

$$\Delta S_{mix} = -\left(\frac{dG_{mix}}{dT}\right)_P = -R[x_A \ln a^\alpha_A + x_B \ln a^\alpha_B].$$

This ideal solution model is an approximation, and deviations from ideal behavior are referred to as *excess properties*, such as the excess free energy of mixing, $\Delta G^E_{mix} = \Delta G_{mix} - \Delta G_{mix}^{ideal}$.
If the curve above is known, it would be possible to evaluate the activity by comparing to the ideal case.

Interestingly, this model is commonly used to describe insertion processes. 
As discussed in [10.1021/ar200329r](https://pubs.acs.org/doi/10.1021/ar200329r) (and others), the free energy of material $A_xBX$ with respect to the "concentration" of inserted species can be approximated as:

$$G_{A_xBX}(x) = G_{BX} + \varepsilon x + RT[x \ln x + (1 - x) \ln (1 - x)],$$

where $\varepsilon$ represents the free energy per added $A$ in the structure, and $x \in [0, 1]$ is the molar fraction of $A$ in $A_xBX$ (and thus $1-x$ is the fraction of holes). 
This range can be more restricted if the phase changes during the charging process, as described below. 
Here, we recognize the entropic term of an ideal solution. 
The corresponding elechtrochemical potential is:

$$\bar\mu_{A_xBX} = \frac{dG_{A_xBX}}{dx} + z\mathcal F\phi  =  \varepsilon + RT\ln\frac{x}{1-x} + z\mathcal F\phi.$$

This result can also be derived from statistical mechanics, using the definition of entropy $S = k_B \ln \Omega$ (see, *e.g.*, [this document](https://dspace.mit.edu/bitstream/handle/1721.1/100188/10-626-spring-2011/contents/lecture-notes/MIT10_626S11_lec07.pdf)).

From this expression, and using the Nernst equation, we can derive the variation of the potential:

$$E^0 = - \frac{1}{\mathcal F}\left[\varepsilon + RT \ln \frac{x}{1 - x} - \mu^0_{A}\right].$$

![](insertion.png)

**Figure:** Evolution of the free energy of $A_xBX$ (top) and of the corresponding potential (bottom) with $x$ (inspired by [10.1021/ar200329r](https://pubs.acs.org/doi/10.1021/ar200329r)).
The potential is proportional to $dG(x)/dx$, as illustrated here for $x=0.55$.

The curve on the bottom is therefore representative of a charging curve for a single phase.

### Thermodynamics of two phases

Now, let's consider a system with two components, $A$ and $B$, but with two possible phases, $\alpha$ and $\beta$. 
For each phase, the free energy equation (as derived previously) is valid. 
By plotting the free energy as a function of composition at different temperatures, we can observe interesting trends in phase stability:

![](twophase_T.svg)

**Figure:** Evolution of $G^\alpha_{AB}$ (red curve) and $G^\beta_{AB}$ (blue curve) with temperature, assuming $T_1 > T_2 > T_3$ (adapted from [10.1016/B978-0-444-53770-6.00003-4](https://doi.org/10.1016/B978-0-444-53770-6.00003-4)). 
The common tangent for $T_2$ and $T_3$ indicates points where the chemical potentials of a given component in $\alpha$ and $\beta$ are equal (discussed below).

**Observations:**

1. **At $T_1$:** The free energy of the $\beta$ phase, $G^\beta$, is lower than that of $\alpha$, $G^\alpha$, at all values of $x_A$. 
   Therefore, the $\beta$ phase is more stable across all compositions, and the system behaves as described for a single phase.
2. **At $T_2 < T_1$:** The curves for $G^\alpha$ and $G^\beta$ intersect. 
   This indicates that between $x_A = 0$ and $x_A = P_1$, the $\beta$ phase is more stable, while between $x_A = Q_1$ and $x_A = 1$, the $\alpha$ phase is more stable. 
   In the region between $P_1$ and $Q_1$, a two-phase equilibrium exists, and the true minimum free energy lies along the common tangent connecting $P_1$ and $Q_1$, since this line is more stable than either phase alone.
3. **At $T_3 < T_2$:** The two-phase equilibrium extends over a different range, from $P_2$ to $Q_2$.

At equilibrium, the (electro)chemical potential of a given component must be equal across the two phases. This condition is expressed as:

$$\mu_A^\alpha = \mu_A^\beta \quad \text{and} \quad \mu_B^\alpha = \mu_B^\beta.$$

This is illustrated by the common tangent between $P_1$ and $Q_1$ (and similarly between $P_2$ and $Q_2$). 
At $x_A = P_1$, the (electro)chemical potential of $A$ in phase $\alpha$ is equal to that in phase $\beta$ at $x_A = Q_1$. 
At a given temperature, for the two-phase region:

$$\mu_A^\alpha = \mu_A^\beta \Leftrightarrow \mu_A^{0,\alpha} + RT\ln a_A^{\alpha} = \mu_A^{0,\beta} + RT\ln a_A^{\beta},$$

which simplifies to:

$$RT \ln \frac{a_A^\alpha}{a_A^\beta} = \Delta G^{0,\alpha \rightarrow \beta}_A,$$

where $\Delta G^{0,\alpha \rightarrow \beta}_A = \mu_A^{0,\beta} - \mu_A^{0,\alpha}$.
This is the equation that provides the $P_1Q_1$ and $P_2Q_2$ tangent in the diagram above.

If we assume an ideal solution (i.e., $a_i = x_i$), then:

$$\frac{x_A^\alpha}{x_A^\beta} = \exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_A}{RT}\right), \quad \frac{1 - x_A^\alpha}{1 - x_A^\beta} = \exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_B}{RT}\right).$$

Here, $x_A^\alpha$ is the molar fraction of $A$ in phase $\alpha$, and $x_A^\beta$ is the molar fraction of $A$ in phase $\beta$. 
Solving this system gives both $x_A^\alpha$ and $x_A^\beta$ at any temperature:

$$
\begin{aligned}
x_A^\alpha &= x_A^\beta \exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_A}{RT}\right), \\\\
x_A^\beta &= \frac{\exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_B}{RT}\right) - 1}{\exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_B}{RT}\right) - \exp\left(\frac{\Delta G^{0,\alpha \rightarrow \beta}_A}{RT}\right)}.
\end{aligned}
$$

It's important to note that $\Delta G^{0,\alpha \rightarrow \beta}_i = \Delta H^{0,\alpha \rightarrow \beta}_i - T \Delta S^{0,\alpha \rightarrow \beta}_i$, meaning that $\Delta G^{0,\alpha \rightarrow \beta}_i$ depends on temperature, while $\Delta H^{0,\alpha \rightarrow \beta}_i$ and $\Delta S^{0,\alpha \rightarrow \beta}_i$ can be approximated as temperature-independent.
Furthermore, $\Delta G^{0,\alpha \rightarrow \beta}_i = 0$ at the temperature of phase transition, $T_i^{\alpha\rightarrow\beta}$, so in first approximation, $\Delta H^{0,\alpha \rightarrow \beta}_i = T_i^{\alpha\rightarrow\beta}\Delta S^{0,\alpha \rightarrow \beta}_i$.

By computing $x_A^\alpha$ and $x_A^\beta$ for a set of temperature, one can mark the limit between the one-phase and the two-phases region.
Plotted in the form of composition versus temperature, this provides a binary phase diagram:

![](twophase_Tx.svg)

**Figure:** Evolution of $G^\alpha_{AB}$ (red curve) and $G^\beta_{AB}$ (blue curve) at a given temperature (top) and corresponding binary phase diagram (bottom).
The dashed line in the phase diagram corresponds to the temperature at which the curves on the top are computed, and the limit of the two-phase region are given by the common tangent.

This is a prototypical phase diagram: at $T>T_A^{\alpha\rightarrow\beta}$, $\beta$ is the only phase for all $x$, while at $T<T_B^{\alpha\rightarrow\beta}$, $\alpha$ is the only phase.
Between these, the number of phase is delimited by the blue and the red line: at $x_A$ before the blue line, $\alpha$ is the only phase, but between the blue and red line, we are in the two-phase region.
In this region, the overall composition is fixed by the intersection between a *tie-line* (horizontal line, constant $T$) and the limit of the two-phase region, $P_2$ and $Q_2$.
In the example given above, the two-phase region will contain, at equlibrium, of a mixture of $\beta$ of composition $x_A=P_2$ and $\alpha$ with a composition $x_A=Q_2$.
The relative concentration of $A$ and $B$ in both phases is given by the [lever rule](https://en.wikipedia.org/wiki/Lever_rule).

Before going back to phase diagrams (which are, in the end, what is relevant), let's have a look at the evolution of the potential with $x_A$.
Taking the derivative, one gets:

![](twophase_E0.svg)

**Figure:** Evolution of the free energy (top) and of the corresponding potential (bottom) with $x_A$.
The two-phase part is characterized by a plateau in the evolution of the potential.

Thus, contrary to the S-shape obtained above for a single phase, a two-phase region is characterized by a **plateau**.
Note that for $\alpha$ and $\beta$ regions, the S-shape curves are still there.

## Phase diagrams and charge/discharge curves

The thermodynamic behavior discussed earlier is encapsulated in the *Gibbs phase rule*, which governs any system in equilibrium. 
Due to the equilibrium of the chemical potential of a compound across all phases, at each point in a phase diagram, the following relation holds:

$$F = C - P + 2,$$

where:

- $C$ is the number of components, i.e., the number of substances required to fully describe the composition of each phase (so far, we have considered $C = 2$),
- $P$ is the number of phases, and
- $F$ is the number of degrees of freedom, or the *variance* of the system.

Consider the temperature-composition (X-T) phase diagram shown below:

![](twophase_phase_rule.svg)

This diagram is valid at a given pressure, which is assumed to be constant. 
In the single-phase regions ($\alpha$ or $\beta$ phases), $P = 1$, so according to the phase rule, $F = 3$. 
This means that pressure, temperature, and composition can all vary independently, giving us three variables.

However, in the two-phase region ($\alpha + \beta$), $P = 2$, so $F = 2$. 
Here, if pressure and temperature are fixed, the composition is constrained by the tie-line, meaning the composition is no longer an independent variable. 
The relative proportions of each phase can be determined using the *lever rule*.

Now, recall that composition and (electro)chemical potential are closely linked. 
In regions where $F < 3$ (and pressure and temperature are fixed), the electrochemical potential becomes constant across all compositions, leading to the appearance of a **plateau** in the potential.
If $F\geq 3$, electrochemical potential is a variable, which results in the **S-shape** curve seen above.

Thus, by analyzing a phase diagram, we can predict the shape of such curves.
For example,

![](example-charge-discharge.jpg)

**Figure:** Example of a charge/discharge curve predicted from a (hypothetical) phase diagram, from [10.1016/j.pecs.2019.01.001](https://doi.org/10.1016/j.pecs.2019.01.001).
[SoC](https://en.wikipedia.org/wiki/State_of_charge) (state of charge) quantifies the remaining capacity available in a battery at a given time.

Conversely, observing the charge/discharge curve provides insight into the phases present in the material during operation:

![](example-charge-discharge2.jpg)

**Figure:** charge/discharge curve of a Lithium-sulfur battery, from [10.1039/C9NA00040B](https://doi.org/10.1039/C9NA00040B).
After 400 mAh g⁻¹, the profile is mostly flat.

Obviously, the changes in the curves are enhanced by taking a derivative of the curve. 
This is referred to as **differential capacity** and **differential voltage**, as discussed in [10.1021/acs.chemmater.2c01976](https://dx.doi.org/10.1021/acs.chemmater.2c01976).

![](example-diff.jpg)

**Figure:** Example of differential capacitance analysis, here for $LiFePO_4$, from [10.3390/en15134520](https://doi.org/10.3390/en15134520).

## Conclusions

In this blog post, we explored the connection between the shape of charge/discharge curves and the underlying phase equilibria. 
We demonstrated that phase diagrams are a valuable tool for understanding battery behavior, offering insight into the electrochemical processes that govern these systems.

However, it's important to remember that kinetics also play a crucial role. 
The rate at which a battery is charged or discharged, often represented by the C-rate, can significantly influence the curve's shape. 
Faster rates may lead to different profiles due to limitations in ion mobility and reaction speed.

But that's a topic for another day! 😉

![](c-rate.png)

**Figure:** Evolution of the charge/discharge curve with the C-rate, from [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5).

## Sources

+ *Advanced Batteries (Materials Science Aspects)* (book), [10.1007/978-0-387-76424-5](https://dx.doi.org/10.1007/978-0-387-76424-5).
+ Chapter 3 of *Physical metallurgy* (book), [10.1016/B978-0-444-53770-6.00003-4](https://doi.org/10.1016/B978-0-444-53770-6.00003-4).
+ *Potentially Confusing: Potentials in Electrochemistry*, [10.1021/acsenergylett.0c02443](https://dx.doi.org/10.1021/acsenergylett.0c02443).
+ *Understanding Li Diffusion in Li-Intercalation Compounds*, [10.1021/ar200329r](https://pubs.acs.org/doi/10.1021/ar200329r).
+ *Differential Analysis of Galvanostatic Cycle Data from Li-Ion Batteries: Interpretative Insights and Graphical Heuristics*, [10.1021/acs.chemmater.2c01976](https://dx.doi.org/10.1021/acs.chemmater.2c01976).

The codes of this post (including images) are available here: <https://github.com/pierre-24/blog.pierrebeaujean.net/tree/master/content/posts/0003-charge-discharge-curves>.
