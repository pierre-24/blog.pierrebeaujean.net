---
title: "On the electrical double layer (EDL)"
date: 2026-07-29T20:46:19+02:00
draft: true
tags: 
- lang:en
- batteries
---

You know what this blog misses? Science.
Let's do science.

Today, I would like to disentangle a concept that often comes up in my work on batteries: the [*electrical double layer*](https://en.wikipedia.org/wiki/Double_layer_%28surface_science%29).
While this concept originates from the study of capacitors, it has many applications, including in battery science. 
Here, it tells us something about the organization of ions and solvents near an electrode, and about how the electric field and potential extend into the bulk.
So, this might come in handy when studying reactions at interfaces, like what I'm trying to do in [my recent publication](https://dx.doi.org/10.1039/D6TA02841A) on the topic.

Note, however, that I will mostly focus on the continuum models developed to describe this behavior, which come with their own issues.
We'll have a look at those at the end.

## But first, (electric) capacitors!

According to [Wikipedia](https://en.wikipedia.org/wiki/Capacitor):

> A capacitor is a device that stores electrical energy by accumulating electric charges on two closely spaced surfaces that are insulated from each other.

When a potential difference is applied across the two plates of a capacitor, charges accumulate on them, creating an electric field between the plates. 
The charges on the two plates have opposite signs and therefore attract each other.
Thanks to the insulator, or, rather, the **dielectric**, between the two plates, no charge can flow from one plate to the other. 
The capacitor can therefore store electrical energy as long as the potential difference is maintained.
When the potential is removed, the stored energy can be released.

Capacitors are characterized by their **capacitance**, given by the ratio of the charge stored on each plate to the applied voltage.
Since the charge can depend on the voltage, a differential form is sometimes more useful, defined as

$$C = \frac{dQ}{dV}.$$

For a simple parallel-plate capacitor, the capacitance can be computed as

$$C = \varepsilon_0\varepsilon_r\frac{A}{d},$$

where $A$ is the area of the plates and $d$ is the distance between them.

But why do we care?

Well, because one of the simplest ways of thinking about what happens at an electrode–electrolyte interface is to see it as something that looks pretty much like a capacitor.
When the potential applied to the electrode differs from its [potential of zero charge](https://en.wikipedia.org/wiki/Point_of_zero_charge#Application_in_electrochemistry) (PZC), charge builds up at its surface. At a potential lower than the PZC, this excess surface charge is negative.
As a consequence, charges of the opposite sign in the electrolyte are attracted towards the surface, creating... two layers of opposite charge.
So, we have a capacitor. 
And a pretty good one, since $d$ is of the order of an Angstrom.
In fact, when well tuned, one might create [supercapacitors](https://en.wikipedia.org/wiki/Supercapacitor) based on this principle.

We will see below that this picture may be a bit too simple.
However, that has not prevented electrochemists from defining a differential capacitance for this situation, written as

$$C_{dl} = \frac{d\sigma_M}{d\psi_0},$$

where $\sigma_M$ is the surface charge density and $\psi_0$ is the electrode potential.

## Continuum models for the EDL

We begin with the **one-dimensional Poisson equation**, which governs the electrostatic potential $\psi(x)$ at a distance $x$ from a flat electrode surface:

$$\frac{d^2\psi(x)}{dx^2}=-\frac{\rho(x)}{\varepsilon},$$

where $\rho(x)$ is the local volume charge density and $\varepsilon=\varepsilon_0\varepsilon_r$ is the dielectric permittivity of the solvent.
In the following, it will be assumed that the potential in the solution is zero, so that $\lim_{x\to\infty} \psi(x) = 0$.

### Derivation of the Helmholtz model

Let's start with the simplest possible model of the EDL: the **Helmholtz model**.

The basic idea is very simple. 
We assume that the counter-ions form a rigid, compact layer at a fixed distance $x=d$ from the electrode surface. 
This distance is generally taken to be the **distance of closest approach**, since short-range repulsive interactions, such as Pauli repulsion, prevent the ions from getting any closer to the electrode.

In the region between the electrode ($x=0$) and the ion layer ($0<x<d$), there are no free charges in the solution. Therefore, $\rho(x)=0$, and the Poisson equation reduces to the one-dimensional **Laplace equation**:

$$\frac{d^2\psi}{dx^2}=0.$$

Integrating once with respect to $x$ gives

$$\frac{d\psi}{dx}=C_1.$$

Since the electric field is related to the potential through

$$E=-\frac{d\psi}{dx},$$

the electric field in the gap is constant.

Integrating once more gives the expected linear potential profile:

$$\psi(x)=C_1x+C_2.$$

We can now apply the appropriate boundary conditions.

1. At the electrode surface ($x=0$), the potential is $\psi_0$:

   $$\psi(0)=\psi_0,$$

   and therefore

   $$\psi_0=C_1(0)+C_2\quad\Longrightarrow\quad C_2=\psi_0.$$

2. At the plane of closest approach ($x=d$), we take the potential of the bulk solution as our reference, so that

   $$\psi(d)=0.$$

   This gives

    $$0=C_1d+\psi_0\quad\Longrightarrow\quad C_1=-\frac{\psi_0}{d}.$$

Substituting these constants into the potential profile gives

$$\psi(x)=\psi_0\left(1-\frac{x}{d}\right).$$

The potential therefore decreases linearly from $\psi_0$ at the electrode to zero at the ion layer.

We can now connect this potential profile to the surface charge density. From Gauss's law, the surface charge density on the electrode, $\sigma_M$, is related to the electric field at the electrode surface:

$$\sigma_M=-\varepsilon\left.\frac{d\psi}{dx}\right|_{x=0}.$$

Using the potential profile above,

$$\sigma_M=\varepsilon\frac{\psi_0}{d}.$$

Finally, differentiating the surface charge density with respect to the electrode potential gives the **Helmholtz differential capacitance**:

$$C_H=\frac{d\sigma_M}{d\psi_0}=\frac{\varepsilon}{d}.$$

In other words, the Helmholtz model is really a parallel-plate capacitor in disguise: the electrode provides one charged plate, while the compact counter-ion layer provides the other.

### Gouy–Chapman Model Derivation

While the Helmholtz model (which dates back to about 1850) assumed that counter-ions form a rigid layer at a fixed distance from the electrode, the **Gouy–Chapman model** (~1910), in light of the devellopment of the Debye-Hückel theory, takes a very different approach: ions are treated as point charges that are free to move in the solution and are assumed to be in thermal equilibrium.
Thus, their local number concentration $n_i(x)$ therefore follows a **Boltzmann distribution** in response to the local electrostatic potential $\psi(x)$:

$$n_i(x)=n^0_i\exp\left[-\frac{z_i e\psi(x)}{k_B T}\right] = n_0\exp[-\beta z_ie\psi(x)],$$

where $z_i$ is the valence of ion $i$, $e$ is the elementary charge, $k_B$ is the Boltzmann constant, $T$ is the temperature, and $\beta = (k_B T)^{-1}$. 
Here, $n_0$ is the bulk number concentration of each ionic species.

For a symmetric $z:z$ electrolyte, with $z_\pm=\pm z$ and $n_+^0 = n_-^0 = n_0$, the two ionic concentrations are therefore

$$n_+(x)=n_0\exp[-\beta ze\psi(x)],\qquad n_-(x)=n_0\exp[\beta ze\psi(x)].$$

The local charge density is then, using the definition of [hyperbolic sine](https://en.wikipedia.org/wiki/Hyperbolic_functions#Exponential_definitions) [$\sinh x = \frac{1}{2}(e^x-e^{-x})$]:

$$\rho(x)=ze\left[n_+(x)-n_-(x)\right]=-2zen_0\sinh[\beta ze\psi(x)].$$

Substituting this expression into the Poisson equation gives the **Poisson–Boltzmann equation**:

$$\frac{d^2\psi}{dx^2}=\frac{2zen_0}{\varepsilon}\sinh[\beta ze\psi(x)].$$

This equation is nonlinear, so we need to do a little more work than in the Helmholtz case.

To obtain a first integral, we multiply both sides by $2,d\psi/dx$:

$$2\frac{d\psi}{dx}\frac{d^2\psi}{dx^2}=\frac{4zen_0}{\varepsilon}\sinh[\beta ze\psi(x)]\frac{d\psi}{dx}.$$

The left-hand side can be recognized as the derivative of $(d\psi/dx)^2$. Integrating both sides with respect to $x$ gives

$$\left(\frac{d\psi}{dx}\right)^2=\frac{4n_0}{\beta\varepsilon}\cosh[\beta ze\psi(x)]+C.$$

We can now use the boundary conditions in the bulk solution, far away from the electrode ($x\rightarrow\infty$):

$$\psi\rightarrow0,\qquad\frac{d\psi}{dx}\rightarrow0.$$

Since $\cosh(0)=1$, we obtain

$$0=\frac{4n_0}{\beta\varepsilon}+C\quad\Longrightarrow\quad C=-\frac{4n_0}{\beta\varepsilon}.$$

Therefore,

$$\left(\frac{d\psi}{dx}\right)^2=\frac{4n_0}{\beta\varepsilon}\left[\cosh(\beta ze\psi)-1\right].$$

Using the identity $\cosh(y)-1=2\sinh^2(y/2)$, this becomes

$$\left(\frac{d\psi}{dx}\right)^2=\frac{8n_0}{\beta\varepsilon}\sinh^2\left(\frac{\beta ze\psi}{2}\right).$$

Taking the negative square root, since the potential decreases with distance from a positively charged electrode, gives

$$\frac{d\psi}{dx}=-\sqrt{\frac{8n_0}{\beta\varepsilon}}\sinh\left(\frac{\beta ze\psi}{2}\right).$$

We can now connect the potential to the charge on the electrode. 
From Gauss's law, the surface charge density on the metal electrode, $\sigma_M$, is related to the electric field at the electrode surface:

$$\sigma_M=-\varepsilon\left.\frac{d\psi}{dx}\right|_{x=0}.$$

At the electrode surface, $x=0$ and $\psi=\psi_0$, so

$$\sigma_M=\sqrt{\frac{8\varepsilon n_0}{\beta}}\sinh\left(\frac{\beta ze\psi_0}{2}\right).$$

The **Gouy–Chapman differential capacitance** is obtained by differentiating the surface charge density with respect to the surface potential:

$$C_{\mathrm{GC}}=\frac{d\sigma_M}{d\psi_0}=\sqrt{2\beta z^2e^2\varepsilon n_0}\cosh\left(\frac{\beta ze\psi_0}{2}\right).$$

We can make this expression a little more intuitive by introducing the inverse **Debye screening length**,

$$\kappa=\sqrt{\frac{2\beta z^2e^2n_0}{\varepsilon}}.$$

The Gouy–Chapman capacitance then takes the particularly simple form

$$C_{\mathrm{GC}}=\varepsilon\kappa\cosh\left(\frac{\beta ze\psi_0}{2}\right).$$

The $\cosh$ term tells us something interesting: unlike the Helmholtz capacitance, the Gouy–Chapman capacitance is not constant. It increases with the magnitude of the electrode potential. In fact, we will see below that it has a characteristic parabola-like shape.

The characteristic length scale is therefore the [**Debye length**](https://en.wikipedia.org/wiki/Debye_length#In_an_electrolyte_solution), $\lambda_D=1/\kappa$.
Note that it is inversely proportional to the square root of the concentration, $n_0$, which means that when $n_0\to\infty$, $\lambda_D\to0$.
This is, obviously, unphysical. There is, however, an interesting limiting case. When the surface potential is small, such that $\beta ze\psi_0\ll1$, we can use the approximation $\sinh(y)\approx y$. The electric-field equation then becomes

$$\frac{d\psi}{dx}\approx-\kappa\psi.$$

Integrating from the electrode surface, where $x=0$ and $\psi=\psi_0$, gives an exponential decay:

$$\psi(x)=\psi_0e^{-\kappa x}.$$

In the linearized Gouy–Chapman model, the electrostatic potential is screened exponentially over roughly $\lambda_D$.

This is already quite different from the Helmholtz picture. Instead of a sharp, compact layer of counter-ions at a fixed distance, the Gouy–Chapman model gives us a **diffuse layer**, in which the ionic concentrations gradually return to their bulk values as we move away from the electrode.


### The Stern Model

The Gouy–Chapman model gives us a much more realistic picture of the EDL than the Helmholtz model: ions are not confined to a single plane but form a **diffuse layer** whose concentration gradually approaches the bulk value.
However, it also comes with an important limitation. 
By treating ions as point charges, the Gouy–Chapman model allows them to approach the electrode arbitrarily closely. 
At sufficiently large surface potentials, this can lead to unrealistically high counter-ion concentrations near the electrode (and, for that matter, capacitance).

The **Stern model** provides a simple way to fix this by simply combining the two models. 
Close to the electrode, ions cannot approach beyond a certain distance because of their finite size and short-range repulsive interactions. 
This creates a compact, charge-free region, just as in the Helmholtz model. 
Beyond this region, ions are free to distribute themselves according to the Gouy–Chapman model.
In other words, the EDL is divided into two regions:

1. A **compact layer**, extending from the electrode surface at $x=0$ to the plane of closest approach at $x=d$.
2. A **diffuse layer**, extending from $x=d$ into the bulk solution.

The potential therefore drops partly across the compact layer and partly across the diffuse layer.
There is no free charge in the solution between $x=0$ and $x=d$, so the potential satisfies the Laplace equation given above.
As in the Helmholtz model, this gives a linear potential profile. 
However, we now need to distinguish between the electrode potential $\psi_0$ and the potential at the beginning of the diffuse layer, which we will call $\psi_d=\psi(d)$.

The potential drop across the compact layer is therefore, starting from the expression found in the Helmotz model, given by:

$$\sigma_M=\varepsilon_c\frac{\psi_0-\psi_d}{d} \Leftrightarrow \psi_0-\psi_d=\frac{\sigma_Md}{\varepsilon_c}.$$

This can be rearranged to give

$$\sigma_M=C_H(\psi_0-\psi_d),$$

where we recognize $C_H$ from above.
So far, this is exactly the same result as before. The difference is that $\psi_d$ is no longer assumed to be zero: the diffuse layer still has a finite potential at $x=d$.
Then, beyond the plane of closest approach, the Gouy–Chapman model applies. 
We can therefore reuse the result derived above, but with $\psi_d$ instead of $\psi_0$ as the potential at the beginning of the diffuse layer.

The charge density associated with the diffuse layer is

$$\sigma_M=\sqrt{\frac{8\varepsilon_{\mathrm{bulk}}n_0}{\beta}}\sinh\left(\frac{\beta ze\psi_d}{2}\right).$$

The Stern model is therefore obtained by combining the compact-layer and diffuse-layer contributions.
We can write the relationship between the electrode potential and the diffuse-layer potential as

$$\psi_0=\psi_d+\frac{\sigma_M}{C_H}.$$

Substituting the Gouy–Chapman expression for $\sigma_M$ gives an implicit equation for $\psi_d$:

$$\psi_0=\psi_d+\frac{1}{C_H}\sqrt{\frac{8\varepsilon_{\mathrm{bulk}}n_0}{\beta}}\sinh\left(\frac{\beta ze\psi_d}{2}\right).$$

This equation generally cannot be inverted analytically to obtain $\psi_d$ as an explicit function of $\psi_0$. But that is not really a problem: it already gives us the full Stern model, and $\psi_d$ can be obtained numerically for any given electrode potential.

We can go one step further and derive the differential capacitance of the complete EDL.
Since the potential in the electrolyte is zero, the total potential drop is equivalent to the sum of the potential drops across the compact and diffuse layers:

$$\psi_0=\frac{\sigma_M}{C_H}+\psi_d.$$

Differentiating with respect to $\sigma_M$ gives

$$\frac{d\psi_0}{d\sigma_M}=\frac{1}{C_H}+\frac{d\psi_d}{d\sigma_M}.$$

The second term is simply the inverse of the Gouy–Chapman differential capacitance evaluated at $\psi_d$:

$$C_{\mathrm{GC}}(\psi_d)=\frac{d\sigma_M}{d\psi_d}.$$

Therefore,

$$\frac{d\psi_0}{d\sigma_M}=\frac{1}{C_H}+\frac{1}{C_{\mathrm{GC}}(\psi_d)} \Leftrightarrow \frac{1}{C_{\mathrm{S}}}=\frac{1}{C_H}+\frac{1}{C_{\mathrm{GC}}(\psi_d)}.$$

This looks exactly like two capacitors connected in series.

The analogy is not accidental. The compact layer and the diffuse layer each sustain part of the total potential drop, and the same surface charge passes through both. Their capacitances therefore combine in series.

For the compact layer, the Helmholtz capacitance is

$$C_H=\frac{\varepsilon_{\mathrm{H}}}{d},$$

where $\varepsilon_{\mathrm{H}}$ is the permittivity of the compact layer and $d$ is its thickness.
For the diffuse layer, the Gouy–Chapman capacitance is

$$C_{\mathrm{GC}}(\psi_d)=\varepsilon_{\mathrm{bulk}}\kappa\cosh\left(\frac{\beta ze\psi_d}{2}\right),$$

where $\varepsilon_{\mathrm{bulk}}$ is the bulk solvent permittivity and $\kappa$ is the inverse Debye length.

We therefore finally obtain

$$\frac{1}{C_{\mathrm{S}}}=\frac{d}{\varepsilon_{\mathrm{H}}}+\frac{1}{\varepsilon_{\mathrm{bulk}}\kappa\cosh\left(\frac{\beta ze\psi_d}{2}\right)}.$$

There is a small but important subtlety here: the Stern capacitance is not simply obtained by putting $\psi_0$ into the Gouy–Chapman expression. The Gouy–Chapman contribution depends on the potential **at the beginning of the diffuse layer**, $\psi_d$, rather than directly on the electrode potential $\psi_0$.
This distinction becomes increasingly important as the compact-layer potential drop becomes significant.

The Stern model therefore gives us a useful compromise between the two limiting pictures.
At low potentials, or when the compact layer is very thin, the diffuse layer can dominate and the model approaches the Gouy–Chapman picture.
Conversely, if the diffuse-layer capacitance becomes very large, the compact layer becomes the limiting contribution:

$$C_{\mathrm{GC}}\gg C_H\quad\Longrightarrow\quad C_{\mathrm{S}}\approx C_H.$$

In this limit, most of the potential drop occurs across the compact layer, and the EDL behaves approximately like a Helmholtz capacitor.

### Grahame's Extension and the Helmholtz Planes (IHP & OHP)

In 1947, David C. Grahame refined the Stern model to account for the fact that some ions can lose part of their solvation shell and adsorb specifically onto the electrode surface, while others remain fully solvated. This led to a distinction between two planes within the compact part of the EDL: the **Inner Helmholtz Plane (IHP)** and the **Outer Helmholtz Plane (OHP)**.

The interface is divided along the spatial coordinate $x$, taken perpendicular to the electrode surface at $x=0$.

1. **Inner Helmholtz Plane (IHP), at $x=x_{\mathrm{IHP}}$**.

   The IHP is the locus of the centers of **specifically adsorbed ions**. These ions are typically partially or completely desolvated and can interact directly with the metal through short-range chemical and non-electrostatic interactions. Because these interactions are not purely electrostatic, specifically adsorbed ions can, in principle, adsorb even when they have the same sign of charge as the electrode.

   The region between the metal surface and the IHP, $0<x<x_{\mathrm{IHP}}$, also contains strongly oriented solvent molecules, particularly water dipoles.

2. **Outer Helmholtz Plane (OHP), at $x=x_{\mathrm{OHP}}$**.

   The OHP is the locus of the centers of **non-specifically adsorbed ions** at their distance of closest approach to the electrode. These ions retain their primary solvation shell and are assumed to interact with the electrode primarily through long-range electrostatic forces.

3. **Diffuse layer, $x>x_{\mathrm{OHP}}$**.

   Beyond the OHP, ions are free to distribute themselves according to the electrostatic potential. This region is described by the Gouy–Chapman Poisson–Boltzmann model and extends continuously into the bulk solution.

The resulting picture is therefore

$$\text{metal}\quad\longrightarrow\quad\text{IHP}\quad\longrightarrow\quad\text{OHP}\quad\longrightarrow\quad\text{diffuse layer}\quad\longrightarrow\quad\text{bulk}.$$

The total charge of the interface must vanish. If $\sigma_M$ is the charge density on the metal, $\sigma_{\mathrm{IHP}}$ is the charge density associated with specifically adsorbed ions, and $\sigma_d$ is the integrated charge density of the diffuse layer, electroneutrality requires

$$\sigma_M+\sigma_{\mathrm{IHP}}+\sigma_d=0.$$

This equation is useful because it separates the charge stored in the different parts of the EDL.

Between the metal and the IHP, there is no volume charge in the continuum description. The potential therefore obeys the Laplace equation. If the dielectric permittivity in this region is $\varepsilon_{\mathrm{IHP}}$, the electric field is constant and determined by the metal surface charge:

$$\frac{d\psi}{dx}=-\frac{\sigma_M}{\varepsilon_{\mathrm{IHP}}}.$$

Integrating from the metal surface, where $\psi(0)=\psi_0$, to the IHP, where $\psi(x_{\mathrm{IHP}})=\psi_{\mathrm{IHP}}$, gives, as in the Stern model above,

$$\psi_0-\psi_{\mathrm{IHP}}=\frac{\sigma_Mx_{\mathrm{IHP}}}{\varepsilon_{\mathrm{IHP}}}.$$

The IHP itself is represented as a sheet of charge, $\sigma_{\mathrm{IHP}}$. The electric field therefore changes discontinuously across this plane according to Gauss's law. If the permittivity immediately on either side is taken to be the same, this gives

$$\left.\frac{d\psi}{dx}\right| _ {x _ {\mathrm{IHP}}^+}-\left.\frac{d\psi}{dx}\right| _ {x _ {\mathrm{IHP}}^-}=-\frac{\sigma_{\mathrm{IHP}}}{\varepsilon_{\mathrm{OHP}}}.$$

Consequently, in the region between the IHP and OHP, the field is determined by the total charge enclosed between the metal and that point:

$$\frac{d\psi}{dx}=-\frac{\sigma_M+\sigma_{\mathrm{IHP}}}{\varepsilon_{\mathrm{OHP}}}.$$

The potential drop between the two Helmholtz planes is therefore

$$\psi_{\mathrm{IHP}}-\psi_{\mathrm{OHP}}=\frac{(\sigma_M+\sigma_{\mathrm{IHP}})(x_{\mathrm{OHP}}-x_{\mathrm{IHP}})}{\varepsilon_{\mathrm{OHP}}},$$

where $\psi_{\mathrm{OHP}}=\psi(x_{\mathrm{OHP}})$ is the potential at the OHP.

Putting the two compact regions together gives the total potential drop between the electrode and the OHP:

$$\psi_0-\psi_{\mathrm{OHP}}=(\psi_0-\psi_{\mathrm{IHP}})+(\psi_{\mathrm{IHP}}-\psi_{\mathrm{OHP}})=\frac{\sigma_Mx_{\mathrm{IHP}}}{\varepsilon_{\mathrm{IHP}}}+\frac{(\sigma_M+\sigma_{\mathrm{IHP}})(x_{\mathrm{OHP}}-x_{\mathrm{IHP}})}{\varepsilon_{\mathrm{OHP}}}.$$

This is essentially the Grahame extension of the compact-layer contribution derived above for the Stern model. The key difference is that the charge in the IHP is no longer simply part of a fixed counter-charge layer: $\sigma_{\mathrm{IHP}}$ can itself vary with the electrode charge.

We can therefore define the **compact-layer capacitance** as

$$C_{\mathrm{compact}}=\frac{\partial\sigma_M}{\partial(\psi_0-\psi_{\mathrm{OHP}})}.$$

Differentiating the potential drop above gives

$$C_{\mathrm{compact}}=\left[\frac{x_{\mathrm{IHP}}}{\varepsilon_{\mathrm{IHP}}}+\frac{x_{\mathrm{OHP}}-x_{\mathrm{IHP}}}{\varepsilon_{\mathrm{OHP}}}\left(1+\frac{d\sigma_{\mathrm{IHP}}}{d\sigma_M}\right)\right]^{-1}.$$

The expression for $C_{\mathrm{compact}}$ contains the key new ingredient introduced by Grahame: **specific adsorption**. If $\sigma_{\mathrm{IHP}}$ changes significantly with the electrode charge, the compact-layer capacitance can deviate strongly from the simple Helmholtz prediction.

This also has an important consequence for the experimentally measured differential capacitance. Specific adsorption can produce pronounced features in the capacitance: if the amount of specifically adsorbed charge changes rapidly with electrode potential, the term $d\sigma_{\mathrm{IHP}}/d\sigma_M$ can become large, strongly modifying the compact-layer response and producing the characteristic **capacitance humps** observed experimentally.

There is another useful distinction here. The compact and diffuse contributions have different physical origins. The compact-layer contribution is determined primarily by the interfacial structure and adsorption properties, whereas the diffuse-layer contribution depends strongly on the electrolyte concentration through the Debye screening length.

This separation is conceptually useful. Changing the salt concentration should primarily modify the diffuse-layer contribution, while changing the electrode material, solvent structure, or specific adsorption chemistry can modify the compact-layer contribution.

So far, however, we have only considered the region up to the OHP. For $x>x_{\mathrm{OHP}}$, the Gouy–Chapman description takes over. The diffuse layer starts at the OHP, so the relevant surface potential in the Poisson–Boltzmann solution is now $\psi_{\mathrm{OHP}}$, rather than the metal potential $\psi_0$.

The charge density in the diffuse layer is

$$\sigma_d=-\sqrt{\frac{8\varepsilon_{\mathrm{bulk}}n_0}{\beta}}\sinh\left(\frac{\beta ze\psi_{\mathrm{OHP}}}{2}\right).$$

The corresponding diffuse-layer capacitance is, once again,

$$C_{\mathrm{GC}}(\psi_{\mathrm{OHP}})=\varepsilon_{\mathrm{bulk}}\kappa\cosh\left(\frac{\beta ze\psi_{\mathrm{OHP}}}{2}\right),$$

where $\kappa$ is the inverse Debye length.

We can now regard the Grahame model as two capacitive contributions in series: the **compact layer**, extending from the metal to the OHP and containing the IHP, and the **diffuse layer**, extending beyond the OHP. Since the same interfacial charge must be balanced across the two regions, their differential capacitances combine in series.

The total differential capacitance of the Grahame model, which we denote by $C_G$, is therefore

$$\frac{1}{C_G}=\frac{1}{C_{\mathrm{compact}}}+\frac{1}{C_{\mathrm{GC}}(\psi_{\mathrm{OHP}})}.$$

Unlike the simple Helmholtz capacitance of the Stern model, however, the compact-layer capacitance is not simply determined by the thickness and dielectric permittivity of the layer. The reason is that the charge associated with specifically adsorbed ions can itself change as the electrode charge changes.

The important point is that the EDL is no longer just a question of electrostatics. Once we distinguish the IHP from the OHP, we are implicitly acknowledging that **ion–surface interactions, solvation, and specific adsorption** can affect the structure of the interface.





