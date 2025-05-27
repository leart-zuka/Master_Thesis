---
id: "Light_Modes"
date: "2025-05-13"
subtitle: 
tags:
    - light
---

# Light_Modes

When some light beam propagates in free space or in a transparent homogenius medium (like a [fibers](general_knowledge/Fibers.md)), its transverse intensity profile generally changes during propagation.

These profiles however are not arbitrary but are self-consistent during propagation, they are known as what we call `modes`.

## Free-space

All light waves must, in order to be considered a wave, satisfy the function known as `wavefunction`:
$$ \nabla^2 u - \frac{1}{c^2 \frac{\partial^2 u}{\partial t^2}} = 0 $$

Plugging a complex wavefunction of the form $U(\vec{r},t) = U(\vec{r})exp(i2\pi\nu t)$ into our wavefunction gives us the so called `Helmholtz equation`:
$$(\nabla^2 + k^2)U(\vec{r}) = 0$$
with 
$$k=\frac{2\pi\nu}{c} = \frac{\omega}{c}$$
which is called the `wavenumber`

There are now different kinds of waves which solve the Helmholtz equation, which shall be covered now

### Plane Waves

The mathematically easiest special wave that solves the Helmholtz equation is the so called `Plane Wave` 
![plane_wave.png](assets/imgs/plane_wave.png)

It has the form:
$$U(\vec{r}) = A exp(-i \vec{k}\cdot\vec{r})$$
where A is a complex constant called the `complex envelope` and k is the so called `wavevector`.
Technically one can also include a phase of the form $-t\cdot\omega$ which is typically included in the exponent, here it's encoded in the complex envelope.
The intensity profile would then be $I(\vec{r}) = |A|^2$, meaning it has constant intensity in space so that it carries infinite power. This wave is clearly an idealization since it exists everywhere and at all times.


### Paraxial Waves

A wave is said to be paraxial if its wavefront normals are paraxial rays. One way of constructing a paraxial wave is to start with a place wave, and modulate its complex envelope A, making it a slowly varying function of position $A(\vec{r})$ so that the complex amplitude of the modulated wave becomes:
$$U(\vec{r}) = A(\vec{r})exp(-ikz)$$
And even though this wave satisfies the Helmholtz equation, the complex envelope $A(\vec{r})$ must satisfy another partial differential equation obtained by plugging the modulated complex amplitude into the Helmholtz equation itself, which reads:
$$\nabla^2_\perp A - i2k \frac{\partial A}{\partial z} = 0$$

which is called the `slowly varying envelope approximation`, as we're assuming that the complex amplitude varies slowly with respect to z. But we shall simply call it the `paraxial Helmholtz equation`.

---

The reason why this detour was important is because we're going to need to use that function in order to work with the Gaussian Beam

### Gaussian Beam

Just like before the `Gaussian Beam` is a paraxial wave. Meaning it's a plane wave $e^{-ikz}$ modulated by a complex envelope $A(\vec{r})$ that is slowly varying function of position. The complex amplitude is:
$$U(\vec{r}) = A(\vec{r}) exp(-ikz)$$

Solving the paraxial Helmholtz equation provides paraboloidal wave for which 
$$A(\vec{r}) = \frac{A_1}{z} exp(-ik\frac{\rho^2}{2z}),\quad \rho^2 = x^2 + y^2$$
where $A_1$ is a constant.

Shifting the envelope by a purely imaginary part (e.g. $\xi = -iz_0$ where $z_0$ is the `Rayleigh range`), again gives us a solution of the paraxial Helmholtz equation. It basically looks exactly the same as before just instead of the solution depending on $q(z) = z+iz_0$ which we then call `Complex Envelope`.

In order to seperate the amplitude and phase of this complex envelope, we write $\frac{1}{q(z)} = \frac{1}{z+iz_0}$ in terms of real and imaginary parts defined by two new real function $R(z)$ and $W(z)$, such that
$$\frac{1}{q(z)} = \frac{1}{R(z)} - i\frac{\lambda}{\pi W^2(z)}$$.

We can plug all of this back into our original complex amplitude giving us now:
![guassian_beam_plus_beam_params.png](assets/imgs/guassian_beam_plus_beam_params.png)

This beam now has the very distinct intensity profile:
![intensity_profile_gaussian_beam.png](assets/imgs/intensity_profile_gaussian_beam.png)

### Hermite - Gaussian Beams

And again the `Gaussian Beams` are not the only solutions to the `paraxial Helmholtz equation`.
Consider a wave whose complex envelope is a cmodulated version of the Gaussian beam,
![hermite_gaussian_beam_equation.png](assets/imgs/hermite_gaussian_beam_equation.png)

where $\mathfrak{X(\cdot)}$, $\mathfrak{Y(\cdot)}$, $\mathfrak{Z(\cdot)}$, are real functions, and $A_G(x,y,z)$ is the gaussian beam equation.
Plugging this equation into the paraxial Helmholtz equation gives us the following:
![hermite_guassian_beam_paraxial_helmholtz_equation.png](assets/imgs/hermite_guassian_beam_paraxial_helmholtz_equation.png)

Where we can then reduce this equation into three ordinary differential equations in accordance with `separation of variables`:
![helmholtz_gaussian_beam_differential_equation.png](assets/imgs/helmholtz_gaussian_beam_differential_equation.png)

The solutions to these differential equations are well known, as they are the `Hermite Polynomials`.
Plugging all of this into the equation for our complex amplitude we get:
![hermite_gaussian_beam_complex_amplitude_equation.png](assets/imgs/hermite_gaussian_beam_complex_amplitude_equation.png)

This equation has the following intensity profile:
![intensity_profile_hermite_gaussian_beam.png](assets/imgs/intensity_profile_hermite_gaussian_beam.png)

## Fiber Modes

Again just like in free space, our electric and magnetic fields must satisfy the Helmholtz equation:

$$\nabla^2U+ n^2k^2_0U=0$$

where $n=n_1$ in the core ($r < a$) and $n=n_2$ in the cladding ($ r > a $) and $k_0 = \frac{2\pi}{\lambda_0}$.
Assuming radius b of the cladding is sufficiently large so we can assume it's infinitely large, in cylinder coordinates the Helmholtz equation is:
![helmholtz_equation_in_fiber.png](assets/imgs/helmholtz_equation_in_fiber.png)

A function for our complex amplitude `U` which solves this equation is of the form:
$U(r,\phi,z) = u(r)e^{-il\phi}e^{-i\beta z},\quad l=0,\pm1, \pm2$

Writing down now the ordinary differential equation for $u(r)$ for both the cladding and core seperately we get:
![ode_core_cladding.png](assets/imgs/ode_core_cladding.png)

We can see that the ODEs get solved by `Bessel functions`, where the parameters $k_T$ and $\gamma$ need to fullfill the following relation:
$$k_T^2 + \gamma^2 = (n_1^2 - n_2^2)k_0^2 = NA^2\cdot k_0^2$$

The two parameters can then be redefined into the parameters $X=k_Ta$ and $Y=\gamma a$, finally giving us the equation:
$$X^2 + Y^2 = V^2$$

with $V=2\pi\frac{a}{\lambda_0}NA$ being the `fiber parameter` or `V parameter`, which is an important parameter that governs the number of modes of the fiber and their propagation constant.

Now using the equation for the complex amplitude for the z components of our electric and magnetic field and inserting the boundary condition $r=a$, allows us to establish a relation between the coefficients of proportionality between the function for the cladding and the one for the core. Using Maxwell's equations we ge the radial components.

For each azimuthal index `l`, the characteristic equation has multiple solutions yielding discrete propagation constant $\beta_{lm}, m=1,2,...$, each solution representing a mode, with each corresponding $k_T$ and $\gamma$ governing the spatial distributions in the core and in the cladding respectively.

Now in a weakly guiding fiber (i.e., $n_1 \approx n_2$ or $\Delta \ll 1$) so that guided rays are paraxial. The longitudinal components of the electric and magnetic fields are then much weaker than the transverse components and the guided waves are approximately transverse electromagnetic (TEM).
Linear polarized modes are usually denoted as $LP_{lm}$.
We can then get the characteristic equation obtained just like before, with the addition that the functions need to have a continuous derivative at $r = a$, thus giving us:
![weakly_guiding_fibers_relation_equation_bessel_equations.png](assets/imgs/weakly_guiding_fibers_relation_equation_bessel_equations.png)

Using the relations for the derivatives of the Bessel functions and using the normalized parameters defined before we get the `characteristic equation`:
![characteristic_equation_weakly_guiding_fiber.png](assets/imgs/characteristic_equation_weakly_guiding_fiber.png)

This equation can then be solved graphically and thus giving us distributions for each mode.


## Cavity/Resonator Modes

Looking at the simplest case (i.e. Planar-Mirrow Resonators), we know that our EM-wave and its complex amplitude both have to satisfy the Helmholtz equation.
Furthermore we note that at the boundaries ($z=0$ and $z=d$), the transverse component of the electric field vanishes, so that $U(\vec{r})=0$, thus giving us:

$$U(\vec{r})=A \cdot sin(kz)$$

From this we know, that in order to satisfy the boundary conditions our k must be equal to:

$$k_q = \frac{q\pi}{d}$$

in order for our wave to gain a phase of $2\pi$ after a full roundtrip, and maintain its shape.
`q` in this case is called the `mode number`.
Since our k is limited to certain values of q, it follows that the frequency is restricted to discrete values:
$$\nu_q = q\frac{c}{2d}$$,
which are the resonance frequencies of the resonator, which are seperated by a constant frequency difference of:
$$\nu_f = \frac{c}{2d}$$

### Resonance Frequencies of Spherical-Mirror Resonators

#### Guassian Modes

Spherical mirror resonators behave a bit more differently. We know that the guassian beam (check [light_Modes](general_knowledge/Light_Modes.md)) has a phase of:

$$\varphi(\rho,z) = kz - \xi(z) + \frac{k\rho^2}{2R(z)}$$

with $\xi(z) = tan^{-1}(\frac{z}{z_0})$ and $\rho^2 = x^2+y^2$.
Looking at the phase difference at mirror 1 ($z=z_1$) and mirror 2 ($z=z_2$) which is:

$$\varphi(0,z_0) - \varphi(0,z_1) = k(z_2-z_1) - [\xi(z_2) - \xi(z_1)] \newline 
kd-\Delta\xi$$

where $\Delta\xi=\xi(z_2)-\xi(z_1)$.

Since our beam does this trip twice before returning to its original position our phase difference times 2, which then needs to be equal to $2\pi q, q=0,\pm1,\pm2,...$.
The condition then becomes:
$$\nu_q=q\nu_F + \frac{\Delta\xi}{\pi}\nu_F$$

#### Hermite-Gaussian Modes

For a Hermite-Gaussian beam, the situation is similar however with a little bit of an extra on top.
The phase of the (l,m) mode on the beam axis is:
$$\varphi(0,z) = kz - (l+m+1)\xi(z)$$

Doing basically the same steps again we get the following resonance frequencies:
$$\nu_{l,m,q} = q\nu_F + (l+m+1)\frac{\Delta\xi}{\pi}\nu_F$$

`NOTE: 
        Modes of different q, but the same (l,m), have identical intensity distributions. They are known as longitudinal or axial modes.
        The indices (l,m) label different spatial dependences on the transverse coordinates x,y; these represent different transverse modes`


##### References
Saleh, Bahaa & Teich, Malvin. (1991). Fundamentals of Photonics. 
