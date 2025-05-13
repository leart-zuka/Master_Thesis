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

The mathemtically easiest special wave that solves the Helmholtz equation is the so called `Plane Wave` 
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

##### References
Saleh, Bahaa & Teich, Malvin. (1991). Fundamentals of Photonics. 
