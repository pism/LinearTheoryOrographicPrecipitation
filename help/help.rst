.. default-role:: math

This plugin implements the linear model of orographic precipitation
described in Smith and Barstad, `"A linear theory of
orographic precipitation"`_\, *Journal of the Atmospheric Sciences*
61, no. 12 (2004): 1377-1391.

The model consists of two steady-state advection equations, one for
the vertically-integrated cloud water density (`q_c`) and one
for the vertically-integrated hydrometeordensity (`q_h`).

Transforming into Fourier space, one gets the following equation
relating the Fourier transform of the terrain `\hat h` to the
transform of the precipitation field `\hat P`.

.. math::

   \hat P(k, l) =
   \frac{C_w i \sigma \hat h(k, l)}{(1 - im H_w)(1 + i \sigma \tau_f)(1 + i \sigma\tau_c)}.

Here the first factor in the denominator describes airflow dynamics,
the second and third factors describe cloud delays and advection, and

.. math::

   m = \left[ \left( \frac{N_m^{2} - \sigma^2}{\sigma^2} \right) (k^2 + l^2) \right]^{1/2}.

.. csv-table:: Notation
   :header-rows: 1

   Symbol, Units, Description
   `P`, mm / hour, Precipitation
   `C_w`, kg / m\ :sup:`3`, Uplift sensitivity factor
   `h`, m, Terrain elevation
   `H_w`, m, Vapor scale height
   `N_m`, 1 / s, Moist stability frequency
   `\tau_f`, s, Fallout time
   `\tau_c`, s, Cloud conversion time
   `\sigma`, , intrinsic frequency

.. _"A linear theory of orographic precipitation": https://doi.org/10.1175/1520-0469(2004)061%3C1377:ALTOOP%3E2.0.CO;2
