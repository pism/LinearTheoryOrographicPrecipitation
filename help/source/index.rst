.. LinearTheoryOrographicPrecipitation documentation master file, created by
   sphinx-quickstart on Sun Feb 12 17:11:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LinearTheoryOrographicPrecipitation's documentation!
============================================

Contents:

.. toctree::
   :maxdepth: 2

The model explicitly represents of a rich set of processes: the linear atmospheric response to flow over terrain; the stability of the atmospheric column to vertical displacement; the vertical depth scale of atmospheric moisture; and the timescales associated with the cloud processes that convert condensation into hydrometeors (snow or rain), and the fallout of those hydrometeors. The model captures the orographic rain shadow as well as patterns at the ridge-valley scale. It is linear and so can be solved efficiently using Fourier transform methods; and it can be used for idealized cross sections across the ice cap as well as for the full three-dimensional terrains. A preliminary calculation for the northern Antarctic peninsula (Fig.~1D) shows that the Smith and Barstad model is capable to emulate the along-wind precipitation predicted by the high-resolution regional climate model RACMO/ANT. The Smith and Barstad model consists of two simple equations:

.. math::

\frac{Dq_c}{Dt} & \approx & {\bf u} \cdot \nabla q_{c} & = &S\left({\bf u}, \nabla z, N^{2}, H_w \right)-\frac{q_{c}}{\tau_{c}}, \\
\frac{Dq_h}{Dt} & \approx & {\bf u} \cdot \nabla q_{h} & = & \frac{q_{c}}{\tau_{c}} - \frac{q_{h}}{\tau_{f}},

where :math:`q_{c}` and :math:`q_{h}` are the vertically integrated cloud water density and hydrometeor density respectively; and :math:`{\bf u}` is the prevailing wind at low levels. :math:`S` is the condensation source function and depends on :math:`{\bf u}`, orographic slopes ($\nabla z$), atmospheric stability (:math:`N^{2})`, and vertical moisture scale height (:math:`H_{w}`). :math:`\tau_{c}` is the time required to convert cloud water into hydrometeors, and :math:`\tau_{f}` is the time required for hydrometeors to reach the ground. Both timescales are :math:`\mathcal{O}(10^{3}`s). Precipitation, :math:`P`, is given by the rate of hydrometeor fallout, :math:`P = \frac{q_{h}}{\tau_{f}}`.
