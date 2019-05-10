
"""
kgwave = "Klein-Gordon Wave"

Package for plotting and visualizing solutions of massive KG wave equation.

Intro
-----

Formalism follows conventions of Birrel and Davies (1982) Chapter 2 with metric signature :math:`(+---)`. Units :math:`\\hbar=c=1`.

:: 

   ......................................................
    Solve KG on interval:

              0             L/4            L/2
              |--------------|--------------|
             y=0                           y=0

    m=2pi
    L>>(2pi/m)
   ......................................................


Consider the flat-space (1+1)d massive KG wave equation

.. math::
   (\\partial_t^2 - \\partial_x^2 + m^2) \\, \\varphi = 0

for a real scalar field :math:`\\varphi`, on the interval (0, L/2) with closed (:math:`\\varphi` =0) boundary conditions.

Choose mass :math:`m=2\\pi` and measure space and time in units of the Compton wavelength :math:`\\lambda_m = 2\\pi/m =1`. The Compton wavenumber is :math:`k_m = 2 \\pi / \\lambda_m = m`.

You can immediately read off the dispersion relation

.. math::
   \\omega^2 = k^2 + m^2

which is made precise later.

The group velocity is

.. math::
   v_g = \\tfrac{d\\omega}{dk} = \\tfrac{k}{\\omega} = \\tfrac{k/m}{\\sqrt{1+(k/m)^2}} < 1,

and the phase velocity is :math:`v_p = \\tfrac{\\omega}{k} = \\tfrac{1}{v_g}`.

In true units, :math:`\\omega/k` has units of :math:`c`, and :math:`v_g = c^2 k/\\omega`.


Energy
------

The equation of motion is obtained from the Lagrangian density

.. math::
   \\mathcal{L} = \\tfrac{1}{2} \\, \\left( (\\partial_t \\varphi)^2 - (\\partial_x \\varphi)^2 - m^2 \\varphi^2 \\right).

The corresponding Hamiltonian density is

.. math::
   \\mathcal{H} = \\tfrac{1}{2} \\, \\left( (\\partial_t \\varphi)^2 + (\\partial_x \\varphi)^2 + m^2 \\varphi^2 \\right)

with canonical momentum density :math:`\\pi = \\partial_t \\varphi`.

The Hamiltonian, which depends on the (1+1) splitting provided by the coordinates, is

.. math::
   H = \\int_{0}^{L/2} \\mathcal{H} \\; dx.

The gravitational stress-energy-momentum (stress) tensor is

.. math::
   T_{ab} = \\partial_a \\varphi \\, \\partial_b \\varphi - \\mathcal{L} \\, \\eta_{ab} ,

so that

.. math::
   T_{tt} &= \\mathcal{H} \\\\
   T_{tx} = T_{xt} &= \\partial_t \\varphi \\; \\partial_x \\varphi \\\\
   T_{xx} &= \\tfrac{1}{2} \\, \\left( (\\partial_t \\varphi)^2 + (\\partial_x \\varphi)^2 - m^2 \\varphi^2 \\right).

The energy flux vector :math:`F^\\mu = {T^\\mu}_t` relative to a stationary observer (sitting at x=const) is

.. math::
   F^t &= \\mathcal{H} \\\\
   F^x &= - \\, \\partial_t \\varphi \\; \\partial_x \\varphi .

More generally, :math:`{T^t}_\\nu = T_{t\\nu}` and :math:`{T^x}_\\nu = -T_{x\\nu}`, so the up-down stress tensor in the Minkowski coordinates is

.. math::
   {T^\\mu}_\\nu = \\begin{pmatrix} T_{tt} & T_{tx} \\\\ -T_{xt} & -T_{xx} \\end{pmatrix}.

The proper density :math:`\\rho` and (minus the) proper pressure :math:`-p` are defined as the timelike and spacelike eigenvalues of :math:`{T^a}_b`, respectively. The energy conditions in this context are: NEC (:math:`\\rho+p\\geq 0`), WEC (:math:`\\rho \\geq 0` +NEC), FEC (:math:`\\rho^2 \\geq p^2`).

Define :math:`\\alpha=- \\, \\partial_t \\varphi/\\partial_x \\varphi`. One generally expects :math:`|\\alpha|>1` (since :math:`|\\omega| > |k|`). Let us suppose, for the time being, this holds true. Then by diagonalizing the stress tensor one obtains

.. math::
   \\rho &= \\mathcal{H} - (\\partial_x \\varphi)^2 \\qquad & &\\textrm{with eigenvector} \\qquad &  &(\\alpha,1) \\\\
   -p    &= \\mathcal{H} - (\\partial_t \\varphi)^2 \\qquad & &\\textrm{with eigenvector} \\qquad & &(1,\\alpha).

It's easy to check that the eigenvectors are orthonormal, and also that both :math:`\\rho>0` and :math:`\\rho+p>0` explicitly, so both NEC and WEC are satisfied. Also, :math:`\\rho^2 - p^2 = m^2 \\varphi^2 \\, ((\\partial_t \\varphi)^2 - (\\partial_x \\varphi)^2) \\geq 0`, so FEC is also satisfied.

Energy conditions are only possible, therefore, if :math:`|\\alpha| \\leq 1`.




"""
