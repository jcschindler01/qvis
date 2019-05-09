
"""
kgwave = "Klein-Gordon Wave"

Package for plotting and visualizing solutions of massive KG wave equation.

Formalism follows conventions of Birrel and Davies (1982) with metric signature :math:`(+---)`. Units :math:`\\hbar=c=1`.

:: 

   ......................................................
    Solve KG on interval:

              0             L/4            L/2
              |--------------|--------------|
             y=0                           y=0

    m=2pi
    L>>(2pi/m)
   ......................................................


Consider the 2d massive KG wave equation

.. math::
   (\\partial_t^2 - \\partial_x^2 + m^2) \\, \\varphi = 0

for a real scalar field :math:`\\varphi`, on the interval (0, L/2) with closed (:math:`\\varphi` =0) boundary conditions.

Choose mass :math:`m=2\\pi` and measure space and time in units of the Compton wavelength :math:`\\lambda_m = 2\\pi/m =1`.

"""
