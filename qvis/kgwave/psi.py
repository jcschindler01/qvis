
"""
Classical solution based on positive-frequency complex solution.
"""

import numpy as np


class psi:

	"""
	Solution of wave equation.
	"""

	def __init__(self, n, cn, mu):
		## attributes
		self.mm = 1. * np.pi * mu
		self.n  = 1. * n.astype(float)
		self.ck = 1. * cn.astype(complex) / np.sqrt(np.sum(np.abs(cn)**2))
		self.k  = 1. * np.pi * self.n
		self.w  = 1. * np.sqrt(self.mm**2 + self.k**2)

	def psi(self, x, t):
		"""
		Positive freq complex solution.
		psi(x,t) = sum_k c_k e^(-iw(k)t) sqrt(2) sin(kx)
		"""
		## initialize
		z = 0.0j * x
		## calculate
		for i in range(len(self.n)):
			z += self.ck[i] * np.exp(-1.0j*self.w[i]*t) * np.sqrt(2) * np.sin(self.k[i]*x)
		## return
		return z

	def psidot(self, x, t):
		"""
		Time derivative.
		psidot(x,t) = sum_k (-iw(k)c_k) e^(-iw(k)t) sqrt(2) sin(kx)
		"""
		## initialize
		z = 0.0j * x
		## calculate
		for i in range(len(self.n)):
			z += -1.0j * self.w[i] * self.ck[i] * np.exp(-1.0j*self.w[i]*t) * np.sqrt(2) * np.sin(self.k[i]*x)
		## return
		return z

	def psiprime(self, x, t):
		"""
		Spatial derivative.
		psi(x,t) = sum_k k c_k e^(-iw(k)t) sqrt(2) cos(kx)
		"""
		## initialize
		z = 0.0j * x
		## calculate
		for i in range(len(self.n)):
			z += self.k[i] * self.ck[i] * np.exp(-1.0j*self.w[i]*t) * np.sqrt(2) * np.cos(self.k[i]*x)
		## return
		return z

	def y(self, x, t):
		"""
		Real part.
		"""
		return self.psi(x,t).real

	def ydot(self, x, t):
		"""
		Real part.
		"""
		return self.psidot(x,t).real

	def yprime(self, x, t):
		"""
		Real part.
		"""
		return self.psiprime(x,t).real

	def energy(self, x, t):
		"""
		H = (1/2) * ( ydot^2 + yprime^2 + m^2 y^2)
		"""
		return 0.5 * ( self.ydot(x,t)**2 + self.yprime(x,t)**2 + self.mm**2 * self.y(x,t)**2)




def psi_from_init(x0, y0, ydot0, mu, nmax=10):
	"""
	Determine cn based on initial data for 1 <= n <= nmax.
	Formula is 
	c_k = xi_k . chi_k
	where
	xi_k = sqrt(2) sin(kx)
	and 
	chi_k = y0 + i ydot0 / w_k .
	"""
	## init
	mm = 1. * np.pi * mu
	n  = np.arange(1., float(nmax), 1.)
	k  = 1. * np.pi * n
	w  = np.sqrt(mm**2 + k**2)
	cn = 0.0j * n
	dx = 1./float(len(x0))
	## find cn
	for i in range(len(n)):
		xi = np.sqrt(2.) * np.sin(k[i]*x0)
		chi = y0 + 1.0j * ydot0 / w[i]
		cn[i] = np.sum(xi*chi*dx)
	## return psi
	return psi(n, cn, mu)



def psi_from_wavefunction(x0, psi0, mu, nmax=10):
	"""
	Determine cn based on initial data for 1 <= n <= nmax.
	Formula is 
	c_k = xi_k . psi0.
	"""
	## init
	mm = 1. * np.pi * mu
	n  = np.arange(1., float(nmax), 1.)
	k  = 1. * np.pi * n
	w  = np.sqrt(mm**2 + k**2)
	cn = 0.0j * n
	dx = 1./float(len(x0))
	## find cn
	for i in range(len(n)):
		xi = 0.0j + np.sqrt(2.) * np.sin(k[i]*x0)
		cn[i] = np.sum(xi*psi0*dx)
	## return psi
	return psi(n, cn, mu)






