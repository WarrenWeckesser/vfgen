#
# test_scipy.py
#

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import odeint
import linearosc


z0 = np.array([1, 0])
t = np.array([0, np.pi])

ysol = odeint(linearosc.vectorfield, z0, t, Dfun=linearosc.jacobian,
              atol=1e-12, rtol=1e-9, tfirst=True)

assert_allclose(ysol[-1], [-1, 0], rtol=1e-6, atol=1e-9)
