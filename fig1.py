# tennis racket theorem on the stability of rigid body rotation
# referece: V. D. Barger and M. G. Olsson
#   "Classical Mechanics: A Modern Perspective"
#    2nd edition, section 7.7

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from EulerEq import I1,I2,I3,w1,w2,w3,KM

label = [r'$\omega_1$', r'$\omega_2$', r'$\omega_3$']
plt.figure(figsize=(5,6))

t = sp.symbols('t')
I = {I1: 1.3, I2: 12.4, I3: 13.7}
f = KM.rhs().subs(I)[3:] # discard KDE
f = sp.lambdify([t,(w1,w2,w3)], f)

t = np.linspace(0,5,200)

init = [8, -1, 2]
s = solve_ivp(f, t[[0,-1]], init, t_eval=t)
plt.subplot(311)
plt.plot(t, s.y.T)
plt.ylabel(r'$\omega$ / rad/sec')
plt.legend(label, loc='upper right')

init = [2, 8, -1]
s = solve_ivp(f, t[[0,-1]], init, t_eval=t)
plt.subplot(312)
plt.plot(t, s.y.T)
plt.ylabel(r'$\omega$ / rad/sec')

init = [-1, 2, 8]
s = solve_ivp(f, t[[0,-1]], init, t_eval=t)
plt.subplot(313)
plt.plot(t, s.y.T)
plt.xlabel(r'$t$ = time / sec')
plt.ylabel(r'$\omega$ / rad/sec')
plt.legend(label, loc='upper right')

plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
