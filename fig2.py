# eigen value analysis of linearized Euler equation

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from EulerEq import w1,w2,w3,I1,I2,I3,KM

I = {I1: 1.3, I3: 13.7}
w0 = {w1:0, w2:1, w3:0}
A,_,_ = KM.linearize(op_point=w0, A_and_B=True)
A = A[3:,3:] # discard KDE

plt.figure(figsize=(5, 3.75))

x = np.linspace(0j, 20, 200)
for l in A.eigenvals():
    if l==0: continue
    f = sp.lambdify(I2, l.subs(I))
    plt.plot(x, np.real(f(x)), 'r')
    plt.plot(x, np.imag(f(x)), 'b')

plt.legend([r'Re $\lambda$', r'Im $\lambda$'])
plt.xlabel(r'$I_2$ / ${\rm g m^2}$')
plt.ylabel(r'Re $\lambda$, Im $\lambda$')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
