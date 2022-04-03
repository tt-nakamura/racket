# Euler equation of rigid body rotation

import sympy as sp
import sympy.physics.mechanics as mec

q1,q2,q3 = mec.dynamicsymbols('q1,q2,q3')
w1,w2,w3 = mec.dynamicsymbols('w1,w2,w3')
I1,I2,I3,m = sp.symbols('I1,I2,I3,m')

n = mec.ReferenceFrame('n')
e = n.orientnew('e', 'Body', (q1,q2,q3), 'ZXY')
W = e.ang_vel_in(n)
w = w1*e.x + w2*e.y + w3*e.z
e.set_ang_vel(n,w)

O = mec.Point('O')
O.set_vel(n,0)

I = mec.inertia(e,I1,I2,I3)
B = mec.RigidBody('B', O, e, m, (I,O))

kde = [(w-W).dot(e) for e in e]

KM = mec.KanesMethod(n, (q1,q2,q3), (w1,w2,w3), kde)
KM.kanes_equations([B])
