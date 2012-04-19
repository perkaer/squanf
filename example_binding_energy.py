import sys
from scipy.sparse import csr_matrix as sm


sys.path.insert(0, '.')

import squanf_class

Ne = 1
Nh = 1

eh = squanf_class.squanf(Ne=Ne, Nh=Nh, use_spin=True)

# short hands

e = eh.e
eD = eh.eD
h = eh.h
hD = eh.hD

# Hamiltonian

H = sm(e(0, 'u').shape)

Eex = 1.
V1 = 20.
V2 = 10.

for spin in ['u', 'd']:
    H = H + Eex / 2 * eD(0, spin) * e(0, spin)
    H = H + Eex / 2 * hD(0, spin) * h(0, spin)

print H.todense().astype(int)
