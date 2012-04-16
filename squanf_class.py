#!/usr/bin/env python

from scipy.sparse import csr_matrix as sm
from scipy.sparse import kron as kron
from scipy.sparse import identity as sp_identity
import numpy as np

# this is based on eq. (2.34) in
# http://e-collection.library.ethz.ch/eserv/eth:30080/eth-30080-02.pdf

# ************ base functions for creating fermion ops


def signed_identity(n):
    """
    returns the kronecker product of n signed identity 2 x 2 matrices
    sI = [[1, 0],
          [0, -1]]
    out = sI \otimes sI \otimes sI \otimes .....
    """
    sI = sm([[1, 0], [0, -1]])
    out = sI.copy()

    for x in xrange(n - 1):
        out = kron(out, sI)
    return out


def unsigned_identity(n):
    """
    returns the kronecker product of n 2 x 2 identity matrices
    """
    return sp_identity(2 * n)


def f_destroy(n, Ntot):
    """
    returns a fermionic descruction operator for the n'th single-particle
    state out of N_tot states.
    assumes that the states are ordered as:
    |n_0, n_1, n_2, .... , n_(Ntot - 2), n_(Ntot - 1)>
    """

    if n < 0 or n >= Ntot:
        raise ValueError('n has to be in [0, 1, ..., Ntot - 1]')

    # destruction operator for 1-state space
    # |basis 0> = |0>, |basis 1> = c^+|0>,
    c = sm([[0, 0], [1, 0]])

    if n == 0:
        return sm(kron(c, unsigned_identity(Ntot - 1)))
    elif n == Ntot - 1:
        return sm(kron(signed_identity(Ntot - 1), c))
    else:
        # signed identity and c
        tmp = kron(signed_identity(n), c)
        # this times unsigned indentity
        return sm(kron(tmp, unsigned_identity(Ntot - 1 - n)))


def f_create(n, Ntot):
    """
    hc of f_destroy
    """
    return f_destroy(n, Ntot).conj().transpose()


# *********** end of base functions

class squanf(object):
    """
    class for creating fermionic second quantization operators
    """

    def __init__(self, Ne=0, Nh=0, use_spin=False):

        # number of electron single-particle states
        self.Ne = Ne

        # number of hole single-particle states
        self.Nh = Nh

        # use spin
        self.use_spin = use_spin

        # total number of single-particle states
        if self.use_spin:
            self.Ntot = 2 * self.Ne + 2 * self.Nh
        else:
            self.Ntot = self.Ne + self.Nh

    def _map_idx(self, eh, state, spin=None):
        """
        maps eh, state, spin to a single number according to the following
        order in occupation number basis:
        |n_h0u, n_h0d, n_h1u, n_h1d, .... , n_e0u, n_e0d, n_e1u, n_e1d, ... > =
        (h0u^+)^n_h0u * (h0d^+)^n_h0d * (h1u^+)^n_h1u * (h1d^+)dn_h1u *
        .... * (e0u^+)^n_e0u * (e0d^+)^n_e0d * (e1u^+)^n_e1u * (e1d^+)^n_e1d *
         ... |0>
        """

        if ((spin is None and self.use_spin)
            or (spin is not None and not self.use_spin)):
            raise ValueError('inconsistent spin var')
        elif spin is not None:
            if spin not in ['u', 'd']:
                raise ValueError('spin should be "u" or "d"')

        if self.use_spin:
            Nspin = 2
            if spin == 'u':
                spin = 0
            else:
                spin = 1
        else:
            Nspin = 1
            spin = 0

        if eh == 'h':
            return spin + Nspin * state
        elif eh == 'e':
            h_max_idx = self.use_spin + Nspin * (self.Nh - 1) + 1
            return h_max_idx + spin + Nspin * state
        else:
            raise ValueError('eh argument must be "e" or "h"')

    def create(self, eh, state, spin=None):
        """
        returns creation operator either electron or hole

        input:
            eh: "h" for hole and "e" for electron
            state: state index \elem [0; N_eh]
            spin (dep. on use_spin): "u" for spin up and "d" for spin down

        output:
            sparse matrix representation of operator
        """
        return f_create(self._map_idx(eh, state, spin=spin), self.Ntot)

    def destroy(self, eh, state, spin=None):
        """
        returns destruction operator either electron or hole

        input:
            eh: "h" for hole and "e" for electron
            state: state index \elem [0; N_eh]
            spin (dep. on use_spin): "u" for spin up and "d" for spin down

        output:
            sparse matrix representation of operator
        """
        return f_destroy(self._map_idx(eh, state, spin=spin), self.Ntot)
