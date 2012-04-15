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

    def _map_idx(self, eh, state, *spin):
        """
        maps eh, state, spin to a single number according to the following
        order in occupation number basis:
        |n_h0u, n_h0d, n_h1u, n_h1d, .... , n_e0u, n_e0d, n_e1u, n_e1d, ... > =
        (h0u^+)^n_h0u * (h0d^+)^n_h0d * (h1u^+)^n_h1u * (h1d^+)dn_h1u *
        .... * (e0u^+)^n_e0u * (e0d^+)^n_e0d * (e1u^+)^n_e1u * (e1d^+)^n_e1d *
         ... |0>
        """

        pass

    def create(self, eh, state, *spin):
        """
        returns creation operator either electron or hole
        """

        # call _map_idx to get uber state number and use base fcts
        pass
