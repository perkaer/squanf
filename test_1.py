
import sys

sys.path.insert(0, '.')

import squanf_class

Ne = 2
Nh = 3

## test without spin
eh = squanf_class.squanf(Ne=Ne, Nh=Nh)

print 'test w/o spin'
print 'list numbers from 0 to %i' % (eh.Ntot - 1)

# holes
for state in xrange(eh.Nh):
    print '*' * 10
    print 'h' + str(state) + '->' + str(eh._map_idx('h', state))
    print eh.create('h', state).shape

# electrons
for state in xrange(eh.Ne):
    print '*' * 10
    print 'e' + str(state) + '->' + str(eh._map_idx('e', state))
    print eh.create('e', state).shape

## test with spin
eh = squanf_class.squanf(Ne=Ne, Nh=Nh, use_spin=True)

print 'test w spin'
print 'list numbers from 0 to %i' % (eh.Ntot - 1)

# holes
for state in xrange(eh.Nh):
    for spin in ['u', 'd']:
        print '*' * 10
        print 'h' + str(state) + spin + '->' \
            + str(eh._map_idx('h', state, spin))
        print eh.create('h', state, spin).shape

# electrons
for state in xrange(eh.Ne):
    for spin in ['u', 'd']:
        print '*' * 10
        print 'e' + str(state) + spin + '->' \
            + str(eh._map_idx('e', state, spin))
        print eh.create('e', state, spin).shape
