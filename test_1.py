
import sys

sys.path.insert(0, '.')

import squanf_class

Ne = 4
Nh = 0

## test without spin
eh = squanf_class.squanf(Ne=Ne, Nh=Nh)

print 'test w/o spin'
print 'list numbers from 0 to %i' % (eh.Ntot - 1)

# holes
for state in xrange(eh.Nh):
    print 'h' + str(state) + '->' + str(eh._map_idx('h', state))

# electrons
for state in xrange(eh.Ne):
    print 'e' + str(state) + '->' + str(eh._map_idx('e', state))

## test with spin
eh = squanf_class.squanf(Ne=Ne, Nh=Nh, use_spin=True)

print 'test w spin'
print 'list numbers from 0 to %i' % (eh.Ntot - 1)

# holes
for state in xrange(eh.Nh):
    for spin in ['u', 'd']:
        print 'h' + str(state) + spin + '->' \
            + str(eh._map_idx('h', state, spin=spin))

# electrons
for state in xrange(eh.Ne):
    for spin in ['u', 'd']:
        print 'e' + str(state) + spin + '->' \
            + str(eh._map_idx('e', state, spin=spin))
