import pickle
import sys

pf = sys.argv[1]

with open(pf, 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    allobjs = pickle.load(f)

print(allobjs)
print(type(allobjs))
print(len(allobjs))
