# Fast, isolated PARI ellrank probe for the U=21 elliptic curve.
from sage.all import *
from sage.env import SAGE_VERSION
import os, traceback

A = ZZ(47643948)
B = ZZ(159472346229504)
C = A*A - 4*B
E = EllipticCurve(QQ, [0,A,0,B,0])
Ep = EllipticCurve(QQ, [0,-2*A,0,C,0])
P = E([ZZ(7784100), ZZ(67822243920)])
Pp = Ep([ZZ(14594052), ZZ(81434810160)])
effort = ZZ(os.environ.get('PARI_EFFORT','0'))

print('SAGE_VERSION=', SAGE_VERSION, flush=True)
print('EFFORT=', effort, flush=True)
print('E=', E, flush=True)
print('EP=', Ep, flush=True)

for name, curve, point in [('E',E,P),('EP',Ep,Pp)]:
    print('\nBEGIN_ELLRANK_'+name, flush=True)
    try:
        raw = curve.pari_curve().ellrank(effort, [[point[0],point[1]]])
        print('RAW=', raw, flush=True)
        print('LOWER=', raw[0], flush=True)
        print('UPPER=', raw[1], flush=True)
        print('CASSELS_S=', raw[2], flush=True)
        print('POINTS=', raw[3], flush=True)
    except BaseException as exc:
        print('FAIL=', type(exc).__name__, str(exc), flush=True)
        traceback.print_exc()
    print('END_ELLRANK_'+name, flush=True)

print('\nBEGIN_SAGE_PROVEN_RANK', flush=True)
try:
    print(E.rank(use_database=False, algorithm='pari', pari_effort=effort, proof=True), flush=True)
except BaseException as exc:
    print('FAIL=', type(exc).__name__, str(exc), flush=True)
    traceback.print_exc()
print('END_SAGE_PROVEN_RANK', flush=True)
