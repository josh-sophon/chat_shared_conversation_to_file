# Exact U=53 global descent/rank/search harness for Structured Unit-Rigidity.
from sage.all import *
from sage.env import SAGE_VERSION
from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
import sys, traceback

proof.all(True)
mode = sys.argv[1]
param = ZZ(sys.argv[2]) if len(sys.argv) > 2 else ZZ(0)
U=ZZ(53)
E0=U*U-12*U+64
C0=23*U*U-304*U+1472
R0=U**4-100*U*U+4096
A=36*E0*C0
B=16128*E0*E0*R0
C=A*A-4*B
E=EllipticCurve(QQ,[0,A,0,B,0])
Ep=EllipticCurve(QQ,[0,-2*A,0,C,0])
D1=ZZ(-17031795449)
D2=ZZ(-15659)

def emit(k,v): print('%s=%s' % (k,repr(v)),flush=True)
def attempt(k,fn):
    try: emit(k,fn())
    except BaseException as exc:
        emit(k+'_ERROR',(type(exc).__name__,str(exc)))
        traceback.print_exc()

emit('SAGE_VERSION',SAGE_VERSION)
emit('MODE',mode); emit('PARAM',param)
emit('U_PARAMS',(U,E0,C0,R0,A,B,C))
emit('E',E); emit('EP',Ep)
emit('TORSION',E.torsion_subgroup())
emit('ROOT_NUMBER',E.root_number())
emit('COVER_ISOMORPHISM_CHECK',bool(D1*D2*48**2==B))

if mode == 'rank-pari':
    effort=int(param)
    attempt('ELLRANK_E',lambda:E.pari_curve().ellrank(effort))
    attempt('ELLRANK_EP',lambda:Ep.pari_curve().ellrank(effort))
    if effort <= 10:
        attempt('PROVEN_RANK_E',lambda:E.rank(use_database=False,algorithm='pari',pari_effort=effort,proof=True))
        attempt('GENS_E',lambda:E.gens(use_database=False,algorithm='pari',pari_effort=effort,proof=True))
        try:
            gens=E.gens(use_database=False,algorithm='pari',pari_effort=effort,proof=True)
            emit('SATURATION_E',E.saturation(gens)); emit('REGULATOR_E',E.regulator_of_points(gens))
            for i,P in enumerate(gens): emit('GEN_E_%s'%i,(P[0],P[1],sign(P[0])))
        except BaseException as exc: emit('SATURATION_E_ERROR',(type(exc).__name__,str(exc)))
elif mode == 'rank-mwrank':
    lim=int(param)
    attempt('RANK_MWRANK',lambda:E.rank(use_database=False,algorithm='mwrank_lib',verbose=True,proof=True))
    attempt('GENS_MWRANK',lambda:E.gens(use_database=False,algorithm='mwrank_lib',proof=False,descent_second_limit=lim,sat_bound=1000))
    try:
        gens=E.gens(use_database=False,algorithm='mwrank_lib',proof=False,descent_second_limit=lim,sat_bound=1000)
        emit('SATURATION_MWRANK',E.saturation(gens)); emit('REGULATOR_MWRANK',E.regulator_of_points(gens))
        for i,P in enumerate(gens): emit('GEN_MWRANK_%s'%i,(P[0],P[1],sign(P[0])))
    except BaseException as exc: emit('SATURATION_MWRANK_ERROR',(type(exc).__name__,str(exc)))
elif mode == 'twoiso':
    H=int(param)
    attempt('TWOISO',lambda:two_descent_by_two_isogeny(E,global_limit_small=min(H,1000),global_limit_large=H,verbosity=4,selmer_only=False,proof=True))
elif mode == 'quartic':
    d=ZZ(param)
    H=ZZ(sys.argv[3])
    f=PolynomialRing(QQ,'x').gen()
    Q=d*f**4+A*f**2+B/d
    emit('D',d);emit('BOUND',H);emit('Q',Q)
    attempt('HYPERELLRATPOINTS',lambda:pari(Q).hyperellratpoints(H,1))
elif mode == 'dual-quartic':
    d=ZZ(param)
    H=ZZ(sys.argv[3])
    f=PolynomialRing(QQ,'x').gen()
    Q=d*f**4-2*A*f**2+C/d
    emit('D',d);emit('BOUND',H);emit('Q',Q)
    attempt('HYPERELLRATPOINTS',lambda:pari(Q).hyperellratpoints(H,1))
else:
    raise ValueError(mode)
emit('DONE',1)
