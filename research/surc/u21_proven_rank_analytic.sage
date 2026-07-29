# Independent rigorous rank probe using Sage's no-conjecture rank path.
from sage.all import *
from sage.env import SAGE_VERSION
import traceback

proof.all(True)
A=ZZ(47643948)
B=ZZ(159472346229504)
E=EllipticCurve(QQ,[0,A,0,B,0])
P=E([ZZ(7784100),ZZ(67822243920)])
print('SAGE_VERSION=',SAGE_VERSION,flush=True)
print('E=',E,flush=True)
print('P=',P,flush=True)

def attempt(name,fn):
    print('\n========== '+name+' ==========',flush=True)
    try:
        v=fn(); print('STATUS=SUCCESS',flush=True); print(repr(v),flush=True); return v
    except BaseException as exc:
        print('STATUS=FAILURE',flush=True); print(type(exc).__name__,str(exc),flush=True); traceback.print_exc(); return None

attempt('CONDUCTOR',lambda:E.conductor())
attempt('ROOT_NUMBER',lambda:E.root_number())
attempt('ANALYTIC_RANK_UPPER_BOUND',lambda:E.analytic_rank_upper_bound())
attempt('RIGOROUS_RANK_ONLY_USE_MWRANK_FALSE',lambda:E.rank(
    use_database=False,only_use_mwrank=False,algorithm='pari',proof=True,pari_effort=4))
attempt('RIGOROUS_SHA_ORDER_IF_RANK_LE1',lambda:E.sha().an(use_database=False,descent_second_limit=20))
attempt('BOUND_KOLYVAGIN',lambda:E.sha().bound_kolyvagin())
print('\nDONE',flush=True)
