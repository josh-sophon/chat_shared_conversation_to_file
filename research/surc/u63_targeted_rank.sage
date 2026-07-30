# Targeted exact U=63 rank search. PARI's first ellrank entry is conjectural;
# theorem-level rank requires saturated points equal to the proven upper bound.
from sage.all import *
from sage.env import SAGE_VERSION
import json, traceback
proof.all(True)
U=ZZ(63)
E0=U*U-12*U+64
C0=23*U*U-304*U+1472
R0=U**4-100*U*U+4096
A=36*E0*C0
B=16128*E0*E0*R0
E=EllipticCurve(QQ,[0,A,0,B,0])
P=E([ZZ(2883045636),ZZ(322196916140208)])
report={'sage_version':str(SAGE_VERSION),'U':63,'E':str(E),'P':[str(P[0]),str(P[1])],'P_on_E':bool(P in E),'root_number':int(E.root_number()),'B_squarefree_part':str(B.squarefree_part())}

def attempt(k,fn):
    try:
        v=fn(); report[k]=repr(v); print(k,repr(v),flush=True); return v
    except BaseException as exc:
        report[k+'_error']='%s: %s'%(type(exc).__name__,exc); report[k+'_traceback']=traceback.format_exc(); print(k+'_ERROR',report[k+'_error'],flush=True); return None

for effort in [0,1,2,4,8,10,12,16,20]:
    attempt('ellrank_seeded_%s'%effort,lambda effort=effort:E.pari_curve().ellrank(effort,[[P[0],P[1]]]))
rank=attempt('proven_rank_pari20',lambda:E.rank(use_database=False,algorithm='pari',pari_effort=20,proof=True))
gens=attempt('gens_pari20',lambda:E.gens(use_database=False,algorithm='pari',pari_effort=20,proof=True))
if rank is not None: report['proven_rank_int']=int(rank)
if gens is not None:
    report['gens_json']=[[str(Q[0]),str(Q[1]),int(sign(Q[0])),bool(QQ(Q[0]).is_square())] for Q in gens]
    report['saturation']=repr(E.saturation(gens))
    report['regulator']=str(E.regulator_of_points(gens))
attempt('rank_mwrank',lambda:E.rank(use_database=False,algorithm='mwrank_lib',proof=True,verbose=True))
report['status']='completed'
open('research/surc/u63_targeted_rank.json','w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
