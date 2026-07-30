# Exact U=56 pair-gate rank/Kummer audit for Structured Unit-Rigidity.
from sage.all import *
from sage.env import SAGE_VERSION
import json, traceback
proof.all(True)
U=ZZ(56)
E0=U*U-12*U+64
C0=23*U*U-304*U+1472
R0=U**4-100*U*U+4096
A=36*E0*C0
B=16128*E0*E0*R0
E=EllipticCurve(QQ,[0,A,0,B,0])
report={'sage_version':str(SAGE_VERSION),'U':int(U),'E0':int(E0),'C0':int(C0),'R0':int(R0),'A':int(A),'B':int(B),'E':str(E),'root_number':int(E.root_number()),'torsion':str(E.torsion_subgroup())}

def attempt(k,fn):
    try:
        v=fn(); report[k]=repr(v); print(k,repr(v),flush=True); return v
    except BaseException as exc:
        report[k+'_error']='%s: %s'%(type(exc).__name__,exc); report[k+'_traceback']=traceback.format_exc(); print(k+'_ERROR',report[k+'_error'],flush=True); return None

for effort in [0,1,2,4,8,12]:
    r=attempt('ellrank_effort_%s'%effort,lambda effort=effort:E.pari_curve().ellrank(effort))
    if r is not None and ZZ(r[0])==ZZ(r[1]):
        report['exact_effort']=effort; break
attempt('proven_rank',lambda:E.rank(use_database=False,algorithm='pari',pari_effort=12,proof=True))
gens=attempt('gens',lambda:E.gens(use_database=False,algorithm='pari',pari_effort=12,proof=True))
if gens is not None:
    report['gens_json']=[[str(P[0]),str(P[1]),int(sign(P[0]))] for P in gens]
    report['saturation']=repr(E.saturation(gens))
    report['regulator']=str(E.regulator_of_points(gens))
report['B_squarefree_part']=str(B.squarefree_part())
report['status']='completed'
open('research/surc/u56_rank_scan.json','w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
