# Exact adversarial audit for the U=53 rank/root-number correction.
from sage.all import *
from sage.env import SAGE_VERSION
import json
proof.all(True)
U=ZZ(53)
E0=U*U-12*U+64
C0=23*U*U-304*U+1472
R0=U**4-100*U*U+4096
A=36*E0*C0
B=16128*E0*E0*R0
E=EllipticCurve(QQ,[0,A,0,B,0])
P=E([ZZ(1282642596),ZZ(97562298101328)])
M=E.global_minimal_model()
bad=list(E.conductor().prime_divisors())
rows=[]
for p in bad:
    ld=E.local_data(p)
    row={
      'p':int(p),
      'root_number':int(E.root_number(p)),
      'kodaira':str(ld.kodaira_symbol()),
      'conductor_valuation':int(ld.conductor_valuation()),
      'discriminant_valuation':int(ld.discriminant_valuation()),
      'tamagawa':int(ld.tamagawa_number()),
      'reduction_type':str(ld.reduction_type()),
    }
    rows.append(row)
local_product=prod(ZZ(r['root_number']) for r in rows)
infinity=int(E.root_number(0))
global_direct=int(E.root_number())
assert global_direct == infinity*local_product
pari_rank=E.pari_curve().ellrank(0,[[P[0],P[1]]])
rank=E.rank(use_database=False,algorithm='pari',pari_effort=0,proof=True)
gens=E.gens(use_database=False,algorithm='pari',pari_effort=0,proof=True)
sat=E.saturation(gens)
report={
 'sage_version':str(SAGE_VERSION),
 'E':str(E),'minimal_model':str(M),'minimal_ainvs':[str(x) for x in M.ainvs()],
 'conductor':str(E.conductor()),'bad_primes':[int(p) for p in bad],
 'infinity_root_number':infinity,'finite_root_product':int(local_product),
 'global_root_number_direct':global_direct,'global_root_number_reconstructed':int(infinity*local_product),
 'local_rows':rows,
 'pari_ellrank_effort0':repr(pari_rank),'proven_rank':int(rank),
 'gens':[[str(Q[0]),str(Q[1])] for Q in gens],
 'saturation':repr(sat),'regulator':str(E.regulator_of_points(gens)),
 'known_point_on_curve':bool(P in E),
 'status':'passed'
}
open('research/surc/u53_root_rank_audit.json','w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
