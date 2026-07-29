from sage.all import *
from sage.env import SAGE_VERSION
from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
import json,os,sys
proof.all(True)
U=ZZ(sys.argv[1]); OUT=sys.argv[2]
E0=U*U-12*U+64; C0=23*U*U-304*U+1472; R0=U**4-100*U*U+4096
A=36*E0*C0; B=16128*E0*E0*R0
E=EllipticCurve(QQ,[0,A,0,B,0])
res=two_descent_by_two_isogeny(E,global_limit_small=10,global_limit_large=100,verbosity=2,selmer_only=True,proof=True)
report={'sage_version':str(SAGE_VERSION),'U':int(U),'result':[int(x) for x in res], 'meaning':'(n1,n2,n1_prime,n2_prime); with selmer_only, n2 and n2_prime are the two isogeny Selmer sizes'}
os.makedirs(os.path.dirname(OUT),exist_ok=True)
open(OUT,'w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
