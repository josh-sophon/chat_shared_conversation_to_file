# Search one U=47 original-side two-isogeny quartic with PARI ratpoints.
from sage.all import *
from sage.env import SAGE_VERSION
import json, os, sys, traceback
proof.all(True)
d=ZZ(sys.argv[1]); H=ZZ(sys.argv[2]); OUT=sys.argv[3]
U=ZZ(47); E0=U*U-12*U+64; C0=23*U*U-304*U+1472; R0=U**4-100*U*U+4096
A=36*E0*C0; B=16128*E0*E0*R0
R.<x>=PolynomialRing(QQ)
Q=d*x^4+A*x^2+B//d
report={'sage_version':str(SAGE_VERSION),'U':47,'class':int(d),'height_bound':str(H),'quartic':str(Q),'status':'started'}
try:
    pts=pari(Q).hyperellratpoints(H)
    report['status']='success'; report['points']=repr(pts); report['point_count']=len(pts)
except BaseException as exc:
    report['status']='failure'; report['type']=type(exc).__name__; report['message']=str(exc); report['traceback']=traceback.format_exc()
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,'w') as fh: json.dump(report,fh,indent=2,sort_keys=True); fh.write('\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
