# Targeted missing-generator search for U=43,47,50.
from sage.all import *
from sage.env import SAGE_VERSION
from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
import json, os, sys, traceback
proof.all(True)
U=ZZ(sys.argv[1]); MODE=sys.argv[2]; PARAM=ZZ(sys.argv[3]); OUT=sys.argv[4]
E0=U*U-12*U+64; C0=23*U*U-304*U+1472; R0=U**4-100*U*U+4096
A=36*E0*C0; B=16128*E0*E0*R0
E=EllipticCurve(QQ,[0,A,0,B,0])
x0=36*(U-6)**2*(3*U-32)**2
y0=144*(U-6)*(3*U-32)*(26*U**4-771*U**3+9130*U**2-49344*U+106496)
P0=E([x0,y0])
report={'sage_version':str(SAGE_VERSION),'U':int(U),'mode':MODE,'parameter':int(PARAM),'A':str(A),'B':str(B),'curve':str(E),'known_point':[str(P0[0]),str(P0[1])],'status':'started'}
def pj(P): return 'O' if P.is_zero() else [str(P[0]),str(P[1])]
def attempt(name,fn):
    try:
        v=fn(); report[name]={'status':'success','repr':repr(v)}; return v
    except BaseException as exc:
        report[name]={'status':'failure','type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}; return None
raw=attempt('pari_ellrank',lambda:E.pari_curve().ellrank(int(PARAM),[[P0[0],P0[1]]]))
gens=None
if MODE=='pari':
    gens=attempt('sage_gens_pari',lambda:E.gens(use_database=False,algorithm='pari',proof=False,pari_effort=int(PARAM)))
elif MODE=='mwrank':
    gens=attempt('sage_gens_mwrank',lambda:E.gens(use_database=False,algorithm='mwrank_lib',proof=False,verbose=True,descent_second_limit=int(PARAM),sat_bound=2000))
elif MODE=='twoiso':
    attempt('two_isogeny_descent',lambda:two_descent_by_two_isogeny(E,global_limit_small=1000,global_limit_large=int(PARAM),verbosity=3,selmer_only=False,proof=True))
    gens=attempt('sage_gens_after_twoiso',lambda:E.gens(use_database=False,algorithm='mwrank_lib',proof=False,verbose=True,descent_second_limit=30,sat_bound=2000))
else: raise ValueError('unknown mode')
if gens is not None:
    report['gens']=[pj(P) for P in gens]; report['gens_count']=len(gens)
    if len(gens)>=2:
        report['height_pairing_det']=str(E.regulator_of_points(gens[:2]))
        sat=attempt('saturation',lambda:E.saturation(gens))
        if sat is not None:
            report['saturation_basis']=[pj(P) for P in sat[0]]; report['saturation_index']=str(sat[1]); report['saturation_regulator']=str(sat[2]); report['status']='two_generators_found'
    else: report['status']='one_generator_only'
else: report['status']='no_generator_result'
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,'w') as fh: json.dump(report,fh,indent=2,sort_keys=True); fh.write('\n')
print(json.dumps(report,indent=2,sort_keys=True),flush=True)
