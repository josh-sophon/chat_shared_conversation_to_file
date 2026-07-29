# Check one exact U=36 Mordell-Weil candidate selected by index.
from sage.all import *
import json, os, sys

U=ZZ(36); E0=ZZ(928); A=ZZ(679385088); B=ZZ(21585334395469824)
E=EllipticCurve(QQ,[0,A,0,B,0])
G1=E([QQ(-624275993600)/1849, QQ(447435009251901440)/79507])
G2=E([ZZ(187142400),ZZ(5864068362240)])
T=E([0,0])
CANDIDATES=[(5,62,0),(31,478,0),(45,125,0),(65,-126,0),(119,38,0),(135,-198,0),(141,-248,0),(169,81,0),(171,343,0),(189,-343,0),(189,161,0),(195,-88,0),(243,433,0),(305,-474,0),(321,94,0),(335,362,0),(425,-412,0),(471,-138,0),(477,-217,0)]
idx=int(sys.argv[1]); n,m,eps=CANDIDATES[idx]
OUT='research/surc/results-u36-candidates/candidate_%02d.json'%idx

def inverse_z(P):
    if P.is_zero(): return []
    x,y=QQ(P[0]),QQ(P[1]); e0=QQ(E0)
    X=x/e0**2; Y=y/e0**3; w=QQ((U-8)**2)/e0
    h=9*(7*w-1)**2; c1=21*w+3; c2=9*(7*w-1)*(147*w*w+61)
    den=64*(X-h)
    if den==0: return ['pole']
    return [(yy-c1*(X-h)-c2)/den for yy in (Y,-Y)]

P=n*G1+m*G2+eps*T
recoveries=[]
for z in inverse_z(P):
    if z=='pole':
        recoveries.append({'kind':'inverse_pole'}); continue
    d2=QQ(E0)*z; t2=d2+4*QQ(U)
    if d2>=0 and t2>=0 and d2.is_square() and t2.is_square():
        d=d2.sqrt(); tt=t2.sqrt()
        for ds in (d,-d):
            for ts in (tt,-tt):
                c=(ts-ds)/2; e=(ts+ds)/2
                if c*e==U:
                    recoveries.append({'z':str(z),'c':str(c),'e':str(e),'target':bool(c>=4 and e>c)})
report={'status':'passed','index':idx,'n':n,'m':m,'epsilon':eps,
        'point':[str(P[0]),str(P[1])],
        'x_numerator_bits':int(abs(P[0].numerator()).nbits()),
        'x_denominator_bits':int(P[0].denominator().nbits()),
        'recoveries':recoveries,'target_hit':bool(any(r.get('target') for r in recoveries)),
        'negation_scope':'The omitted negative coefficient tuple gives -P; inverse_z tests both Y signs, so it has the identical z-candidate set.'}
os.makedirs(os.path.dirname(OUT),exist_ok=True)
open(OUT,'w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'index':idx,'n':n,'m':m,'target_hit':report['target_hit'],'recoveries':len(recoveries)},sort_keys=True),flush=True)
