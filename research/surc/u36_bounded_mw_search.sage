# Exact bounded Mordell-Weil search for the U=36 pair-recovery fiber.
from sage.all import *
from sage.env import SAGE_VERSION
import json, os

proof.all(True)
U=ZZ(36); E0=ZZ(928); A=ZZ(679385088); B=ZZ(21585334395469824)
E=EllipticCurve(QQ,[0,A,0,B,0])
G1=E([QQ(-624275993600)/1849, QQ(447435009251901440)/79507])
G2=E([ZZ(187142400),ZZ(5864068362240)])
T=E([0,0])
BOUND=500
PRIMES=[11,13,17,23,31,37,41,47,53,59,67,71,79,83,89,97]
OUT='research/surc/results-u36-bounded/u36_bounded_search.json'


def point_key(P):
    return ('O',) if P.is_zero() else (int(P[0]),int(P[1]))


def pair_compatible_ff(P,F):
    if P.is_zero(): return True
    e0=F(E0); uu=F(U); x=P[0]; y=P[1]
    X=x/e0**2; Y=y/e0**3; w=F((U-8)**2)/e0
    h=F(9)*(F(7)*w-1)**2
    c1=F(21)*w+3
    c2=F(9)*(F(7)*w-1)*(F(147)*w*w+61)
    den=F(64)*(X-h)
    if den==0: return True
    for yy in (Y,-Y):
        z=(yy-c1*(X-h)-c2)/den
        if (e0*z).is_square() and (e0*z+F(4)*uu).is_square(): return True
    return False


def allowed_for_prime(p):
    F=GF(p); Ep=EllipticCurve(F,[0,F(A),0,F(B),0])
    g1=Ep([F(G1[0]),F(G1[1])]); g2=Ep([F(G2[0]),F(G2[1])]); t=Ep([0,0])
    o1=ZZ(g1.order()); o2=ZZ(g2.order()); nm=lcm(o1,2)
    allowed=set()
    for n in range(nm):
        if n%2==0: continue
        for m in range(o2):
            for eps in (0,1):
                P=n*g1+m*g2+eps*t
                if pair_compatible_ff(P,F): allowed.add((int(n),int(m),int(eps)))
    return int(nm),int(o2),allowed

constraints=[]
for p in PRIMES:
    nm,mm,S=allowed_for_prime(p)
    constraints.append((p,nm,mm,S))
constraints.sort(key=lambda row: QQ(len(row[3]))/(row[1]*row[2]*2))

candidates=[]
for n in range(-BOUND,BOUND+1):
    if n%2==0: continue
    for m in range(-BOUND,BOUND+1):
        for eps in (0,1):
            if all((n%nm,m%mm,eps) in S for _,nm,mm,S in constraints):
                candidates.append((n,m,eps))


def inverse_z(P):
    if P.is_zero(): return []
    x,y=QQ(P[0]),QQ(P[1]); e0=QQ(E0); uu=QQ(U)
    X=x/e0**2; Y=y/e0**3; w=QQ((U-8)**2)/e0
    h=9*(7*w-1)**2; c1=21*w+3; c2=9*(7*w-1)*(147*w*w+61)
    den=64*(X-h)
    if den==0: return ['pole']
    return [(yy-c1*(X-h)-c2)/den for yy in (Y,-Y)]

hits=[]; exact_rows=[]
for n,m,eps in candidates:
    P=n*G1+m*G2+eps*T
    recovered=[]
    for z in inverse_z(P):
        if z=='pole':
            recovered.append({'kind':'inverse_pole'})
            continue
        d2=QQ(E0)*z; t2=d2+4*QQ(U)
        if d2>=0 and t2>=0 and d2.is_square() and t2.is_square():
            d=d2.sqrt(); tt=t2.sqrt()
            for ds in (d,-d):
                for ts in (tt,-tt):
                    c=(ts-ds)/2; e=(ts+ds)/2
                    if c*e==U:
                        recovered.append({'z':str(z),'c':str(c),'e':str(e),'target':bool(c>=4 and e>c)})
    row={'n':n,'m':m,'epsilon':eps,'point':'O' if P.is_zero() else [str(P[0]),str(P[1])],
         'x_numerator_bits':0 if P.is_zero() else int(abs(P[0].numerator()).nbits()),
         'x_denominator_bits':0 if P.is_zero() else int(P[0].denominator().nbits()),
         'recoveries':recovered}
    exact_rows.append(row)
    if any(r.get('target') for r in recovered): hits.append(row)

report={'status':'passed','sage_version':str(SAGE_VERSION),'U':36,'bound':BOUND,
        'basis':[[str(G1[0]),str(G1[1])],[str(G2[0]),str(G2[1])]],
        'primes':PRIMES,'constraint_summary':[{'p':p,'n_modulus':nm,'m_modulus':mm,'allowed_count':len(S)} for p,nm,mm,S in constraints],
        'modular_candidate_count':len(candidates),'modular_candidates':[list(x) for x in candidates],
        'exact_rows':exact_rows,'target_hit_count':len(hits),'target_hits':hits,
        'evidence_scope':'Exhaustive only for |n|,|m|<=500 in the displayed saturated basis and epsilon in {0,1}. It is not a global Mordell-Weil sieve.'}
os.makedirs(os.path.dirname(OUT),exist_ok=True)
open(OUT,'w').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'passed','candidate_count':len(candidates),'target_hit_count':len(hits)},sort_keys=True),flush=True)
