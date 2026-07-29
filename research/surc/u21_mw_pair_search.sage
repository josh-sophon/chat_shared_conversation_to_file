# Exact Mordell-Weil basis audit and bounded pair-recovery search for U=21.
from sage.all import *
from sage.env import SAGE_VERSION
import json, traceback

proof.all(True)
A=ZZ(47643948); B=ZZ(159472346229504); U=ZZ(21); E0=ZZ(253)
E=EllipticCurve(QQ,[0,A,0,B,0])
T=E([0,0])
raw1=E([ZZ(7784100),ZZ(67822243920)])
raw2=E([
 QQ(-35341982703018681458668692158016)/QQ(5010410469853324398643225),
 QQ(335463976347290836626457566990041825275711400448)/QQ(11215275834124745849151273434962304125)
])
w=QQ(169)/253
h=9*(7*w-1)^2
c1=21*w+3
c2=9*(7*w-1)*(147*w^2+61)

def squareclass(q):
    if q==0: return ZZ(0)
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

def zvals(P):
    if P.is_zero() or P[0]==0: return []
    X=P[0]/E0^2; Y=P[1]/E0^3
    den=64*(X-h)
    if den==0: return []
    return [(Y-c1*(X-h)-c2)/den,(-Y-c1*(X-h)-c2)/den]

def recover(z):
    V=E0*z; S2=V+4*U
    if V<=0 or S2<=0 or not V.is_square() or not S2.is_square(): return None
    d=V.sqrt(); t=S2.sqrt()
    for dd in [d,-d]:
      for tt in [t,-t]:
        c=(tt-dd)/2; e=(tt+dd)/2
        if c>e: c,e=e,c
        if c>=4 and c*e==U:
            return (c,e,V,S2)
    return None

print('SAGE_VERSION=',SAGE_VERSION,flush=True)
print('BEGIN_RANK',flush=True)
r=E.rank(use_database=False,algorithm='pari',pari_effort=10,proof=True)
print('RANK=',r,flush=True)
gens=E.gens(use_database=False,algorithm='pari',pari_effort=10,proof=True)
print('GENS=',gens,flush=True)
print('GEN_SQUARECLASSES=',[squareclass(P[0]) for P in gens],flush=True)
print('RAW_SATURATION=',E.saturation([raw1,raw2]),flush=True)
print('END_RANK',flush=True)

bound=ZZ(80)
mults=[]
for G in gens:
    table={0:E(0)}
    Q=E(0)
    for n in range(1,int(bound)+1):
        Q=Q+G; table[n]=Q; table[-n]=-Q
    mults.append(table)

stats={'points':0,'inverse_z':0,'target_interval':0,'V_square':0,'S_square':0,'both_square':0,'hits':[]}
seen=set()
for n in range(-int(bound),int(bound)+1):
  for m in range(-int(bound),int(bound)+1):
    P=mults[0][n]+mults[1][m]
    for eps in [0,1]:
      Q=P+(T if eps else E(0))
      key=(Q[0],Q[1]) if not Q.is_zero() else ('O',)
      if key in seen: continue
      seen.add(key); stats['points']+=1
      for z in zvals(Q):
        stats['inverse_z']+=1
        if 0<z<QQ(1)/16: stats['target_interval']+=1
        V=E0*z; S2=V+4*U
        vs=(V>=0 and V.is_square()); ss=(S2>=0 and S2.is_square())
        if vs: stats['V_square']+=1
        if ss: stats['S_square']+=1
        if vs and ss:
          stats['both_square']+=1
          rec=recover(z)
          if rec is not None:
            c,e,V0,S20=rec
            hit={'n':n,'m':m,'eps':eps,'point':str(Q),'z':str(z),'c':str(c),'e':str(e),'V':str(V0),'S2':str(S20),'alpha':str(squareclass(Q[0]))}
            stats['hits'].append(hit)
            print('HIT=',json.dumps(hit,sort_keys=True),flush=True)
print('SEARCH_BOUND=',bound,flush=True)
print('SEARCH_STATS=',json.dumps(stats,sort_keys=True),flush=True)
print('DONE',flush=True)
