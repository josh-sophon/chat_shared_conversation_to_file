from sage.all import *
from sage.env import SAGE_VERSION
proof.all(True)
A=ZZ(47643948); B=ZZ(159472346229504)
E=EllipticCurve(QQ,[0,A,0,B,0])
P1=E([ZZ(7784100),ZZ(67822243920)])
P2=E([QQ(-35341982703018681458668692158016)/QQ(5010410469853324398643225),QQ(335463976347290836626457566990041825275711400448)/QQ(11215275834124745849151273434962304125)])
print('SAGE_VERSION=',SAGE_VERSION,flush=True)
print('RANK=',E.rank(use_database=False,algorithm='pari',pari_effort=10,proof=True),flush=True)
print('GENS=',E.gens(use_database=False,algorithm='pari',pari_effort=10,proof=True),flush=True)
print('SATURATION_RAW=',E.saturation([P1,P2]),flush=True)
print('REG_RAW=',E.regulator_of_points([P1,P2]),flush=True)
print('REG_GENS=',E.regulator_of_points(E.gens(use_database=False,algorithm='pari',pari_effort=10,proof=True)),flush=True)
