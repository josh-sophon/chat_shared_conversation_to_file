# U=21 global descent harness for the Structured Unit-Rigidity research case.
# This script reports exact/proven computations separately from heuristic or
# analytic diagnostics. It deliberately catches failures so all completed
# evidence is retained in the log.

from sage.all import *
from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
import json
import traceback

proof.all(True)

A = ZZ(47643948)
B = ZZ(159472346229504)
C = A*A - 4*B
E = EllipticCurve(QQ, [0, A, 0, B, 0])
Ep = EllipticCurve(QQ, [0, -2*A, 0, C, 0])
P = E([ZZ(7784100), ZZ(67822243920)])
Pp = Ep([ZZ(14594052), ZZ(81434810160)])

def section(name):
    print("\n========== %s ==========" % name, flush=True)

def attempt(name, fn):
    section(name)
    try:
        value = fn()
        print("STATUS=SUCCESS", flush=True)
        print(repr(value), flush=True)
        return value
    except BaseException as exc:
        print("STATUS=FAILURE", flush=True)
        print("TYPE=%s" % type(exc).__name__, flush=True)
        print("MESSAGE=%s" % exc, flush=True)
        traceback.print_exc()
        return None

section("ENVIRONMENT")
print("SAGE_VERSION=%s" % SAGE_VERSION, flush=True)
print("A=%s" % A, flush=True)
print("B=%s" % B, flush=True)
print("C=%s" % C, flush=True)
print("E=%s" % E, flush=True)
print("EP=%s" % Ep, flush=True)
print("P=%s" % P, flush=True)
print("PP=%s" % Pp, flush=True)
print("P_ON_E=%s" % (P in E), flush=True)
print("PP_ON_EP=%s" % (Pp in Ep), flush=True)
print("E_MINIMAL_MODEL=%s" % E.global_minimal_model(), flush=True)
print("EP_MINIMAL_MODEL=%s" % Ep.global_minimal_model(), flush=True)
print("E_DISCRIMINANT=%s" % E.discriminant(), flush=True)
print("EP_DISCRIMINANT=%s" % Ep.discriminant(), flush=True)
print("E_CONDUCTOR=%s" % E.conductor(), flush=True)
print("EP_CONDUCTOR=%s" % Ep.conductor(), flush=True)
print("E_TORSION=%s" % E.torsion_subgroup(), flush=True)
print("EP_TORSION=%s" % Ep.torsion_subgroup(), flush=True)
print("E_ROOT_NUMBER=%s" % E.root_number(), flush=True)
print("EP_ROOT_NUMBER=%s" % Ep.root_number(), flush=True)

attempt("TWO_ISOGENY_SELMER_ONLY", lambda: two_descent_by_two_isogeny(
    E, global_limit_small=10, global_limit_large=1000,
    verbosity=2, selmer_only=True, proof=True))

# Raw PARI ellrank outputs contain proven upper bounds and Cassels-pairing
# information. Supply the known point so it is not rediscovered.
for effort in [0, 1, 2, 4, 6, 8, 10]:
    attempt("PARI_ELLRANK_E_EFFORT_%s" % effort,
            lambda effort=effort: E.pari_curve().ellrank(effort, [[P[0], P[1]]]))
    attempt("PARI_ELLRANK_EP_EFFORT_%s" % effort,
            lambda effort=effort: Ep.pari_curve().ellrank(effort, [[Pp[0], Pp[1]]]))

attempt("SAGE_PROVEN_RANK_PARI_EFFORT10", lambda: E.rank(
    use_database=False, algorithm='pari', pari_effort=10, proof=True))
attempt("SAGE_PROVEN_RANK_MWRANK", lambda: E.rank(
    use_database=False, algorithm='mwrank_lib', verbose=True, proof=True))
attempt("SAGE_GENS_PARI", lambda: E.gens(
    use_database=False, algorithm='pari', proof=False, pari_effort=10))
attempt("SAGE_GENS_MWRANK", lambda: E.gens(
    use_database=False, algorithm='mwrank_lib', proof=False,
    descent_second_limit=20, sat_bound=1000))

attempt("TWO_ISOGENY_GLOBAL_SEARCH_1E4", lambda: two_descent_by_two_isogeny(
    E, global_limit_small=100, global_limit_large=10000,
    verbosity=2, selmer_only=False, proof=True))
attempt("TWO_ISOGENY_GLOBAL_SEARCH_1E6", lambda: two_descent_by_two_isogeny(
    E, global_limit_small=1000, global_limit_large=1000000,
    verbosity=2, selmer_only=False, proof=True))

attempt("PARI_ELL2COVER_E", lambda: E.pari_curve().ell2cover())
attempt("PARI_ELL2COVER_EP", lambda: Ep.pari_curve().ell2cover())

# Analytic data are diagnostics only. They are explicitly not promoted to a
# theorem unless an independent algebraic proof closes the rank.
attempt("ANALYTIC_RANK_E", lambda: E.analytic_rank())
attempt("ANALYTIC_RANK_EP", lambda: Ep.analytic_rank())
attempt("ANALYTIC_SHA_E", lambda: E.sha().an_numerical(prec=100, use_database=False, proof=False))

section("DONE")
print("COMPLETED=1", flush=True)
