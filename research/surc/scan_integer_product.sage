# Exact low-height product-fiber scan for the Structured Unit-Rigidity case.
# Computes a proven Mordell-Weil basis when Sage/PARI can certify it, then
# searches for a conservative one-prime pair-recovery sieve.

from sage.all import *
from sage.env import SAGE_VERSION
import json, sys, traceback

proof.all(True)
U = ZZ(sys.argv[1])
OUT = sys.argv[2]

E0 = U*U - 12*U + 64
C0 = 23*U*U - 304*U + 1472
R0 = U**4 - 100*U*U + 4096
A = 36*E0*C0
B = 16128*E0*E0*R0
E = EllipticCurve(QQ, [0, A, 0, B, 0])
T = E([0,0])

report = {
    'status': 'started',
    'sage_version': str(SAGE_VERSION),
    'U': int(U), 'E0': int(E0), 'C0': int(C0), 'R0': int(R0),
    'A': int(A), 'B': int(B),
    'curve': str(E),
    'discriminant': str(E.discriminant()),
    'target_rule': 'A rational target pair has negative two-isogeny Kummer sign; finite-field inverse-map poles are retained conservatively.'
}

def point_json(P):
    if P.is_zero():
        return 'O'
    return [str(P[0]), str(P[1])]

def write_report():
    with open(OUT, 'w') as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
        fh.write('\n')

try:
    raw_rank = E.pari_curve().ellrank(10)
    report['pari_ellrank_effort10'] = repr(raw_rank)
except BaseException as exc:
    report['pari_ellrank_effort10_error'] = '%s: %s' % (type(exc).__name__, exc)

try:
    rank = E.rank(use_database=False, algorithm='pari', pari_effort=10, proof=True)
    gens = E.gens(use_database=False, algorithm='pari', pari_effort=10, proof=True)
    report['rank'] = int(rank)
    report['gens'] = [point_json(P) for P in gens]
    report['gens_count'] = len(gens)
    if len(gens) != rank:
        raise RuntimeError('generator count does not equal proven rank')
    if gens:
        sat = E.saturation(gens)
        report['saturation'] = repr(sat)
        report['regulator'] = str(E.regulator_of_points(gens))
        # Sage returns (basis, index, regulator); require exact index one.
        if len(sat) < 2 or ZZ(sat[1]) != 1:
            raise RuntimeError('returned generators are not certified saturated')
    else:
        report['saturation'] = 'rank zero: empty basis'
    report['torsion'] = str(E.torsion_subgroup())
except BaseException as exc:
    report['status'] = 'basis_not_certified'
    report['basis_error'] = '%s: %s' % (type(exc).__name__, exc)
    report['basis_traceback'] = traceback.format_exc()
    write_report()
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    sys.exit(0)

# The torsion Kummer class B is positive for U>16.  The sign of a point's
# Kummer class is therefore the parity sum of basis generators with x<0.
sign_bits = [1 if P[0] < 0 else 0 for P in gens]
report['generator_kummer_sign_bits'] = sign_bits
report['torsion_kummer_sign_bit'] = 0
if rank == 0 or not any(sign_bits):
    report['status'] = 'excluded_globally_positive_kummer_image'
    report['conclusion'] = 'All Mordell-Weil Kummer classes have positive real sign, whereas every target pair has negative sign.'
    write_report()
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    sys.exit(0)


def point_key(P):
    if P.is_zero():
        return ('O',)
    return (int(P[0]), int(P[1]))


def reduce_point(P, F, Ebar, p):
    if P.is_zero():
        return Ebar(0)
    x, y = QQ(P[0]), QQ(P[1])
    if x.denominator() % p == 0 or y.denominator() % p == 0:
        return None
    return Ebar([F(x), F(y)])


def pair_compatible(P, F):
    # Retain the point at infinity and every inverse-map pole.
    if P.is_zero():
        return True
    e0 = F(E0); uu = F(U)
    x, y = P[0], P[1]
    X = x / (e0**2)
    Y = y / (e0**3)
    w = F((U-8)**2) / e0
    h = F(9) * (F(7)*w - 1)**2
    c1 = F(21)*w + 3
    c2 = F(9) * (F(7)*w - 1) * (F(147)*w*w + 61)
    den = F(64) * (X-h)
    if den == 0:
        return True
    for yy in (Y, -Y):
        z = (yy - c1*(X-h) - c2) / den
        if (e0*z).is_square() and (e0*z + F(4)*uu).is_square():
            return True
    return False


def augmented_subgroup(Ebar, reduced_gens, bits, Tbar):
    generators = list(zip(reduced_gens, bits)) + [(Tbar, 0)]
    O = Ebar(0)
    states = {(point_key(O), 0): O}
    queue = [(O,0)]
    pos = 0
    while pos < len(queue):
        P, bit = queue[pos]; pos += 1
        for G, gbit in generators:
            Q = P + G
            state = (point_key(Q), bit ^^ gbit)
            if state not in states:
                states[state] = Q
                queue.append((Q, bit ^^ gbit))
    return states

prime_rows = []
excluding = None
for p in prime_range(5, 251):
    if E0 % p == 0 or E.discriminant() % p == 0:
        continue
    F = GF(p)
    Ebar = EllipticCurve(F, [0, F(A), 0, F(B), 0])
    red_gens = [reduce_point(P,F,Ebar,p) for P in gens]
    Tbar = reduce_point(T,F,Ebar,p)
    if any(P is None for P in red_gens) or Tbar is None:
        continue
    states = augmented_subgroup(Ebar, red_gens, sign_bits, Tbar)
    compatible_keys = {point_key(P) for P in Ebar.points() if pair_compatible(P,F)}
    negative_states = []
    for (k, bit), P in states.items():
        if bit == 1 and k in compatible_keys:
            negative_states.append(k)
    row = {
        'p': int(p),
        'curve_order': int(Ebar.cardinality()),
        'generator_orders': [int(P.order()) for P in red_gens],
        'torsion_order': int(Tbar.order()),
        'augmented_subgroup_size': len(states),
        'compatible_point_count': len(compatible_keys),
        'negative_compatible_state_count': len(negative_states),
        'negative_compatible_points': [list(k) if k != ('O',) else 'O' for k in sorted(set(negative_states), key=str)[:20]],
    }
    prime_rows.append(row)
    if len(negative_states) == 0:
        excluding = row
        break

report['prime_sieve_rows'] = prime_rows
if excluding is not None:
    report['status'] = 'excluded_by_single_prime_mw_sieve'
    report['excluding_prime'] = excluding['p']
    report['conclusion'] = 'No negative-Kummer Mordell-Weil state survives the conservative pair-recovery test at the displayed good prime.'
else:
    report['status'] = 'survives_primes_below_251'
    report['conclusion'] = 'At least one negative-Kummer state survives every audited good prime below 251; this is reconnaissance, not a rational point.'

write_report()
print(json.dumps(report, indent=2, sort_keys=True), flush=True)
