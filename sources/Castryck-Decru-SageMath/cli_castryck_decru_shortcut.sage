from Crypto.Util.number import *
from public_values_aux import *
import hashlib
import sys

load("castryck_decru_shortcut.sage")


# fmt: off
(
    E_start0, E_start1, E_start2, E_start3, E_start4,
    a, b,
    EB0, EB0i, EB1, EB1i, EB2, EB2i, EB3, EB3i, EB4, EB4i,
    Pax0, Pax1, Pay0, Pay1,
    Pbx0, Pbx1, Pby0, Pby1,
    Qax0, Qax1, Qay0, Qay1,
    Qbx0, Qbx1, Qby0, Qby1,
    psiPax, psiPaxi, psiPay, psiPayi,
    psiQax, psiQaxi, psiQay, psiQayi,
) = map(int, sys.argv[1:])
# fmt: on

p = 2 ^ a * 3 ^ b - 1
assert p in Primes()

K, i = GF(p ^ 2, names="i", modulus=x ^ 2 + 1).objgen()

curve_start = (E_start0, E_start1, E_start2, E_start3, E_start4)

assert curve_start in ((0, 6, 0, 1, 0), (0, 0, 0, 1, 0))

EE = EllipticCurve(K, curve_start)


if curve_start == (0, 6, 0, 1, 0):
    two_i = generate_distortion_map(EE)

else:
    # y ^ 2 = x ^ 3 + x
    # https://github.com/pcw109550/write-up/tree/master/2022/RCTF/S2DH
    EE.set_order((p + 1) ** 2)
    PR, x = PolynomialRing(K, name="x").objgen()
    phi = EllipticCurveIsogeny(EE, x)
    E1728 = phi.codomain()
    for iota in E1728.automorphisms():
        P = E1728.random_point()
        if iota(iota(P)) == -P:
            two_i = phi.post_compose(iota).post_compose(phi.dual())
            break
    # def two_i(A):
    #     x0, y0 = A[0], A[1]
    #     x1 = -x0
    #     y1 = i * y0
    #     return 2 * E(x1, y1)


EEb = EllipticCurve(
    K, [EB0 + EB0i * i, EB1 + EB1i * i, EB2 + EB2i * i, EB3 + EB3i * i, EB4 + EB4i * i]
)


PPa = EE(Pax0 + Pax1 * i, Pay0 + Pay1 * i)
PPb = EE(Pbx0 + Pbx1 * i, Pby0 + Pby1 * i)
QQa = EE(Qax0 + Qax1 * i, Qay0 + Qay1 * i)
QQb = EE(Qbx0 + Qbx1 * i, Qby0 + Qby1 * i)

psiPPa = EEb(psiPax + psiPaxi * i, psiPay + psiPayi * i)
psiQQa = EEb(psiQax + psiQaxi * i, psiQay + psiQayi * i)

# P3 and Q3 are passed by global variables...
P3 = PPb
Q3 = QQb

print("Start CastryckDecruAttack")

try:
    result = CastryckDecruAttack(EE, PPa, QQa, EEb, psiPPa, psiQQa, two_i, a, b, num_cores=1)
except Exception as e:
    print("Failed", e)
    result = None
print("Result", result)
