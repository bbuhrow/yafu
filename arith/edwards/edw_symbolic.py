"""
Phase 0, step 1: symbolic checks (sympy) of the formulas the Edwards plan
relies on (arith/avxecm_edwards_plan.md).  The paper formulas were extracted
from PDFs by machine, so each one is rechecked here as an identity.

  python edw_symbolic.py

Checks:
  1. Edwards -> Montgomery map: (x, y) on a x^2 + y^2 = 1 + d x^2 y^2 gives
     (u, v) = ((1+y)/(1-y), u/x) on B v^2 = u^3 + A u^2 + u with
     A = 2(a+d)/(a-d), B = 4/(a-d); and (A+2)/4 = a/(a-d).
  2. Gelin-Kleinjung-Lenstra Theorem 3.5, j = 1 and j = 4: (x_j, y_j) lies on
     -x^2 + y^2 = 1 - e_j^4 x^2 y^2 for every t.
  3. Corollary 3.6, family (iv): with (x, y) on y^2 = x^3 - x^2 - 9x + 9 and
     t = (4x + 4)/(y - 4), e_4^2 + 1 is the square given in the paper.
     Also family (i) (j = 1, y^2 = x^3 - 36x, t = (x+6)/(x-6)): e_1 = (3x/y)^2.
  4. The first curves: m = 1 gives e = 15/8 and (28/195, 3152/3495) (iv),
     e = 16/9 and (12/91, 27/29) (i).
"""
import sys
import sympy as sp

x, y, t, a, d, X, Y = sp.symbols('x y t a d X Y')


def is_zero(expr):
    return sp.simplify(sp.together(expr)) == 0


def check(name, ok):
    print('%-60s %s' % (name, 'OK' if ok else 'FAIL'))
    if not ok:
        check.fails += 1
check.fails = 0


# ---- 1. Edwards -> Montgomery ----
A = 2 * (a + d) / (a - d)
B = 4 / (a - d)
u = (1 + y) / (1 - y)
xsq = (1 - y ** 2) / (a - d * y ** 2)          # x^2 from the curve equation
# B v^2 = B u^2 / x^2
lhs = B * u ** 2 / xsq
rhs = u ** 3 + A * u ** 2 + u
check('Edwards (x,y) -> Montgomery (u, u/x) on B v^2 = u^3+Au^2+u', is_zero(lhs - rhs))
check('(A+2)/4 = a/(a-d)', is_zero((A + 2) / 4 - a / (a - d)))
check('a = -1, d = -e^4: (A+2)/4 = 1/(1-e^4)',
      is_zero(((A + 2) / 4).subs({a: -1, d: -sp.Symbol('e') ** 4})
              - 1 / (1 - sp.Symbol('e') ** 4)))

# ---- 2. Theorem 3.5, j = 1 and j = 4 ----
par = {
    1: (3 * (t ** 2 - 1) / (8 * t),
        128 * t ** 3 / (27 * t ** 6 + 63 * t ** 4 - 63 * t ** 2 - 27),
        (9 * t ** 4 - 2 * t ** 2 + 9) / (9 * t ** 4 - 9)),
    4: ((t ** 2 + 4 * t) / (t ** 2 - 4),
        (2 * t ** 3 + 2 * t ** 2 - 8 * t - 8) / (t ** 4 + 6 * t ** 3 + 12 * t ** 2 + 16 * t),
        (t ** 6 + 6 * t ** 5 + 10 * t ** 4 - 16 * t ** 3 - 48 * t ** 2 - 32 * t - 32)
        / (t ** 6 + 6 * t ** 5 + 10 * t ** 4 + 16 * t ** 3 + 48 * t ** 2 + 64 * t)),
}
for j, (e, x1, y1) in par.items():
    check('GKL Thm 3.5 j=%d: point on -x^2+y^2 = 1 - e^4 x^2 y^2' % j,
          is_zero(-x1 ** 2 + y1 ** 2 - 1 + e ** 4 * x1 ** 2 * y1 ** 2))


def mod_curve(expr, cubic):
    """numerator of expr reduced modulo y^2 = cubic(x); returns 0 iff identity"""
    num = sp.numer(sp.together(expr))
    num = sp.expand(num)
    r = sp.rem(sp.Poly(num, Y), sp.Poly(Y ** 2 - cubic, Y))
    return sp.simplify(r.as_expr())


# ---- 3. families ----
cub4 = X ** 3 - X ** 2 - 9 * X + 9
t4 = (4 * X + 4) / (Y - 4)
e4 = par[4][0].subs(t, t4)
sq4 = (X ** 4 + 4 * X ** 3 + 14 * X ** 2 - 108 * X + 153) / \
      (X ** 4 - 4 * X ** 3 - 18 * X ** 2 - 16 * X * Y + 12 * X + 48 * Y + 9)
check('GKL Cor 3.6 (iv): e4^2 + 1 = (paper square)^2 on E\'', mod_curve(e4 ** 2 + 1 - sq4 ** 2, cub4) == 0)

cub1 = X ** 3 - 36 * X
e1 = par[1][0].subs(t, (X + 6) / (X - 6))
check('GKL Cor 3.2 (i): e1 = (3x/y)^2 on y^2 = x^3 - 36x', mod_curve(e1 - (3 * X / Y) ** 2, cub1) == 0)

# ---- 4. first curves ----
def at(j, tv):
    return [sp.nsimplify(f.subs(t, tv)) for f in par[j]]

e, x1, y1 = at(4, sp.Rational(4 * 5 + 4, 8 - 4))       # Q = (5, 8): t = 6
check('(iv) m=1: e = 15/8, point (28/195, 3152/3495)',
      (e, x1, y1) == (sp.Rational(15, 8), sp.Rational(28, 195), sp.Rational(3152, 3495)))
check('(iv) m=1: x1 e != +-1 (non-torsion test)', abs(x1 * e) != 1)
e, x1, y1 = at(1, sp.Rational(-3 + 6, -3 - 6))           # Q = (-3, 9): t = -1/3
# t = -1/3 is excluded in Thm 3.1 (it is +-3^-1) -- the paper's first curve
# uses another sign/point; report what it is
print('    (i) m=1 via t=(x+6)/(x-6) = -1/3 is in the excluded set; '
      'paper lists e = 16/9, (12/91, 27/29)')
e16 = sp.Rational(16, 9)
pt = (sp.Rational(12, 91), sp.Rational(27, 29))
check('(i) paper curve: (12/91, 27/29) on -x^2+y^2 = 1 - (16/9)^4 x^2 y^2',
      -pt[0] ** 2 + pt[1] ** 2 - 1 + e16 ** 4 * pt[0] ** 2 * pt[1] ** 2 == 0)

print('%d failure(s)' % check.fails)
sys.exit(1 if check.fails else 0)
