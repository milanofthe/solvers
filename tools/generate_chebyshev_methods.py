"""Construct the Runge-Kutta-Chebyshev families.

These are the one family here whose point is not accuracy. An explicit method of
`s` stages can reach at most `-2 s^2` along the negative real axis, and the
Chebyshev polynomials are what attains it: the first order family is stable out
to about `1.9 s^2` and the second order one to about `0.65 s^2`, against `-2.78`
for classical Runge-Kutta whatever its cost. On a parabolic problem, where the
eigenvalues are real and spread out and the accuracy wanted is modest, that is
the whole game.

They are also not tableaux anybody writes down. The method is a three term
recursion whose coefficients come from the Chebyshev polynomials evaluated just
outside their interval, which is what keeps its storage independent of the
stage count:

    Y_0 = y
    Y_1 = Y_0 + mu~_1 h f(Y_0)
    Y_j = (1 - mu_j - nu_j) Y_0 + mu_j Y_{j-1} + nu_j Y_{j-2}
          + mu~_j h f(Y_{j-1}) + gamma~_j h f(Y_0)

The recursion is unrolled into a Butcher tableau here, which is the honest
translation and costs the storage the recursion was designed to save. What it
buys is that every piece of analysis in this project then applies unchanged.

The damping `eps` moves the polynomial off the interval endpoints so that the
stability boundary is an interval with a margin rather than a chain of points
touching one. Two thirteenths is the usual choice and the one the RKC solver
ships with.

References
----------
* P. J. van der Houwen, B. P. Sommeijer, "On the internal stability of explicit
  m-stage Runge-Kutta methods for large m-values", ZAMM 60, 1980,
  doi:10.1002/zamm.19800601005
* B. P. Sommeijer, L. F. Shampine, J. G. Verwer, "RKC: an explicit solver for
  parabolic PDEs", J. Comput. Appl. Math. 88, 1998,
  doi:10.1016/S0377-0427(97)00219-7
"""

import io
import json
import os
import re

from mpmath import mp, mpf

mp.dps = 60

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "methods", "rkc")
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

DAMPING = mpf(2) / 13

HOUWEN80 = {
    "authors": "van der Houwen, P. J., & Sommeijer, B. P.",
    "title": "On the internal stability of explicit m-stage Runge-Kutta methods for large m-values",
    "year": 1980,
    "source": "ZAMM - Journal of Applied Mathematics and Mechanics, 60(10), 479-485",
    "doi": "10.1002/zamm.19800601005",
}
SOMMEIJER98 = {
    "authors": "Sommeijer, B. P., Shampine, L. F., & Verwer, J. G.",
    "title": "RKC: an explicit solver for parabolic PDEs",
    "year": 1998,
    "source": "Journal of Computational and Applied Mathematics, 88(2), 315-326",
    "doi": "10.1016/S0377-0427(97)00219-7",
}


def chebyshev(s, x):
    """`T_j(x)` and its first two derivatives, for `j = 0..s`."""
    t = [mpf(1), x]
    d1 = [mpf(0), mpf(1)]
    d2 = [mpf(0), mpf(0)]
    for j in range(2, s + 1):
        t.append(2 * x * t[j - 1] - t[j - 2])
        d1.append(2 * t[j - 1] + 2 * x * d1[j - 1] - d1[j - 2])
        d2.append(4 * d1[j - 1] + 2 * x * d2[j - 1] - d2[j - 2])
    return t, d1, d2


def unroll(s, mu, nu, mu_tilde, gamma_tilde):
    """The three term recursion as a Butcher tableau.

    `rows[j][k]` is what stage value `Y_j` owes to `h f(Y_k)`. The method has
    `s` stages, the last row is the propagating solution.
    """
    rows = [[mpf(0)] * s for _ in range(s + 1)]
    rows[1][0] = mu_tilde[1]
    for j in range(2, s + 1):
        for k in range(s):
            rows[j][k] = mu[j] * rows[j - 1][k] + nu[j] * rows[j - 2][k]
        rows[j][j - 1] += mu_tilde[j]
        rows[j][0] += gamma_tilde[j]
    return rows


def first_order(s):
    """`R(z) = T_s(w0 + w1 z) / T_s(w0)`, which is consistent and no more."""
    w0 = 1 + DAMPING / s**2
    t, d1, _ = chebyshev(s, w0)
    w1 = t[s] / d1[s]

    mu = [mpf(0)] * (s + 1)
    nu = [mpf(0)] * (s + 1)
    mu_tilde = [mpf(0)] * (s + 1)
    gamma_tilde = [mpf(0)] * (s + 1)
    mu_tilde[1] = w1 / w0
    for j in range(2, s + 1):
        mu[j] = 2 * w0 * t[j - 1] / t[j]
        nu[j] = -t[j - 2] / t[j]
        mu_tilde[j] = 2 * w1 * t[j - 1] / t[j]
    return unroll(s, mu, nu, mu_tilde, gamma_tilde)


def second_order(s):
    """The second order member, which needs the `a_j` correction to the past."""
    w0 = 1 + DAMPING / s**2
    t, d1, d2 = chebyshev(s, w0)
    w1 = d1[s] / d2[s]

    b = [mpf(0)] * (s + 1)
    for j in range(2, s + 1):
        b[j] = d2[j] / (d1[j] ** 2)
    b[0] = b[2]
    b[1] = b[2]
    a = [1 - b[j] * t[j] for j in range(s + 1)]

    mu = [mpf(0)] * (s + 1)
    nu = [mpf(0)] * (s + 1)
    mu_tilde = [mpf(0)] * (s + 1)
    gamma_tilde = [mpf(0)] * (s + 1)
    mu_tilde[1] = b[1] * w1
    for j in range(2, s + 1):
        mu[j] = 2 * b[j] * w0 / b[j - 1]
        nu[j] = -b[j] / b[j - 2]
        mu_tilde[j] = 2 * b[j] * w1 / b[j - 1]
        gamma_tilde[j] = -a[j - 1] * mu_tilde[j]
    return unroll(s, mu, nu, mu_tilde, gamma_tilde)


def order_of(rows, s):
    """The quadrature order of the unrolled tableau, as a check on the algebra."""
    weights = rows[s]
    nodes = [sum(rows[i]) for i in range(s)]
    reached = 0
    for k in range(1, 6):
        total = sum(weights[i] * nodes[i] ** (k - 1) for i in range(s))
        if abs(total - mpf(1) / k) > mpf(10) ** -30:
            break
        reached = k
    return reached


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def write(method):
    os.makedirs(OUT, exist_ok=True)
    with io.open(os.path.join(OUT, method["id"] + ".json"), "w", encoding="utf-8") as handle:
        handle.write(render(method))


def main():
    written = []
    for order, build in ((1, first_order), (2, second_order)):
        for s in (5, 10, 20):
            rows = build(s)
            reached = order_of(rows, s)
            if reached < order:
                raise SystemExit(f"rkc{order} with {s} stages integrates only to order {reached}")

            method = {
                "id": f"rkc{order}_{s}",
                "name": f"Runge-Kutta-Chebyshev {order}, {s} stages",
                "class": "runge_kutta",
                "family": "rkc",
                "order": order,
                "properties": {"a_stable": False, "l_stable": False},
                "tableau": {
                    "a": [[round(float(v), 16) for v in rows[i][:i]] for i in range(s)],
                    "b": [round(float(v), 16) for v in rows[s]],
                },
                "references": [HOUWEN80, SOMMEIJER98],
            }
            write(method)
            written.append((method["id"], s, reached))

    print(f"{len(written)} Chebyshev method files written")
    for identifier, s, reached in written:
        print(f"  {identifier:<12} {s:>2} stages   quadrature order {reached}")


if __name__ == "__main__":
    main()
