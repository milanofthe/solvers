"""Construct the collocation families from their nodes.

A Gauss, Radau or Lobatto method is not a table somebody wrote down, it is what
follows from a choice of nodes. Transcribing the tableaux of the higher stage
members would be error prone and would also hide that; building them from the
node polynomials instead makes the construction the source of truth, and the
Rust side then verifies the order that construction is supposed to give.

  nodes           roots of a node polynomial on [0, 1]
  A               sum_j a_ij c_j^(k-1) = c_i^k / k,  k = 1..s   (collocation)
  b               sum_j b_j c_j^(k-1) = 1 / k,       k = 1..s

All three node polynomials are the same object differentiated a different number
of times,

  Gauss       d^s     / dx^s     [ x^s     (x-1)^s     ]
  Radau IIA   d^(s-1) / dx^(s-1) [ x^(s-1) (x-1)^s     ]
  Lobatto     d^(s-2) / dx^(s-2) [ x^(s-1) (x-1)^(s-1) ]

which is built here with integer coefficients and then rooted in fifty digits.
That precision is not decoration. The linear systems that turn nodes into a
tableau are Vandermonde systems, and at eight stages a double precision solve
loses six digits: the coefficients would no longer satisfy the order conditions
they were derived from, and the analysis would correctly report a method of
lower order than the construction gives.

Lobatto IIIB and IIIC are not collocation methods and are built from their own
simplifying assumptions instead, which is noted where it happens.
"""

import io
import json
import os
import re
from fractions import Fraction

from mpmath import mp, mpf, matrix, lu_solve, polyroots

mp.dps = 50

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "methods")
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

BUTCHER64 = {
    "authors": "Butcher, J. C.",
    "title": "Implicit Runge-Kutta processes",
    "year": 1964,
    "source": "Mathematics of Computation, 18(85), 50-64",
    "doi": "10.1090/S0025-5718-1964-0159424-9",
}
EHLE69 = {
    "authors": "Ehle, B. L.",
    "title": "On Pade approximations to the exponential function and A-stable methods for the numerical solution of initial value problems",
    "year": 1969,
    "source": "Research Report CSRR 2010, University of Waterloo",
}
HW2 = {
    "authors": "Hairer, E., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
    "year": 1996,
    "source": "Springer Series in Computational Mathematics, Vol. 14, 2nd edition",
    "doi": "10.1007/978-3-642-05221-7",
}
BUTCHER16 = {
    "authors": "Butcher, J. C.",
    "title": "Numerical Methods for Ordinary Differential Equations",
    "year": 2016,
    "source": "Wiley, 3rd edition",
    "doi": "10.1002/9781119121534",
}


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def write(directory, method):
    folder = os.path.join(OUT, directory)
    os.makedirs(folder, exist_ok=True)
    with io.open(os.path.join(folder, method["id"] + ".json"), "w", encoding="utf-8") as handle:
        handle.write(render(method))
    return method["id"]


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------


def node_polynomial(left, right, derivatives):
    """`x^left (x-1)^right` differentiated `derivatives` times.

    Integer coefficients throughout, lowest power first, so nothing is lost
    before the roots are taken.
    """
    poly = [0] * (left + right + 1)
    for k in range(right + 1):
        sign = -1 if (right - k) % 2 else 1
        binomial = Fraction(1)
        for i in range(k):
            binomial = binomial * (right - i) / (i + 1)
        poly[left + k] = sign * int(binomial)
    for _ in range(derivatives):
        poly = [poly[i] * i for i in range(1, len(poly))]
    return poly


def roots_on_unit_interval(poly):
    """The real roots, sorted, at working precision."""
    highest_first = list(reversed(poly))
    while highest_first and highest_first[0] == 0:
        highest_first.pop(0)
    found = polyroots([mpf(c) for c in highest_first], maxsteps=200, extraprec=400)
    return sorted(mp.re(root) for root in found)


def gauss_nodes(s):
    """Gauss-Legendre points: the shifted Legendre polynomial by Rodrigues."""
    return roots_on_unit_interval(node_polynomial(s, s, s))


def radau_iia_nodes(s):
    """Radau IIA points, which include the right endpoint and leave the left free."""
    return roots_on_unit_interval(node_polynomial(s - 1, s, s - 1))


def lobatto_nodes(s):
    """Lobatto points, which pin both endpoints."""
    return roots_on_unit_interval(node_polynomial(s - 1, s - 1, s - 2))


# ---------------------------------------------------------------------------
# Tableaux
# ---------------------------------------------------------------------------


def vandermonde(nodes):
    """Rows `c_j^(k-1)`, the system every one of these solves against."""
    s = len(nodes)
    return matrix([[nodes[j] ** k for j in range(s)] for k in range(s)])


def quadrature_weights(c):
    """The weights that integrate a polynomial of degree s - 1 exactly."""
    s = len(c)
    return lu_solve(vandermonde(c), matrix([mpf(1) / (k + 1) for k in range(s)]))


def collocation_matrix(c):
    """`A` from the collocation conditions, one row at a time."""
    s = len(c)
    system = vandermonde(c)
    rows = []
    for i in range(s):
        targets = matrix([c[i] ** (k + 1) / (k + 1) for k in range(s)])
        rows.append(list(lu_solve(system, targets)))
    return rows


def lobatto_iiic_matrix(c, b):
    """Lobatto IIIC pins the first column to `b_1` and asks the rest to satisfy
    the collocation conditions one degree short. That first column is what buys
    the damping at infinity that IIIA lacks."""
    s = len(c)
    rest = c[1:]
    system = matrix([[rest[j] ** k for j in range(s - 1)] for k in range(s - 1)])
    rows = []
    for i in range(s):
        targets = matrix(
            [
                c[i] ** (k + 1) / (k + 1) - b[0] * (c[0] ** k if k > 0 else mpf(1))
                for k in range(s - 1)
            ]
        )
        rows.append([b[0]] + list(lu_solve(system, targets)))
    return rows


def lobatto_iiib_matrix(c, b):
    """Lobatto IIIB is defined by the adjoint conditions, one column at a time:
    sum_i b_i c_i^(k-1) a_ij = b_j (1 - c_j^k) / k."""
    s = len(c)
    system = matrix([[b[i] * c[i] ** k for i in range(s)] for k in range(s)])
    columns = []
    for j in range(s):
        targets = matrix([b[j] * (1 - c[j] ** (k + 1)) / (k + 1) for k in range(s)])
        columns.append(list(lu_solve(system, targets)))
    return [[columns[j][i] for j in range(s)] for i in range(s)]


def tableau(a, b, c):
    return {
        "a": [[round(float(v), 16) for v in row] for row in a],
        "b": [round(float(v), 16) for v in b],
        "c": [round(float(v), 16) for v in c],
    }


def method(identifier, name, order, a, b, c, properties, references, aliases=None):
    entry = {
        "id": identifier,
        "name": name,
        "class": "runge_kutta",
        "family": "irk",
        "order": order,
    }
    if aliases:
        entry["aliases"] = aliases
    entry["properties"] = properties
    entry["tableau"] = tableau(a, b, c)
    entry["references"] = references
    return entry


# ---------------------------------------------------------------------------
# The families
# ---------------------------------------------------------------------------

written = []

# Gauss-Legendre: the maximum order a Runge-Kutta method of s stages can reach,
# and symplectic with it.
for s in range(2, 9):
    c = gauss_nodes(s)
    b = quadrature_weights(c)
    a = collocation_matrix(c)
    written.append(
        write(
            "irk",
            method(
                f"gauss_legendre_{2 * s}",
                f"Gauss-Legendre {2 * s}",
                2 * s,
                a,
                b,
                c,
                {
                    "a_stable": True,
                    "l_stable": False,
                    "symplectic": True,
                    "symmetric": True,
                    "stage_order": s,
                },
                [BUTCHER64, BUTCHER16, HW2],
                aliases=[f"gauss{2 * s}"],
            ),
        )
    )

# Radau IIA: one order below Gauss, but L-stable and stiffly accurate, which is
# what matters on a stiff problem.
for s in range(2, 8):
    c = radau_iia_nodes(s)
    b = quadrature_weights(c)
    a = collocation_matrix(c)
    written.append(
        write(
            "irk",
            method(
                f"radau_iia_{2 * s - 1}",
                f"Radau IIA {2 * s - 1}",
                2 * s - 1,
                a,
                b,
                c,
                {
                    "a_stable": True,
                    "l_stable": True,
                    "stiffly_accurate": True,
                    "stage_order": s,
                },
                [BUTCHER64, EHLE69, HW2],
                aliases=[f"radau{2 * s - 1}"],
            ),
        )
    )

# Lobatto, all three families on the same nodes.
for s in range(2, 7):
    c = lobatto_nodes(s)
    b = quadrature_weights(c)
    if s >= 3:
        written.append(
            write(
                "irk",
                method(
                    f"lobatto_iiia_{2 * s - 2}",
                    f"Lobatto IIIA {2 * s - 2}",
                    2 * s - 2,
                    collocation_matrix(c),
                    b,
                    c,
                    {
                        "a_stable": True,
                        "l_stable": False,
                        "stiffly_accurate": True,
                        "symmetric": True,
                        "stage_order": s,
                    },
                    [BUTCHER64, EHLE69, HW2],
                ),
            )
        )
        written.append(
            write(
                "irk",
                method(
                    f"lobatto_iiib_{2 * s - 2}",
                    f"Lobatto IIIB {2 * s - 2}",
                    2 * s - 2,
                    lobatto_iiib_matrix(c, b),
                    b,
                    c,
                    {"a_stable": True, "l_stable": False, "symmetric": True},
                    [BUTCHER64, EHLE69, HW2],
                ),
            )
        )
    written.append(
        write(
            "irk",
            method(
                f"lobatto_iiic_{2 * s - 2}",
                f"Lobatto IIIC {2 * s - 2}",
                2 * s - 2,
                lobatto_iiic_matrix(c, b),
                b,
                c,
                {
                    "a_stable": True,
                    "l_stable": True,
                    "stiffly_accurate": True,
                    "stage_order": s - 1,
                },
                [BUTCHER64, EHLE69, HW2],
            ),
        )
    )

print(f"{len(written)} collocation methods written")
for identifier in written:
    print(" ", identifier)
