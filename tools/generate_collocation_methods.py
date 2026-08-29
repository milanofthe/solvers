"""Construct the collocation families from their nodes.

A Gauss, Radau or Lobatto method is not a table somebody wrote down, it is what
follows from a choice of nodes. Transcribing the tableaux of the higher stage
members would be error prone and would also hide that; building them from the
node polynomials instead makes the construction the source of truth, and the
Rust side then verifies the order that construction is supposed to give.

  nodes           roots of a node polynomial on [0, 1]
  A               sum_j a_ij c_j^(k-1) = c_i^k / k,  k = 1..s   (collocation)
  b               sum_j b_j c_j^(k-1) = 1 / k,       k = 1..s

Lobatto IIIB and IIIC are not collocation methods and are built from their own
simplifying assumptions instead, which is noted where it happens.
"""

import io
import json
import os
import re

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "methods")
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

BUTCHER64 = {
    "authors": "Butcher, J. C.",
    "title": "Implicit Runge-Kutta processes",
    "year": 1964,
    "source": "Mathematics of Computation, 18(85), 50-64",
    "doi": "10.1090/S0025-5718-1964-0159424-9",
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


def gauss_nodes(s):
    """Gauss-Legendre points on [0, 1]: the roots of the shifted Legendre polynomial."""
    nodes, _ = np.polynomial.legendre.leggauss(s)
    return np.sort((nodes + 1.0) / 2.0)


def radau_iia_nodes(s):
    """Radau IIA points: the roots of d^(s-1)/dx^(s-1) [ x^(s-1) (x-1)^s ], which
    includes the right endpoint and leaves the left one free."""
    poly = np.polynomial.polynomial.polymul(
        np.polynomial.polynomial.polypow([0, 1], s - 1),
        np.polynomial.polynomial.polypow([-1, 1], s),
    )
    for _ in range(s - 1):
        poly = np.polynomial.polynomial.polyder(poly)
    return np.sort(np.real(np.polynomial.polynomial.polyroots(poly)))


def lobatto_nodes(s):
    """Lobatto points: the roots of d^(s-2)/dx^(s-2) [ x^(s-1) (x-1)^(s-1) ],
    which pins both endpoints."""
    poly = np.polynomial.polynomial.polymul(
        np.polynomial.polynomial.polypow([0, 1], s - 1),
        np.polynomial.polynomial.polypow([-1, 1], s - 1),
    )
    for _ in range(s - 2):
        poly = np.polynomial.polynomial.polyder(poly)
    return np.sort(np.real(np.polynomial.polynomial.polyroots(poly)))


# ---------------------------------------------------------------------------
# Tableaux
# ---------------------------------------------------------------------------


def quadrature_weights(c):
    """The weights that integrate a polynomial of degree s - 1 exactly."""
    s = len(c)
    vandermonde = np.vander(c, s, increasing=True).T
    targets = np.array([1.0 / (k + 1) for k in range(s)])
    return np.linalg.solve(vandermonde, targets)


def collocation_matrix(c):
    """`A` from the collocation conditions, one row at a time."""
    s = len(c)
    vandermonde = np.vander(c, s, increasing=True).T
    a = np.zeros((s, s))
    for i in range(s):
        targets = np.array([c[i] ** (k + 1) / (k + 1) for k in range(s)])
        a[i] = np.linalg.solve(vandermonde, targets)
    return a


def lobatto_iiic_matrix(c, b):
    """Lobatto IIIC pins the first column to `b_1` and asks the rest to satisfy
    the collocation conditions one degree short. That first column is what buys
    the damping at infinity that IIIA lacks."""
    s = len(c)
    a = np.zeros((s, s))
    a[:, 0] = b[0]
    rest = c[1:]
    vandermonde = np.vander(rest, s - 1, increasing=True).T
    for i in range(s):
        targets = np.array(
            [c[i] ** (k + 1) / (k + 1) - b[0] * (c[0] ** k if k > 0 else 1.0) for k in range(s - 1)]
        )
        a[i, 1:] = np.linalg.solve(vandermonde, targets)
    return a


def lobatto_iiib_matrix(c, b):
    """Lobatto IIIB is defined by the adjoint conditions, one column at a time:
    sum_i b_i c_i^(k-1) a_ij = b_j (1 - c_j^k) / k."""
    s = len(c)
    a = np.zeros((s, s))
    system = np.array([[b[i] * c[i] ** k for i in range(s)] for k in range(s)])
    for j in range(s):
        targets = np.array([b[j] * (1.0 - c[j] ** (k + 1)) / (k + 1) for k in range(s)])
        a[:, j] = np.linalg.solve(system, targets)
    return a


def tableau(a, b, c):
    return {
        "a": [[round(float(v), 15) for v in row] for row in a],
        "b": [round(float(v), 15) for v in b],
        "c": [round(float(v), 15) for v in c],
    }


def method(identifier, name, order, a, b, c, properties, aliases=None):
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
    entry["references"] = [BUTCHER64, BUTCHER16, HW2]
    return entry


written = []

# Gauss-Legendre: the maximum order a Runge-Kutta method of s stages can reach,
# and symplectic with it.
for s in (4, 5):
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
                aliases=[f"gauss{2 * s}"],
            ),
        )
    )

# Radau IIA: one order below Gauss, but L-stable and stiffly accurate, which is
# what matters on a stiff problem.
for s in (4, 5):
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
                aliases=[f"radau{2 * s - 1}"],
            ),
        )
    )

# Lobatto, all three families on the same nodes.
for s in (4,):
    c = lobatto_nodes(s)
    b = quadrature_weights(c)
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
            ),
        )
    )

print(f"{len(written)} collocation methods written")
for identifier in written:
    print(" ", identifier)
