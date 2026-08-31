"""Construct the strong stability preserving families from their Shu-Osher form.

An SSP method is not chosen for its order. It is chosen because the step can be
written as a convex combination of forward Euler steps, so whatever bound holds
for Euler (positivity, a maximum principle, total variation) survives the whole
step, up to a step size limit `C h_FE`. That coefficient `C` is the thing being
optimized, and the optimal methods of each family are known in closed form.

The families are published in Shu-Osher form, which is where the convexity is
visible:

  u^(i) = (1 - sum_j alpha_ij) u^n + sum_j alpha_ij u^(j) + h sum_j beta_ij F(u^(j))

with all alpha and beta non-negative. The Butcher arrays follow from

  A = [(I - alpha)^-1 beta]_{1..s, 1..s},   b = [(I - alpha)^-1 beta]_{s+1, 1..s}

which is a linear solve, done here in exact arithmetic. Writing the construction
down instead of the tableaux keeps the reason the coefficients look the way they
do, and the radius of absolute monotonicity computed on the Rust side then has
to come back as the `C` the paper claims.

Usage:
    python tools/generate_ssp_methods.py
"""

import io
import json
import os
import re

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods")
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

KETCHESON08 = {
    "authors": "Ketcheson, D. I.",
    "title": "Highly efficient strong stability preserving Runge-Kutta methods with low-storage implementations",
    "year": 2008,
    "source": "SIAM Journal on Scientific Computing, 30(4), 2113-2136",
    "doi": "10.1137/07070485X",
}
KRAAIJEVANGER91 = {
    "authors": "Kraaijevanger, J. F. B. M.",
    "title": "Contractivity of Runge-Kutta methods",
    "year": 1991,
    "source": "BIT Numerical Mathematics, 31(3), 482-528",
    "doi": "10.1007/BF01933264",
}
GOTTLIEB01 = {
    "authors": "Gottlieb, S., Shu, C.-W., & Tadmor, E.",
    "title": "Strong stability-preserving high-order time discretization methods",
    "year": 2001,
    "source": "SIAM Review, 43(1), 89-112",
    "doi": "10.1137/S003614450036757X",
}
GKS11 = {
    "authors": "Gottlieb, S., Ketcheson, D. I., & Shu, C.-W.",
    "title": "Strong Stability Preserving Runge-Kutta and Multistep Time Discretizations",
    "year": 2011,
    "source": "World Scientific",
    "doi": "10.1142/7498",
}
SHU_OSHER88 = {
    "authors": "Shu, C.-W., & Osher, S.",
    "title": "Efficient implementation of essentially non-oscillatory shock-capturing schemes",
    "year": 1988,
    "source": "Journal of Computational Physics, 77(2), 439-471",
    "doi": "10.1016/0021-9991(88)90177-5",
}
KMG09 = {
    "authors": "Ketcheson, D. I., Macdonald, C. B., & Gottlieb, S.",
    "title": "Optimal implicit strong stability preserving Runge-Kutta methods",
    "year": 2009,
    "source": "Applied Numerical Mathematics, 59(2), 373-392",
    "doi": "10.1016/j.apnum.2008.03.034",
}


# ---------------------------------------------------------------------------
# Shu-Osher to Butcher
# ---------------------------------------------------------------------------


def euler_chain(s):
    """The skeleton every one of these starts from: stage i is a forward Euler
    step off stage i-1, and the last row is the update."""
    alpha = sp.zeros(s + 1, s)
    for i in range(1, s + 1):
        alpha[i, i - 1] = 1
    return alpha


def butcher(alpha, beta):
    """`A` and `b` from the Shu-Osher pair, exactly."""
    s = alpha.shape[1]
    a_pad = sp.zeros(s + 1, s + 1)
    b_pad = sp.zeros(s + 1, s + 1)
    a_pad[:, :s] = alpha
    b_pad[:, :s] = beta
    k = sp.simplify((sp.eye(s + 1) - a_pad).inv() * b_pad)
    a = [[k[i, j] for j in range(s)] for i in range(s)]
    b = [k[s, j] for j in range(s)]
    return a, b


def coefficient(value):
    """A JSON coefficient: an integer where it can be, a string otherwise, so the
    order conditions stay exact."""
    value = sp.nsimplify(sp.simplify(value))
    if value.is_Integer:
        return int(value)
    return str(value).replace(" ", " ")


def tableau(a, b):
    """Explicit tableaux keep their triangle, the rest is written in full."""
    s = len(b)
    explicit = all(a[i][j] == 0 for i in range(s) for j in range(i, s))
    rows = []
    for i in range(s):
        row = a[i][: i + 1] if explicit else a[i]
        rows.append([coefficient(v) for v in row])
    c = [coefficient(sum(a[i])) for i in range(s)]
    return {"a": rows, "b": [coefficient(v) for v in b], "c": c}


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------


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


def method(identifier, name, family, order, a, b, references, aliases=None, properties=None):
    entry = {
        "id": identifier,
        "name": name,
        "class": "runge_kutta",
        "family": family,
        "order": order,
    }
    if aliases:
        entry["aliases"] = aliases
    if properties:
        entry["properties"] = properties
    entry["tableau"] = tableau(a, b)
    entry["references"] = references
    return entry


def identifier(prefix, stages, order):
    """Two single digits run together the way the field writes them, anything
    wider gets a separator so `ssprk10_2` cannot be read as `ssprk1_02`."""
    if stages < 10:
        return f"{prefix}{stages}{order}"
    return f"{prefix}{stages}_{order}"


written = []

# ---------------------------------------------------------------------------
# Explicit: second order, s-1 is the largest SSP coefficient s stages can give
# ---------------------------------------------------------------------------
for s in range(2, 11):
    r = sp.Rational(s - 1)
    alpha = euler_chain(s)
    alpha[s, s - 1] = sp.Rational(s - 1, s)
    beta = alpha / r
    alpha[s, 0] = sp.Rational(1, s)
    a, b = butcher(alpha, beta)
    written.append(
        write(
            "ssprk",
            method(
                identifier("ssprk", s, 2),
                f"SSPRK({s},2)",
                "ssprk",
                2,
                a,
                b,
                [SHU_OSHER88, GOTTLIEB01, GKS11]
                if s == 2
                else [KRAAIJEVANGER91, GOTTLIEB01, KETCHESON08, GKS11],
            ),
        )
    )

# ---------------------------------------------------------------------------
# Explicit: third order, available whenever the stage count is a square
# ---------------------------------------------------------------------------
for n in range(2, 6):
    s = n * n
    r = sp.Rational(s - n)
    k = n * (n + 1) // 2
    alpha = euler_chain(s)
    alpha[k, k - 1] = sp.Rational(n - 1, 2 * n - 1)
    beta = alpha / r
    alpha[k, (n - 1) * (n - 2) // 2] = sp.Rational(n, 2 * n - 1)
    a, b = butcher(alpha, beta)
    written.append(
        write(
            "ssprk",
            method(
                identifier("ssprk", s, 3),
                f"SSPRK({s},3)",
                "ssprk",
                3,
                a,
                b,
                [KETCHESON08, GKS11],
                aliases=["ssprk34"] if s == 4 else None,
            ),
        )
    )

# ---------------------------------------------------------------------------
# Implicit: unconditionally SSP for the first order family, and the optimal
# second and third order families, all singly diagonally implicit
# ---------------------------------------------------------------------------
for s in range(2, 6):
    a = [[sp.Rational(1, s) if j <= i else sp.Integer(0) for j in range(s)] for i in range(s)]
    b = [sp.Rational(1, s) for _ in range(s)]
    written.append(
        write(
            "dirk",
            method(
                identifier("sspirk", s, 1),
                f"SSPIRK({s},1)",
                "dirk",
                1,
                a,
                b,
                [KMG09, GKS11],
            ),
        )
    )

for s in range(2, 6):
    r = sp.Integer(2 * s)
    alpha = euler_chain(s)
    beta = alpha / r
    for i in range(s):
        beta[i, i] = sp.Rational(1, r)
    a, b = butcher(alpha, beta)
    written.append(
        write(
            "dirk",
            method(
                identifier("sspirk", s, 2),
                f"SSPIRK({s},2)",
                "dirk",
                2,
                a,
                b,
                [KMG09, GKS11],
            ),
        )
    )

for s in range(2, 6):
    r = s - 1 + sp.sqrt(s * s - 1)
    alpha = euler_chain(s)
    alpha[s, s - 1] = ((s + 1) * r) / (s * (r + 2))
    beta = alpha / r
    for i in range(s):
        beta[i, i] = (1 - sp.sqrt(sp.Rational(s - 1, s + 1))) / 2
    a, b = butcher(alpha, beta)
    written.append(
        write(
            "dirk",
            method(
                identifier("sspirk", s, 3),
                f"SSPIRK({s},3)",
                "dirk",
                3,
                a,
                b,
                [KMG09, GKS11],
            ),
        )
    )

print(f"{len(written)} SSP methods written")
for entry in written:
    print(" ", entry)
