"""Import named tableaux from NodePy's catalogue.

NodePy carries a curated set of published Runge-Kutta methods as exact sympy
arrays, and five of the very high order explicit pairs as sixty digit tables in
their authors' own notation. Both are transcriptions somebody else already made
carefully, which is a better starting point than making the transcription again.

Nothing here is trusted. Every method that comes out of this is read back by the
Rust side, which derives the order from the coefficients and refuses the file if
it does not match the order claimed for it. A tableau that arrives wrong is
therefore dropped rather than published, which has happened before.

Two things are checked here rather than there:

  * a method already in the library under another name is reported instead of
    written twice, because the same tableau under two ids would make the counts
    lie.
  * a rational stays a rational only while its parts are small enough that the
    products inside an order condition cannot overflow.

Requires the NodePy checkout:
    git clone https://github.com/ketch/nodepy C:/Repositories/TEMP/nodepy

Usage:
    python tools/import_nodepy_tableaux.py
"""

import io
import json
import os
import re
import sys
import warnings
from fractions import Fraction

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods")
NODEPY = "C:/Repositories/TEMP/nodepy"
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

# A rational is only written as one while its parts stay small enough that the
# products inside an order condition cannot overflow a 128 bit rational.
EXACT_LIMIT = 10**18

sys.path.insert(0, NODEPY)
import sympy as sp  # noqa: E402
from nodepy import runge_kutta_method as rk  # noqa: E402
from nodepy.loadmethod import load_rk_from_file  # noqa: E402

# ---------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------

HNW93 = {
    "authors": "Hairer, E., Norsett, S. P., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations I: Nonstiff Problems",
    "year": 1993,
    "source": "Springer Series in Computational Mathematics, Vol. 8, 2nd edition",
    "doi": "10.1007/978-3-540-78862-1",
}
HW96 = {
    "authors": "Hairer, E., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
    "year": 1996,
    "source": "Springer Series in Computational Mathematics, Vol. 14, 2nd edition",
    "doi": "10.1007/978-3-642-05221-7",
}
MERSON57 = {
    "authors": "Merson, R. H.",
    "title": "An operational method for the study of integration processes",
    "year": 1957,
    "source": "Proceedings of the Symposium on Data Processing, Weapons Research Establishment, Salisbury",
}
ZONNEVELD63 = {
    "authors": "Zonneveld, J. A.",
    "title": "Automatic Numerical Integration",
    "year": 1963,
    "source": "Mathematical Centre Tracts 8, Mathematisch Centrum, Amsterdam",
}
FEHLBERG69 = {
    "authors": "Fehlberg, E.",
    "title": "Low-order classical Runge-Kutta formulas with stepsize control and their application to some heat transfer problems",
    "year": 1969,
    "source": "NASA Technical Report R-315",
}
BUTCHER64 = {
    "authors": "Butcher, J. C.",
    "title": "On Runge-Kutta processes of high order",
    "year": 1964,
    "source": "Journal of the Australian Mathematical Society, 4(2), 179-194",
    "doi": "10.1017/S1446788700023387",
}
CALVO90 = {
    "authors": "Calvo, M., Montijano, J. I., & Randez, L.",
    "title": "A new embedded pair of Runge-Kutta formulas of orders 5 and 6",
    "year": 1990,
    "source": "Computers & Mathematics with Applications, 20(1), 15-24",
    "doi": "10.1016/0898-1221(90)90064-Q",
}
HIGHAM_HALL90 = {
    "authors": "Higham, D. J., & Hall, G.",
    "title": "Embedded Runge-Kutta formulae with stable equilibrium states",
    "year": 1990,
    "source": "Journal of Computational and Applied Mathematics, 29(1), 25-33",
    "doi": "10.1016/0377-0427(90)90192-3",
}
SHARP_SMART93 = {
    "authors": "Sharp, P. W., & Smart, E.",
    "title": "Explicit Runge-Kutta pairs with one more derivative evaluation than the minimum",
    "year": 1993,
    "source": "SIAM Journal on Scientific Computing, 14(2), 338-348",
    "doi": "10.1137/0914021",
}
PRINCE_DORMAND81 = {
    "authors": "Prince, P. J., & Dormand, J. R.",
    "title": "High order embedded Runge-Kutta formulae",
    "year": 1981,
    "source": "Journal of Computational and Applied Mathematics, 7(1), 67-75",
    "doi": "10.1016/0771-050X(81)90010-3",
}
TSITOURAS11 = {
    "authors": "Tsitouras, C.",
    "title": "Runge-Kutta pairs of order 5(4) satisfying only the first column simplifying assumption",
    "year": 2011,
    "source": "Computers & Mathematics with Applications, 62(2), 770-775",
    "doi": "10.1016/j.camwa.2011.06.002",
}
WANG_SPITERI07 = {
    "authors": "Wang, R., & Spiteri, R. J.",
    "title": "Linear instability of the fifth-order WENO method",
    "year": 2007,
    "source": "SIAM Journal on Numerical Analysis, 45(5), 1871-1901",
    "doi": "10.1137/050637868",
}
SPITERI_RUUTH02 = {
    "authors": "Spiteri, R. J., & Ruuth, S. J.",
    "title": "A new class of optimal high-order strong-stability-preserving time discretization methods",
    "year": 2002,
    "source": "SIAM Journal on Numerical Analysis, 40(2), 469-491",
    "doi": "10.1137/S0036142901389025",
}
RUUTH06 = {
    "authors": "Ruuth, S. J.",
    "title": "Global optimization of explicit strong-stability-preserving Runge-Kutta methods",
    "year": 2006,
    "source": "Mathematics of Computation, 75(253), 183-207",
    "doi": "10.1090/S0025-5718-05-01772-2",
}
KETCHESON08 = {
    "authors": "Ketcheson, D. I.",
    "title": "Highly efficient strong stability preserving Runge-Kutta methods with low-storage implementations",
    "year": 2008,
    "source": "SIAM Journal on Scientific Computing, 30(4), 2113-2136",
    "doi": "10.1137/07070485X",
}
NORSETT74 = {
    "authors": "Norsett, S. P.",
    "title": "Semi-explicit Runge-Kutta methods",
    "year": 1974,
    "source": "Report Mathematics and Computation No. 6/74, University of Trondheim",
}
FEAGIN12 = {
    "authors": "Feagin, T.",
    "title": "High-order explicit Runge-Kutta methods using m-symmetry",
    "year": 2012,
    "source": "Neural, Parallel & Scientific Computations, 20(3-4), 437-458",
}
CURTIS75 = {
    "authors": "Curtis, A. R.",
    "title": "High-order explicit Runge-Kutta formulae, their uses, and limitations",
    "year": 1975,
    "source": "IMA Journal of Applied Mathematics, 16(1), 35-55",
    "doi": "10.1093/imamat/16.1.35",
}

# ---------------------------------------------------------------------------
# What to take
# ---------------------------------------------------------------------------

# key in NodePy's catalogue -> what the library calls it
CATALOGUE = [
    ("Merson43", "merson43", "Merson 4(3)", "erk", 4, 3, [MERSON57, HNW93]),
    ("Zonneveld43", "zonneveld43", "Zonneveld 4(3)", "erk", 4, 3, [ZONNEVELD63, HNW93]),
    ("Fehlberg43", "rkf43", "Fehlberg 4(3)", "erk", 4, 3, [FEHLBERG69, HNW93]),
    ("BuRK65", "butcher5", "Butcher 5", "erk", 5, None, [BUTCHER64, HNW93]),
    ("CMR6", "cmr65", "Calvo-Montijano-Randez 6(5)", "erk", 6, 5, [CALVO90]),
    ("HH5", "hh54", "Higham-Hall 5(4)", "erk", 5, 4, [HIGHAM_HALL90]),
    ("HH5S", "hh54s", "Higham-Hall 5(4) stable", "erk", 5, 4, [HIGHAM_HALL90]),
    ("SS3", "ss32", "Sharp-Smart 3(2)", "erk", 3, 2, [SHARP_SMART93]),
    ("PD8", "pd87", "Prince-Dormand 8(7)", "erk", 8, 7, [PRINCE_DORMAND81]),
    ("Tsit5", "tsit54", "Tsitouras 5(4)", "erk", 5, 4, [TSITOURAS11]),
    ("NSSP32", "nssp32", "NSSP(3,2)", "ssprk", 2, None, [WANG_SPITERI07]),
    ("NSSP33", "nssp33", "NSSP(3,3)", "ssprk", 3, None, [WANG_SPITERI07]),
    ("SSP54", "ssprk54", "SSPRK(5,4)", "ssprk", 4, None, [SPITERI_RUUTH02]),
    ("SSP63", "ssprk63", "SSPRK(6,3)", "ssprk", 3, None, [RUUTH06]),
    ("SSP75", "ssprk75", "SSPRK(7,5)", "ssprk", 5, None, [RUUTH06]),
    ("SSP85", "ssprk85", "SSPRK(8,5)", "ssprk", 5, None, [RUUTH06]),
    ("SSP95", "ssprk95", "SSPRK(9,5)", "ssprk", 5, None, [RUUTH06]),
    ("SSP104", "ssprk10_4", "SSPRK(10,4)", "ssprk", 4, None, [KETCHESON08]),
    ("SDIRK23", "sdirk23", "SDIRK 2 stage 3rd order", "dirk", 3, None, [NORSETT74, HW96]),
    ("SDIRK34", "sdirk34", "SDIRK 3 stage 4th order", "dirk", 4, None, [NORSETT74, HW96]),
    ("SDIRK54", "sdirk54", "SDIRK 5 stage 4th order", "dirk", 4, None, [HW96]),
]

# file in NodePy's coefficient directory -> what the library calls it
# Feagin's tables are carried without their error estimates, so those three are
# the single formula and are named that way.
FILES = [
    ("rk108.txt", "feagin10", "Feagin 10", 10, None, [FEAGIN12]),
    ("rk1210.txt", "feagin12", "Feagin 12", 12, None, [FEAGIN12]),
    ("rk1412.txt", "feagin14", "Feagin 14", 14, None, [FEAGIN12]),
    ("rk108curtis.txt", "curtis108", "Curtis 10(8)", 10, 8, [CURTIS75]),
]


# ---------------------------------------------------------------------------
# Coefficients
# ---------------------------------------------------------------------------


def coefficient(value):
    """One entry, exact where the source is exact and a double where it is not.

    The source decides, not the appearance: a rational stays a rational, and a
    float stays a float. Trying to recognise a fraction inside a decimal is
    exactly the mistake this library exists to avoid, and on a sixty digit table
    truncated to a double it invents coefficients that were never published.
    """
    expression = sp.sympify(value)
    if expression.is_Integer:
        return int(expression)
    if expression.is_Rational:
        fraction = Fraction(int(expression.p), int(expression.q))
        if abs(fraction.numerator) < EXACT_LIMIT and fraction.denominator < EXACT_LIMIT:
            return f"{fraction.numerator}/{fraction.denominator}"
    return float(expression)


def tableau(a, b, b_embedded=None):
    """Explicit tableaux keep their triangle; `c` is left out so the loader
    derives it from the row sums and the analysis checks it."""
    s = len(b)
    explicit = all(a[i][j] == 0 for i in range(s) for j in range(i, s))
    rows = []
    for i in range(s):
        row = a[i][: i + 1] if explicit else a[i]
        rows.append([coefficient(v) for v in row])
    out = {"a": rows, "b": [coefficient(v) for v in b]}
    if b_embedded is not None:
        out["b_embedded"] = [coefficient(v) for v in b_embedded]
    return out


def embedded_weights(b, b_embedded):
    """The embedded weights, whichever of the two conventions the source used.

    Some tables state the second formula, others state its difference from the
    first. The two are told apart by what they sum to: a set of weights sums to
    one, a difference sums to zero.
    """
    if b_embedded is None or all(float(v) == 0.0 for v in b_embedded):
        return None
    total = float(sum(float(v) for v in b_embedded))
    if abs(total - 1.0) < 1e-8:
        return b_embedded
    if abs(total) < 1e-8:
        return [x - y for x, y in zip(b, b_embedded)]
    raise SystemExit(f"embedded weights sum to {total}, which is neither convention")


def numeric(a, b):
    """The tableau as plain doubles, for comparing against what is already there."""
    return [[float(v) for v in row] for row in a], [float(v) for v in b]


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    return INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text) + "\n"


def existing():
    """Every tableau already in the library, as doubles, by id."""
    out = {}
    for folder in sorted(os.listdir(OUT)):
        directory = os.path.join(OUT, folder)
        if not os.path.isdir(directory):
            continue
        for name in sorted(os.listdir(directory)):
            entry = json.load(io.open(os.path.join(directory, name), encoding="utf-8"))
            table = entry.get("tableau")
            if not table or "a" not in table:
                continue
            try:
                size = len(table["b"])
                matrix = [[0.0] * size for _ in range(size)]
                for i, row in enumerate(table["a"]):
                    for j, value in enumerate(row):
                        matrix[i][j] = float(sp.sympify(str(value).replace("sqrt", "sqrt")))
                weights = [float(sp.sympify(str(v))) for v in table["b"]]
            except Exception:
                continue
            out[entry["id"]] = (matrix, weights)
    return out


def same(left, right, tol=1e-12):
    (a1, b1), (a2, b2) = left, right
    if len(b1) != len(b2):
        return False
    for row1, row2 in zip(a1, a2):
        for x, y in zip(row1, row2):
            if abs(x - y) > tol:
                return False
    return all(abs(x - y) <= tol for x, y in zip(b1, b2))


def write(family, entry):
    folder = os.path.join(OUT, family)
    os.makedirs(folder, exist_ok=True)
    with io.open(os.path.join(folder, entry["id"] + ".json"), "w", encoding="utf-8") as handle:
        handle.write(render(entry))


def emit(identifier, name, family, order, embedded_order, a, b, b_embedded, references, library):
    duplicate = next(
        (k for k, v in library.items() if k != identifier and same(v, numeric(a, b))), None
    )
    if duplicate:
        print(f"  {identifier}: already in the library as {duplicate}, skipped")
        return None
    entry = {
        "id": identifier,
        "name": name,
        "class": "runge_kutta",
        "family": family,
        "order": order,
    }
    if embedded_order is not None:
        entry["embedded_order"] = embedded_order
    entry["tableau"] = tableau(a, b, b_embedded)
    entry["references"] = references
    write(family, entry)
    library[identifier] = numeric(a, b)
    return identifier


library = existing()
written = []

print("catalogue:")
methods = rk.loadRKM("All")
for key, identifier, name, family, order, embedded_order, references in CATALOGUE:
    method = methods[key]
    a = [[method.A[i, j] for j in range(len(method.b))] for i in range(len(method.b))]
    b = list(method.b)
    b_embedded = list(method.bhat) if getattr(method, "bhat", None) is not None else None
    b_embedded = embedded_weights(b, b_embedded)
    if b_embedded is None:
        embedded_order = None
    elif getattr(method, "embedded_method", None) is not None and method.embedded_method.order(
        tol=1e-11
    ) > method.order(tol=1e-11):
        # The reference stores this pair the other way round, with the higher
        # order formula as the estimate.
        b, b_embedded = b_embedded, b
    result = emit(
        identifier, name, family, order, embedded_order, a, b, b_embedded, references, library
    )
    if result:
        written.append(result)

print("high order pairs:")
for filename, identifier, name, order, embedded_order, references in FILES:
    method = load_rk_from_file(filename)
    size = len(method.b)
    a = [[method.A[i, j] for j in range(size)] for i in range(size)]
    weights = list(method.b)
    estimate = embedded_weights(weights, list(method.bhat))
    result = emit(
        identifier,
        name,
        "erk",
        order,
        embedded_order if estimate is not None else None,
        a,
        weights,
        estimate,
        references,
        library,
    )
    if result:
        written.append(result)

print(f"{len(written)} methods written")
for identifier in written:
    print(" ", identifier)
