"""Import named explicit Runge-Kutta tableaux from a reference implementation.

The reference stores a tableau as a constructor full of named assignments,
`a21`, `a32`, `b1`, `btilde1`, one per nonzero entry. That is a convention, not
a format, but it is a consistent one, so the constructor is translated into
Python and evaluated and the names are read back into a matrix. Nothing is
retyped, which is the part of this job a person gets wrong.

Three conventions have to be undone:

  * `b` is often absent because the method is FSAL. Then the last row of `a` is
    the weight vector, which is what FSAL means.
  * `btilde` is the difference `b - b_hat`, not the embedded weights.
  * `c` is stated but not stored here. It is the row sums of `a` by definition,
    so leaving it out lets the loader derive it, and the analysis then checks
    the abscissae rather than trusting them.

Provenance follows the source. Where the reference writes a rational, the
coefficient is exact and is written as a fraction; where it writes a decimal, it
stays a double. Fractions past what a hundred and twenty eight bit rational can
multiply out are written as doubles too, because an exact coefficient that
overflows during the order conditions is worse than an honest approximation.

Usage:
    python tools/import_runge_kutta_tableaux.py
"""

import io
import json
import math
import os
import re
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods")
LIB = "C:/Repositories/TEMP/OrdinaryDiffEq.jl/lib"
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

# A rational is only written as one while its parts stay small enough that the
# products inside an order condition cannot overflow.
EXACT_LIMIT = 10**18

SOURCES = {
    "low": "OrdinaryDiffEqLowOrderRK/src/low_order_rk_tableaus.jl",
    "high": "OrdinaryDiffEqHighOrderRK/src/high_order_rk_tableaus.jl",
    "verner": "OrdinaryDiffEqVerner/src/verner_tableaus.jl",
    "tsit": "OrdinaryDiffEqTsit5/src/tsit_tableaus.jl",
}

OWREN92 = {
    "authors": "Owren, B., & Zennaro, M.",
    "title": "Derivation of efficient, continuous, explicit Runge-Kutta methods",
    "year": 1992,
    "source": "SIAM Journal on Scientific and Statistical Computing, 13(6), 1488-1501",
    "doi": "10.1137/0913072",
}
VERNER10 = {
    "authors": "Verner, J. H.",
    "title": "Numerically optimal Runge-Kutta pairs with interpolants",
    "year": 2010,
    "source": "Numerical Algorithms, 53(2-3), 383-396",
    "doi": "10.1007/s11075-009-9290-3",
}

WANTED = [
    dict(
        source="low",
        constructor="OwrenZen3ConstantCache",
        id="owrenzen3",
        name="Owren-Zennaro 3(2)",
        family="erk",
        order=3,
        embedded_order=2,
        references=[OWREN92],
    ),
    dict(
        source="low",
        constructor="OwrenZen4ConstantCache",
        id="owrenzen4",
        name="Owren-Zennaro 4(3)",
        family="erk",
        order=4,
        embedded_order=3,
        references=[OWREN92],
    ),
    dict(
        source="low",
        constructor="OwrenZen5ConstantCache",
        id="owrenzen5",
        name="Owren-Zennaro 5(4)",
        family="erk",
        order=5,
        embedded_order=4,
        references=[OWREN92],
    ),
    dict(
        source="low",
        constructor="BS5ConstantCache",
        id="rkbs54",
        name="Bogacki-Shampine 5(4)",
        family="erk",
        order=5,
        embedded_order=4,
        aliases=["bs5"],
        references=[
            {
                "authors": "Bogacki, P., & Shampine, L. F.",
                "title": "An efficient Runge-Kutta (4,5) pair",
                "year": 1996,
                "source": "Computers & Mathematics with Applications, 32(6), 15-28",
                "doi": "10.1016/0898-1221(96)00141-1",
            }
        ],
    ),
    dict(
        source="low",
        constructor="Ralston4ConstantCache",
        id="ralston4",
        name="Ralston 4",
        family="erk",
        order=4,
        references=[
            {
                "authors": "Ralston, A.",
                "title": "Runge-Kutta methods with minimum error bounds",
                "year": 1962,
                "source": "Mathematics of Computation, 16(80), 431-437",
                "doi": "10.1090/S0025-5718-1962-0150954-0",
            }
        ],
    ),
    dict(
        source="tsit",
        constructor="Tsit5ConstantCacheActual",
        id="tsit5",
        name="Tsitouras 5(4)",
        family="erk",
        order=5,
        embedded_order=4,
        references=[
            {
                "authors": "Tsitouras, C.",
                "title": "Runge-Kutta pairs of order 5(4) satisfying only the first column simplifying assumption",
                "year": 2011,
                "source": "Computers & Mathematics with Applications, 62(2), 770-775",
                "doi": "10.1016/j.camwa.2011.06.002",
            }
        ],
    ),
    dict(
        source="high",
        constructor="TanYam7ConstantCache",
        id="tanyam7",
        name="Tanaka-Yamashita 7(6)",
        family="erk",
        order=7,
        embedded_order=6,
        references=[
            {
                "authors": "Tanaka, M., Muramatsu, S., & Yamashita, S.",
                "title": "On the optimization of some nine-stage seventh-order Runge-Kutta method",
                "year": 1992,
                "source": "Information Processing Society of Japan, 33(12), 1512-1526",
            }
        ],
    ),
    dict(
        source="high",
        constructor="TsitPap8ConstantCache",
        id="tsitpap8",
        name="Tsitouras-Papakostas 8(7)",
        family="erk",
        order=8,
        embedded_order=7,
        references=[
            {
                "authors": "Tsitouras, C., & Papakostas, S. N.",
                "title": "Cheap error estimation for Runge-Kutta pairs",
                "year": 1999,
                "source": "SIAM Journal on Scientific Computing, 20(6), 2067-2088",
                "doi": "10.1137/S1064827597315509",
            }
        ],
    ),
    dict(
        source="verner",
        constructor="Vern6Tableau",
        id="rkv65",
        name="Verner 6(5)",
        family="erk",
        order=6,
        embedded_order=5,
        references=[VERNER10],
    ),
    dict(
        source="verner",
        constructor="Vern8Tableau",
        id="rkv87",
        name="Verner 8(7)",
        family="erk",
        order=8,
        embedded_order=7,
        references=[VERNER10],
    ),
]


# ---------------------------------------------------------------------------
# Reading one constructor
# ---------------------------------------------------------------------------


def balanced(text, start):
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    raise ValueError("unbalanced parenthesis")


def strip_converts(text):
    while True:
        match = re.search(r"convert\(\s*T[12]?\s*,", text)
        if not match:
            return text
        close = balanced(text, match.start() + len("convert"))
        text = text[: match.start()] + "(" + text[match.end() : close] + ")" + text[close + 1 :]


def bodies(text, name):
    """Every definition of a constructor, in the order they appear.

    The reference gives two: one in decimals for compiled floats and one in
    rationals for everything else. The exact one is preferred and comes second.
    """
    found = []
    pattern = re.compile(r"^function " + re.escape(name) + r"\(.*?\n(.*?)^end$", re.M | re.S)
    for match in pattern.finditer(text):
        found.append(match.group(1))
    return found


def translate(body):
    """Julia assignments into Python, with rationals kept as rationals."""
    body = strip_converts(body)
    body = re.sub(r"#.*", "", body)
    # A Julia line continues while its brackets are open, so join first.
    joined = []
    depth = 0
    current = ""
    for line in body.splitlines():
        current += " " + line.strip()
        depth += line.count("(") - line.count(")")
        if depth <= 0:
            joined.append(current.strip())
            current = ""
            depth = 0
    lines = []
    for line in joined:
        if not line or line.startswith("return") or "=" not in line:
            continue
        name, _, value = line.partition("=")
        name = name.strip()
        if not re.fullmatch(r"[A-Za-z_]\w*", name):
            continue
        value = re.sub(r"(\d+)\s*//\s*(-?\d+)", r"Fraction(\1, \2)", value)
        value = value.replace("^", "**")
        if "big" in value or "parse" in value:
            return None
        lines.append(f"{name} = {value.strip()}")
    return "\n".join(lines)


def evaluate(body):
    """Run the assignments, one at a time.

    A constructor ends by building things this does not care about, out of
    names this does not have. Running line by line takes the coefficients and
    lets the rest fail where it stands.
    """
    source = translate(body)
    if source is None:
        return None
    scope = {"Fraction": Fraction, "sqrt": math.sqrt, "float": float}
    for line in source.splitlines():
        try:
            exec(compile(line, "<tableau>", "exec"), scope)
        except Exception:  # noqa: BLE001 - a line this does not need
            continue
    return scope


# ---------------------------------------------------------------------------
# Names into a tableau
# ---------------------------------------------------------------------------

ENTRY = re.compile(r"^a(\d\d)$|^a(\d+)_(\d+)$")
WEIGHT = re.compile(r"^b(\d+)$")
ERROR = re.compile(r"^btilde(\d+)$")


def indices(name):
    match = ENTRY.match(name)
    if not match:
        return None
    if match.group(1):
        return int(match.group(1)[0]), int(match.group(1)[1])
    return int(match.group(2)), int(match.group(3))


def assemble(scope):
    """The matrix and the two weight vectors, from the names alone."""
    entries = {}
    weights = {}
    errors = {}
    for name, value in scope.items():
        if name in ("Fraction", "__builtins__"):
            continue
        if not isinstance(value, (int, float, Fraction)):
            continue
        pair = indices(name)
        if pair:
            entries[pair] = value
            continue
        match = WEIGHT.match(name)
        if match:
            weights[int(match.group(1))] = value
            continue
        match = ERROR.match(name)
        if match:
            errors[int(match.group(1))] = value

    if not entries:
        return None
    stages = max(max(i for i, _ in entries), max(weights or [0]), max(errors or [0]))

    a = [[0] * stages for _ in range(stages)]
    for (i, j), value in entries.items():
        if j >= i or i > stages:
            return None
        a[i - 1][j - 1] = value

    if weights:
        b = [weights.get(i + 1, 0) for i in range(stages)]
    else:
        # No weights stated means the last stage is the solution: FSAL.
        b = list(a[stages - 1])
        if not any(b):
            return None

    b_embedded = None
    if errors:
        b_embedded = [b[i] - errors.get(i + 1, 0) for i in range(stages)]
    return a, b, b_embedded


def coefficient(value):
    """A rational stays one while it is small enough to be multiplied out."""
    if isinstance(value, Fraction):
        if value.denominator == 1 and abs(value.numerator) < EXACT_LIMIT:
            return value.numerator
        if abs(value.numerator) < EXACT_LIMIT and value.denominator < EXACT_LIMIT:
            return f"{value.numerator}/{value.denominator}"
        return float(value)
    if isinstance(value, int):
        return value
    if value == int(value):
        return int(value)
    return round(float(value), 17)


def integrates_to(a, b, order):
    """Whether the weights integrate a polynomial to the claimed order.

    `sum_i b_i c_i^(k-1) = 1/k` for `k = 1..p` is necessary for order `p` and
    costs nothing to check. It is here as a gate rather than as an analysis: a
    tableau read out of a naming convention can be read wrong, and one that
    cannot even do quadrature has been. The Rust side does the real work on
    what survives.
    """
    c = [sum(row) for row in a]
    for k in range(1, order + 1):
        total = sum(float(b[i]) * float(c[i]) ** (k - 1) for i in range(len(b)))
        if abs(total - 1.0 / k) > 1e-10:
            return False
    return True


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def main():
    text = {key: io.open(os.path.join(LIB, path), encoding="utf-8").read() for key, path in SOURCES.items()}
    written = []
    for entry in WANTED:
        candidates = bodies(text[entry["source"]], entry["constructor"])
        if not candidates:
            print(entry["id"], ": constructor not found")
            continue
        # The exact definition comes last where there are two of them.
        table = None
        for body in reversed(candidates):
            scope = evaluate(body)
            if scope is None:
                continue
            table = assemble(scope)
            if table:
                break
        if not table:
            print(entry["id"], ": could not be read")
            continue
        a, b, b_embedded = table
        if not integrates_to(a, b, entry["order"]):
            print(entry["id"], ": read out of the source but does not integrate to order", entry["order"])
            continue

        method = {
            "id": entry["id"],
            "name": entry["name"],
            "class": "runge_kutta",
            "family": entry["family"],
            "order": entry["order"],
        }
        if entry.get("aliases"):
            method["aliases"] = entry["aliases"]
        if entry.get("embedded_order") and b_embedded:
            method["embedded_order"] = entry["embedded_order"]
        method["tableau"] = {
            "a": [[coefficient(v) for v in row[:i]] for i, row in enumerate(a)],
            "b": [coefficient(v) for v in b],
        }
        if b_embedded:
            method["tableau"]["b_embedded"] = [coefficient(v) for v in b_embedded]
        method["references"] = entry["references"]

        folder = os.path.join(OUT, entry["family"])
        os.makedirs(folder, exist_ok=True)
        io.open(os.path.join(folder, entry["id"] + ".json"), "w", encoding="utf-8").write(
            render(method)
        )
        exact = sum(1 for row in a for v in row if isinstance(v, Fraction))
        written.append((entry["id"], len(b), exact))

    print(f"{len(written)} Runge-Kutta method files written")
    for identifier, stages, exact in written:
        print(f"  {identifier:<12} {stages:>2} stages   {exact} rational entries")


if __name__ == "__main__":
    main()
