"""Import Rosenbrock tableaux from a reference implementation.

The published ROW coefficients live in the substituted form every code runs,

    a = alpha Gamma^-1,  C = diag(1/gamma_ii) - Gamma^-1,  m^T = b^T Gamma^-1

while a method file here stores `alpha`, `gamma` and `b`, which is the form the
order conditions are written in. This undoes the substitution.

The source is read rather than retyped. A constructor in the reference is
straight line arithmetic, so it is translated into Python and evaluated instead
of being transcribed by hand, which is the part of this job a person gets wrong.
Every conversion is then checked against two quantities the source states
independently of the matrices: the abscissae `c`, which have to come out as the
row sums of `alpha`, and the time derivative weights `d`, which have to come out
as the row sums of `gamma`. A mistake anywhere in the algebra breaks one of
those, and the Rust side afterwards derives the order from nothing but the
result.

Usage:
    python tools/import_rosenbrock_tableaux.py [path to rosenbrock_tableaus.jl]
"""

import io
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods", "rosenbrock")
LIB = "C:/Repositories/TEMP/OrdinaryDiffEq.jl/lib"
DEFAULT_SOURCES = [
    os.path.join(LIB, "OrdinaryDiffEqRosenbrockTableaus/src/rosenbrock_tableaus.jl"),
    os.path.join(LIB, "OrdinaryDiffEqRosenbrock/src/rosenbrock_tableaus.jl"),
]
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

# Which constructors to take, and what the entry should say about itself. The
# order is not recorded here: it is derived from the coefficients by the Rust
# side, and this file only claims what the publication claims so the two can be
# held against each other.
RANG05 = {
    "authors": "Rang, J., & Angermann, L.",
    "title": "New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1",
    "year": 2005,
    "source": "BIT Numerical Mathematics, 45(4), 761-787",
    "doi": "10.1007/s10543-005-0035-y",
}
RANG15 = {
    "authors": "Rang, J.",
    "title": "Improved traditional Rosenbrock-Wanner methods for stiff ODEs and DAEs",
    "year": 2015,
    "source": "Journal of Computational and Applied Mathematics, 286, 128-144",
    "doi": "10.1016/j.cam.2015.03.010",
}
RANG16 = {
    "authors": "Rang, J.",
    "title": "The Prothero and Robinson example: Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods",
    "year": 2016,
    "source": "Applied Numerical Mathematics, 108, 37-56",
    "doi": "10.1016/j.apnum.2016.04.012",
}
STEINEBACH20 = {
    "authors": "Steinebach, G.",
    "title": "Improvement of Rosenbrock-Wanner method RODASP",
    "year": 2020,
    "source": "Progress in Differential-Algebraic Equations II, Springer, 165-184",
    "doi": "10.1007/978-3-030-53905-4_6",
}
STEINEBACH23 = {
    "authors": "Steinebach, G.",
    "title": "Construction of Rosenbrock-Wanner method Rodas5P and numerical benchmarks within the Julia Differential Equations package",
    "year": 2023,
    "source": "BIT Numerical Mathematics, 63(2), 27",
    "doi": "10.1007/s10543-023-00967-x",
}
STEINEBACH24 = {
    "authors": "Steinebach, G.",
    "title": "Rosenbrock methods within OrdinaryDiffEq.jl: overview, recent developments and applications",
    "year": 2024,
    "source": "Proceedings of the JuliaCon Conferences",
}
STEINEBACH25 = {
    "authors": "Steinebach, G.",
    "title": "Rodas6P and Tsit5DA: two new Rosenbrock-type methods for DAEs",
    "year": 2025,
    "source": "arXiv:2511.21252",
    "url": "https://arxiv.org/abs/2511.21252",
}

WANTED = {
    "ROS2": dict(
        id="ros2",
        name="ROS2",
        order=2,
        embedded_order=1,
        properties={"a_stable": True, "l_stable": True},
        references=[dict(
            authors="Verwer, J. G., Spee, E. J., Blom, J. G., & Hundsdorfer, W.",
            title="A second-order Rosenbrock method applied to photochemical dispersion problems",
            year=1999,
            source="SIAM Journal on Scientific Computing, 20(4), 1456-1480",
            doi="10.1137/S1064827597326651",
        )],
    ),
    "ROS3P": dict(
        id="ros3p",
        name="ROS3P",
        order=3,
        embedded_order=2,
        properties={"a_stable": True, "l_stable": False},
        references=[dict(
            authors="Lang, J., & Verwer, J.",
            title="ROS3P - an accurate third-order Rosenbrock solver designed for parabolic problems",
            year=2001,
            source="BIT Numerical Mathematics, 41(4), 731-738",
            doi="10.1023/A:1021900219772",
        )],
    ),
    "Rodas3": dict(
        id="rodas3",
        name="RODAS3",
        order=3,
        embedded_order=2,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Sandu, A., Verwer, J. G., Blom, J. G., Spee, E. J., Carmichael, G. R., & Potra, F. A.",
            title="Benchmarking stiff ODE solvers for atmospheric chemistry problems II: Rosenbrock solvers",
            year=1997,
            source="Atmospheric Environment, 31(20), 3459-3472",
            doi="10.1016/S1352-2310(97)83212-8",
        )],
    ),
    "GRK4T": dict(
        id="grk4t",
        name="GRK4T",
        order=4,
        embedded_order=3,
        properties={"a_stable": False},
        references=[dict(
            authors="Kaps, P., & Rentrop, P.",
            title="Generalized Runge-Kutta methods of order four with stepsize control for stiff ordinary differential equations",
            year=1979,
            source="Numerische Mathematik, 33(1), 55-68",
            doi="10.1007/BF01396495",
        )],
    ),
    "GRK4A": dict(
        id="grk4a",
        name="GRK4A",
        order=4,
        embedded_order=3,
        properties={"a_stable": True},
        references=[dict(
            authors="Kaps, P., & Rentrop, P.",
            title="Generalized Runge-Kutta methods of order four with stepsize control for stiff ordinary differential equations",
            year=1979,
            source="Numerische Mathematik, 33(1), 55-68",
            doi="10.1007/BF01396495",
        )],
    ),
    "Ros4LStab": dict(
        id="ros4lstab",
        name="RODAS4 L-stable",
        order=4,
        embedded_order=3,
        properties={"a_stable": True, "l_stable": True},
        references=[dict(
            authors="Hairer, E., & Wanner, G.",
            title="Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
            year=1996,
            source="Springer Series in Computational Mathematics, Vol. 14, 2nd edition, IV.7",
            doi="10.1007/978-3-642-05221-7",
        )],
    ),
    "Rodas4": dict(
        id="rodas4",
        name="RODAS4",
        order=4,
        embedded_order=3,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Hairer, E., & Wanner, G.",
            title="Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
            year=1996,
            source="Springer Series in Computational Mathematics, Vol. 14, 2nd edition, IV.7",
            doi="10.1007/978-3-642-05221-7",
        )],
    ),
    "Rodas4P": dict(
        id="rodas4p",
        name="RODAS4P",
        order=4,
        embedded_order=3,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Steinebach, G.",
            title="Order reduction of ROW methods for DAEs and method of lines applications",
            year=1995,
            source="Technical Report 1741, Technische Hochschule Darmstadt",
        )],
    ),
    "ROS34PW2": dict(
        id="ros34pw2",
        name="ROS34PW2",
        order=3,
        embedded_order=2,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Rang, J., & Angermann, L.",
            title="New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1",
            year=2005,
            source="BIT Numerical Mathematics, 45(4), 761-787",
            doi="10.1007/s10543-005-0035-y",
        )],
    ),
    "Rodas5": dict(
        id="rodas5",
        name="RODAS5",
        order=5,
        embedded_order=4,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Di Marzo, G.",
            title="RODAS5(4) - Methodes de Rosenbrock d'ordre 5(4) adaptees aux problemes differentiels-algebriques",
            year=1993,
            source="MSc mathematics thesis, Faculty of Science, University of Geneva",
        )],
    ),
    "Rodas42": dict(
        id="rodas42",
        name="RODAS4.2",
        order=4,
        embedded_order=3,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[dict(
            authors="Hairer, E., & Wanner, G.",
            title="Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
            year=1996,
            source="Springer Series in Computational Mathematics, Vol. 14, 2nd edition, IV.7",
            doi="10.1007/978-3-642-05221-7",
        )],
    ),
    "RosShamp4": dict(
        id="rosshamp4",
        name="Shampine 4",
        order=4,
        embedded_order=3,
        properties={"a_stable": True},
        references=[dict(
            authors="Shampine, L. F.",
            title="Implementation of Rosenbrock methods",
            year=1982,
            source="ACM Transactions on Mathematical Software, 8(2), 93-113",
            doi="10.1145/355993.355994",
        )],
    ),
    "Veldd4": dict(
        id="veldd4",
        name="van Veldhuizen D 4",
        order=4,
        embedded_order=3,
        # D-stable rather than A-stable, which is what the paper is about.
        properties={},
        references=[dict(
            authors="van Veldhuizen, M.",
            title="D-stability and Kaps-Rentrop methods",
            year=1984,
            source="Computing, 32(3), 229-237",
            doi="10.1007/BF02243574",
        )],
    ),
    "Velds4": dict(
        id="velds4",
        name="van Veldhuizen S 4",
        order=4,
        embedded_order=3,
        properties={"a_stable": True},
        references=[dict(
            authors="van Veldhuizen, M.",
            title="D-stability and Kaps-Rentrop methods",
            year=1984,
            source="Computing, 32(3), 229-237",
            doi="10.1007/BF02243574",
        )],
    ),
    "ROS34PW1a": dict(
        id="ros34pw1a",
        name="ROS34PW1a",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG05],
    ),
    "ROS34PW1b": dict(
        id="ros34pw1b",
        name="ROS34PW1b",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG05],
    ),
    "ROS34PW3": dict(
        id="ros34pw3",
        name="ROS34PW3",
        order=4,
        # A W-method of order three, which with an exact Jacobian reaches four.
        # Its embedded pair is second order, not third.
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG05],
    ),
    "ROS2S": dict(
        id="ros2s",
        name="ROS2S",
        order=2,
        embedded_order=1,
        properties={"a_stable": True},
        references=[RANG15],
    ),
    "ROS3PRL": dict(
        id="ros3prl",
        name="ROS3PRL",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG15],
    ),
    "ROS2PR": dict(
        id="ros2pr",
        name="ROS2PR",
        order=2,
        embedded_order=1,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[RANG16],
    ),
    "ROS3PR": dict(
        id="ros3pr",
        name="ROS3PR",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG16],
    ),
    "ROS3PRL2": dict(
        id="ros3prl2",
        name="ROS3PRL2",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG16],
    ),
    "ROS34PRw": dict(
        id="ros34prw",
        name="ROS34PRw",
        order=3,
        embedded_order=2,
        properties={"a_stable": True},
        references=[RANG15],
    ),
    "Rodas4P2": dict(
        id="rodas4p2",
        name="RODAS4P2",
        order=4,
        embedded_order=3,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[STEINEBACH20],
    ),
    "Rodas3P": dict(
        id="rodas3p",
        name="RODAS3P",
        order=3,
        embedded_order=2,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[STEINEBACH24],
    ),
    "Rodas5P": dict(
        id="rodas5p",
        name="RODAS5P",
        order=5,
        embedded_order=4,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[STEINEBACH23],
    ),
    "Rodas5Pe": dict(
        id="rodas5pe",
        name="RODAS5Pe",
        order=5,
        embedded_order=4,
        properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
        references=[STEINEBACH24],
    ),
    "Rodas6P": dict(
        id="rodas6p",
        name="RODAS6P",
        order=6,
        embedded_order=5,
        properties={"a_stable": True, "l_stable": True},
        references=[STEINEBACH25],
    ),
}



# ---------------------------------------------------------------------------
# Translating one constructor
# ---------------------------------------------------------------------------


def balanced(text, start):
    """Index just past the parenthesis opened at `start`."""
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
    """`convert(T, expr)` carries a type, not a value. Drop the wrapper."""
    while True:
        match = re.search(r"convert\(\s*T2?\s*,", text)
        if not match:
            return text
        close = balanced(text, match.start() + len("convert"))
        text = text[: match.start()] + "(" + text[match.end() : close] + ")" + text[close + 1 :]


class Table:
    """An array that indexes from one, so the source can be run as it is
    written rather than rewritten index by index."""

    def __init__(self, data):
        self.data = np.asarray(data, dtype=float)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self.data[key[0] - 1, key[1] - 1]
        return self.data[key - 1]

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self.data[key[0] - 1, key[1] - 1] = value
        else:
            self.data[key - 1] = value

    def push(self, value):
        self.data = np.append(self.data, float(value))

    @property
    def shape(self):
        return self.data.shape


def as_array(value):
    return value.data if isinstance(value, Table) else np.asarray(value, dtype=float)


def translate(body):
    """Julia straight line arithmetic into the equivalent Python."""
    # A statement runs on while its brackets are open, so join before parsing.
    joined = []
    current = ""
    depth = 0
    for raw in body.splitlines():
        raw = raw.split("#")[0]
        current += " " + raw.strip()
        depth += raw.count("(") + raw.count("[") - raw.count(")") - raw.count("]")
        if depth <= 0:
            joined.append(current.strip())
            current = ""
            depth = 0
    if current.strip():
        joined.append(current.strip())

    lines = []
    for line in joined:
        if not line or line.startswith("return"):
            continue
        line = strip_converts(line)
        line = re.sub(r"\bzeros\(\s*T2?\s*,\s*", "Table(np.zeros(_shape(", line)
        line = re.sub(r"\bT2?\[", "Table([", line)
        line = re.sub(r"\bsize\(\s*(\w+)\s*,\s*1\s*\)", r"\1.shape[0]", line)
        line = re.sub(r"\bpush!\(\s*(\w+)\s*,", r"\1.push(", line)
        line = re.sub(r"\bzero\(\s*T2?\s*\)", "0.0", line)
        line = re.sub(r"\bone\(\s*T2?\s*\)", "1.0", line)
        line = re.sub(r"\binv\(", "1.0/(", line)
        # A Julia range in a comprehension, and its rational literals.
        line = re.sub(r"\b1:\((.+?)\)", r"range(1, (\1) + 1)", line)
        line = re.sub(r"//", "/", line)
        line = line.replace("^", "**")
        # Close the two constructors that were opened with an extra bracket.
        line = _close(line, "Table(np.zeros(_shape(", 3)
        line = _close(line, "Table([", 2)
        lines.append(line)
    return "\n".join(lines)


def _close(line, opener, closers):
    """`zeros(T, 3, 3)` and `T[...]` were opened with more brackets than they
    were written with; give back the ones that are missing."""
    if opener not in line:
        return line
    start = line.index(opener) + len(opener) - 1
    if not opener.endswith("(") and "]" not in line:
        return line
    end = balanced(line, start) if opener.endswith("(") else line.rindex("]")
    if opener.endswith("("):
        return line[:end] + ")" * (closers - 1) + line[end:]
    return line[: end + 1] + ")" + line[end + 1 :]


RETURN = re.compile(r"return\s+RodasTableau(?:\{[^(]*\})?\(([^)]*)\)")


def returned(body):
    """The eight coefficient slots, by the name each constructor passes them
    under. Several of them hand over a module level table directly rather than
    building a local one, so the names cannot be assumed."""
    match = RETURN.search(body)
    if not match:
        raise SystemExit("no RodasTableau return found")
    names = [part.strip() for part in match.group(1).split(",")]
    return dict(zip(["A", "C", "gamma", "c", "d", "H", "b", "btilde"], names))


def evaluate(body, constants):
    scope = dict(constants)
    scope.update({
        "np": np,
        "sqrt": np.sqrt,
        "Table": Table,
        "_shape": lambda *dims: tuple(dims),
    })
    exec(compile(translate(body), "<tableau>", "exec"), scope)
    return scope


def constants_of(text):
    """The module level coefficient tables, which several constructors share."""
    out = {}
    pattern = re.compile(r"^const (\w+) = (\[.*?\])\s*$", re.M | re.S)
    for match in pattern.finditer(text):
        inner = match.group(2).strip()[1:-1]
        if "," in inner:
            # A vector, written with separators.
            data = np.array([float(v) for v in inner.replace(",", " ").split()])
        else:
            # A matrix, whose rows are lines.
            data = np.array([
                [float(v) for v in row.split()] for row in inner.splitlines() if row.strip()
            ])
        out[match.group(1)] = Table(data)
    return out


def constructors(text):
    """Every `function NameRodasTableau(...) ... end` in the file."""
    out = {}
    pattern = re.compile(r"^function (\w+?)(?:Rodas)?Tableau\(::Type\{T\}[^\n]*\n(.*?)^end$", re.M | re.S)
    for match in pattern.finditer(text):
        out[match.group(1)] = match.group(2)
    return out


# ---------------------------------------------------------------------------
# Undoing the substitution
# ---------------------------------------------------------------------------


def invert_lower(m):
    """Invert a lower triangular matrix by forward substitution.

    A general inverse would do, but not exactly: the matrix here has a constant
    diagonal, so its inverse has to have one too, and a general factorization
    returns the same number six times with the last bits different. That is
    enough to make a stepper believe every stage needs its own factorization.
    """
    n = m.shape[0]
    inverse = np.zeros((n, n))
    for column in range(n):
        for row in range(column, n):
            acc = 1.0 if row == column else 0.0
            for k in range(column, row):
                acc -= m[row, k] * inverse[k, column]
            inverse[row, column] = acc / m[row, row]
    return inverse


def to_alpha_gamma(a, c_matrix, gamma, m, error_weights):
    """From the running coefficients back to `alpha`, `gamma` and `b`."""
    s = a.shape[0]
    inverse = np.diag(np.full(s, 1.0 / gamma)) - c_matrix
    big_gamma = invert_lower(inverse)
    alpha = a @ big_gamma
    b = big_gamma.T @ m
    b_embedded = None
    if error_weights is not None:
        # The reference stores `m - mhat`, so the embedded weights come back by
        # subtraction before the same transformation.
        b_embedded = big_gamma.T @ (m - error_weights)
    return alpha, big_gamma, b, b_embedded


def clean(value, tolerance=1e-13):
    """Snap the values that are exactly zero or exactly one in the source but
    arrived through a matrix inverse."""
    rounded = np.round(value)
    return np.where(np.abs(value - rounded) < tolerance, rounded, value)


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def number(value):
    if abs(value - round(value)) < 1e-15:
        return int(round(value))
    return round(float(value), 15)


def main():
    paths = sys.argv[1:] or DEFAULT_SOURCES
    text = "\n".join(io.open(path, encoding="utf-8").read() for path in paths)
    found = constructors(text)
    constants = constants_of(text)

    os.makedirs(OUT, exist_ok=True)
    written = []
    for name, entry in WANTED.items():
        if name not in found:
            print("  missing in the source:", name)
            continue
        scope = evaluate(found[name], constants)
        slots = returned(found[name])
        a = as_array(scope[slots["A"]])
        c_matrix = as_array(scope[slots["C"]])
        gamma = float(scope[slots["gamma"]])
        m = as_array(scope[slots["b"]])
        error_weights = scope.get(slots["btilde"])
        error_weights = None if error_weights is None else as_array(error_weights)

        # Some of the tables leave off the last column of `C`, which is zero by
        # construction. Square it up before the inverse.
        if c_matrix.shape[1] < a.shape[0]:
            padding = a.shape[0] - c_matrix.shape[1]
            c_matrix = np.hstack([c_matrix, np.zeros((c_matrix.shape[0], padding))])

        alpha, big_gamma, b, b_embedded = to_alpha_gamma(a, c_matrix, gamma, m, error_weights)
        alpha = clean(alpha)
        big_gamma = clean(big_gamma)

        # The two independent checks.
        stated_c = as_array(scope[slots["c"]])
        stated_d = as_array(scope[slots["d"]])
        residual_c = np.max(np.abs(alpha.sum(axis=1) - stated_c))
        residual_d = np.max(np.abs(big_gamma.sum(axis=1) - stated_d))
        if residual_c > 1e-10 or residual_d > 1e-10:
            raise SystemExit(
                "%s: the conversion disagrees with the source, c off by %.2e, d off by %.2e"
                % (name, residual_c, residual_d)
            )

        s = alpha.shape[0]
        method = {
            "id": entry["id"],
            "name": entry["name"],
            "class": "rosenbrock",
            "family": "rosenbrock",
            "order": entry["order"],
        }
        if entry.get("embedded_order") is not None and b_embedded is not None:
            method["embedded_order"] = entry["embedded_order"]
        method["properties"] = entry["properties"]
        method["rosenbrock"] = {
            "alpha": [[number(alpha[i, j]) for j in range(i)] for i in range(s)],
            "gamma": [[number(big_gamma[i, j]) for j in range(i + 1)] for i in range(s)],
            "b": [number(v) for v in b],
        }
        if b_embedded is not None:
            method["rosenbrock"]["b_embedded"] = [number(v) for v in b_embedded]
        method["references"] = entry["references"]

        path = os.path.join(OUT, entry["id"] + ".json")
        io.open(path, "w", encoding="utf-8").write(render(method))
        written.append((entry["id"], s, residual_c, residual_d))

    print("%d Rosenbrock method files written" % len(written))
    for identifier, stages, residual_c, residual_d in written:
        print("  %-12s %d stages   c %.1e   d %.1e" % (identifier, stages, residual_c, residual_d))


if __name__ == "__main__":
    main()
