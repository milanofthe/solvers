"""Import IMEX pairs from a reference implementation.

An additive method is two tableaux, and the reference writes them as two
matrices filled entry by entry, which is the one convention in that file that
can be read without knowing anything about the method: `Ae` and `be` are the
explicit half, `Ai` and `bi` the implicit one, on shared abscissae `c`.

Nothing is assumed about how the halves fit together. The loader refuses a pair
whose halves disagree on their abscissae, and the Rust side derives the order of
the pair from the two coloured trees, which is a stronger statement than either
half satisfying its own conditions. A pair that does not reach the order claimed
for it is reported here and no file is written.

Requires the reference checkout:
    git clone https://github.com/SciML/OrdinaryDiffEq.jl C:/Repositories/TEMP/OrdinaryDiffEq.jl

Usage:
    python tools/import_imex_tableaux.py
"""

import importlib.util
import io
import json
import os
import re
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods")
SOURCE = "OrdinaryDiffEqSDIRK/src/imex_tableaus.jl"
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

_spec = importlib.util.spec_from_file_location(
    "rk_import", os.path.join(HERE, "import_runge_kutta_tableaux.py")
)
rk_import = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rk_import)

ARS97 = {
    "authors": "Ascher, U. M., Ruuth, S. J., & Spiteri, R. J.",
    "title": "Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations",
    "year": 1997,
    "source": "Applied Numerical Mathematics, 25(2-3), 151-167",
    "doi": "10.1016/S0168-9274(97)00056-1",
}
PARESCHI_RUSSO05 = {
    "authors": "Pareschi, L., & Russo, G.",
    "title": "Implicit-explicit Runge-Kutta schemes and applications to hyperbolic systems with relaxation",
    "year": 2005,
    "source": "Journal of Scientific Computing, 25(1), 129-155",
    "doi": "10.1007/s10915-004-4636-4",
}
CALVO01 = {
    "authors": "Calvo, M. P., de Frutos, J., & Novo, J.",
    "title": "Linearly implicit Runge-Kutta methods for advection-reaction-diffusion equations",
    "year": 2001,
    "source": "Applied Numerical Mathematics, 37(4), 535-549",
    "doi": "10.1016/S0168-9274(00)00061-1",
}

WANTED = [
    dict(constructor="ARS222Tableau", id="ars222", name="ARS(2,2,2)", order=2, references=[ARS97]),
    dict(constructor="ARS232Tableau", id="ars232", name="ARS(2,3,2)", order=2, references=[ARS97]),
    dict(constructor="ARS343Tableau", id="ars343", name="ARS(3,4,3)", order=3, references=[ARS97]),
    dict(constructor="ARS443Tableau", id="ars443", name="ARS(4,4,3)", order=3, references=[ARS97]),
    dict(constructor="IMEXSSP222Tableau", id="imexssp222", name="IMEX-SSP2(2,2,2)", order=2,
         references=[PARESCHI_RUSSO05]),
    dict(constructor="IMEXSSP2322Tableau", id="imexssp2322", name="IMEX-SSP2(3,2,2)", order=2,
         references=[PARESCHI_RUSSO05]),
    dict(constructor="IMEXSSP3332Tableau", id="imexssp3332", name="IMEX-SSP3(3,3,2)", order=2,
         references=[PARESCHI_RUSSO05]),
    dict(constructor="IMEXSSP3433Tableau", id="imexssp3433", name="IMEX-SSP3(4,3,3)", order=3,
         references=[PARESCHI_RUSSO05]),
    dict(constructor="CFNLIRK3Tableau", id="cfnlirk3", name="CFNLIRK3", order=3,
         references=[CALVO01]),
]

MATRIX = re.compile(r"^\s*(Ai|Ae)\[\s*(\d+)\s*,\s*(\d+)\s*\]\s*=\s*(.+?)\s*$")
VECTOR = re.compile(r"^\s*(bi|be)\[\s*(\d+)\s*\]\s*=\s*(.+?)\s*$")
# The weights are written either entry by entry or as one array literal.
LITERAL = re.compile(r"^\s*(bi|be)\s*=\s*T\[(.+)\]\s*$")
SIZE = re.compile(r"^\s*s\s*=\s*(\d+)\s*$")


def split_arguments(text):
    """The entries of an array literal, respecting the brackets inside them."""
    parts, depth, current = [], 0, ""
    for character in text:
        if character in "([":
            depth += 1
        elif character in ")]":
            depth -= 1
        if character == "," and depth == 0:
            parts.append(current)
            current = ""
            continue
        current += character
    if current.strip():
        parts.append(current)
    return parts


def value(expression, scope):
    """One right hand side, exactly where the source is exact."""
    text = rk_import.strip_converts(expression)
    text = re.sub(r"\b(?:one|zero)\(\w+\)", lambda m: "1" if "one" in m.group(0) else "0", text)
    text = re.sub(r"\bT2?\((\d+)\)", r"\1", text)
    text = re.sub(r"BigInt\((-?\d+)\)\s*//\s*BigInt\((-?\d+)\)", r"Fraction(\1, \2)", text)
    text = re.sub(r"(-?\d+)\s*//\s*(-?\d+)", r"Fraction(\1, \2)", text)
    text = text.replace("^", "**")
    if "//" in text:
        return None
    try:
        return eval(text, dict(scope))  # noqa: S307 - a coefficient, not a program
    except Exception:  # noqa: BLE001
        return None


def read(body):
    """The two matrices and the two weight vectors, or nothing."""
    import math

    scope = {"Fraction": Fraction, "sqrt": math.sqrt, "__builtins__": {}}
    size = None
    for line in body.splitlines():
        match = SIZE.match(line)
        if match:
            size = int(match.group(1))
    if size is None:
        return None

    a = {"Ai": [[0] * size for _ in range(size)], "Ae": [[0] * size for _ in range(size)]}
    b = {"bi": [0] * size, "be": [0] * size}

    # Scalars first, so the entries that refer to them can be evaluated.
    for line in body.splitlines():
        stripped = re.sub(r"#.*", "", line).strip()
        if not stripped or "=" not in stripped:
            continue
        name, _, rest = stripped.partition("=")
        name = name.strip()
        if name.isidentifier() and name not in ("s",):
            got = value(rest, scope)
            if got is not None:
                scope[name] = got

    for line in body.splitlines():
        stripped = re.sub(r"#.*", "", line)
        match = MATRIX.match(stripped)
        if match:
            got = value(match.group(4), scope)
            if got is None:
                return None
            a[match.group(1)][int(match.group(2)) - 1][int(match.group(3)) - 1] = got
            continue
        match = VECTOR.match(stripped)
        if match:
            got = value(match.group(3), scope)
            if got is None:
                return None
            b[match.group(1)][int(match.group(2)) - 1] = got
            continue
        match = LITERAL.match(stripped)
        if match:
            entries = split_arguments(match.group(2))
            if len(entries) != size:
                return None
            for index, entry in enumerate(entries):
                got = value(entry, scope)
                if got is None:
                    return None
                b[match.group(1)][index] = got

    if not any(any(row) for row in a["Ae"]) or not any(b["be"]):
        return None
    return a["Ae"], b["be"], a["Ai"], b["bi"]


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    return INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text) + "\n"


def half(a, b, explicit):
    """One tableau. `c` is left out so the loader derives it and the two halves
    are then checked against each other rather than against what is written."""
    size = len(b)
    rows = []
    for i in range(size):
        row = a[i][:i] if explicit else a[i]
        rows.append([rk_import.coefficient(v) for v in row])
    return {"a": rows, "b": [rk_import.coefficient(v) for v in b]}


def main():
    text = io.open(os.path.join(rk_import.LIB, SOURCE), encoding="utf-8").read()
    written, refused = [], []
    for entry in WANTED:
        table = None
        for body in rk_import.bodies(text, entry["constructor"]):
            table = read(body)
            if table:
                break
        if not table:
            refused.append(entry["id"])
            continue
        ae, be, ai, bi = table

        method = {
            "id": entry["id"],
            "name": entry["name"],
            "class": "additive_runge_kutta",
            "family": "imex",
            "order": entry["order"],
            "additive": {
                "explicit": half(ae, be, True),
                "implicit": half(ai, bi, False),
            },
            "references": entry["references"],
        }
        folder = os.path.join(OUT, "imex")
        os.makedirs(folder, exist_ok=True)
        io.open(os.path.join(folder, entry["id"] + ".json"), "w", encoding="utf-8").write(
            render(method)
        )
        written.append((entry["id"], len(be)))

    print(f"{len(written)} IMEX pairs written")
    for identifier, stages in written:
        print(f"  {identifier:<14} {stages:>2} stages")
    if refused:
        print("not read out of the source:")
        for identifier in refused:
            print(" ", identifier)


if __name__ == "__main__":
    main()
