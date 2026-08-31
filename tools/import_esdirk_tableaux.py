"""Import the ESDIRK and SDIRK tableaux from a reference implementation.

These are stored in a different convention from the explicit ones, because the
family shares a structure worth exploiting: one diagonal entry for every
implicit stage, so one factorization serves them all. The reference therefore
states the diagonal once, as a single number, and only the strictly lower part
row by row:

  * `a_ii = gamma` for every implicit stage, and `a_11 = 0` where the first
    stage is explicit, which is what the E in ESDIRK means.
  * `a_21` is not stated. It follows from the abscissa: `c_2 = a_21 + gamma`.
  * `b` is not stated either. These methods are stiffly accurate, so the last
    row of `A` is the weight vector, which is the whole point of the design.
  * `btilde` is the difference `b - bhat`, so the embedded weights are recovered
    by subtracting it.

Every one of those is an assumption about somebody else's file, so none of them
is trusted. What comes out is read back by the Rust side, which derives the
order and the stage structure from the coefficients; a tableau that does not
satisfy the order claimed for it is reported here and no file is written.

The reference states the rest of its diagonally implicit methods under other
conventions again: KenCarp4, KenCarp5, KenCarp47, KenCarp58, the ESDIRK L[2]SA
family, Cash4 and Hairer42 come out of this reader with a last row that is not
their weight vector, so they are not imported here rather than imported wrong.
Their weights would have to be read from the struct each one builds.

Requires the reference checkout:
    git clone https://github.com/SciML/OrdinaryDiffEq.jl C:/Repositories/TEMP/OrdinaryDiffEq.jl

Usage:
    python tools/import_esdirk_tableaux.py
"""

import importlib.util
import io
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "methods")
SOURCE = "OrdinaryDiffEqSDIRK/src/sdirk_tableaus.jl"
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

# The Julia reader lives in the explicit importer; only the assembly differs.
_spec = importlib.util.spec_from_file_location(
    "rk_import", os.path.join(HERE, "import_runge_kutta_tableaux.py")
)
rk_import = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rk_import)

KVAERNO04 = {
    "authors": "Kvaerno, A.",
    "title": "Singly diagonally implicit Runge-Kutta methods with an explicit first stage",
    "year": 2004,
    "source": "BIT Numerical Mathematics, 44(3), 489-502",
    "doi": "10.1023/B:BITN.0000046811.70614.38",
}
KENNEDY_CARPENTER03 = {
    "authors": "Kennedy, C. A., & Carpenter, M. H.",
    "title": "Additive Runge-Kutta schemes for convection-diffusion-reaction equations",
    "year": 2003,
    "source": "Applied Numerical Mathematics, 44(1-2), 139-181",
    "doi": "10.1016/S0168-9274(02)00138-1",
}
KENNEDY_CARPENTER19 = {
    "authors": "Kennedy, C. A., & Carpenter, M. H.",
    "title": "Higher-order additive Runge-Kutta schemes for ordinary differential equations",
    "year": 2019,
    "source": "Applied Numerical Mathematics, 136, 183-205",
    "doi": "10.1016/j.apnum.2018.10.007",
}
KENNEDY_CARPENTER16 = {
    "authors": "Kennedy, C. A., & Carpenter, M. H.",
    "title": "Diagonally implicit Runge-Kutta methods for ordinary differential equations. A review",
    "year": 2016,
    "source": "NASA/TM-2016-219173",
}
CASH79 = {
    "authors": "Cash, J. R.",
    "title": "Diagonally implicit Runge-Kutta formulae with error estimates",
    "year": 1979,
    "source": "IMA Journal of Applied Mathematics, 24(3), 293-301",
    "doi": "10.1093/imamat/24.3.293",
}
HW96 = {
    "authors": "Hairer, E., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
    "year": 1996,
    "source": "Springer Series in Computational Mathematics, Vol. 14, 2nd edition",
    "doi": "10.1007/978-3-642-05221-7",
}

WANTED = [
    dict(constructor="Kvaerno3Tableau", id="kvaerno3", name="Kvaerno 3(2)", family="esdirk",
         order=3, embedded_order=2, references=[KVAERNO04]),
    dict(constructor="Kvaerno4Tableau", id="kvaerno4", name="Kvaerno 4(3)", family="esdirk",
         order=4, embedded_order=3, references=[KVAERNO04]),
    dict(constructor="Kvaerno5Tableau", id="kvaerno5", name="Kvaerno 5(4)", family="esdirk",
         order=5, embedded_order=4, references=[KVAERNO04]),
    dict(constructor="KenCarp3Tableau", id="kencarp3", name="Kennedy-Carpenter 3(2)",
         family="esdirk", order=3, embedded_order=2, references=[KENNEDY_CARPENTER03]),
]

ROW = re.compile(r"^a(\d)(\d)$")
ERROR = re.compile(r"^btilde(\d+)$")
NODE = re.compile(r"^c(\d+)$")


def assemble(scope):
    """The tableau, from the diagonal and the strictly lower part."""
    gamma = scope.get("γ")
    if gamma is None:
        return None
    lower, errors, nodes = {}, {}, {}
    for name, value in scope.items():
        if not isinstance(value, (int, float)):
            continue
        match = ROW.match(name)
        if match:
            i, j = int(match.group(1)), int(match.group(2))
            if j < i:
                lower[(i, j)] = value
            continue
        match = ERROR.match(name)
        if match:
            errors[int(match.group(1))] = value
            continue
        match = NODE.match(name)
        if match:
            nodes[int(match.group(1))] = value
    if not lower:
        return None

    stages = max(max(i for i, _ in lower), max(errors or [0]), max(nodes or [0]))
    a = [[0.0] * stages for _ in range(stages)]
    for (i, j), value in lower.items():
        if i > stages:
            return None
        a[i - 1][j - 1] = value
    # The first stage is explicit and the rest share the diagonal. The second
    # row is not stated because its abscissa fixes it.
    for i in range(1, stages):
        a[i][i] = gamma
    if stages >= 2 and a[1][0] == 0.0:
        a[1][0] = nodes.get(2, 2 * gamma) - gamma

    b = list(a[stages - 1])
    if not any(b):
        return None
    b_embedded = None
    if errors:
        difference = [errors.get(i + 1, 0.0) for i in range(stages)]
        b_embedded = [b[i] - difference[i] for i in range(stages)]
    return a, b, b_embedded


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    return INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text) + "\n"


def main():
    text = io.open(os.path.join(rk_import.LIB, SOURCE), encoding="utf-8").read()
    written, refused = [], []
    for entry in WANTED:
        table = None
        for body in reversed(rk_import.bodies(text, entry["constructor"])):
            scope = rk_import.evaluate(body)
            if scope is None:
                continue
            candidate = assemble(scope)
            if candidate and rk_import.integrates_to(candidate[0], candidate[1], entry["order"]):
                table = candidate
                break
        if not table:
            refused.append(entry["id"])
            continue
        a, b, b_embedded = table

        method = {
            "id": entry["id"],
            "name": entry["name"],
            "class": "runge_kutta",
            "family": entry["family"],
            "order": entry["order"],
        }
        if b_embedded and entry.get("embedded_order"):
            method["embedded_order"] = entry["embedded_order"]
        method["tableau"] = {
            "a": [[rk_import.coefficient(v) for v in row] for row in a],
            "b": [rk_import.coefficient(v) for v in b],
        }
        if b_embedded:
            method["tableau"]["b_embedded"] = [rk_import.coefficient(v) for v in b_embedded]
        method["references"] = entry["references"]

        folder = os.path.join(OUT, entry["family"])
        os.makedirs(folder, exist_ok=True)
        io.open(os.path.join(folder, entry["id"] + ".json"), "w", encoding="utf-8").write(
            render(method)
        )
        written.append((entry["id"], len(b)))

    print(f"{len(written)} diagonally implicit method files written")
    for identifier, stages in written:
        print(f"  {identifier:<16} {stages:>2} stages")
    if refused:
        print("no definition read out of the source integrates to the claimed order:")
        for identifier in refused:
            print(" ", identifier)


if __name__ == "__main__":
    main()
