"""Emit solvers method files from the tableaux already implemented in pathsim.

Transcribing Butcher tableaux by hand is the single most typo prone part of a
project like this. pathsim already carries verified tableaux, so they are read
programmatically and re-emitted in the solvers JSON schema. The order conditions
are then checked independently on the Rust side, which catches any mistake.
"""

import io
import json
import os
import re
import sys
from fractions import Fraction

sys.path.insert(0, r"C:\Repositories\pathsim\src")

import pathsim.solvers as ps
from pathsim.solvers._rungekutta import ExplicitRungeKutta, DiagonallyImplicitRungeKutta

OUT = r"C:\Repositories\solvers\methods"

FAMILIES = [
    ("esdirk", "esdirk"),
    ("dirk", "dirk"),
    ("ssprk", "ssprk"),
    ("rk", "erk"),
]

FAMILY_NAMES = {"esdirk": "esdirk", "dirk": "dirk", "ssprk": "ssprk", "rk": "erk"}

# id -> (family, directory, human name override)
NAME_OVERRIDES = {
    "rk4": "Classical Runge-Kutta",
    "rkbs32": "Bogacki-Shampine 3(2)",
    "rkck54": "Cash-Karp 5(4)",
    "rkdp54": "Dormand-Prince 5(4)",
    "rkdp87": "Dormand-Prince 8(7)",
    "rkf21": "Fehlberg 2(1)",
    "rkf45": "Fehlberg 4(5)",
    "rkf78": "Fehlberg 7(8)",
    "rkv65": "Verner 6(5)",
    "ssprk22": "SSPRK(2,2) Heun",
    "ssprk33": "SSPRK(3,3) Shu-Osher",
    "ssprk34": "SSPRK(3,4)",
    "dirk2": "DIRK 2(1)",
    "dirk3": "DIRK 3(2)",
    "esdirk32": "ESDIRK 3(2)",
    "esdirk4": "ESDIRK 4",
    "esdirk43": "ESDIRK 4(3)",
    "esdirk54": "ESDIRK 5(4)",
    "esdirk85": "ESDIRK 8(5)",
}


# Denominators above this are not plausibly the intended fraction, so a float
# that only rationalizes past it is left as a float. Below it, recovering the
# fraction makes the tableau exact and the order conditions become identities
# rather than a comparison against a tolerance.
MAX_DENOMINATOR = 10 ** 6


def rationalize(value):
    """The fraction a float came from, or None when it is not recoverable."""
    value = float(value)
    if value == 0.0:
        return Fraction(0)
    if value == int(value) and abs(value) < 2 ** 53:
        return Fraction(int(value))
    candidate = Fraction(value).limit_denominator(MAX_DENOMINATOR)
    return candidate if float(candidate) == value else None


def coefficient(value):
    """Emit an exact fraction when the float unambiguously came from one."""
    if isinstance(value, Fraction):
        if value.denominator == 1:
            return int(value)
        return f"{value.numerator}/{value.denominator}"
    exact = rationalize(value)
    if exact is None:
        return float(value)
    return coefficient(exact)


def coefficients(values):
    return [coefficient(v) for v in values]


def difference(b, tr):
    """`b - tr` computed exactly whenever both sides rationalize.

    Doing the subtraction in floating point first would leave the embedded
    weights one unit in the last place off their true fractions, which is
    enough to lose exactness for the whole tableau.
    """
    out = []
    for left, right in zip(b, tr):
        a = rationalize(left)
        c = rationalize(right)
        if a is not None and c is not None:
            out.append(a - c)
        else:
            out.append(float(left) - float(right))
    return out


def family_of(name):
    lowered = name.lower()
    for prefix, directory in FAMILIES:
        if lowered.startswith(prefix):
            return prefix, directory
    return "other", "other"


def parse_references(doc):
    """Pull the numpydoc References block apart into structured entries."""
    if not doc:
        return []
    match = re.search(r"References\s*\n\s*-+\s*\n(.*)", doc, re.S)
    if not match:
        return []
    body = match.group(1)
    entries = re.split(r"\n\s*\.\.\s*\[\d+\]\s*", "\n" + body.strip())
    references = []
    for entry in entries:
        entry = " ".join(entry.split())
        if not entry:
            continue
        reference = {}
        doi = re.search(r":doi:`([^`]+)`", entry)
        if doi:
            reference["doi"] = doi.group(1).strip()
            entry = entry[: doi.start()].strip()
        year = re.search(r"\((\d{4})[a-z]?\)", entry)
        if year:
            reference["year"] = int(year.group(1))
            reference["authors"] = entry[: year.start()].strip().rstrip(".,")
            rest = entry[year.end():].strip().lstrip(".").strip()
        else:
            rest = entry
        title = re.search(r'"([^"]+)"', rest)
        if title:
            reference["title"] = title.group(1).strip().rstrip(".")
            rest = (rest[: title.start()] + rest[title.end():]).strip()
        source = rest.strip().strip(".").strip(",").strip()
        if source:
            reference["source"] = source
        if reference:
            references.append(reference)
    return references


def summary(doc):
    if not doc:
        return None
    lines = []
    for line in doc.strip().splitlines():
        stripped = line.strip()
        if not stripped:
            break
        lines.append(stripped)
    text = " ".join(lines)
    return text or None


def tableau(instance):
    """Return (a, b, b_hat, c) as plain float lists."""
    s = instance.s
    c = [float(v) for v in instance.eval_stages]
    a = [[0.0] * s for _ in range(s)]

    if isinstance(instance, ExplicitRungeKutta):
        # BT[i] holds the combination applied after stage i, so it is row i+1
        # of A, and the last entry is b.
        for i in range(s - 1):
            row = instance.BT[i] or []
            for j, value in enumerate(row):
                a[i + 1][j] = float(value)
        b = [0.0] * s
        for j, value in enumerate(instance.BT[s - 1]):
            b[j] = float(value)
    elif isinstance(instance, DiagonallyImplicitRungeKutta):
        # BT[i] is row i of A including the diagonal.
        for i in range(s):
            row = instance.BT.get(i) or []
            for j, value in enumerate(row):
                a[i][j] = float(value)
        final = getattr(instance, "A", None)
        if final is not None:
            b = [float(v) for v in final]
        else:
            b = list(a[s - 1])
    else:
        raise TypeError(f"unsupported solver type {type(instance)}")

    b_hat = None
    if getattr(instance, "TR", None) is not None:
        b_hat = difference(b, [float(v) for v in instance.TR])
    return a, b, b_hat, c


def trim(rows):
    """Drop trailing zeros so the file shows the triangular shape."""
    out = []
    for i, row in enumerate(rows):
        last = 0
        for j, value in enumerate(row):
            if value != 0.0:
                last = j + 1
        out.append(row[:max(last, min(i + 1, len(row)))])
    return out


INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)


def render(method):
    """Pretty print with each coefficient row kept on a single line."""
    text = json.dumps(method, indent=2, ensure_ascii=False)
    # Collapse the innermost arrays, which are exactly the coefficient rows.
    text = INNER_ARRAY.sub(
        lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]",
        text,
    )
    return text + "\n"


def main():
    written = 0
    for name in dir(ps):
        cls = getattr(ps, name)
        if not isinstance(cls, type):
            continue
        if not issubclass(cls, (ExplicitRungeKutta, DiagonallyImplicitRungeKutta)):
            continue
        if cls in (ExplicitRungeKutta, DiagonallyImplicitRungeKutta):
            continue
        instance = cls(0.0)
        if not instance.BT:
            continue

        a, b, b_hat, c = tableau(instance)
        identifier = name.lower()
        prefix, directory = family_of(name)

        method = {
            "id": identifier,
            "name": NAME_OVERRIDES.get(identifier, name),
            "class": "runge_kutta",
            "family": FAMILY_NAMES.get(prefix, prefix),
            "order": int(instance.n),
        }
        if b_hat is not None and int(instance.m) > 0:
            method["embedded_order"] = int(instance.m)
        description = summary(cls.__doc__)
        if description:
            method["description"] = description
        method["tableau"] = {
            "a": [coefficients(row) for row in trim(a)],
            "b": coefficients(b),
            "c": coefficients(c),
        }
        if b_hat is not None:
            method["tableau"]["b_embedded"] = coefficients(b_hat)
        references = parse_references(cls.__doc__)
        if references:
            method["references"] = references

        folder = os.path.join(OUT, directory)
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, identifier + ".json")
        with io.open(path, "w", encoding="utf-8") as handle:
            handle.write(render(method))
        written += 1
        print(f"{identifier:<10} s={instance.s} order={instance.n} -> {directory}/")

    print(f"\n{written} method files written")


if __name__ == "__main__":
    main()
