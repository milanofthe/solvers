"""Write the classical method families that are stated in closed form.

Everything here is either a textbook tableau small enough to write out exactly,
or a multistep family that is fully described by which of its coefficients are
free. The order of every one of them is checked independently by the Rust side,
so a slip shows up as a failing verification rather than as a silent bug.
"""

import io
import json
import os
import re

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "methods")

INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

HNW1 = {
    "authors": "Hairer, E., Noersett, S. P., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations I: Nonstiff Problems",
    "year": 1993,
    "source": "Springer Series in Computational Mathematics, Vol. 8, 2nd edition",
    "doi": "10.1007/978-3-540-78862-1",
}
HW2 = {
    "authors": "Hairer, E., & Wanner, G.",
    "title": "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems",
    "year": 1996,
    "source": "Springer Series in Computational Mathematics, Vol. 14, 2nd edition",
    "doi": "10.1007/978-3-642-05221-7",
}
BUTCHER64 = {
    "authors": "Butcher, J. C.",
    "title": "Implicit Runge-Kutta processes",
    "year": 1964,
    "source": "Mathematics of Computation, 18(85), 50-64",
    "doi": "10.1090/S0025-5718-1964-0159424-9",
}
BUTCHER16 = {
    "authors": "Butcher, J. C.",
    "title": "Numerical Methods for Ordinary Differential Equations",
    "year": 2016,
    "source": "Wiley, 3rd edition",
    "doi": "10.1002/9781119121534",
}
CURTISS52 = {
    "authors": "Curtiss, C. F., & Hirschfelder, J. O.",
    "title": "Integration of stiff equations",
    "year": 1952,
    "source": "Proceedings of the National Academy of Sciences, 38(3), 235-243",
    "doi": "10.1073/pnas.38.3.235",
}
DAHLQUIST63 = {
    "authors": "Dahlquist, G.",
    "title": "A special stability problem for linear multistep methods",
    "year": 1963,
    "source": "BIT Numerical Mathematics, 3(1), 27-43",
    "doi": "10.1007/BF01963532",
}
DAHLQUIST56 = {
    "authors": "Dahlquist, G.",
    "title": "Convergence and stability in the numerical integration of ordinary differential equations",
    "year": 1956,
    "source": "Mathematica Scandinavica, 4, 33-53",
    "doi": "10.7146/math.scand.a-10454",
}
NORSETT74 = {
    "authors": "Noersett, S. P.",
    "title": "Semi-explicit Runge-Kutta methods",
    "year": 1974,
    "source": "Report Mathematics and Computation No. 6/74, University of Trondheim",
}
CROUZEIX75 = {
    "authors": "Crouzeix, M.",
    "title": "Sur l'approximation des equations differentielles operationnelles lineaires par des methodes de Runge-Kutta",
    "year": 1975,
    "source": "These d'Etat, Universite Paris VI",
}
EHLE69 = {
    "authors": "Ehle, B. L.",
    "title": "On Pade approximations to the exponential function and A-stable methods for the numerical solution of initial value problems",
    "year": 1969,
    "source": "Research Report CSRR 2010, University of Waterloo",
}


def _power(theta, degree):
    """`theta^degree` with the `0^0 = 1` convention the conditions need."""
    if degree < 0:
        return 0.0
    return 1.0 if degree == 0 else theta**degree


def lmm_coefficients(alpha_slots, beta_slots):
    """Solve the exactness conditions on a uniform grid, the same way the core
    solves them on a real step size history."""
    k = len(alpha_slots) - 1
    theta = [0.0] + [-float(j) for j in range(1, k + 1)]
    unknowns = [("a", j) for j, slot in enumerate(alpha_slots) if slot == FREE]
    unknowns += [("b", j) for j, slot in enumerate(beta_slots) if slot == FREE]

    rows, rhs, degree = [], [], 0
    while len(unknowns) and len(rows) < len(unknowns) and degree <= 2 * (k + 2):
        row = [
            _power(theta[j], degree) if kind == "a" else -degree * _power(theta[j], degree - 1)
            for kind, j in unknowns
        ]
        fixed = sum(
            float(slot) * _power(theta[j], degree)
            for j, slot in enumerate(alpha_slots)
            if slot != FREE
        )
        fixed -= sum(
            float(slot) * degree * _power(theta[j], degree - 1)
            for j, slot in enumerate(beta_slots)
            if slot != FREE
        )
        if max(abs(v) for v in row) > 1e-14:
            rows.append(row)
            rhs.append(-fixed)
        degree += 1

    alpha = [0.0 if slot == FREE else float(slot) for slot in alpha_slots]
    beta = [0.0 if slot == FREE else float(slot) for slot in beta_slots]
    if unknowns:
        solution = np.linalg.solve(np.array(rows), np.array(rhs))
        for value, (kind, j) in zip(solution, unknowns):
            if kind == "a":
                alpha[j] = value
            else:
                beta[j] = value

    if abs(alpha[0]) > 1e-300:
        scale = alpha[0]
        alpha = [v / scale for v in alpha]
        beta = [v / scale for v in beta]
    return alpha, beta


def lmm_order(alpha_slots, beta_slots):
    """The highest degree the formula integrates exactly."""
    alpha, beta = lmm_coefficients(alpha_slots, beta_slots)
    k = len(alpha) - 1
    theta = [0.0] + [-float(j) for j in range(1, k + 1)]
    magnitude = max(max(abs(v) for v in alpha), max(abs(v) for v in beta), 1.0)
    order = 0
    for degree in range(0, 2 * k + 5):
        residual = sum(alpha[j] * _power(theta[j], degree) for j in range(k + 1))
        residual -= sum(beta[j] * degree * _power(theta[j], degree - 1) for j in range(k + 1))
        if abs(residual) > 1e-9 * magnitude:
            return order
        order = degree
    return order


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def write(directory, method):
    folder = os.path.join(OUT, directory)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, method["id"] + ".json")
    with io.open(path, "w", encoding="utf-8") as handle:
        handle.write(render(method))
    return method["id"]


def rk(identifier, name, family, order, a, b, c, **extra):
    method = {
        "id": identifier,
        "name": name,
        "class": "runge_kutta",
        "family": family,
        "order": order,
    }
    for key in ("aliases", "tags", "embedded_order", "properties", "defaults"):
        if key in extra:
            method[key] = extra.pop(key)
    tableau = {"a": a, "b": b, "c": c}
    if "b_embedded" in extra:
        tableau["b_embedded"] = extra.pop("b_embedded")
    method["tableau"] = tableau
    if "references" in extra:
        method["references"] = extra.pop("references")
    assert not extra, extra
    return method


def lmm(identifier, name, family, order, steps, alpha, beta, **extra):
    derived = lmm_order(alpha, beta)
    if order is not None and order != derived:
        raise SystemExit(f"{identifier}: claims order {order} but the pattern gives {derived}")
    method = {
        "id": identifier,
        "name": name,
        "class": "linear_multistep",
        "family": family,
        "order": derived,
    }
    for key in ("aliases", "tags", "properties", "defaults"):
        if key in extra:
            method[key] = extra.pop(key)
    method["multistep"] = {
        "steps": steps,
        "alpha": alpha,
        "beta": beta,
        "normalization": "alpha0",
    }
    for key in ("startup", "min_steps"):
        if key in extra:
            method["multistep"][key] = extra.pop(key)
    if "references" in extra:
        method["references"] = extra.pop("references")
    assert not extra, extra
    return method


FREE = "free"

written = []

# ---------------------------------------------------------------------------
# Explicit Runge-Kutta
# ---------------------------------------------------------------------------

written.append(write("erk", rk(
    "euler_explicit", "Forward Euler", "erk", 1,
    a=[[0]], b=[1], c=[0],
    aliases=["euler", "forward_euler", "explicit_euler"],
    references=[HNW1],
)))

written.append(write("erk", rk(
    "midpoint_explicit", "Explicit midpoint", "erk", 2,
    a=[[0], ["1/2", 0]], b=[0, 1], c=[0, "1/2"],
    aliases=["explicit_midpoint", "rk2_midpoint"],
    references=[HNW1],
)))

written.append(write("erk", rk(
    "ralston2", "Ralston 2", "erk", 2,
    a=[[0], ["2/3", 0]], b=["1/4", "3/4"], c=[0, "2/3"],
    references=[{
        "authors": "Ralston, A.",
        "title": "Runge-Kutta methods with minimum error bounds",
        "year": 1962,
        "source": "Mathematics of Computation, 16(80), 431-437",
        "doi": "10.1090/S0025-5718-1962-0150954-0",
    }],
)))

written.append(write("erk", rk(
    "heun3", "Heun 3", "erk", 3,
    a=[[0], ["1/3", 0], [0, "2/3", 0]], b=["1/4", 0, "3/4"], c=[0, "1/3", "2/3"],
    references=[HNW1],
)))

written.append(write("erk", rk(
    "rk38", "Runge-Kutta 3/8 rule", "erk", 4,
    a=[[0], ["1/3", 0], ["-1/3", 1, 0], [1, -1, 1, 0]],
    b=["1/8", "3/8", "3/8", "1/8"], c=[0, "1/3", "2/3", 1],
    aliases=["rk4_38"],
    references=[HNW1],
)))

# ---------------------------------------------------------------------------
# Diagonally implicit
# ---------------------------------------------------------------------------

written.append(write("dirk", rk(
    "euler_implicit", "Backward Euler", "dirk", 1,
    a=[[1]], b=[1], c=[1],
    aliases=["backward_euler", "implicit_euler"],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True, "stage_order": 1},
    references=[HW2],
)))

written.append(write("dirk", rk(
    "implicit_midpoint", "Implicit midpoint", "dirk", 2,
    a=[["1/2"]], b=[1], c=["1/2"],
    aliases=["midpoint_implicit", "gauss_legendre_2"],
    properties={"a_stable": True, "l_stable": False, "stiffly_accurate": False,
                "symplectic": True, "symmetric": True, "stage_order": 1},
    references=[HW2],
)))

written.append(write("dirk", rk(
    "trapezoidal", "Trapezoidal rule", "dirk", 2,
    a=[[0, 0], ["1/2", "1/2"]], b=["1/2", "1/2"], c=[0, 1],
    aliases=["crank_nicolson", "trapezoid"],
    properties={"a_stable": True, "l_stable": False, "stiffly_accurate": True, "symmetric": True},
    references=[HW2],
)))

written.append(write("dirk", rk(
    "sdirk2", "SDIRK 2", "dirk", 2,
    a=[["1 - sqrt(2)/2", 0], ["sqrt(2)/2", "1 - sqrt(2)/2"]],
    b=["sqrt(2)/2", "1 - sqrt(2)/2"],
    c=["1 - sqrt(2)/2", 1],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
    references=[EHLE69, HW2],
)))

# Noersett's L-stable three stage SDIRK of order three. gamma is the middle root
# of x^3 - 3x^2 + 3x/2 - 1/6, given here to full double precision.
GAMMA3 = 0.43586652150845899941601945
written.append(write("dirk", rk(
    "sdirk3", "SDIRK 3", "dirk", 3,
    a=[[GAMMA3, 0, 0],
       [(1 - GAMMA3) / 2, GAMMA3, 0],
       [-(6 * GAMMA3 ** 2 - 16 * GAMMA3 + 1) / 4, (6 * GAMMA3 ** 2 - 20 * GAMMA3 + 5) / 4, GAMMA3]],
    b=[-(6 * GAMMA3 ** 2 - 16 * GAMMA3 + 1) / 4, (6 * GAMMA3 ** 2 - 20 * GAMMA3 + 5) / 4, GAMMA3],
    c=[GAMMA3, (1 + GAMMA3) / 2, 1],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
    references=[NORSETT74, CROUZEIX75, HW2],
)))

# ---------------------------------------------------------------------------
# Fully implicit: Radau, Gauss, Lobatto
# ---------------------------------------------------------------------------

written.append(write("irk", rk(
    "radau_iia_3", "Radau IIA 3", "irk", 3,
    a=[["5/12", "-1/12"], ["3/4", "1/4"]],
    b=["3/4", "1/4"], c=["1/3", 1],
    aliases=["radau3", "radau_iia_2s"],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True, "stage_order": 2},
    references=[BUTCHER64, HW2],
)))

written.append(write("irk", rk(
    "radau_iia_5", "Radau IIA 5", "irk", 5,
    a=[["(88 - 7*sqrt(6))/360", "(296 - 169*sqrt(6))/1800", "(-2 + 3*sqrt(6))/225"],
       ["(296 + 169*sqrt(6))/1800", "(88 + 7*sqrt(6))/360", "(-2 - 3*sqrt(6))/225"],
       ["(16 - sqrt(6))/36", "(16 + sqrt(6))/36", "1/9"]],
    b=["(16 - sqrt(6))/36", "(16 + sqrt(6))/36", "1/9"],
    c=["(4 - sqrt(6))/10", "(4 + sqrt(6))/10", 1],
    aliases=["radau5", "radau_iia_3s"],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True, "stage_order": 3},
    references=[BUTCHER64, HW2],
)))

written.append(write("irk", rk(
    "gauss_legendre_4", "Gauss-Legendre 4", "irk", 4,
    a=[["1/4", "1/4 - sqrt(3)/6"], ["1/4 + sqrt(3)/6", "1/4"]],
    b=["1/2", "1/2"],
    c=["1/2 - sqrt(3)/6", "1/2 + sqrt(3)/6"],
    aliases=["gauss4"],
    properties={"a_stable": True, "l_stable": False, "symplectic": True, "symmetric": True, "stage_order": 2},
    references=[BUTCHER64, HW2],
)))

written.append(write("irk", rk(
    "gauss_legendre_6", "Gauss-Legendre 6", "irk", 6,
    a=[["5/36", "2/9 - sqrt(15)/15", "5/36 - sqrt(15)/30"],
       ["5/36 + sqrt(15)/24", "2/9", "5/36 - sqrt(15)/24"],
       ["5/36 + sqrt(15)/30", "2/9 + sqrt(15)/15", "5/36"]],
    b=["5/18", "4/9", "5/18"],
    c=["1/2 - sqrt(15)/10", "1/2", "1/2 + sqrt(15)/10"],
    aliases=["gauss6"],
    properties={"a_stable": True, "l_stable": False, "symplectic": True, "symmetric": True, "stage_order": 3},
    references=[BUTCHER64, HW2],
)))

written.append(write("irk", rk(
    "lobatto_iiia_4", "Lobatto IIIA 4", "irk", 4,
    a=[[0, 0, 0], ["5/24", "1/3", "-1/24"], ["1/6", "2/3", "1/6"]],
    b=["1/6", "2/3", "1/6"], c=[0, "1/2", 1],
    properties={"a_stable": True, "l_stable": False, "stiffly_accurate": True, "symmetric": True, "stage_order": 3},
    references=[BUTCHER16, HW2],
)))

written.append(write("irk", rk(
    "lobatto_iiib_4", "Lobatto IIIB 4", "irk", 4,
    a=[["1/6", "-1/6", 0], ["1/6", "1/3", 0], ["1/6", "5/6", 0]],
    b=["1/6", "2/3", "1/6"], c=[0, "1/2", 1],
    properties={"a_stable": True, "l_stable": False, "symmetric": True},
    references=[BUTCHER16, HW2],
)))

written.append(write("irk", rk(
    "lobatto_iiic_2", "Lobatto IIIC 2", "irk", 2,
    a=[["1/2", "-1/2"], ["1/2", "1/2"]],
    b=["1/2", "1/2"], c=[0, 1],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True, "stage_order": 1},
    references=[BUTCHER16, HW2],
)))

written.append(write("irk", rk(
    "lobatto_iiic_4", "Lobatto IIIC 4", "irk", 4,
    a=[["1/6", "-1/3", "1/6"], ["1/6", "5/12", "-1/12"], ["1/6", "2/3", "1/6"]],
    b=["1/6", "2/3", "1/6"], c=[0, "1/2", 1],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True, "stage_order": 2},
    references=[BUTCHER16, HW2],
)))

# ---------------------------------------------------------------------------
# Backward differentiation formulas
# ---------------------------------------------------------------------------

BDF_NOTES = {
    1: "Identical to backward Euler, and the order the variable order driver falls back to.",
    2: "A-stable and the highest order BDF that is. The workhorse for stiff problems where robustness matters more than order.",
    3: "A(86.03) stable. The usual default of a variable order BDF code.",
    4: "A(73.35) stable.",
    5: "A(51.84) stable. Beyond this the stability wedge closes quickly.",
    6: "A(17.84) stable, and the last BDF that is zero stable at all.",
}

# A k step formula needs k history points, and they have to come from a one
# step method of at least the same order and on the same step size. Generating
# them by ramping the step size instead leaves the history so unevenly spaced
# that the variable step coefficients lose all conditioning at high order.
# Stiff families get a stiff starter.
BDF_STARTUP = {1: "esdirk32", 2: "esdirk32", 3: "esdirk43", 4: "esdirk43", 5: "esdirk54", 6: "esdirk85"}

for k in range(1, 7):
    written.append(write("bdf", lmm(
        f"bdf{k}", f"BDF{k}", "bdf", k, k,
        alpha=[FREE] * (k + 1),
        beta=[1] + [0] * k,
        startup=BDF_STARTUP[k],
        aliases=[f"gear{k}"] if k > 1 else ["gear1"],
        properties={"a_stable": k <= 2},
        references=[CURTISS52, DAHLQUIST63, HW2],
    )))

# ---------------------------------------------------------------------------
# Adams families
# ---------------------------------------------------------------------------

for k in range(1, 9):
    written.append(write("adams", lmm(
        f"adams_bashforth_{k}", f"Adams-Bashforth {k}", "adams_bashforth", k, k,
        alpha=[1, -1] + [0] * (k - 1),
        beta=[0] + [FREE] * k,
        startup="rkdp54" if k <= 5 else "rkdp87",
        aliases=[f"ab{k}"],
        references=[HNW1, DAHLQUIST56],
    )))

for k in range(1, 7):
    written.append(write("adams", lmm(
        f"adams_moulton_{k + 1}", f"Adams-Moulton {k + 1}", "adams_moulton", k + 1, k,
        alpha=[1, -1] + [0] * (k - 1),
        beta=[FREE] * (k + 1),
        startup="esdirk54" if k <= 4 else "esdirk85",
        aliases=[f"am{k + 1}"],
        references=[HNW1, DAHLQUIST56],
    )))

# ---------------------------------------------------------------------------
# Nystrom and Milne-Simpson: the same machinery with a different free pattern
# ---------------------------------------------------------------------------

written.append(write("adams", lmm(
    "nystrom2", "Nystrom 2", "nystrom", 2, 2,
    alpha=[1, 0, -1], beta=[0, FREE, FREE],
    startup="rkdp54",
    min_steps=2,
    aliases=["explicit_midpoint_lmm", "leapfrog"],
    references=[HNW1],
)))

written.append(write("adams", lmm(
    "milne_simpson", "Milne-Simpson", "milne_simpson", 4, 2,
    alpha=[1, 0, -1], beta=[FREE, FREE, FREE],
    startup="rkdp54",
    min_steps=2,
    aliases=["simpson_lmm"],
    references=[HNW1, DAHLQUIST56],
)))

# ---------------------------------------------------------------------------
# Further explicit methods
# ---------------------------------------------------------------------------

written.append(write("erk", rk(
    "kutta3", "Kutta 3", "erk", 3,
    a=[[0], ["1/2", 0], [-1, 2, 0]], b=["1/6", "2/3", "1/6"], c=[0, "1/2", 1],
    aliases=["rk3"],
    references=[HNW1],
)))

written.append(write("erk", rk(
    "ralston3", "Ralston 3", "erk", 3,
    a=[[0], ["1/2", 0], [0, "3/4", 0]], b=["2/9", "1/3", "4/9"], c=[0, "1/2", "3/4"],
    references=[{
        "authors": "Ralston, A.",
        "title": "Runge-Kutta methods with minimum error bounds",
        "year": 1962,
        "source": "Mathematics of Computation, 16(80), 431-437",
        "doi": "10.1090/S0025-5718-1962-0150954-0",
    }],
)))

written.append(write("erk", rk(
    "heun_euler21", "Heun-Euler 2(1)", "erk", 2,
    a=[[0], [1, 0]], b=["1/2", "1/2"], c=[0, 1],
    embedded_order=1,
    b_embedded=[1, 0],
    references=[HNW1],
)))

# ---------------------------------------------------------------------------
# Further diagonally implicit methods
# ---------------------------------------------------------------------------

written.append(write("dirk", rk(
    "crouzeix3", "Crouzeix 3", "dirk", 3,
    a=[["1/2 + sqrt(3)/6", 0], ["-sqrt(3)/3", "1/2 + sqrt(3)/6"]],
    b=["1/2", "1/2"],
    c=["1/2 + sqrt(3)/6", "1/2 - sqrt(3)/6"],
    properties={"a_stable": True, "l_stable": False},
    references=[CROUZEIX75, HW2],
)))

written.append(write("dirk", rk(
    "qin_zhang2", "Qin-Zhang 2", "dirk", 2,
    a=[["1/4", 0], ["1/2", "1/4"]], b=["1/2", "1/2"], c=["1/4", "3/4"],
    aliases=["symplectic_dirk2"],
    properties={"a_stable": True, "l_stable": False, "symplectic": True},
    references=[{
        "authors": "Qin, M. Z., & Zhang, M. Q.",
        "title": "Symplectic Runge-Kutta algorithms for Hamiltonian systems",
        "year": 1992,
        "source": "Journal of Computational Mathematics, Supplementary Issue, 205-215",
    }],
)))

written.append(write("esdirk", rk(
    "trbdf2", "TR-BDF2", "esdirk", 2,
    a=[
        [0, 0, 0],
        ["1 - sqrt(2)/2", "1 - sqrt(2)/2", 0],
        ["sqrt(2)/4", "sqrt(2)/4", "1 - sqrt(2)/2"],
    ],
    b=["sqrt(2)/4", "sqrt(2)/4", "1 - sqrt(2)/2"],
    c=[0, "2 - sqrt(2)", 1],
    aliases=["ode23tb"],
    properties={"a_stable": True, "l_stable": True, "stiffly_accurate": True},
    references=[{
        "authors": "Hosea, M. E., & Shampine, L. F.",
        "title": "Analysis and implementation of TR-BDF2",
        "year": 1996,
        "source": "Applied Numerical Mathematics, 20(1-2), 21-37",
        "doi": "10.1016/0168-9274(95)00115-8",
    }],
)))

# ---------------------------------------------------------------------------
# Further multistep methods
# ---------------------------------------------------------------------------

for k in (3, 4):
    written.append(write("adams", lmm(
        f"nystrom{k}", f"Nystrom {k}", "nystrom", None, k,
        alpha=[1, 0, -1] + [0] * (k - 2),
        beta=[0] + [FREE] * k,
        startup="rkdp54",
        min_steps=2,
        references=[HNW1],
    )))

print(f"{len(written)} classical method files written")
for identifier in written:
    print(" ", identifier)
