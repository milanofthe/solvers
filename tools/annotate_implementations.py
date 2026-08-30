"""Record where each method can be found under some other name.

A reader who arrives knowing `ode45` and a reader who arrives knowing
Dormand-Prince have the same method in mind. Nothing on the way from one to the
other is written down in one place, and this is the one thing in the project
that cannot be derived from a tableau: it is a fact about other people's
software.

So it is entered by hand, and only where it is certain. A wrong entry here is
worse than a missing one, because it is the kind of claim nobody checks.

Usage:
    python tools/annotate_implementations.py
"""

import io
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = os.path.join(HERE, "..", "methods")
INNER_ARRAY = re.compile(r"\[[^\[\]{}]*?\]", re.S)

# library -> { method id: the name it goes by there }
KNOWN = {
    "SciPy": {
        "rkdp54": "RK45",
        "rkbs32": "RK23",
        "rkdp87": "DOP853",
        "radau_iia_5": "Radau",
    },
    "MATLAB": {
        "rkdp54": "ode45",
        "rkbs32": "ode23",
        "trapezoidal": "ode23t",
        "trbdf2": "ode23tb",
    },
    "OrdinaryDiffEq.jl": {
        "euler_explicit": "Euler",
        "euler_implicit": "ImplicitEuler",
        "midpoint_explicit": "Midpoint",
        "implicit_midpoint": "ImplicitMidpoint",
        "trapezoidal": "Trapezoid",
        "heun3": "Heun",
        "rk4": "RK4",
        "rkbs32": "BS3",
        "rkbs54": "BS5",
        "owrenzen3": "OwrenZen3",
        "owrenzen4": "OwrenZen4",
        "owrenzen5": "OwrenZen5",
        "ralston4": "Ralston4",
        "rkdp54": "DP5",
        "rkdp87": "DP8",
        "radau_iia_5": "RadauIIA5",
        "radau_iia_9": "RadauIIA9",
        "trbdf2": "TRBDF2",
        "ssprk22": "SSPRK22",
        "ssprk33": "SSPRK33",
        "ssprk34": "SSPRK43",
        "ros2": "ROS2",
        "ros3p": "ROS3P",
        "ros34pw1a": "ROS34PW1a",
        "ros34pw1b": "ROS34PW1b",
        "ros34pw2": "ROS34PW2",
        "ros34pw3": "ROS34PW3",
        "ros34prw": "ROS34PRw",
        "ros2pr": "ROS2PR",
        "ros3pr": "ROS3PR",
        "ros3prl2": "ROS3PRL2",
        "grk4a": "GRK4A",
        "grk4t": "GRK4T",
        "ros4lstab": "Ros4LStab",
        "rosshamp4": "RosShamp4",
        "veldd4": "Veldd4",
        "velds4": "Velds4",
        "rodas3": "Rodas3",
        "rodas3p": "Rodas3P",
        "rodas4": "Rodas4",
        "rodas42": "Rodas42",
        "rodas4p": "Rodas4P",
        "rodas4p2": "Rodas4P2",
        "rodas5": "Rodas5",
        "rodas5p": "Rodas5P",
        "rodas5pe": "Rodas5Pe",
        "rodas6p": "Rodas6P",
    },
    "GSL": {
        "rk4": "rk4",
        "rkf45": "rkf45",
        "rkck54": "rkck",
        "rkdp87": "rk8pd",
        "euler_implicit": "rk1imp",
        "implicit_midpoint": "rk2imp",
        "gauss_legendre_4": "rk4imp",
    },
    "Boost.Odeint": {
        "rk4": "runge_kutta4",
        "rkck54": "runge_kutta_cash_karp54",
        "rkdp54": "runge_kutta_dopri5",
        "rkf78": "runge_kutta_fehlberg78",
    },
}


def render(method):
    text = json.dumps(method, indent=2, ensure_ascii=False)
    text = INNER_ARRAY.sub(lambda m: "[" + " ".join(m.group(0)[1:-1].split()) + "]", text)
    return text + "\n"


def main():
    by_method = {}
    for library, mapping in KNOWN.items():
        for identifier, name in mapping.items():
            by_method.setdefault(identifier, []).append({"library": library, "name": name})

    seen = set()
    touched = 0
    for root, _, files in os.walk(METHODS):
        for filename in sorted(files):
            if not filename.endswith(".json"):
                continue
            path = os.path.join(root, filename)
            method = json.load(io.open(path, encoding="utf-8"))
            entries = by_method.get(method["id"])
            if entries is None:
                if "implementations" in method:
                    del method["implementations"]
                    io.open(path, "w", encoding="utf-8").write(render(method))
                continue
            seen.add(method["id"])

            # Rebuilt from the table every time, so removing a line here removes
            # it from the library rather than leaving it behind.
            ordered = {}
            for key, value in method.items():
                ordered[key] = value
                if key == "references":
                    ordered["implementations"] = entries
            if "implementations" not in ordered:
                ordered["implementations"] = entries
            io.open(path, "w", encoding="utf-8").write(render(ordered))
            touched += 1

    missing = sorted(set(by_method) - seen)
    if missing:
        raise SystemExit("no method file for: %s" % ", ".join(missing))
    print(f"{touched} methods annotated across {len(KNOWN)} libraries")


if __name__ == "__main__":
    main()
