# solvers

A unified framework for ordinary differential equation solvers, in Rust, with a
browser interface that derives every property of a method from its coefficients.

The organizing idea is that **a method is data**. A Butcher tableau or a
multistep coefficient pattern lives in a JSON file next to the DOI of the paper
it was published in; the code is generic over that data. What follows from it:

- one Runge-Kutta stepper covering explicit, diagonally implicit and fully
  implicit tableaux, dispatching on the sparsity of `A`;
- one multistep engine that solves for its variable step coefficients from the
  order conditions, so BDF, Adams-Bashforth, Adams-Moulton, Nystrom and
  Milne-Simpson are the same code with a different pattern of free coefficients;
- one set of nonlinear solvers and one parametrized error controller, shared by
  every implicit method;
- analyses that read the same files, so a plot cannot drift away from the method
  it claims to describe.

Because the properties are derived rather than declared, the library checks
itself: `solvers verify` compares every claim in every file against what the
coefficients actually satisfy.

## Contents

| | |
|---|---|
| `crates/solvers-core` | the framework: coefficients, linear algebra, steppers, controllers, analyses |
| `crates/solvers-cli` | command line access, `solvers` |
| `crates/solvers-wasm` | browser bindings |
| `methods/` | the method library, one JSON file per method |
| `web/` | the interface |

## The method library

54 methods across the classical families.

| family | members |
|---|---|
| explicit Runge-Kutta | Euler, midpoint, Ralston, Heun, classical RK4, 3/8 rule, Fehlberg 2(1)/4(5)/7(8), Bogacki-Shampine 3(2), Cash-Karp 5(4), Dormand-Prince 5(4)/8(7), Verner 6(5) |
| SSP Runge-Kutta | SSPRK(2,2), SSPRK(3,3), SSPRK(3,4) |
| diagonally implicit | backward Euler, implicit midpoint, trapezoidal, SDIRK 2 and 3, DIRK 2(1) and 3(2) |
| ESDIRK | 3(2), 4, 4(3), 5(4), 8(5) |
| fully implicit | Radau IIA 3 and 5, Gauss-Legendre 4 and 6, Lobatto IIIA/IIIB/IIIC |
| BDF | orders 1 to 6 |
| Adams | Bashforth 1 to 5, Moulton 2 to 5 |
| other multistep | Nystrom, Milne-Simpson |

A method file states the coefficients and, optionally, what the publication
claims about them. Coefficients may be exact fractions, arithmetic expressions,
or plain floats:

```json
{
  "id": "sdirk2",
  "class": "runge_kutta",
  "order": 2,
  "properties": { "a_stable": true, "l_stable": true },
  "tableau": {
    "a": [["1 - sqrt(2)/2", 0], ["sqrt(2)/2", "1 - sqrt(2)/2"]],
    "b": ["sqrt(2)/2", "1 - sqrt(2)/2"],
    "c": ["1 - sqrt(2)/2", 1]
  },
  "references": [{ "authors": "Ehle, B. L.", "year": 1969, "doi": "..." }]
}
```

Provenance decides exactness: a fraction written as a string is exact and the
order conditions are then verified as identities, while a float stays a float
and is checked numerically. The report says which of the two happened.

A multistep file does not store coefficients at all, only which of them are
free:

```json
{
  "id": "bdf3",
  "class": "linear_multistep",
  "multistep": {
    "steps": 3,
    "alpha": ["free", "free", "free", "free"],
    "beta": [1, 0, 0, 0],
    "startup": "esdirk43"
  }
}
```

The values are solved for at every step from the actual step size history, which
is what makes the variable step form a property of the framework rather than of
the individual method.

## What is derived, not declared

- **Order**, from the rooted tree conditions, in exact arithmetic where the
  tableau allows it. Stage order and abscissa consistency alongside it.
- **Stability**, from `R(z) = det(I - zA + z e b^T) / det(I - zA)`, obtained as a
  characteristic polynomial rather than sampled. That gives the poles, the value
  at infinity, an exact A-stability test through the E-polynomial, and the real
  and imaginary axis limits. For multistep methods the same questions are
  answered through the root condition, including the A(alpha) wedge angle.
- **Cost per step**, with FSAL reuse deducted.

The derived numbers reproduce the published ones: the A(alpha) angles of BDF3 to
BDF6 come out as 86.11, 73.64, 52.09 and 17.85 degrees against the tabulated
86.03, 73.35, 51.84 and 17.84, and RK4's stability interval as -2.785 on the
real axis and 2.828 on the imaginary axis.

## Notes on the implementation

A few decisions that are not obvious and that matter:

**The stiff error estimate is filtered.** An implicit method's raw embedded
difference inherits the stiff eigenvalues of the problem, so it reports an error
for modes the method is in fact damping perfectly. Passing it through
`(I - h*gamma*J)^{-1}` first is what lets the step size grow past the explicit
stability limit, which is the whole point of using an L-stable method.

**Fully implicit stage systems are diagonalized.** Writing `A^{-1} = T L T^{-1}`
turns one factorization of size `s*n` into one real and `(s-1)/2` complex ones of
size `n`. Measured on a semi discretized diffusion problem with `n = 80`, that is
3.9 times faster for Radau IIA 5 and 4.3 times for Gauss-Legendre 6, with results
identical to 1e-14. What remains is a small dense transform applied across long
vectors, which is the shape the vectorized kernels are built for.

**A method with no embedded pair still runs adaptively**, through step doubling.
It costs three steps instead of one, but it is what makes Gauss, Radau and
Lobatto usable with an error controller at all.

**Every error controller is the same digital filter** in Soederlind's form, with
different gains. I, PI, PID, the H211 and H312 families and Gustafsson's
predictive controller are parameter sets, not implementations.

## Verification

```
cargo test
solvers verify
```

- Every one of the 54 methods converges at its stated order, measured on a
  nonlinear problem with a closed form solution. A linear problem is not enough:
  it only measures the order at which the stability function agrees with the
  exponential, which for several high order methods is higher than their true
  order.
- Every method file agrees with the analysis of its own coefficients.
- The diagonalized stage solve reproduces the dense one.
- Every controller and every nonlinear solver drives an implicit method to
  tolerance.
- Radau IIA 5, ESDIRK 8(5) and BDF5 agree to seven digits on Robertson.

## Command line

```
solvers list                              every method with its derived properties
solvers problems                          the test problems
solvers analyze rkdp54                    order and stability as JSON
solvers verify                            check every file against its analysis
solvers convergence bdf4 nonlinear_decay  fixed step convergence study
solvers cost esdirk43 robertson           adaptive work precision diagram
solvers solve radau_iia_5 robertson 1e-8  integrate and report the statistics
```

## The interface

```
wasm-pack build crates/solvers-wasm --target web --out-dir ../../web/pkg --release
python -m http.server 8099 --directory web
```

Then open <http://127.0.0.1:8099>. The method files are compiled into the
module, so there is no server side: stability regions, convergence studies and
work precision diagrams are all computed in the page by the same code the
command line runs.

The views are the library with its derived properties, a gallery of stability
regions, one page per method with its tableau and references, convergence,
work precision, and a comparison of the step size controllers.

## References

The framework follows the standard references, and every method file carries the
DOI of its own source.

- E. Hairer, S. P. Noersett, G. Wanner, *Solving Ordinary Differential Equations
  I: Nonstiff Problems*, 2nd ed., Springer 1993,
  [10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
- E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II: Stiff and
  Differential-Algebraic Problems*, 2nd ed., Springer 1996,
  [10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
- J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd
  ed., Wiley 2016,
  [10.1002/9781119121534](https://doi.org/10.1002/9781119121534)
- G. Soederlind, *Digital filters in adaptive time-stepping*, ACM TOMS 29(1),
  2003, [10.1145/641876.641877](https://doi.org/10.1145/641876.641877)
- H. F. Walker, P. Ni, *Anderson acceleration for fixed-point iterations*, SIAM
  J. Numer. Anal. 49(4), 2011,
  [10.1137/10078356X](https://doi.org/10.1137/10078356X)

## License

MIT.
