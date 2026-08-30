# solvers

ODE integration methods as data. A Butcher tableau or a multistep coefficient
pattern lives in a JSON file with the DOI it was published under; order, stage
order, stability and cost are derived from those coefficients rather than
declared.

74 methods: explicit and SSP Runge-Kutta, DIRK and ESDIRK, Gauss-Legendre,
Radau IIA, Lobatto IIIA/IIIB/IIIC, BDF, Adams-Bashforth and Adams-Moulton,
Nystrom, Milne-Simpson.

## Layout

    crates/solvers-core   framework
    crates/solvers-cli    command line
    crates/solvers-wasm   browser bindings
    methods/              method files
    tools/                method file generators
    web/                  interface

## Build

    cargo test
    ./build-web.sh        builds the wasm module and serves the interface

## Command line

    solvers list
    solvers problems
    solvers analyze rkdp54
    solvers verify
    solvers convergence bdf4 nonlinear_decay
    solvers cost esdirk43 robertson
    solvers solve radau_iia_5 robertson 1e-8

## Method file

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
  "references": [{ "authors": "Ehle, B. L.", "year": 1969 }]
}
```

A fraction written as a string is exact and the order conditions are verified as
identities; a floating point number stays one and is checked numerically.

A multistep file declares which coefficients are free instead of their values:

```json
{
  "id": "bdf3",
  "class": "linear_multistep",
  "order": 3,
  "multistep": {
    "steps": 3,
    "alpha": ["free", "free", "free", "free"],
    "beta": [1, 0, 0, 0],
    "startup": "esdirk43"
  }
}
```

They are solved for at every step from the actual step size history.

## The library as data

The whole library, coefficients as the files state them together with everything
derived from them:

    https://solvers.milanrother.com/api/methods.json
    https://solvers.milanrother.com/api/methods/<id>.json
    https://solvers.milanrother.com/api/problems.json
    https://solvers.milanrother.com/api/index.json

or `solvers export <directory>` from a checkout.

## Verification

`cargo test` and `solvers verify`.

Every method converges at its stated order on a nonlinear problem with a closed
form solution, every file agrees with the analysis of its own coefficients, the
diagonalized stage solve reproduces the dense one, and Radau IIA 5, ESDIRK 8(5)
and BDF5 agree to seven digits on Robertson.

Adams-Bashforth 8 is exempt from the empirical order test: its stability region
forces a step size at which it is already exact to double precision. Its order
is still checked from its coefficients.

## References

- Hairer, Noersett, Wanner, *Solving Ordinary Differential Equations I*, 2nd ed.,
  Springer 1993, [10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
- Hairer, Wanner, *Solving Ordinary Differential Equations II*, 2nd ed.,
  Springer 1996, [10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
- Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed.,
  Wiley 2016, [10.1002/9781119121534](https://doi.org/10.1002/9781119121534)
- Soederlind, *Digital filters in adaptive time-stepping*, ACM TOMS 29(1), 2003,
  [10.1145/641876.641877](https://doi.org/10.1145/641876.641877)

Every method file carries the DOI of its own source.

## Citing

`CITATION.cff` at the root, and GitHub renders it as a citation on the
repository page.

## License

MIT, by [Milan Rother](https://milanrother.com).
