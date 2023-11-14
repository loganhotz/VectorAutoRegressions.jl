# VectorAutoRegressions.jl

Welcome to the documentation for the VectorAutoRegressions.jl package.





## Creating a VAR Model

Like many statistical modeling approaches, `VectorAutoRegressions.jl` separates the model
itself from the data. The typical model one sees in textbooks is
```math
y_t = A_1 y_{t-1} + \cdots + A_p y_{t-p} + A_c + \varepsilon_t
```
With many variables in the `y` vector and many lags `p`, expressing these models in the
usual way via a `@formula` macro would be incredibly tedious.

To simplify this representation problem, VectorAutoRegressions.jl introduces the `@var`
macro. In the simplest case it accepts two arguments:
1. The model definition. This encodes information about constant terms and lags.
2. A definition of the model's vector, whose elements correspond to individual data series.



## Examples

To create a reduced form VAR with `p = 4` lags and the three variables
- `inflation`
- `unemployment`
- `interest_rate`

one should define it using the `@var` macro, the lag operator `L`, and the vector definition
operator `\equiv`:

```julia-repl
julia> mod = @var(y ~ 1 + L(y, 4),
                  y ≡ [inflation, unemployment, interest_rate])
```

Importantly, the vector symbol in the model defintion (`y`) is the same as that in the
vector regressor definition.

If one wanted to impose a Taylor rule-type policy on the interest rate, a la Stock and
Watson (2001), the `@var` macro accepts an additional argument. These authors implement the
monetary policy rule by assuming it's set according to
```math
R_t = φ_π\bar{π}_t + φ_u\bar{u}_t + L(R_t, 4) + L(π_t, 4) + L(u_t, 4)
```
Where `R` is the interest rate, `π` is inflation, `u` is the unemployment rate, symbols with
bars over them indicate four-quarter moving averages. The `φ` parameters are exogenously set
to `φπ = 1.5` and `φu = -1.25`. To implement this using the `@var` macro, just use
```julia-repl
julia> smod = @var(y ~ 1 + L(y, 4),
                  y ≡ [r, π, u],
                  r => r - 1.5*sma(π, 4) + 1.25*sma(u, 4) = 1 + L(y))
```
The `sma` function represents simple moving averages, and is defined by the
`VectorAutoRegressions.jl` package. 






## References
Much of the documentation, focus, and coverage of this package is based on
[Helmut Lutkepohl (2005)](https://link.springer.com/book/10.1007/978-3-540-27752-1).
Internally, much of the notation follows it too, but deviates where needed.





## Other VAR packages
There are a couple existing Julia packages for estimating and interacting with VARs:

- [VectorAutoregressions](https://github.com/lucabrugnolini/VectorAutoregressions.jl)
- [VecAutoReg](https://github.com/lucabrugnolini/VectorAutoregressions.jl)
- [VARModels](https://github.com/tomaskrehlik/VARmodels.jl)
- [VectorAR](https://github.com/justinjoliver/julia-VectorAR.jl)

These packages are all listed [here](https://discourse.julialang.org/t/time-series-in-julia-working-list/62539).
In my opinion, none of these are written in a particularly Julian paradigm, and excepting
[VectorAutoregressions](https://github.com/lucabrugnolini/VectorAutoregressions.jl), seem
to be mostly written with a particular paper, or set of papers, in mind. They do not appear
to be in active development, but that particular affliction weighs on many projects. Perhaps
this one may persist. As the the old saying goes, one must imagine the open-source developer
happy.
