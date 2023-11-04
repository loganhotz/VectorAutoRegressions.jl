# VectorAutoRegressions

[![Build Status](https://github.com/loganhotz/VectorAutoRegressions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/loganhotz/VectorAutoRegressions.jl/actions/workflows/CI.yml?query=branch%3Amain)



### Syntax
Vector autoregressions have several equations, by definition, with each equation containing
many regressors. Using [StatsModels](https://juliastats.org/StatsModels.jl/stable/)
`@formula` macro to write them all out is incredibly tedious, and so to simplify the
construction of VAR models, this package introduces the vector lag operator `L`, and an
accompanying vector definition operator $\equiv$.

As an example, suppose one wanted to estimate a VAR with the vector of variables `[inflation,
unemployment, interest_rate]`, with four lags of each variable and a constant. To create the
VAR model that represents this system, the `VectorAutoRegressions.jl` syntax is
```julia
julia> var = @var(
           y ~ 1 + L(y, 4),
           y ≡ [inflation, unemployment, interest_rate]
       )
```
If `df` is a `DataFrame` with columns `inflation`, `unemployment`, and `interest_rate`, then
the model can be estimated by calling `fit!(var, df)`.


### Notes
The main branch of `LazilyInitializedFields` doesn't allow for supertypes. My commit, with
SHA1 `#489be6b` allows for that. To load the right version of the `LazilyInitializedFields`
package, be sure to install as in
```
>>> using Pkg; Pkg.add("LazilyInitializedFields#489be6b")
```
