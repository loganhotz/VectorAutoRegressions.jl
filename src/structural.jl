"""
Enforcing structural restrictions
"""


"""
    AbstractStructuralVectorAutoRegression{P} <: AbstractVectorAutoRegression{P}

Abstract supertype for structural VARs. Aliased as `AbstractSVAR`
"""
abstract type AbstractStructuralVectorAutoRegression{P} <: AbstractVectorAutoRegression{P} end
const AbstractSVAR = AbstractStructuralVectorAutoRegression



@lazy struct StructuralVectorAutoRegression{P} <: AbstractStructuralVectorAutoRegression{P}
    formula::FormulaTerm
    vars::VectorRegressor

    rule_exorig::Rule # original expression of the Rule
    @lazy rule::Rule
    @lazy H::Matrix{<:Union{Missing, Float64}}

    @lazy data::Any # need type for `@lazy` macro
    @lazy schemas::OrderedDict{Symbol, FormulaTerm}

    @lazy intercept::OrderedDict{Symbol, Bool}
    @lazy coef::OrderedDict{Symbol, Vector{Float64}}
    @lazy fitted::OrderedDict{Symbol, Vector{Float64}}
    
    @lazy ma::MovingAverage{P}
end
const SVAR = StructuralVectorAutoRegression

"""
    StructuralVectorAutoRegression{P} <: AbstractStructuralVectorAutoRegression{P}

A representation of a structural VAR. To create a SVAR, use the `@var` macro.

# Examples
To create the backward-looking Taylor-rule VAR from Stock and Watson (2001), one should
run the following code:

```julia-repl
julia> include("data/manager.jl") # exposes some data retrieval codes
julia> df = load_data(:sw2001)
julia> svar = @var(y ~ 1 + L(y, 4),
                   y ≡ [gdpd, urate, ffr],
                   ffr => ffr - sma(gdpd, 4) - sma(urate, 4) = 1 + L(y, 4))
```
"""
function StructuralVectorAutoRegression(f::FormulaTerm, v::VectorRegressor, r::Rule)
    return StructuralVectorAutoRegression{maximum_lag(f)}(f, v, r,
                                                          uninit, uninit, uninit, uninit,
                                                          uninit, uninit, uninit, uninit)
end



"""
    assign_data!(m::StructuralVectorAutoRegression, data)

Equip the VAR with concrete data, as well as the data needed for the structural estimation
"""
function assign_data!(m::StructuralVectorAutoRegression, data)
    validate_structural_equation!(m.rule_exorig.eqn, m)
    
    # get a vector of the formulas representing individual equations in the VAR
    fmls = interpolate(m.formula, m.vars)

    schemas    = OrderedDict{Symbol, FormulaTerm}()
    intercepts = OrderedDict{Symbol, Bool}()
    for fml in fmls
        fml = apply_schema(fml, continuous_schema(data))

        schemas[fml.lhs.sym] = fml
        intercepts[fml.lhs.sym] = InterceptTerm{true} ∈ typeof.(fml.rhs.terms) # hate this
    end

    rule_interp = interpolate(m.rule_exorig, m.vars)
    rule_sch    = apply_schema(rule_interp, continuous_schema(data))

    @init! m.data      = data
    @init! m.schemas   = schemas
    @init! m.intercept = intercepts
    @init! m.rule      = rule_sch    

    return m
end


function StatsAPI.fit!(m::StructuralVectorAutoRegression{P}) where {P}

    coef, fitted = fit_schemas(m.schemas, m.data, m.vars)

    @init! m.coef   = coef
    @init! m.fitted = fitted
    @init! m.ma     = MovingAverage(m)

    H = structural_effect(m)
    @init! m.H = H

    return m
end

function validate_structural_equation!(t::Term, m::StructuralVectorAutoRegression)
    var_name = t.sym
    (var_name ∈ m.vars) || throw(ArgumentError("structural eqn not found in model: $t"))
end
function validate_structural_equation!(t::ConstantTerm, m::StructuralVectorAutoRegression)
    (t.n ≤ nvars(m)) || throw(ArgumentError("model only has $(nvars(m)) equations: $t"))
end



"""
    structural_effect(m::StructuralVectorAutoRegression)

Compute the structural effect implied by each of the model's rules
"""
function structural_effect(m::StructuralVectorAutoRegression)

    # locate which endogenous variable the rule corresponds with
    validate_structural_equation!(m.rule.eqn, m)
    idx = structural_equation(m.rule.eqn, m)

    # estimate the structural shock
    Y_full = reduce(+, modelcols(m.rule.lhs, m.data))
    X_full = hcat(modelcols(m.rule.rhs, m.data)...)
    Y, X   = drop_missing_observations(Y_full, X_full)

    β = X \ Y
    ε = Y - X*β

    u = residual_block(m)
    H = missings(Float64, nvars(m), nvars(m))
    H[:, idx] = vec(ε\u) # assumed only one structural shock => vec is fine

    return H
end



"""
    structural_equation(m::AbstractStructuralVectorAutoRegression)

Return the index of the equation that the structural constraint operates on
"""
structural_equation(m::AbstractSVAR)          = structural_equation(m.rule.eqn, m)
structural_equation(t::Term, m::SVAR)         = findfirst(t.sym .== m.vars)
structural_equation(t::ConstantTerm, m::SVAR) = t.n



"""
    structural_variable(m::AbstractStructuralVectorAutoRegression)

Return the index of the variable that the structural constraint operates on
"""
structural_variable(m::AbstractSVAR) = m.vars[structural_equation(m.rule.eqn, m)]



function Base.show(io::IO, ::MIME"text/plain", m::AbstractSVAR{P}) where {P}
    labels   = ["model:", "rule:", "vars:", "lag:"]
    features = [m.formula, m.rule_exorig, m.vars, P]

    if isfitted(m)
        push!(labels, "nobs:")
        push!(features, nobs(m))
    end

    # left-align all the rows before joining
    width = maximum(length, labels) + 1
    rows  = join.(zip(rpad.(labels, width), features), "")

    s = typeof(m).name.wrapper
    print(io, join([s, rows...], "\n    "))
end
