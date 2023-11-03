"""
common interface for VAR models
"""

function Base.show(io::IO, m::AbstractVectorAutoRegression{P}) where {P}
    s = typeof(m).name.wrapper
    print(io, "$s{lags=$P}")
end

function Base.show(io::IO, ::MIME"text/plain", m::AbstractVectorAutoRegression{P}) where {P}
    labels   = ["model:", "vars:", "lag:"]
    features = [m.formula, m.vars, P]

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


@lazy mutable struct VectorAutoRegression{P} <: AbstractVectorAutoRegression{P}
    formula::FormulaTerm
    vars::VectorRegressor

    @lazy data::Any # need type for `@lazy` macro
    @lazy schemas::OrderedDict{Symbol, FormulaTerm}

    @lazy intercept::OrderedDict{Symbol, Bool}
    @lazy coef::OrderedDict{Symbol, Vector{Float64}}
    @lazy fitted::OrderedDict{Symbol, Vector{Float64}}

    @lazy ma::MovingAverage{P}
end
const VAR = VectorAutoRegression

function VectorAutoRegression(f::FormulaTerm, v::VectorRegressor)
    return VectorAutoRegression{maximum_lag(f)}(f, v,
                                                uninit, uninit, uninit, uninit, uninit,
                                                uninit)
end


function MovingAverage(m::AbstractVectorAutoRegression)
    !allequal([sch.rhs for sch ∈ values(m.schemas)]) && error("all RHS need to match")

    Z       = modelmatrix(m, m.vars[1])
    ΓY0_inv = inv( Symmetric(Z' * Z) )[2:end, 2:end] # remove the intercept
    Σα_hat  = Symmetric( kron(ΓY0_inv, model_cov(m)) )

    return MovingAverage(coef_matrices(m), Σα_hat)
end


"""
    assign_data!(m::AbstractVectorAutoRegression, data)

Equip the VAR with concrete data
"""
function assign_data!(m::AbstractVectorAutoRegression, data)
    # get a vector of the formulas representing individual equations in the VAR
    fmls = interpolate(m.formula, m.vars)

    schemas    = OrderedDict{Symbol, FormulaTerm}()
    intercepts = OrderedDict{Symbol, Bool}()
    for fml in fmls
        fml = apply_schema(fml, continuous_schema(data))

        schemas[fml.lhs.sym] = fml
        intercepts[fml.lhs.sym] = InterceptTerm{true} ∈ typeof.(fml.rhs.terms) # hate this
    end

    @init! m.data = data
    @init! m.schemas = schemas
    @init! m.intercept = intercepts

    return m
end



"""
    cholesky_decomp(m::AbstractVectorAutoRegression)

Return the lower-triangular matrix of the Cholesky decomposition of the model error cov
"""
cholesky_decomp(m::AbstractVectorAutoRegression) = cholesky(model_cov(m)).L



"""
    coef(m::AbstractVectorAutoRegression)
    coef(m::AbstractVectorAutoRegression, v::Symbol)

Return the coefficients in each of the equations of the VAR. Optionally, the coefnames of a
single equation can be returned
"""
StatsBase.coef(m::AbstractVectorAutoRegression)            = m.coef
StatsBase.coef(m::AbstractVectorAutoRegression, v::Symbol) = m.coef[v]



"""
    coefnames(m::AbstractVectorAutoRegression)
    coefnames(m::AbstractVectorAutoRegression, v::Symbol)

Return the coefficient names in each of the equations of the VAR. Optionally, the coefnames
of a single equation can be returned
"""
function StatsModels.coefnames(m::AbstractVectorAutoRegression)
    return Dict(s => sorted_coefnames(sch.rhs, m.vars) for (s, sch) ∈ pairs(m.schemas))
end
function StatsModels.coefnames(m::AbstractVectorAutoRegression, v::Symbol)
    return sorted_coefnames(m.schemas[v].rhs, m.vars)
end



"""
    coef_block(m::AbstractVectorAutoRegression)

Return a matrix of the coefficients, organized into blocks indexed by the lag of their
corresponding endogenous variable
"""
function coef_block(m::AbstractVectorAutoRegression)
    cn       = [c for (_, c) ∈ coefnames(m)]
    all_vars = sorted_coefnames(unique(vcat(cn...)), m.vars)

    A = zeros(Float64, nvars(m), length(all_vars))
    for (n, v) ∈ enumerate(m.vars)
        c   = cn[n]
        idx = indexin(all_vars, c)

        A[n, idx] = coef(m, v)
    end

    return A
end



"""
    coef_matrices(m::AbstractVectorAutoRegression{P})

Return a `P`-tuple of the `N × N` matrices of the coefficient
"""
function coef_matrices(m::AbstractVectorAutoRegression{P}) where {P}
    # remove intercept if necessary
    A = coef_block(m)
    hasintercept(m) && (A = A[:, 2:end])

    N = nvars(m)
    return Tuple( [A[:, 1+N*(i-1):N*i] for i ∈ 1:P] )
end



"""
    function companion_form(m::AbstractVectorAutoRegression)

Return the companion form of the estimated model. If the model is
```math
Y_t = A_c + A_1*Y_{t-1} + ... + A_p*Y_{t-p} + ε_t
```
this function returns
```math
[
    A_1 A_2 ... A_{p-1} A_p
    I   0   ... 0       0
    0   I   ... 0       0
    .   .       .       .
    0   0   ... I       0
]
```
"""
companion_form(m::AbstractVectorAutoRegression) = companion_form(coef_matrices(m))
function companion_form(As::NTuple{P, Matrix{Float64}}) where {P}
    N = size(As[1], 1)
    C = zeros(Float64, N*P, N*P)

    # first block row
    for (i, A) ∈ enumerate(As)
        C[1:N, (i-1)*N+1:i*N] = A
    end

    # subsequent identity matrices
    eye = Matrix(1.0*I, N, N)
    for i ∈ 1:P-1
        C[i*N+1:(i+1)*N, (i-1)*N+1:i*N] = eye
    end

    return C
end
function companion_form(m::AbstractVectorAutoRegression{P}, s::Symbol) where {P}
    # constructs the \bm{B} matrix in (3.5.10)
    if s === :intercept
        K    = nvars(m)
        KPm1 = K*(P-1)

        top       = zeros(Float64, 1, K*P+1)
        top[1, 1] = 1

        bottom = zeros(Float64, KPm1, K*P+1)
        bottom[:, 2:KPm1+1] = Matrix(1.0*I, KPm1, KPm1)

        return vcat(top, coef_block(m), bottom)
    else
        error("unrecognized companion form symbol: $s")
    end
end



"""
    data(m::AbstractVectorAutoRegression)

Returns the associated data of a model
"""
data(m::AbstractVectorAutoRegression) = m.data



"""
    dof(m::AbstractVectorAutoRegression)

Returns `nobs(m) - nvars(m)*nlags(m) - 1`
"""
dof(m::AbstractVectorAutoRegression) = nobs(m) - nvars(m)*nlags(m) - 1



"""
    fitted(m::AbstractVectorAutoRegression)
    fitted(m::AbstractVectorAutoRegression, v::Symbol)

Return the fitted values of each equation in the VAR.
"""
fitted(m::AbstractVectorAutoRegression)            = m.fitted
fitted(m::AbstractVectorAutoRegression, v::Symbol) = m.fitted[v]



"""
    hasdata(m::AbstractVectorAutoRegression)

Check if the model has had data assigned to it. See also [`assign_data!`](@ref)
"""
hasdata(m::AbstractVectorAutoRegression) = @isinit m.data



"""
    hasintercept(m::AbstractVectorAutoRegression)
    hasintercept(m::AbstractVectorAutoRegression, v::Symbol)

Check if any of the equations in the model has a constant.
"""
hasintercept(m::AbstractVectorAutoRegression)            = any( values(m.intercept) )
hasintercept(m::AbstractVectorAutoRegression, v::Symbol) = m.intercept[v]



"""
    isfitted(m::AbstractVectorAutoRegression)

Check if the model has been fitted
"""
isfitted(m::AbstractVectorAutoRegression) = @isinit m.coef



"""
    fit!(m::AbstractVectorAutoRegression, data::AbstractDataFrame)
    fit!(m::AbstractVectorAutoRegression)

Estimate a VAR. If the `data` argument isn't provided, it's assumed the model is already
equipped with a dataset via `assign_data!`
"""
function fit!(m::AbstractVectorAutoRegression{P}, data; kwargs...) where {P}
    hasdata(m) && throw(ArgumentError("to-be-estimated model already has data"))

    # write schemas `L` -> `VARTerm` and attach the data to the model
    m = assign_data!(m, data)
    return StatsAPI.fit!(m; kwargs...)
end
function fit!(m::AbstractVectorAutoRegression{P}; kwargs...) where {P}
    !hasdata(m) && throw(ArgumentError("to-be-estimated model does not have data"))
    return StatsAPI.fit!(m; kwargs...)
end

function StatsAPI.fit!(m::AbstractVectorAutoRegression{P}) where {P}
    
    coef, fitted = fit_schemas(m.schemas, m.data, m.vars)

    @init! m.coef   = coef
    @init! m.fitted = fitted
    @init! m.ma     = MovingAverage(m)

    return m
end



function fit_schemas(schemas::OrderedDict, data, vars::VectorRegressor)
    coef   = OrderedDict{Symbol, Vector{Float64}}()
    fitted = OrderedDict{Symbol, Vector{Float64}}()
    for (sym, sch) in pairs(schemas)
        lhs_full, rhs_full = modelcols(sch, data)
        lhs, rhs           = drop_missing_observations(lhs_full, rhs_full)

        # sorts from [var.lag] ordering to [lag.var]
        c      = coefnames(sch.rhs)
        sc     = sorted_coefnames(c, vars)
        sc_idx = indexin(sc, c)

        β, Yhat = _ols_output(lhs, rhs[:, sc_idx])
        coef[sym]   = β
        fitted[sym] = Yhat
    end

    return coef, fitted
end




"""
    modelmatrix(m::AbstractVectorAutoRegression)
    modelmatrix(m::AbstractVectorAutoRegression, v::Symbol)

Return the RHS data in each of the equations of the VAR. Optionally, the RHS data of a
single equation can be returned
"""
function modelmatrix(m::AbstractVectorAutoRegression)
    mm = Dict{Symbol, Matrix{Float64}}()
    for (v, sch) ∈ pairs(m.schemas)
        lhs_full, rhs_full = modelcols(sch, m.data)
        lhs, rhs           = drop_missing_observations(lhs_full, rhs_full)

        mm[v] = rhs
    end

    return mm
end
function modelmatrix(m::AbstractVectorAutoRegression, v::Symbol)
    lhs_full, rhs_full = modelcols(m.schemas[v], m.data)
    lhs, rhs           = drop_missing_observations(lhs_full, rhs_full)

    return rhs
end



"""
    model_cov(m::AbstractVectorAutoRegression)

Returns Σ ≡ ( ε * ε' ) / dof(m)
"""
function model_cov(m::AbstractVectorAutoRegression)
    ε = residual_block(m)
    return (ε' * ε) / dof(m)
end



"""
    moving_average(m::AbstractVectorAutoRegression)

Returns the moving-average (Wold) representation of a model
"""
moving_average(m::AbstractVectorAutoRegression) = m.ma



"""
    mse(m::AbstractVectorAutoRegression, h::Integer = 1, kind::Symbol = :exact)

computes the MSE matrix of `y_{t+h}`. If the VAR(P) process has the MA representation
```math
    y_t = \\mu + \\sum_{i=0}^{\\infty} \\Phi_i u_i, \\quad Var(u) = \\Sigma_u
```
then the "exact" MSE matrix at horizon `h` is
```math
    \\Sigma_y(h) = \\sum_{i=0}^{h-1} \\Phi_i \\Sigma_u \\Phi_i'
```

The "estimated" MSE matrix at horizon `h` is `\\Sigma_y(h) + \\Omega(h)`, where Ω(h)
is given by (3.5.12) in Lutkepohl
"""
function mse(m::AbstractVectorAutoRegression, h::Integer = 1, kind::Symbol = :exact)
    if kind === :exact
        return mse_exact(m, h)
    elseif kind === :estimate
        return mse_estimate(m, h)
    else
        error("unrecognized MSE matrix kind: $kind")
    end
end
function mse_exact(m::AbstractVectorAutoRegression, h::Integer)
    Σu = model_cov(m)
    Σy = zeros(Float64, nvars(m), nvars(m))

    for i ∈ 0:h-1
        Σy += m.ma[i] * Σu * m.ma[i]'
    end

    return Σy
end
function mse_estimate(m::AbstractVectorAutoRegression, h::Integer)
    Σu = model_cov(m)
    B  = companion_form(m, :intercept)

    Z     = modelmatrix(m, varnames(m)[1]) # hacky
    Γ     = Symmetric( (Z' * Z) / nobs(m) )
    Γ_inv = inv(Γ)

    # equation (3.5.12)
    Ω = zeros(Float64, nvars(m), nvars(m))
    for i ∈ 0:h-1
        for j ∈ 0:h-1
            prefix = tr((B')^(h-1-i) * Γ_inv * B^(h-1-j) * Γ)
            Ω     += prefix * m.ma[i] * Σu * m.ma[j]'
        end
    end

    return mse_exact(m, h) + Ω / nobs(m)
end



"""
    nlags(m::AbstractVectorAutoRegression)

Returns the number of lags
"""
nlags(m::AbstractVectorAutoRegression{P}) where {P} = P



"""
    nobs(m::AbstractVectorAutoRegression)

Returns the number of observations used to estimate the model
"""
nobs(m::AbstractVectorAutoRegression) = length( m.fitted[m.vars[1]] )



"""
    nvars(m::AbstractVectorAutoRegression)

Returns the number of variables in the LHS vector
"""
nvars(m::AbstractVectorAutoRegression) = length(m.vars)



"""
    residuals(m::AbstractVectorAutoRegression)
    residuals(m::AbstractVectorAutoRegression, v::Symbol)

Return the residuals from the estimated model. Optionally, the residuals of a single
equation can be returned
"""
function residuals(m::AbstractVectorAutoRegression)
    return Dict(v => response(m, v) - fitted(m, v) for v ∈ keys(m.schemas))
end
residuals(m::AbstractVectorAutoRegression, v::Symbol) = response(m, v) - fitted(m, v)



"""
    residual_block(m::AbstractVectorAutoRegression)

Return a matrix of the residuals for the estimated model. The order is determined by the
declaration order of the endogenous variables in `@var`
"""
residual_block(m::AbstractVectorAutoRegression) = hcat([residuals(m, v) for v ∈ m.vars]...)



"""
    response(m::AbstractVectorAutoRegression)
    response(m::AbstractVectorAutoRegression, v::Symbol)

Return the LHS data in each of the equations of the VAR. Optionally, the LHS data of a
single equation can be returned
"""
function response(m::AbstractVectorAutoRegression)
    rs = Dict{Symbol, Vector{Float64}}()
    for (v, sch) ∈ pairs(m.schemas)
        lhs_full, rhs_full = modelcols(sch, m.data)
        lhs, rhs           = drop_missing_observations(lhs_full, rhs_full)

        rs[v] = lhs
    end

    return rs
end
function response(m::AbstractVectorAutoRegression, v::Symbol)
    lhs_full, rhs_full = modelcols(m.schemas[v], m.data)
    lhs, rhs           = drop_missing_observations(lhs_full, rhs_full)

    return lhs
end



"""
    varnames(m::AbstractVectorAutoRegression)

Return a vector of `Symbol`s of the endogenous variable names
"""
varnames(m::AbstractVectorAutoRegression) = m.vars.args



function _ols_output(Y::AbstractVector, X::AbstractMatrix)
    β    = X \ Y
    Yhat = X * β

    return β, Yhat
end
