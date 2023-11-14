


abstract type ResponseKind end
struct ReducedForm <: ResponseKind end
struct Orthogonal  <: ResponseKind end
struct Structural  <: ResponseKind end


"""
    ImpulseResponse{R<:ResponseKind, T<:AbstractVectorAutoRegression}

representation of an impulse reponse for a VAR of type `T`. `R` is either
- `Orthogonal`, or
- `ReducedForm`
"""
struct ImpulseResponse{R<:ResponseKind, T<:AbstractVectorAutoRegression}
    model::T
    shock::VariableID
    paths::Matrix{Float64}
    scale::Real
end

Base.show(io::IO, ir::ImpulseResponse) = print(io, "ImpulseResponse($(ir.shock))")
function Base.show(io::IO, ::MIME"text/plain", ir::ImpulseResponse{R, T}) where {R, T}
    labels   = ["shock:", "length:", "model:", "response kind:"]
    features = [ir.shock, length(ir), T.name.wrapper, R]

    # left-align all the rows before joining
    width = maximum(length, labels) + 1
    rows  = join.(zip(rpad.(labels, width), features), "")

    s = "Impulse Response"
    print(io, join([s, rows...], "\n    "))
end



@recipe function f(ir::ImpulseResponse{R, T}, responses = nothing; q = 0.66) where {R, T}

    m        = model(ir)
    all_vars = varnames(m)
    
    # cannot figure out how to reference this recipe in other recipes; settle for this
    if isnothing(responses)
        resp_vars = all_vars
    elseif isa(responses, AbstractVector)
        resp_vars = responses
    elseif isa(responses, VariableID)
        resp_vars = [responses]
    elseif isa(responses, Integer)
        resp_vars = [all_vars[responses]]
    end
    response_idx = indexin(resp_vars, all_vars)

    print_bands = !isnothing(q)

    # only select the data we need to use
    paths = ir.paths[:, response_idx]
    if print_bands
        ci    = confidence_intervals(ir, q)
        ci    = ci[:, vcat(response_idx, response_idx .+ nvars(m))]
    end

    nr = length(resp_vars)
    x  = 0:length(ir)

    # set up the subplots' shared attributes
    layout         --> (nr, 1)
    plot_title     --> "Response to $(ir.shock) Shock"
    tick_direction --> :out
    xlim           --> (0, length(ir))
    framestyle     --> :origin

    for k ∈ 1:nr
        if print_bands
            @series begin
                seriestype := :path
                subplot    := k

                # ignore series and legend and color cycling
                primary   := false
                linecolor := nothing
                fillalpha := 0.2
                fillrange := ci[:, k]

                # ensure no markers are used for error band
                markershape := :none

                x, ci[:, k+nr]
            end
        end

        @series begin
            seriestype := :line
            subplot    := k
            label      --> string(resp_vars[k])
            linewidth  --> 3

            x, paths[:, k]
        end
    end
end



"""
    length(ir::ImpulseResponse)

the number of periods a response is computed for. excludes the initial "shocked" period
"""
Base.length(ir::ImpulseResponse) = size(ir.paths, 1) - 1



"""
    model(ir::ImpulseResponse)

retrieve the underlying model
"""
model(ir::ImpulseResponse) = ir.model



"""
    nvars(ir::ImpulseResponse)

return the number of variables whose paths are recorded in the impulse response
"""
nvars(ir::ImpulseResponse) = size(ir.paths, 2)



"""
    confidence_intervals(ir::ImpulseResponse[, q::Real = 0.66])

return the lower and upper bounds of the confidence intervals of an impulse response.
`q` denotes the quantile. the lower and upper bounds of the CI for variable `k` are in
columns `k` and `k+nvars(ir)` of the returned array
"""
function confidence_intervals(ir::ImpulseResponse, q::Real = 0.66)
    se = standard_errors(ir)
    T, N = size(se)

    z  = quantile(Normal(), q)
    ci = zeros(Float64, T, 2*N)

    ci[:, 1:N]     = ir.paths - z*se
    ci[:, N+1:end] = ir.paths + z*se

    return ci
end



"""
    shock(m::AbstractVectorAutoRegression, s::VariableID; kwargs)
    shock(m::AbstractVectorAutoRegression, s::Integer; kwargs)

compute the impulse response of the model to shock `s`.

### Keyword Arguments
- `periods::Integer = 24`: the number of periods for which to compute responses. the
    computed IRF will contain `periods+1` observations: the initial shocked period, and
    the subsequent `periods` observations of responses
- `kind::Symbol = :orthogonal` determines whether the impulse responses are orthogonalized
    (the default) or reduced-form (`kind = :reduced`). if `m` is an `SVAR`, `kind` can also
    be set to `structural`
- `normalize::Bool = true` in the case of orthogonalized shocks, this normalizes all the
    impulse responses so that the initial observation of the shocked variable is 1
"""
function shock(
    m::AbstractVectorAutoRegression,
    s::VariableID;
    periods::Integer = 24,
    kind::Symbol     = :orthogonal,
    normalize::Bool  = true
)
    isfitted(m) || error("model must be fit prior to computing IRFs")

    # find which shock we're concerned with
    sidx = findfirst(Symbol(s) .== Symbol.(m.vars))
    isnothing(sidx) ? error("$s is not the name of a variable") : nothing

    if kind === :orthogonal
        k = Orthogonal
    elseif kind === :reduced
        k = ReducedForm
    elseif kind === :structural
        (m isa AbstractSVAR) || throw(ArgumentError("structural shocks only for SVARs"))
        k = Structural
    else
        error("unrecognized impulse kind: $kind")
    end

    irf = _shock_paths(k(), m, sidx, periods)
    if normalize & (kind === :orthogonal)
        scale = irf[1, sidx]
        irf   = irf / scale
    else
        scale = 1.0
    end
    return ImpulseResponse{k, typeof(m)}(m, s, irf, scale)
end
function shock(m::AbstractVectorAutoRegression, s::Integer; kwargs...)
    (s ≤ 0) ? error("shock = $s. shock index must be positive") : nothing
    (s > nvars(m)) ? error("shock = $s. model only has $(nvars(m)) variables") : nothing
    return shock(m, varnames(m)[s]; kwargs...)
end

function shock(
    m::AbstractStructuralVectorAutoRegression;
    periods::Integer = 24,
    kind::Symbol     = :structural,
    normalize::Bool  = true
)
    isfitted(m) || error("model must be fit prior to computing IRFs")

    sidx = structural_equation(m)
    if kind === :orthogonal
        k = Orthogonal
    elseif kind === :reduced
        k = ReducedForm
    elseif kind === :structural
        k = Structural
    else
        error("unrecognized impulse kind: $kind")
    end

    irf = _shock_paths(k(), m, sidx, periods)
    if normalize & (kind === :orthogonal)
        scale = irf[1, sidx]
        irf   = irf / scale
    else
        scale = 1.0
    end
    return ImpulseResponse{k, typeof(m)}(m, varnames(m)[sidx], irf, scale)
end





"""
    standard_errors(ir::ImpulseResponse)

compute the standard errors of an impulse response

### References
Lutkepohl (3.7)
"""
function standard_errors(ir::ImpulseResponse{R, T}) where {R, T}
    shock = ir.shock
    sidx  = findfirst(shock .== varnames(model(ir)))

    unscaled_se = _shock_paths_standard_errors(R(), model(ir), sidx, length(ir))
    return unscaled_se / ir.scale
end



"""
    _shock_paths(m::AbstractVectorAutoRegression, args...)
    _shock_paths(ma::MovingAverage, P::AbstractMatrix, args...)

compute the impulse response paths of all the variables in a VAR
"""
function _shock_paths(kind::ResponseKind, m::AbstractVectorAutoRegression, args...)
    return _shock_paths(m.ma, UniformScaling(1.0), args...)
end
function _shock_paths(kind::Orthogonal, m::AbstractVectorAutoRegression, args...)
    L = cholesky_decomp(m)
    return _shock_paths(m.ma, L, args...)
end
function _shock_paths(kind::Structural, m::AbstractVectorAutoRegression, args...)
    t = typeof(m).name.wrapper
    error("cannot compute structural shock for $t type")
end
function _shock_paths(kind::Structural, m::AbstractSVAR, args...)
    return _shock_paths(m.ma, m.H, args...)
end
function _shock_paths(
    ma::MovingAverage,
    H::Union{AbstractMatrix, UniformScaling},
    sidx::Integer,
    periods::Integer
)
    # account for the shocked period
    periods_and_shock = periods + 1

    # preallocate the time series array
    n   = nvars(ma)
    irf = zeros(Float64, (periods_and_shock, n))

    # structural shock
    ε       = H[:, sidx]
    # ε       = zeros(Float64, n)
    # ε[sidx] = 1

    # compute shocks
    for t ∈ 1:periods_and_shock
        # irf[t, :] = (ma[t-1] * H) * ε
        irf[t, :] = ma[t-1]*ε
    end

    return irf
end
    


function _shock_paths_standard_errors(
    kind::ReducedForm,
    m::AbstractVectorAutoRegression{P},
    sidx::Integer,
    periods::Integer,
) where {P}

    periods_and_shock = periods + 1
    K                 = nvars(m)

    Σ = zeros(Float64, periods_and_shock, K)
    for i ∈ 1:periods
        Φi_cov    = covariance(m.ma, i)
        Φi_se     = reshape(sqrt.(diag(Φi_cov)), K, K)
        Σ[i+1, :] = Φi_se[:, sidx]
    end

    return Σ
end
function _shock_paths_standard_errors(
    kind::Orthogonal,
    m::AbstractVectorAutoRegression{P},
    sidx::Integer,
    periods::Integer,
) where {P}

    periods_and_shock = periods + 1
    K                 = nvars(m)

    # see Lutkepohl appendix A.12
    Lₖ      = elimination_matrix(K)
    Kₖₖ     = commutation_matrix(K, K)
    Dₖ_plus = pinv(duplication_matrix(K))

    H = cholesky_decomp(m)

    # in Lutkepohl, this is the H matrix in (3.7.9)
    Ik   = Matrix(1.0*I, K, K)
    aech = Lₖ' * inv( Lₖ*(I + Kₖₖ)*kron(H, Ik)*Lₖ')

    # covariance matrices of the VAR coefficients, residuals, and of the innovations
    Σα = m.ma.Σα
    Σu = model_cov(m)
    Σσ = 2 * Dₖ_plus * kron(Σu, Σu) * Dₖ_plus' # in Remark 4 of (3.7); pg 113

    HIk = kron(H', Ik) # save a few computations
    Σ   = zeros(Float64, periods_and_shock, K)
    for i ∈ 0:periods
        if i == 0
            Ci_part = zeros(Float64, K^2, K^2)
        else
            Ci      = HIk * get_G(m.ma, i)
            Ci_part = nobs(m) * Ci * Σα * Ci' # 'nobs' is taken care of below
        end

        Ci_bar = kron(Ik, m.ma[i]) * aech
        Σi_cov = (Ci_part + Ci_bar * Σσ * Ci_bar') / nobs(m)
        Σi_se  = reshape(sqrt.(diag(Σi_cov)), K, K)
        Σ[i+1, :] = Σi_se[:, sidx]
    end

    return Σ
end
function _shock_paths_standard_errors(kind::Structural, args...)
    error("structural standard errors not available yet")
end
