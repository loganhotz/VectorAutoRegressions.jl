
"""
    MovingAverage{P}

a moving-average representation of a VAR with `P` lags
"""
mutable struct MovingAverage{P}
    As::NTuple{P, Matrix{Float64}}
    Σα::Symmetric{Float64, Matrix{Float64}} # \hat{Σα} in (3.7.4) and (3.7.5)
    Φs::Vector{Matrix{Float64}}
    Gs::Dict{Int64, Matrix{Float64}}
end
function MovingAverage(As::NTuple{P, MT}, Σ::Symmetric{Float64, MT}) where {P, MT}
    return MovingAverage{P}(As, Σ, Vector{MT}(), Dict{Int64, MT}())
end

Base.show(io::IO, ma::MovingAverage{P}) where {P} = print(io, "MovingAverage(lags=$P)")
function Base.getindex(ma::MovingAverage{P}, k::Integer) where {P}
    (k == 0) ? (return Matrix(1.0*I, size(ma.As[1]))) : nothing
    (k < 0)  ? throw(DomainError("$k. MA indices must be non-negative integers")) : nothing
    (k ≤ length(ma.Φs)) ? (return ma.Φs[k]) : nothing

    Φk = zeros(Float64, size(ma.As[1]))
    for j ∈ 1:k
        if j ≤ P
            Φk += ma[k-j] * ma.As[j]
        end
    end
    push!(ma.Φs, Φk)

    return Φk
end

function get_G(ma::MovingAverage{P}, h::Integer) where {P}
    # Lutkepohl (3.7.6)
    (h == 0) ? (return zeros(Float64, nvars(ma)^2, P*nvars(ma)^2)) : nothing
    (h < 0)  ? throw(DomainError("$h. indices must be non-negative integers")) : nothing

    if h ∉ keys(ma.Gs)
        K  = nvars(ma)
        Ap = companion_form(ma.As)'

        Gh = zeros(Float64, K^2, P*K^2)
        for i = 0:h-1
            Ap_ = Ap^(h-1-i)
            Gh += kron(Ap_[1:K, :], ma[i])
        end
        ma.Gs[h] = Gh
    end

    return ma.Gs[h]
end





"""
    covariance(ma::MovingAverage, h::Integer)

Return the variance-covariance matrix of the `h`-th horizon MA coefficient matrix of a model
"""
function covariance(ma::MovingAverage, h::Integer)
    Gh = get_G(ma, h)
    Σα = ma.Σα

    return Gh * Σα * G'
end



"""
    nvars(ma::MovingAverage)

Returns the number of variables in the MA representation
"""
nvars(ma::MovingAverage) = size(ma.As[1], 1)



"""
    nlags(ma::MovingAverage)

Returns the number of lags in the reduced-form representation of a model
"""
nlags(ma::MovingAverage{P}) where {P} = P
