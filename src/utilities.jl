


"""
    commutation_matrix(m::Integer, n::Integer)

generates the mn × mn commutation matrix Kₘₙ which satisfies
```math
Kₘₙvec(G) = vec(G')
```
for any (m × n) matrix
"""
function commutation_matrix(m::Integer, n::Integer)
    mn  = m*n
    idx = reshape(1:mn, m, n)'
    K   = Matrix(Int8(1)*I, mn, mn)
    return BitMatrix(K[:, idx[:]])
end



"""
    duplication_matrix(m::Integer)

generates the m² × m*(m+1)÷2 matrix Dₘ which satisfies
```math
vec(A) = Dₘ vech(H)
```
for any symmetric m × m matrix A
"""
function duplication_matrix(m::Integer)
    ncols = m*(m+1) ÷ 2
    eye   = Matrix(Int8(1)*I, ncols, ncols)

    D_rows = [vec(permutedims(unvech(idx))) for idx in eachrow(eye)]
    D      = permutedims(mapreduce(permutedims, vcat, D_rows))

    return BitMatrix(D)
end



"""
    elimination_matrix(m::Integer)

generates the (m/2)*(m+1) × m² elimination matrix Lₘ which satisfies
```math
vech(F) = Lₘ vec(F)
```
for any (m × m) matrix F
"""
function elimination_matrix(m::Integer)
    vech_idx = BitVector(vec(LowerTriangular(ones(Int8, m, m))))
    eye      = Matrix(Int8(1)*I, m^2, m^2)
    return BitMatrix(eye[vech_idx, :])
end



"""
    unvech(x::AbstractVector)

undoes the `x = vech(X)` operation; i.e. returns the `X` matrix that satisfies this identity

### Notes
- `x` must be a vector with a triangle number length, or else an assignment error is thrown
"""
function unvech(x::AbstractVector)
    nx    = length(x)
    nrows = round(Int, (-1 + sqrt(1 + 8*nx))/2) # quadratic formula

    X                        = zeros(Float64, nrows, nrows)
    X[tril!(trues(size(X)))] = x # assign `x` to the lower triangular portion of `X`
    X                        = X + X'

    X[diagind(X)] /= 2
    return X
end



"""
    drop_missing_observations(As...)

Given a collection of arrays with the same number of rows, drop all observations that have
at least one missing observation
"""
function drop_missing_observations(As...)
    idx = .!missing_observations(As...)
    return map(
               x -> (ndims(x) == 2) ? disallowmissing(x[idx, :]) : disallowmissing(x[idx]),
               As
    )
end



"""
    missing_observations(As...)

Return a BitVector recording the observations with a missing value in any of the matrices in
`As`
"""
function missing_observations(As...)
    !all(x -> ndims(x) ≤ 2, As) && throw(ArgumentError("all arrays must have 1 or 2 dims"))
    !allequal( size.(As, 1) ) && throw(ArgumentError("all arrays must have same # of rows"))

    return vec( .|( map(x -> any(ismissing.(x), dims=2), As)... ) )
end
