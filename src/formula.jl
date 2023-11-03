


const VAR_CONTEXT = Any



"""
    L(v, l₀ = 1, lₜ = 1)

Lag operator for constructing VAR formulas.

### Examples
`@formula(y ~ 1 + L(y))` estimates
```math
y_t = A_1 y_{t-1} + A_c + ε_t
```
`@formula(y ~ 1 + L(y, 2) + L(x, 0, 1)` estimates
```math
y_t = A_1 y_{t-1} + A_2 y_{t-2} + B_0 x_t + B_1 x_{t-1} + A_c + ε_t
```
"""
function L end

# removes Expr wrapping when printing to REPL
Base.show(io::IO, t::FunctionTerm{typeof(L)}) = print(io, "L($(join(t.args, ", ")))")



"""
    sma(v, k)

Simple moving average construction for structural Rules. `k` is the number of lags to
include
"""
function sma end





# VARTerm parsing for the VAR formula definitions

struct VARTerm{T<:AbstractTerm} <: AbstractTerm
    term::T
    l0::Int # initial lag
    lt::Int # last lag
end
VARTerm(t::AbstractTerm, l0::ConstantTerm, lt::ConstantTerm) = VARTerm(t, l0.n, lt.n)

# prettier schema printing
Base.show(io::IO, v::VARTerm) = print(io, "L($(v.term), $(v.l0), $(v.lt))")



struct VectorRegressor
    term::Symbol
    args::Vector{Symbol}
end
function VectorRegressor(ex::Expr)
    is_call(ex, :≡) || throw(ArgumentError("expected definition operator ≡, got $(ex.head)"))
    lhs, rhs = ex.args[2], ex.args[3]

    (lhs isa Symbol) || throw(ArgumentError("regressor lhs should be a symbol, got $(lhs)"))
    if !Meta.isexpr(rhs, :vect)
        throw(ArgumentError("regressor rhs should be a vector, got $(rhs)"))
    end

    for v in rhs.args
        (v isa Symbol) || throw(ArgumentError("elements should be symbols, got $v"))
    end

    return VectorRegressor(lhs, [rhs.args...])
end

# iterator interface
Base.IteratorSize(::VectorRegressor)   = Base.HasLength()
Base.IteratorEltype(::VectorRegressor) = Base.HasEltype()
Base.eltype(::VectorRegressor)         = Symbol
Base.show(io::IO, v::VectorRegressor)  = print(io, "$(v.term) ≡ $(v.args)")
Base.length(v::VectorRegressor)        = Base.length(v.args)
Base.size(v::VectorRegressor)          = (Base.length(v), )
Base.iterate(v::VectorRegressor, s::Integer = 1) = s > length(v) ? nothing : (v.args[s], s+1)

# vector interface
Base.getindex(v::VectorRegressor, i) = Base.getindex(v.args, i)
Base.firstindex(::VectorRegressor)   = 1
Base.lastindex(v::VectorRegressor)   = Base.length(v)





"""
    interpolate(fml::FormulaTerm, vr::VectorRegressor)

Interpolate the elements of `vr` into the formula and return `length(vr)`-vector of formulas
"""
function interpolate(fml::FormulaTerm, vr::VectorRegressor)
    # unpack the two sides of the defining formula and verify it matches the VectorRegressor
    # symbol
    lhs, rhs = fml.lhs, fml.rhs
    matches(lhs, vr) || error("regressor vector does not match lhs variable (got $lhs)")

    fmls = Vector{FormulaTerm}(undef, length(vr))
    for (i, v) in enumerate(vr)
        lhs_term  = Term(v)
        rhs_tuple = Tuple(interpolate(rhs, vr))
        fmls[i]   = FormulaTerm(lhs_term, rhs_tuple)
    end

    return fmls
end
function interpolate(rhs::NTuple{N, AbstractTerm}, vr::VectorRegressor) where {N}
    return reduce(vcat, [interpolate(t, vr) for t ∈ rhs])
end
interpolate(t::ConstantTerm, ::VectorRegressor) = t
function interpolate(t::FunctionTerm{typeof(L)}, vr::VectorRegressor)
    name = t.args[1]
    matches(name, vr) || error("regressor vector does not match rhs variable (got $name)")

    if length(t.args) == 1 # L(term)
        l0, lt = 1, 1
    elseif length(t.args) == 2 # L(term, lt)
        l0, lt = 1, t.args[2]
    elseif length(t.args) == 3 # L(term, l0, lt)
        l0, lt = t.args[2], t.args[3]
    else
        throw(ArgumentError("`L` accepts at most 3 arguments"))
    end

    # original function & expression in the @var macro, e.g. `L` and `L(y, 1, 4)`
    f      = t.f
    exorig = t.exorig

    # at this point, the exorig and args representation of the original `L` call don't match
    args = [[Term(v), ensure_constant_term(l0), ensure_constant_term(lt)] for v ∈ vr]
    fs   = [FunctionTerm(f, arg, exorig) for arg in args]

    return fs
end

matches(t::AbstractTerm, args...) = error("lhs of VAR model should be a symbol (got $t)")
matches(t::Term, v::VectorRegressor)   = t.sym == v.term
matches(s::Symbol, v::VectorRegressor) = s == v.term

ensure_constant_term(x::ConstantTerm) = x
ensure_constant_term(x::Integer)      = ConstantTerm(x)





function StatsModels.apply_schema(
    t::FunctionTerm{typeof(L)},
    sch::StatsModels.Schema,
    ctx::Type{<:VAR_CONTEXT}
)
    if length(t.args) == 1 # L(term)
        term   = first(t.args)
        l0, lt = 1, 1

    elseif length(t.args) == 2 # L(terms, lt)
        term, lt = t.args
        checklag(lt)
        l0, lt   = 1, lt.n

    elseif length(t.args) == 3 # L(terms, l0, lt)
        term, l0, lt = t.args
        checklag(l0)
        checklag(lt)
        l0, lt = l0.n, lt.n

    else
        throw(ArgumentError("`L` accepts at most 3 arguments"))

    end

    # check for no leads, and make sure indices are ordered
    checklag(l0, lt)

    term = apply_schema(term, sch, ctx)
    return apply_schema(VARTerm(term, l0, lt), sch, ctx)
end

checklag(l) = (l isa ConstantTerm) || throw(ArgumentError("lag must be a number (got $l)"))
checklag(l::Integer) = (l < 0) && throw(ArgumentError("lags must be positive (got $l)"))
function checklag(l0::Integer, lt::Integer)
    checklag(l0)
    checklag(lt)
    (lt > l0) || throw(ArgumentError("l0 must be less than lt (got l0 = $l0 > $lt = lt)"))
end



function StatsModels.modelcols(t::VARTerm, d::NamedTuple)
    original = modelcols(t.term, d)
    return reduce(hcat, [lag(original, l) for l ∈ t.l0:t.lt])
end

StatsModels.terms(t::VARTerm)     = terms(t.term)
StatsModels.termvars(t::VARTerm)  = StatsModels.termvars(t.term)
StatsModels.width(t::VARTerm)     = t.lt - t.l0
StatsModels.coefnames(t::VARTerm) = StatsModels.coefnames(t.term).*".".*string.(t.l0:t.lt)



# SMATerm parsing for a structural rule
struct SMATerm{T<:AbstractTerm} <: AbstractTerm
    term::T
    k::Int
end

# prettier schema printing
Base.show(io::IO, s::SMATerm) = print(io, "sma($(s.term), $(s.k))")
Base.show(io::IO, f::FunctionTerm{typeof(sma)}) = print(io, "sma($(f.args[1]), $(f.args[2]))")

function StatsModels.modelcols(t::SMATerm, d::NamedTuple)
    original = modelcols(t.term, d)
    return average_window(original, t.k)
end
function average_window(v, k::Integer) # simple moving average
    short = [sum(@view v[i-k+1:i])/k for i ∈ k:length(v)]
    return vcat(missings(k-1), short)
end


function StatsModels.apply_schema(
    t::FunctionTerm{typeof(sma)},
    sch::StatsModels.Schema,
    ctx::Type{<:VAR_CONTEXT}
)
    (length(t.args) != 2) && throw(ArgumentError("`sma` accepts exactly two arguments"))

    # unpack, and verify
    term, k = t.args
    checklag(k)
    (k.n < 2) && throw(ArgumentError("number of periods to average must be ≥ 2 (got $k)"))

    term = apply_schema(term, sch, ctx)
    return apply_schema(SMATerm(term, k.n), sch, ctx)
end


"""
    continuous_schema(df)
    continuous_schema(d::NamedTuple)

a wrapper for `StatsModels.schema` that ensures any Float-valued column is a ContinuousTerm
"""
function continuous_schema(df)
    hints = Dict{Symbol, Any}()
    for v ∈ Symbol.(names(df))
        et = eltype(df[!, v])
        if any( union_types(et) .<: AbstractFloat )
            hints[v] = ContinuousTerm
        end
    end

    return schema(df, hints)
end
function continuous_schema(d::NamedTuple)
    hints = Dict{Symbol, Any}()
    for v ∈ Symbol.(keys(df))
        et = eltype(df[!, v])
        if any( union_types(et) .<: AbstractFloat )
            hints[v] = ContinuousTerm
        end
    end

    return schema(df, hints)
end

# used in continuous_schema
union_types(x::Union) = (x.a, union_types(x.b)...)
union_types(x::Type)  = (x, )



"""
    maximum_lag(f::Formulaterm)

Return the maximum lag
"""
# maximum_lag(t::FormulaTerm) = Integer(maximum([maximum_lag(x) for x ∈ t.rhs]))
function maximum_lag(t::FormulaTerm)
    if t.rhs isa Tuple # lagged term and a constant, for example
        return Integer(maximum([maximum_lag(x) for x ∈ t.rhs]))
    else # e.g. just a lagged term
        return maximum_lag(t.rhs)
    end
end
maximum_lag(::AbstractTerm) = -Inf
maximum_lag(t::FunctionTerm{typeof(L)}) = (length(t.args) > 1) ? t.args[end].n : 1



"""
    sorted_coefnames(t::MatrixTerm, vr::VectorRegressor)
    sorted_coefnames(v::Vector{String}, vr::VectorRegressor)

Sorts the variables to be ordered [(Intercept) (t-1) (t-2) ... (t-p)]
"""
sorted_coefnames(t::MatrixTerm, vr::VectorRegressor) = sorted_coefnames(coefnames(t), vr)
function sorted_coefnames(vs::Vector{<:AbstractString}, vr::VectorRegressor)
    var_locs = Dict{String, Int64}(string(v) => i for (i, v) in enumerate(vr))
    return sort(vs; by = x -> coefname_sorter(x, var_locs))
end

function coefname_sorter(x, all_vars::Dict{String, Int64})
    parts = split(x, ".")
    (length(parts) == 1) ? (return x) : (return join((parts[2], all_vars[parts[1]]), "."))
end





# Rules parsing for structural identities

struct Rule{T, L, R} # not sure if this should be a subtype of `StatsModels.AbstractTerm`
    eqn::T
    lhs::L
    rhs::R
end

Base.show(io::IO, r::Rule) = print(io, "$(r.lhs) = $(r.rhs)")
function Base.show(io::IO, mime::MIME"text/plain", r::Rule; prefix = "")
    println(io, "Structural Rule")
    print(io, "Left-hand side")
    show(io, mime, r.lhs, prefix = "\n    ")
    println(io)
    print(io, "Right-hand side")
    show(io, mime, r.rhs, prefix = "\n    ")
end



const SPECIALS = (:+, :-, :*, :/)

function Rule(ex::Expr)
    # syntax check before any other parsing
    h = ex.head
    Meta.isexpr(ex, :(=))  || throw(ArgumentError("expected Rule separator =, got $(h)"))
    (length(ex.args) == 2) || throw(ArgumentError("malformed definition in Rule $ex"))

    # eqn_and_lhs has form `eqn => lhs`, and rhs has form `rhs`
    eqn_and_lhs, rhs = ex.args

    if Meta.isexpr(rhs, :block)
        # equals sign causes the RHS to be read in as block; first arg is a LineNumberNode
        rhs = rhs.args[2]
    else
        throw(ArgumentError("malformed rhs in Rule, $ex"))

    end

    if is_call(eqn_and_lhs, :(=>))
        eqn = term(eqn_and_lhs.args[2])
        lhs = eqn_and_lhs.args[3]

    elseif is_call(eqn_and_lhs)
        eqn = term(1) # don't like this
        lhs = eqn_and_lhs

    else
        throw(ArgumentError("malformed lhs in Rule, $ex"))

    end

    lhs_rule = parse_rule!(lhs, false) # explicitly pass `protected` in case `lhs <: Symbol`
    rhs_rule = parse_rule!(rhs)

    return :(Rule($eqn, $lhs_rule, $rhs_rule))
end

function parse_rule!(ex::Expr, protected::Bool = false)
    catch_dollar(ex)
    check_call(ex)

    if (ex.args[1] ∈ SPECIALS) && !protected
        # ex.args[1]      = esc(ex.args[1])
        ex.args[2:end] .= parse_rule!.(ex.args[2:end], false)

    else
        # capture non-special call, or special call inside a non-special
        exorig  = deepcopy(ex)
        f       = ex.args[1] # in StatsModels, this is `f = esc(ex.args[1])`
        args    = parse_rule!.(ex.args[2:end], true)
        ex.args = [:FunctionTerm,
                   f,
                   :[$(args...)],
                   Meta.quot(exorig)]

        ex
    end

    return ex
end

parse_rule!(::Nothing, protected) = :(nothing)
parse_rule!(s::Symbol, protected) = :(Term($(Meta.quot(s))))
parse_rule!(n::Number, protected) = :(ConstantTerm($n))



function catch_dollar(ex::Expr)
    # respect StatsModels' syntax for the moment
    if Meta.isexpr(ex, :$)
        throw(ArgumentError("interpolation with \$ not supported in @formula. "*
                            "Use @eval @formula(...) instead."))
    end
end
function check_call(ex::Expr)
    # respect StatsModels' syntax for the moment
    if !is_call(ex)
        throw(ArgumentError("non-call expression encountered: $ex"))
    end
end



"""
    interpolate(r::Rule, vr::VectorRegressor)

Interpolate the elements of `vr` into the RHS of the rule (since the LHS is assumed to be
a variable in a provided dataset [and not necessarily an endogenous variable in the VAR]).
This returns a single Rule
"""
interpolate(r::Rule, vr::VectorRegressor) = Rule(r.eqn,r.lhs,Tuple(interpolate(r.rhs,vr)))
interpolate(t::FunctionTerm{typeof(sma)}, vr::VectorRegressor) = t


function StatsModels.apply_schema(
    r::Rule,
    sch::StatsModels.Schema,
    ctx::Type{<:VAR_CONTEXT}
)
    lhs = apply_schema(r.lhs, sch, ctx)
    rhs = apply_schema(r.rhs, sch, ctx)

    return Rule(r.eqn, lhs, rhs)
end





# algebra for rule definition
struct NegativeTerm{T} <: AbstractTerm; term::T; end
Base.:-(a::AbstractTerm, b::AbstractTerm)           = a + NegativeTerm(b)
Base.:-(as::StatsModels.TupleTerm, b::AbstractTerm) = (as..., NegativeTerm(b))

function StatsModels.apply_schema(
    t::NegativeTerm,
    sch::StatsModels.Schema,
    ctx::Type{<:VAR_CONTEXT}
)
    term = apply_schema(t.term, sch, ctx)
    return NegativeTerm(term)
end

function StatsModels.modelcols(t::NegativeTerm, d::NamedTuple)
    original = modelcols(t.term, d)
    return -original
end

Base.show(io::IO, t::NegativeTerm) = print(io, "(-1)*$(t.term)")


"""
    @var(exprs)

Capture and parse a formula and vector definition into a VAR object
"""
macro var(exprs...)
    # (length(exprs) ≥ 3) && throw(ArgumentError("too many arguments (got $exprs)"))
    (length(exprs) ≤ 1) && throw(ArgumentError("too few arguments (got $exprs)"))

    # all invocations of `@var` require at least a model formula and variable definition.
    # using `eval` and `Meta.parse` feels kind of hacky
    fml, def = exprs[1], exprs[2]
    fml      = eval(Meta.parse("@formula($fml)"))
    vr       = VectorRegressor(def)

    # no structural rules
    (length(exprs) == 2) && (return VectorAutoRegression(fml, vr))

    (length(exprs) > 3) && error("not implemented yet")
    rule = eval( Rule(exprs[3]) )
    return StructuralVectorAutoRegression(fml, vr, rule)
end

is_call(ex::Expr)             = Meta.isexpr(ex, :call)
is_call(ex::Expr, op::Symbol) = Meta.isexpr(ex, :call) && ex.args[1] == op
is_call(ex::Any)              = false
is_call(ex::Any, op::Any)     = false
