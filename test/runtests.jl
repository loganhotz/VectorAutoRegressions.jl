using VectorAutoRegressions
using Test

# needed for loading data
using CSV
using DataFrames
using Dates

@testset "VectorAutoRegressions.jl" begin
    # load data utility functions
    include("data/utilities.jl")

    # very simple VAR
    mod = @var(y ~ 1 + L(y, 4),
               y ≡ [infl, urate, ffr])
    df = load_data()
    fit!(mod, df)
end
