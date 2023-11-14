module VectorAutoRegressions



using Distributions,
      LazilyInitializedFields,
      LinearAlgebra,
      Missings,
      OrderedCollections,
      RecipesBase,
      Reexport,
      StatsAPI,
      StatsBase

@reexport using StatsModels

import ShiftedArrays: lag



# VectorAutoRegressions.jl
export AbstractVectorAutoRegression,
       # miscellany
       VariableID

# model.jl
export VectorAutoRegression,
       VAR,
       # interface
       assign_data!,
       cholesky_decomp,
       coef,
       coef_block,
       coef_matrices,
       companion_form,
       data,
       dof,
       fitted,
       hasdata,
       hasintercept,
       isfitted,
       modelmatrix,
       model_cov,
       moving_average,
       mse,
       nlags,
       nobs,
       nvars,
       residuals,
       residual_block,
       response,
       varnames,
       # estimation
       fit! # i'm still unsure why StatsAPI doesn't export this function automatically

# structural.jl
export AbstractStructuralVectorAutoRegression,
       StructuralVectorAutoRegression,
       SVAR,
       # interface
       structural_equation,
       structural_variable

# formula.jl
export L,
       Rule,
       sma,
       # utility functions
       @var,
       continuous_schema

# impulse.jl
export ImpulseResponse,
       # impulse response interface
       model,
       # responses
       ResponseKind,
       ReducedForm,
       Orthogonal,
       Structural,
       # methods
       confidence_intervals,
       shock,
       standard_errors

# wold.jl
export MovingAverage,
       # methods
       covariance

# utilities.jl
export drop_missing_observations,
       missing_observations


"""
    AbstractVectorAutoRegression{P} <: StatisticalModel

Abstract supertype for vector autoregressions with `P` lags
"""
abstract type AbstractVectorAutoRegression{P} <: StatisticalModel end



"""
    VariableID = Union{<:AbstractString, Symbol}
"""
const VariableID = Union{<:AbstractString, Symbol}



include("formula.jl")
include("wold.jl")
include("model.jl")
include("structural.jl")
include("impulse.jl")
include("utilities.jl")





end # module
