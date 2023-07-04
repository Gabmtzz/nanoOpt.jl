using LazyGrids,Interpolations,DifferentialEquations,CSV,DataFrames,FFTW,SpecialFunctions

abstract type Material end
abstract type Field end
abstract type Structure end


include("Material/Material.jl")
include("Optics/Optics.jl")
include("Laser/Laser.jl")
include("stratified/stratified.jl")