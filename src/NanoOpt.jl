using LazyGrids,Interpolations,LinearAlgebra,DifferentialEquations,CSV,DataFrames,FFTW,SpecialFunctions

using Base.Threads

abstract type Material end
abstract type Field end
abstract type Structure end
abstract type Line end

include("Material/Material.jl")
include("Optics/Optics.jl")
include("Laser/Laser.jl")
include("Mie/Mie.jl")
include("stratified/stratified.jl")
include("surfaceG2D/GFIEM2Dsc.jl")
include("plasmModes/plasmModes.jl")