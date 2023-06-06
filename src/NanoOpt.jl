using LazyGrids,Interpolations,CSV,DataFrames,FFTW,SpecialFunctions

abstract type Material end
abstract type Field end
abstract type Structure end


include("Material.jl")
include("Optics/Optics.jl")
include("Laser/Laser.jl")