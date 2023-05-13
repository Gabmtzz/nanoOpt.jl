using LazyGrids,FFTW,SpecialFunctions

abstract type Material end
abstract type Field end


include("Material.jl")
include("Optics/Optics.jl")
include("Laser/Laser.jl")