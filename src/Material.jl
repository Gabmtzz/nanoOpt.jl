using LazyGrids
abstract type Material end

struct MaterialFunc <: Material
    eps::Function
    mu::Function
    
    function MaterialFunc(eps::Function,mu::Function)
        eV2nm = 1 / 8.0655477e-4
        epsfun = x -> eps(eV2nm/(2π)*x)
        mufun = x -> mu(eV2nm/(2π)*x)
        
        new(epsfun,mufun)
    end
    
    function MaterialFunc(eps::Number,mu::Number)
        epsfun = x -> eps
        mufun = x -> mu
        
        new(epsfun,mufun)
    end
end

struct MaterialParams <: Material
    ε::Function
    μ::Function
    n::Function
    Z::Function
    k::Function
end

function material(eps,mu)
    fM = MaterialFunc(eps,mu)
    
    ε = x -> fM.eps(x); μ = x ->  fM.mu(x); n = x -> sqrt(fM.eps(x)*fM.mu(x))
    Z = x -> sqrt(fM.mu(x) / fM.eps(x)); k = x -> x*sqrt(fM.eps(x)*fM.mu(x))
    
    MaterialParams(ε,μ,n,Z,k)
end

abstract type Field end

struct dipole <: Field
    mat::MaterialParams
    k₀::Number
    pos::Vector{Number}
    dip::Vector{Number}
    
    function dipole(material::MaterialParams,k₀::Number,pos,dip)
        new(material,k₀,pos,dip)
    end
end

struct Angular <: Field
    mat::MaterialParams
    iefield::Array{ComplexF64, 3}
    k0::Number
    dir::Number
    x::Vector{Number}
    y::Vector{Number}
    kvec::Vector{Matrix{ComplexF64}}
    kx::Vector{Number}
    ky::Vector{Number}
    
    function Angular(x,y,efield, mat, k0,dir,npad)
        k = mat.k(k0)
        Mx,My = npad*length(x),npad*length(y)
        iefield,kx,ky = fourier2(x,y,efield,Mx,My)
        kxg,kyg= ndgrid(kx,ky)
        kz = dir*sqrt.(Complex.(k^2 .- kxg.^2 .- kyg.^2))
        kvec = [Complex.(kxg),Complex.(kyg),kz]
        
        new(mat,iefield,k0,dir,collect(x),collect(y),kvec,kx,ky)
    end
end