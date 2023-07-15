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
    
     function MaterialFunc(eps::Any,mu::Number)
        eV2nm = 1 / 8.0655477e-4
        epsfun = x -> eps(eV2nm/(2π)*x)
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

function epstable(name::String)
    bucket = "/home/martinez/Documents/repos/nanoOpt.jl/src/Material/"
    
    data = CSV.read(bucket*name*".dat",DataFrame, header=false, delim='\t')
    ene,n,k = data[:,1],data[:,2],data[:,3]
    eneG = ene[1]: 0.12416666666666654:ene[end]
    method = BSpline(Cubic(Natural(OnGrid())))
    eps = scale(interpolate((n+im.*k).^2 , method),eneG)
    epsf = x -> eps(x)
    
    epsf
end

function epsdrude(name::String)

    if name == "Ag"
        εb,ωₚ,τ = 3.3,10.0,30.
    elseif name == "Au"
        εb,ωₚ,τ = 10.0,10.0,10.0
    elseif name == "Al"
        εb,ωₚ,τ = 1.0,15.0,1.0
    end

    γ = 0.66/τ

    epsd = x ->  εb - ωₚ^2/(x*(x+im*γ))

    epsd
end