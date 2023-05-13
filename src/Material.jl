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