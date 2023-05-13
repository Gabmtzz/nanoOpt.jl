function efocusHG00(lens::lensfocus,w₀::Number,x::Matrix{Float64},y::Matrix{Float64},z::Number)
    nquad = 101
    
    z = z*ones(size(x))
    k0,n = lens.k0,lens.mat.n(lens.k0)
    tmax = asin(lens.NA/lens.mat.n(k0))
    t,w,_ = lglnodes(nquad)
    t,w = 0.5*(t.+1)*tmax, 0.5*w*tmax
    
    ϕa,ρa = zeros(size(z)),zeros(size(z))

    for i in axes(ϕa,1)
        for j in axes(ϕa,2)
            ϕa[i,j],ρa[i,j] = cart2pol(x[i,j],y[i,j])
        end
    end
    
    I₀ = fun(w₀,z,ρa,0,n*k0,1 .+cos.(t),t,w)
    I₁ = fun(w₀,z,ρa,1,n*k0,sin.(t),t,w)
    I₂ = fun(w₀,z,ρa,2,n*k0,1 .- cos.(t),t,w)
    
    e₁ = 0.5*(I₀- I₂.*cos.(2*ϕa))
    e₂ = 0.5*(-I₂.*sin.(2*ϕa))
    e₃ = 0.5*(-I₁.*cos.(ϕa))

    im*sqrt(n)*k₀*reshape(cat(e₁,e₂,e₃,dims=2),(size(e₁,1),size(e₁,2),3))
end

function fun(w0::Number,z::Matrix{Float64},ρ::Matrix{Float64},n::Int,k::Number,g::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    siz = size(z)
    ρ,z = reshape(ρ,(1,:)),reshape(z,(1,:))
    
    f = im^n* @. exp(-im*k*cos(t)*z)*besselj(n,k*sin(t)*ρ) * exp(-sin(t)^2/w0^2)*sin(t)*sqrt(cos(t))*g
    
    reshape(w'*f,(siz))
end