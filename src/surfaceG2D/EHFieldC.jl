function getIntergalp(i::Int64,k₀::Number,n::Number,m::Int64,v::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    _,SArr = getSvec(m,Str)
    I1,_ = quadgk(x -> (im/4)*hankelh1(0,k₀*n*norm(pos-SQuad(x,i,Str)))*Δ(i,SArr)*ϕ[i,v+1]*fpol(m,v,x),0,1, rtol=1e-8)
    I2,_ = quadgk(x -> -(im/4)*hankelh1(1,k₀*n*norm(pos-SQuad(x,i,Str)))k₀*n*(((SQuad(x,i,Str)-pos)/(norm(SQuad(x,i,Str)-pos))) ⋅ nQuad(x,i,Str))*Δ(i,SArr)*H[i,v+1]*fpol(m,v,x),
            0,1, rtol=1e-8)
    I1-I2
end

function getIntegralScatter(k₀::Number,n::Number,m::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    N = numEl(Str)
    IScatter = 0.0*im

    for i ∈ 1:N
        for v in 0:m
            IScatter += getIntergalp(i,k₀,n,m,v,pos,Str,H,ϕ)
        end
    end
    IScatter
end

function getHfielOutsie(k₀::Number,n₁::Number,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})

    H₀Arr = Get_H0(k₀,n₁,Xout,Yout,(π*α)/180)

    posOut = [Xout Yout]

    HOut = zeros(size(posOut,1))*im

    for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        HOut[p] = H₀Arr[p] - getIntegralScatter(k₀,n₁,m,pos,Str,H,ϕ)
    end

    HOut
end

function getHfieldInside(k₀::Number,n₂::Number,m::Int64,dThr::Number,Str::Structure,Xin::Vector{Float64},Yin::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
     
    posIn = [Xin Yin]

    HIn = zeros(size(posIn,1))*im

    for p ∈ axes(posIn,1)
        pos = posIn[p,:]
        HIn[p] = getIntegralScatter(k₀,n₂,m,pos,Str,H,ϕ)
    end

    HIn
end


function B_Efield(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64})
    -(im/4)*hankelh1(1,k₀*n*norm(r-rs))*k₀*n
end

function A_Efield(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},ns::Vector{Float64})
    B = B_Efield(k₀,n,r,rs)

    (-im/8*(hankelh1(0,k₀*n*norm(r-rs)) - hankelh1(2,k₀*n*norm(r-rs)))*(k₀*n)^2 - (B/norm(r-rs))) * ((rs-r)/norm(rs-r) ⋅ ns)
    
end


function get_zxdndH(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},ns::Vector{Float64})
    A = A_Efield(k₀,n,r,rs,ns)
    B = B_Efield(k₀,n,r,rs)

    xn,yn = [1; 0], [0; 1]

    vx = A*(-(r[2]-rs[2])/norm(r-rs))+B*((yn ⋅ ns)/norm(r-rs))
    vy = A*((r[1]-rs[1])/norm(r-rs))-B*((xn ⋅ ns)/norm(r-rs))

    [vx; vy]
end

function get_zxdH(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64})
    B = B_Efield(k₀,n,r,rs)

    (1/norm(r-rs))*[-(r[2]-rs[2]); r[1]-rs[1]]*B
end

function getIntegralEfp(i::Int64,k₀::Number,n::Number,m::Int64,v::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    _,SArr = getSvec(m,Str)
    
    I1,_ = quadgk(x -> get_zxdH(k₀,n,pos,SQuad(x,i,Str))*ϕ[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-8)
    I2,_ = quadgk(x -> get_zxdndH(k₀,n,pos,SQuad(x,i,Str),nQuad(x,i,Str))*H[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-8)
    I1-I2
end

function getIntegralScatterEf(k₀::Number,n::Number,m::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    N = numEl(Str)
    IScv = [0.0; 0.0]

    for i ∈ 1:N
        for v ∈ 1:m
            IScv += getIntegralEfp(i,k₀,n,m,v,pos,Str,H,ϕ)
        end
    end
    IScv
end

function getEfieldOutside(k₀::Number,n₁::Number,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12

    E₀Arr = Get_E0(k₀,n₁,Xout,Yout,α)
    #E₀Arr = E₀Arr./maximum(norm.(E₀Arr))
    posOut = [Xout Yout]

    ω = k₀/√(μ₀*ε₀)

    Eout = Array{Array{ComplexF64,1}}(undef,size(posOut,1))

    for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        Eout[p] = E₀Arr[p].+ (im/(ω*ε₀*n₁^2))*getIntegralScatterEf(k₀,n₁,m,pos,Str,H,ϕ)
    end
    Eout
end

function getEfieldInside(k₀::Number,n₂::Number,m::Int64,dThr::Number,Str::Structure,Xin::Vector{Float64},Yin::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12

    posIn = [Xin Yin]

    ω = k₀/√(μ₀*ε₀)

    Ein = Array{Array{ComplexF64,1}}(undef,size(posIn,1))

    for p ∈ axes(posIn,1)
        pos = posIn[p,:]
        Ein[p] = -(im/(ω*ε₀*n₂^2))*getIntegralScatterEf(k₀,n₂,m,pos,Str,H,ϕ)
    end
    Ein
end