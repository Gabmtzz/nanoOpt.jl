function getIntergalp(i::Int64,k₀::Number,GA::greenFunct,m::Int64,v::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    Gf = GA.gTot
    dGf = GA.DgTot

    _,SArr = getSvec(m,Str)

    I1,_ = quadgk(x -> Gf(pos,SQuad(x,i,Str),k₀)*Δ(i,SArr)*ϕ[i,v+1]*fpol(m,v,x),0,1, rtol=1e-8)
    I2,_ = quadgk(x -> (conj(dGf(pos,SQuad(x,i,Str),k₀)) ⋅ nQuad(x,i,Str))*Δ(i,SArr)*H[i,v+1]*fpol(m,v,x),
            0,1, rtol=1e-8)
    I1-I2
end

function getIntegralScatter(k₀::Number,mat::MaterialParams,m::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    GH = greenHomo(mat.k)
    
    N = numEl(Str)
    IScatter = 0.0*im

    for i ∈ 1:N
        for v in 0:m
            IScatter += getIntergalp(i,k₀,GH,m,v,pos,Str,H,ϕ)
        end
    end
    IScatter
end

function getIntegralScatter(k₀::Number,layer::layerstructure,m::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},opt::String)
    GA = greenAnalytic(k₀,layer,opt)
    
    N = numEl(Str)
    IScatter = 0.0*im

    for i ∈ 1:N
        for v in 0:m
            IScatter += getIntergalp(i,k₀,GA,m,v,pos,Str,H,ϕ)
        end
    end
    IScatter
end

function kAngFA(α::Float64,Eₗ::Float64,Eₕ::Float64)
    sign(α)*(1-cos(α))*Eₗ/2 - im*Eₕ*sin(α)
end

function dkAngFA(α::Float64,Eₗ::Float64,Eₕ::Float64)
    abs(sin(α))*Eₗ/2 - im*Eₕ*cos(α)#cos
end

function getIntegralHk(kx::Number,k₁::Number,ky1::Number,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    _,SArr = getSvec(m,Str)
    N = numEl(Str)
    IScatter = 0.0*im

    for i ∈ 1:N
        for v ∈ 0:m
            Iv,_ = quadgk(x->((exp(-im*kx*SQuad(x,i,Str)[1])*exp(im*ky1*SQuad(x,i,Str)[2]))/ky1)*(ϕ[i,v+1]-H[i,v+1]*(nQuad(x,i,Str)⋅[-im*kx,im*ky1]))*Δ(i,SArr)*fpol(m,v,x),0,1,rtol=1e-3)
            IScatter += Iv
        end
    end

    (-im/(4π))*IScatter
end

function getIntegralHkFF(kx::Number,k₁::Number,ky1::Number,rp::Number,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},yT::Number)
    _,SArr = getSvec(m,Str)
    N = numEl(Str)
    IScatter = 0.0*im

    for i ∈ 1:N
        for v ∈ 0:m
            Iv,_ = quadgk(x->((exp(-im*kx*SQuad(x,i,Str)[1])*(exp(im*ky1*(yT-SQuad(x,i,Str)[2]))+
                                rp*exp(im*ky1*(yT+SQuad(x,i,Str)[2]))))/ky1)*(ϕ[i,v+1]-H[i,v+1]*(nQuad(x,i,Str)⋅[-im*kx,im*ky1]))*Δ(i,SArr)*fpol(m,v,x),0,1,rtol=1e-3)
            IScatter += Iv
        end
    end

    (-im/(4π))*IScatter
end

function IntegUp(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    k₁ = layer.mat[1].k(k₀)
    ky1 = kyi(kx,k₁)
    
    εₗₑ = layer.mat[end].ε(k₀)
    ε₁ = layer.mat[1].ε(k₀)
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    rp = rtC.r.TM[1]
    
    getIntegralHk(kx,k₁,ky1,m,Str,H,ϕ)*(rp-(εₗₑ - ε₁)/(εₗₑ + ε₁))*exp(im*kx*pos[1])*exp(im*ky1*pos[2])
end

function IntegUp(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},yT::Number)
    k₁ = layer.mat[1].k(k₀)
    ky1 = kyi(kx,k₁)
    
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    rp = rtC.r.TM[1]
    
    getIntegralHkFF(kx,k₁,ky1,rp,m,Str,H,ϕ,yT)*exp(im*kx*pos[1])*exp(im*ky1*(pos[2]-yT))
end

function IntegDw(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    k₁ = layer.mat[1].k(k₀); kₑ = layer.mat[end].k(k₀) 
    ε1,εL = layer.mat[1].ε(k₀),layer.mat[end].ε(k₀)
    
    ky1 = kyi(kx,k₁); kyₑ = kyi(kx,kₑ); ky1 = kyi(kx,k₁); kye = kyi(kx,kₑ) 
    
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    tp = rtC.t.TM[1]
    
    getIntegralHk(kx,k₁,ky1,m,Str,H,ϕ)*(tp*exp(-im*kye*pos[2])-(2εL/(εL+ε1))*exp(-im*ky1*pos[2]))*exp(im*kx*pos[1])
end

function Hins(k₀::Number,pos::Vector{Float64},SParms::SommerfieldParams,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},ps::String)
    layer = SParms.layer
    EL = SParms.EL
    EH = SParms.EH
    
    if ps == "up"
        IvNeg,_ = quadgk(x -> IntegUp(x,k₀,pos,layer,m,Str,H,ϕ),-0.5,-EL,rtol=1e-3)
        IAng,_ = quadgk(x -> IntegUp(kAngFA(x,EL,EH),k₀,pos,layer,m,Str,H,ϕ)*dkAngFA(x,EL,EH),-π,π,rtol=1e-3)
        IvPos,_ = quadgk(x -> IntegUp(x,k₀,pos,layer,m,Str,H,ϕ),EL,0.5,rtol=1e-3)
    elseif ps == "down"
        IvNeg,_ = quadgk(x -> IntegDw(x,k₀,pos,layer,m,Str,H,ϕ),-0.5,-EL,rtol=1e-3)
        IAng,_ = quadgk(x -> IntegDw(kAngFA(x,EL,EH),k₀,pos,layer,m,Str,H,ϕ)*dkAngFA(x,EL,EH),-π,π,rtol=1e-3)
        IvPos,_ = quadgk(x -> IntegDw(x,k₀,pos,layer,m,Str,H,ϕ),EL,0.5,rtol=1e-3)
    end
    
    IvNeg+IAng+IvPos
end

function HFFup(k₀::Number,pos::Vector{Float64},SParms::SommerfieldParams,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},yT::Number)
    layer = SParms.layer
    EL = SParms.EL
    EH = SParms.EH
    
    IvNeg,_ = quadgk(x -> IntegUp(x,k₀,pos,layer,m,Str,H,ϕ,yT),-0.5,-EL)
    IAng,_ = quadgk(x -> IntegUp(kAngFA(x,EL,EH),k₀,pos,layer,m,Str,H,ϕ,yT)*dkAngFA(x,EL,EH),-π,π)
    IvPos,_ = quadgk(x -> IntegUp(x,k₀,pos,layer,m,Str,H,ϕ,yT),EL,0.5)
    
    IvNeg+IAng+IvPos
end

function getHfielOutsie(k₀::Number,mat::MaterialParams,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})

    H₀Arr = Get_H0(k₀,mat.n(k₀),Xout,Yout,(π*α)/180)

    posOut = [Xout Yout]

    HOut = zeros(size(posOut,1))*im

    @threads for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        HOut[p] = H₀Arr[p] - getIntegralScatter(k₀,mat,m,pos,Str,H,ϕ)
    end

    HOut
end

function getHfielOutsie(k₀::Number,SParms::SommerfieldParams,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})

    layer = SParms.layer

    H₀Arr = Get_H0(k₀,Xout,Yout,layer,(π*α)/180)

    posOut = [Xout Yout]

    HOut = zeros(size(posOut,1))*im

    @threads for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        opt = pos[2] ≤ 0 ? "down" : "up"
        HOut[p] = H₀Arr[p] - getIntegralScatter(k₀,layer,m,pos,Str,H,ϕ,opt) + Hins(k₀,pos,SParms,m,Str,H,ϕ,opt)
    end

    HOut
end

function getHfieldInside(k₀::Number,mat::MaterialParams,m::Int64,dThr::Number,Str::Structure,Xin::Vector{Float64},Yin::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
     
    posIn = [Xin Yin]

    HIn = zeros(size(posIn,1))*im

    @threads for p ∈ axes(posIn,1)
        pos = posIn[p,:]
        HIn[p] = getIntegralScatter(k₀,mat,m,pos,Str,H,ϕ)
    end

    HIn
end

# =====================================================================================================================
# =====================================================================================================================
# =====================================================================================================================

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
#-----------------------------------------------------------------
function B_Efield_2(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},c::Number)
    -(im/4)*c*hankelh1(1,k₀*n*norm(r-rs))*k₀*n
end

function A_Efield_2(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},ns::Vector{Float64},c::Number)
    B = B_Efield_2(k₀,n,r,rs,c)

    (-im/8*c*(hankelh1(0,k₀*n*norm(r-rs)) - hankelh1(2,k₀*n*norm(r-rs)))*(k₀*n)^2 - (B/norm(r-rs))) * ((rs-r)/norm(rs-r) ⋅ ns)
    
end

function get_zxdndH_2(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},ns::Vector{Float64},c::Number)
    A = A_Efield_2(k₀,n,r,rs,ns,c)
    B = B_Efield_2(k₀,n,r,rs,c)

    xn,yn = [1; 0], [0; 1]

    vx = A*(-(r[2]-rs[2])/norm(r-rs))+B*((yn ⋅ ns)/norm(r-rs))
    vy = A*((r[1]-rs[1])/norm(r-rs))-B*((xn ⋅ ns)/norm(r-rs))

    [vx; vy]
end

function get_zxdH_2(k₀::Number,n::Number,r::Vector{Float64},rs::Vector{Float64},c::Number)
    B = B_Efield_2(k₀,n,r,rs,c)

    (1/norm(r-rs))*[-(r[2]-rs[2]); r[1]-rs[1]]*B
end

function get_zxdndH_L(k₀::Number,layer::layerstructure,r::Vector{Float64},rs::Vector{Float64},ns::Vector{Float64},op::String)
    mVec(r::Vector{Float64}) = r.*[1,-1]

    n₁ = layer.mat[1].n(k₀)
    ε₁,εLf = layer.mat[1].ε(k₀),layer.mat[end].ε(k₀)

    if op == "up"
        zxdndHT = get_zxdndH(k₀,n₁,r,rs,ns) + get_zxdndH_2(k₀,n₁,mVec(r),rs,ns,(εLf-ε₁) / (εLf+ε₁))
    else

        zxdndHT =  get_zxdndH_2(k₀,n₁,r,rs,ns,(2εLf) / (εLf+ε₁))
    end
    #c = (εLf-ε₁) / (εLf+ε₁)

    zxdndHT
end

function get_zxdH_L(k₀::Number,layer::layerstructure,r::Vector{Float64},rs::Vector{Float64},op::String)
    mVec(r::Vector{Float64}) = r.*[1,-1]
    n₁ = layer.mat[1].n(k₀)
    ε₁,εLf = layer.mat[1].ε(k₀),layer.mat[end].ε(k₀)
    c = (εLf-ε₁) / (εLf+ε₁)

    if op == "up"
        zxdHT = get_zxdH(k₀,n₁,r,rs) + get_zxdH_2(k₀,n₁,mVec(r),rs,(εLf-ε₁) / (εLf+ε₁))
    else

        zxdHT =  get_zxdH_2(k₀,n₁,mVec(r),rs,(2εLf) / (εLf+ε₁))
    end

    zxdHT
end


#-----------------------------------------------------------------


function getIntegralEfp(i::Int64,k₀::Number,n::Number,m::Int64,v::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    _,SArr = getSvec(m,Str)
    
    I1,_ = quadgk(x -> get_zxdH(k₀,n,pos,SQuad(x,i,Str))*ϕ[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-3)
    I2,_ = quadgk(x -> get_zxdndH(k₀,n,pos,SQuad(x,i,Str),nQuad(x,i,Str))*H[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-3)
    I1-I2
end

function getIntegralEfp(i::Int64,k₀::Number,layer::layerstructure,m::Int64,v::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},opt::String)
    _,SArr = getSvec(m,Str)
    
    I1,_ = quadgk(x -> get_zxdH_L(k₀,layer,pos,SQuad(x,i,Str),opt)*ϕ[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-3)
    I2,_ = quadgk(x -> get_zxdndH_L(k₀,layer,pos,SQuad(x,i,Str),nQuad(x,i,Str),opt)*H[i,v+1]*Δ(i,SArr)*fpol(m,v,x),0,1, rtol=1e-3)
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

function getIntegralScatterEf(k₀::Number,layer::layerstructure,m::Int64,pos::Vector{Float64},Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},opt::String)
    N = numEl(Str)
    IScv = [0.0; 0.0]

    for i ∈ 1:N
        for v ∈ 1:m
            IScv += getIntegralEfp(i,k₀,layer,m,v,pos,Str,H,ϕ,opt)
        end
    end
    IScv
end

function IntegEUp(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    k₁ = layer.mat[1].k(k₀)
    ky1 = kyi(kx,k₁)
    
    εₗₑ = layer.mat[end].ε(k₀)
    ε₁ = layer.mat[1].ε(k₀)
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    rp = rtC.r.TM[1]
    
    getIntegralHk(kx,k₁,ky1,m,Str,H,ϕ)*(rp-(εₗₑ - ε₁)/(εₗₑ + ε₁))*exp(im*kx*pos[1])*exp(im*ky1*pos[2])*[im*ky1,-im*kx]
end

function IntegEUp(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},yT::Number)
    k₁ = layer.mat[1].k(k₀)
    ky1 = kyi(kx,k₁)
    
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    rp = rtC.r.TM[1]
    
    getIntegralHkFF(kx,k₁,ky1,rp,m,Str,H,ϕ,yT)*exp(im*kx*pos[1])*exp(im*ky1*(pos[2]-yT))*[im*ky1,-im*kx]
end

function IntegEDw(kx::Number,k₀::Number,pos::Vector{Float64},layer::layerstructure,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    k₁ = layer.mat[1].k(k₀); kₑ = layer.mat[end].k(k₀) 
    ε1,εL = layer.mat[1].ε(k₀),layer.mat[end].ε(k₀)
    
    ky1 = kyi(kx,k₁); kye = kyi(kx,kₑ) 
    
    rtC = rtcoeffs(layer,k₀,[kx,],"up")
    tp = rtC.t.TM[1]
    
    getIntegralHk(kx,k₁,ky1,m,Str,H,ϕ)*exp(im*kx*pos[1])*[-im*(kye*tp*exp(-im*kye*pos[2])-ky1*(2εL/(εL+ε1))*exp(-im*ky1*pos[2])),-im*kx*(tp*exp(-im*kye*pos[2])-(2εL/(εL+ε1))*exp(-im*ky1*pos[2]))]
end

function Eins(k₀::Number,pos::Vector{Float64},SParms::SommerfieldParams,m::Int64,Str::Structure,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64},ps::String)
    layer = SParms.layer
    EL = SParms.EL
    EH = SParms.EH
    
    if ps == "up"
        IvNeg,_ = quadgk(x -> IntegEUp(x,k₀,pos,layer,m,Str,H,ϕ),-0.5,-EL,rtol=1e-3)
        IAng,_ = quadgk(x -> IntegEUp(kAngFA(x,EL,EH),k₀,pos,layer,m,Str,H,ϕ)*dkAngFA(x,EL,EH),-π,π,rtol=1e-3)
        IvPos,_ = quadgk(x -> IntegEUp(x,k₀,pos,layer,m,Str,H,ϕ),EL,0.5,rtol=1e-3)
    elseif ps == "down"
        IvNeg,_ = quadgk(x -> IntegEDw(x,k₀,pos,layer,m,Str,H,ϕ),-0.5,-EL,rtol=1e-3)
        IAng,_ = quadgk(x -> IntegEDw(kAngFA(x,EL,EH),k₀,pos,layer,m,Str,H,ϕ)*dkAngFA(x,EL,EH),-π,π,rtol=1e-3)
        IvPos,_ = quadgk(x -> IntegEDw(x,k₀,pos,layer,m,Str,H,ϕ),EL,0.5,rtol=1e-3)
    end
    
    IvNeg+IAng+IvPos
end

function getEfieldOutside(k₀::Number,n₁::Number,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12

    E₀Arr = Get_E0(k₀,n₁,Xout,Yout,α)
    #E₀Arr = E₀Arr./maximum(norm.(E₀Arr))
    posOut = [Xout Yout]

    ω = k₀/√(μ₀*ε₀)

    Eout = Array{Array{ComplexF64,1}}(undef,size(posOut,1))

    @threads for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        Eout[p] = E₀Arr[p].+ (im/(ω*ε₀*n₁^2))*getIntegralScatterEf(k₀,n₁,m,pos,Str,H,ϕ)
    end
    Eout
end

function getEfieldOutside(k₀::Number,SParms::SommerfieldParams,m::Int64,dThr::Number,Str::Structure,Xout::Vector{Float64},Yout::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12
    layer = SParms.layer
    
    n₁ = layer.mat[1].n(k₀)
    E₀Arr = Get_E0(k₀,layer,Xout,Yout,α)
    #E₀Arr = E₀Arr./maximum(norm.(E₀Arr))
    posOut = [Xout Yout]

    ω = k₀/√(μ₀*ε₀)

    Eout = Array{Array{ComplexF64,1}}(undef,size(posOut,1))

    @threads for p ∈ axes(posOut,1)
        pos = posOut[p,:]
        opt = pos[2] ≤ 0 ? "down" : "up"
        #Eout[p] = E₀Arr[p].+ (im/(ω*ε₀*n₁^2))*(getIntegralScatterEf(k₀,layer,m,pos,Str,H,ϕ,opt) + Eins(k₀,pos,SParms,m,Str,H,ϕ,opt))
        Eout[p] = E₀Arr[p].+ (im/(ω*ε₀*n₁^2))*getIntegralScatterEf(k₀,layer,m,pos,Str,H,ϕ,opt) .+ (im/(ω*ε₀*n₁^2))*Eins(k₀,pos,SParms,m,Str,H,ϕ,opt)
        #Eout[p] = E₀Arr[p].+ (im/(ω*ε₀*n₁^2))*getIntegralScatterEf(k₀,layer,m,pos,Str,H,ϕ)
    end
    Eout
end

function getEfieldInside(k₀::Number,n₂::Number,m::Int64,dThr::Number,Str::Structure,Xin::Vector{Float64},Yin::Vector{Float64},α::Number,H::Matrix{ComplexF64},ϕ::Matrix{ComplexF64})
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12

    posIn = [Xin Yin]

    ω = k₀/√(μ₀*ε₀)

    Ein = Array{Array{ComplexF64,1}}(undef,size(posIn,1))

    @threads for p ∈ axes(posIn,1)
        pos = posIn[p,:]
        Ein[p] = -(im/(ω*ε₀*n₂^2))*getIntegralScatterEf(k₀,n₂,m,pos,Str,H,ϕ)
    end
    Ein
end