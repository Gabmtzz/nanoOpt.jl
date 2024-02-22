function fCoef(kx::Number,ky1::Number,k0::Number,layer::layerstructure)
    rtC = rtcoeffs(layer,k0,[kx,],"up")
    rp = rtC.r.TM[1]
    
    f1 = (r::Vector{Float64}) ->  exp(-im*kx*r[1])*(exp(-im*ky1*r[2]) + rp*exp(im*ky1*r[2]))
    f2 = (r::Vector{Float64}) ->  exp(-im*kx*r[1])*(exp(-im*ky1*r[2]) - rp*exp(im*ky1*r[2]))
    
    f1,f2
end

function getHₛ(θ::Number,Rf::Number,ϕ::Matrix{ComplexF64},H::Matrix{ComplexF64},k0::Number,m::Int64,str::Structure,layer::layerstructure)
    _,SArr= getSvec(m,str)
    
    k₁ = layer.mat[1].k(k0)
    n₁ = layer.mat[1].n(k0)
    
    
    hₛF = 0.0*im
    if θ ≥ 0. && θ < π
        kx = k₁*cos(θ)
        ky1 = kyi(kx,k₁)
        
        f1,f2 = fCoef(kx,ky1,k0,layer)
    
        sum1 = 0.0*im
        for i ∈ 1:numEl(str)
            for v ∈ 0:m
                elI,_ =   quadgk(x -> (ϕ[i,v+1]*f1(SQuad(x,i,str)) - H[i,v+1]*(nQuad(x,i,str)⋅[-im*kx*f1(SQuad(x,i,str)),-im*ky1*f2(SQuad(x,i,str))]))*fpol(m,v,x),0.0,1.0)
                sum1 += Δ(i,SArr)*elI
            end
        end

        hₛF =  (-exp(im*k₁*Rf)/√(Rf))*(exp(im*(π/4))/4)*√(2/(π*k₁))*sum1
        
    elseif θ < 0 && θ > -π
        
        kₙ = layer.mat[end].k(k0)
        nₙ = layer.mat[end].n(k0)
        
        kx = kₙ*cos(θ)
        ky1 = kyi(kx,k₁)
        kyn = kyi(kx,kₙ)
        rtC = rtcoeffs(layer,k₁,[kx,],"up")
        tp = rtC.t.TM[1]
        
        
        sum1 = 0.0*im
        for i ∈ 1:numEl(str)
            for v ∈ 0:m
                elI,_ = quadgk(x -> exp(-im*kx*SQuad(x,i,str)[1])*exp(im*ky1*SQuad(x,i,str)[2])*(ϕ[i,v+1]-H[i,v+1]*(nQuad(x,i,str)⋅[-im*kx,im*ky1]))*fpol(m,v,x),0.0,1.0)
                sum1 += Δ(i,SArr)*elI
            end
        end
        
         hₛF = (-exp(im*kₙ*Rf)/√(Rf))*(exp(im*(π/4))/4)*√(2/(π*kₙ))*(kyn/ky1)*tp*sum1
    end
    
    hₛF
end

function getσOUP(m::Int,k₀,Rf::Number,matScatter::Vector{MaterialParams},str::Structure,dThr::Float64,α::Float64,layer::layerstructure,opt::String="up")
    Δθ = 0.01
    
    n₁ = matScatter[1].n
    n₂ = matScatter[2].n

    X,Y = getSurfPoints(str)
    xm = (maximum(X)-minimum(X))+1
    yp = 2maximum(Y)+2
    ym = 2minimum(Y)
    
    σOUP = Array{Float64}(undef,length(k₀))

    @threads for i ∈ eachindex(k₀)
        k0 = k₀[i]
        SParms = SommerfieldParams(layer,3k0,9e-3)
        Grf = GreenFunctions(matScatter,"layer";xP= (xm , 10), yP = (ym,yp,10), SParms = SParms)
        H,ϕ = getHϕ(m,str,k0,Grf,n₂,n₁,dThr,α,layer,"layer")
    
        θArr  =  opt == "up"  ? (0.1:Δθ:π) : (-0.1:-Δθ:-π)
        c = opt == "up" ? 1.0 : layer.mat[end].n(k0)/layer.mat[1].n(k0)
    
        HFFa = [getHₛ(θArr[i],Rf,ϕ,H,k0,m,str,layer) for i in eachindex(θArr)]

        σOUP[i] =   c*sum(norm(HFFa)^2 *Rf*Δθ)
    
    end
    
    σOUP
end

function getσEXT(m::Int,k₀,Rf::Number,matScatter::Vector{MaterialParams},str::Structure,dThr::Float64,α::Float64,layer::layerstructure)
    #n₁ = matScatter[1].n
    #n₂ = matScatter[2].n
    αᵣ = (π/180)*α
    σEXTu = Array{Float64}(undef,length(k₀))
    σEXTd = Array{Float64}(undef,length(k₀))
    
    X,Y = getSurfPoints(str)
    xm = (maximum(X)-minimum(X))+1
    yp = 2maximum(Y)+2
    ym = 2minimum(Y)

    @threads for i ∈ eachindex(k₀)
        k0 = k₀[i]
        
        k₁ = layer.mat[1].k(k0)
        n₁ = layer.mat[1].n(k0)

        kₙ = layer.mat[end].k(k0)
        nₙ = layer.mat[end].n(k0)

        rtC = rtcoeffs(layer,k0,[k₁*cos(αᵣ),],"up")
        rp = rtC.r.TM[1]; tp = rtC.t.TM[1] 

        θRadR = π/2 - αᵣ 
        θRadRf =  asin((n₁/nₙ)*sin(θRadR))

        αRadT = -θRadRf - π/2
        
        SParms = SommerfieldParams(layer,3k0,9e-3)
        Grf = GreenFunctions(matScatter,"layer";xP= (xm , 10), yP = (ym,yp,10), SParms = SParms)
        H,ϕ = getHϕ(m,str,k0,Grf,matScatter[2].n,matScatter[1].n,dThr,α,layer,"layer")
        
        c = layer.mat[end].n(k0)/layer.mat[1].n(k0)
        σEXTu[i] = -2sqrt(2π/(k₁))*real(rp*conj(getHₛ(αᵣ,Rf,ϕ,H,k0,m,str,layer))*sqrt(Rf)*exp(im*k₁*Rf)*exp(-im*π/4))
        σEXTd[i] = -2*c*sqrt(2π/(kₙ))*real(tp*conj(getHₛ(αRadT,Rf,ϕ,H,k0,m,str,layer))*sqrt(Rf)*exp(im*kₙ*Rf)*exp(-im*π/4))
    end

    σEXTu,σEXTd
end