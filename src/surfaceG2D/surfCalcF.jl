function Get_H0(k₀::Number,n₁::Number,X::Vector{Float64},Y::Vector{Float64},α::Number)
    return exp.(im*k₀*n₁.*(cos(α)*X + sin(α)*Y))
end

function Get_H0(k₀::Number,n₁::Number,X::Vector{Float64},Y::Vector{Float64},layer::layerstructure,α::Number)
    
    rt = rtcoeffs(layer,k₀,[k₀*n₁*sin(π/2 - α),],"up")
    r = rt.r.TM[1]#; t = rt.t.TM[1]
    
    return @. exp(im*k₀*n₁*(cos(α)*X+sin(α)*Y)) + r*exp(-im*k₀*n₁*(cos(α)*X+sin(α)*Y))
end

function GetH0Arr(N::Number,SArr::Matrix{Vector{Float64}},k₀::Number,n₁::Number,α::Number)
    αR =(π*α)/180 
    H₀Arr = zeros(N,size(SArr,2))*im

    for j in axes(SArr,2) 
        SVec = SArr[:,j]
        X = [SVec[i][1] for i in 1:N]; Y = [SVec[i][2] for i in 1:N]
    
        H₀Arr[:,j] = Get_H0(k₀,n₁,X,Y,αR)
    end
    return H₀Arr
end

function GetH0Arr(N::Number,SArr::Matrix{Vector{Float64}},k₀::Number,n₁::Number,layer::layerstructure,α::Number)
    αR =(π*α)/180 
    H₀Arr = zeros(N,size(SArr,2))*im

    for j in axes(SArr,2) 
        SVec = SArr[:,j]
        X = [SVec[i][1] for i in 1:N]; Y = [SVec[i][2] for i in 1:N]
    
        H₀Arr[:,j] = Get_H0(k₀,n₁,X,Y,layer,αR)
    end
    return H₀Arr
end

function Get_E0(k₀::Number,n₁::Number,X::Vector{Float64},Y::Vector{Float64},α::Number)
    μ₀ = 1.25663706144e-6
    ε₀ = 8.85418781762e-12

    αᵣ = (π*α)/180

    H₀Arr = Get_H0(k₀,n₁,X,Y,αᵣ)

    E₀Arr = Array{Array{ComplexF64,1}}(undef,length(H₀Arr))

    for i ∈ eachindex(E₀Arr)
        E₀Arr[i] = √(μ₀/ε₀)*(1/n₁)*H₀Arr[i]*[-sin(αᵣ); cos(αᵣ)]
    end

    E₀Arr
end

function getD(N::Int64)
    D = Matrix(Bidiagonal(zeros(N),ones(N-1),:U))
    D[N,1] = 1.0
    
    Dₕ = D; Dϕ = D
    Dₕ,Dϕ
end

function Δ(i::Int64,SArr::Matrix{Vector{Float64}})
    ny = size(SArr,2)
    
    S1a = SArr[:,1]
    S1e = SArr[:,ny]
    
    return norm(S1e[i] - S1a[i])
end

function fpol(m::Int64,v::Int64,x::Number)
    xst,xend=0.0,1.0
    xₖ = collect(LinRange(xst,xend,m+1))

    mult1,mult2 = 1.0,1.0

    for k in 0:m
        if k ≠ v
            mult1 *= (xₖ[k+1] -x); mult2 *= (xₖ[k+1] - xₖ[v+1])
        end
    end
    mult1/mult2
end

# =====================================================================================
# =====================================================================================

function AelementQuad(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,Grf::GreenFunctions,
    sArr::Matrix{Vector{Float64}},rpos::Vector{Float64},struc::Structure)

    vMax = size(sArr,2)

    gFun = Grf.GrFuncs[u].gTot

    if i ≠ j
        A,_ = quadgk(x -> gFun(sArr[i,v+1],SQuad(x,j,struc),k₀)*Δ(j,sArr)*fpol(m,vˡ,x), 
            0.0, 1.0, rtol=1e-8)
    else 
        if v == 0 || v == vMax-1
            A,_ = quadgk(x -> gFun(sArr[i,v+1],SQuad(x,j,struc),k₀)*Δ(j,sArr)*fpol(m,vˡ,x), 
                0.0, 1.0, rtol=1e-8)
        else
            A,_ = quadgk(x -> gFun(sArr[i,v+1],SQuad(x,j,struc),k₀)*Δ(j,sArr)*fpol(m,vˡ,x), 
                0.0, rpos[v+1],1.0, rtol=1e-8)
        end
    end

A
end

function AelementSum(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,Grf::GreenFunctions,
    sArr::Matrix{Vector{Float64}},struc::Structure)

    gFun = Grf.GrFuncs[u].gTot
    tGr = 0:0.05:1.0

    s = 0.0*im

    for t ∈ tGr
        s += gFun(sArr[i,v+1],SQuad(t,j,struc),k₀)*fpol(m,vˡ,t)
    end

    Δ(j,sArr)*s/length(tGr)

end

function constBel(i::Int64,j::Int64,v::Int64,vˡ::Int64,u::Int64,sV::Vector{Float64},layer::layerstructure,k0::Number,Opt::String)
    
    δ(i::Int64,j::Int64) = i == j ? 1 : 0
    δ(Sᵢ::Vector{Float64}) = [0.;1.]⋅Sᵢ == 0 ? 1 : 0
    
    if Opt == "Homo"
        return δ(i,j)*δ(v,vˡ)*0.5*(δ(u,1)-δ(u,2))
    else
        ε1 = layer.mat[1].ε
        εL = layer.mat[end].ε
        
        return δ(i,j)*δ(v,vˡ)*0.5*(δ(u,1)-δ(u,2)) + 0.5*(δ(i,j)*δ(v,vˡ)*δ(sV)*δ(u,1)*((εL(k0)-ε1(k0))/(εL(k0)+ε1(k0))))
    end
end

function BelementQuad(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,Grf::GreenFunctions,
    sArr::Matrix{Vector{Float64}},rpos::Vector{Float64},struc::Structure,layer::layerstructure,Opt::String)
    es = 1e-6
    vMax = size(sArr,2)

    #δ(i,j) = i == j ? 1 : 0

    cte = constBel(i,j,v,vˡ,u,sArr[i,v+1],layer,k₀,Opt)

    DgFun = Grf.GrFuncs[u].DgTot

    if i ≠ j
        B,_ = quadgk(x -> (conj(DgFun(sArr[i,v+1],SQuad(x,j,struc),k₀)) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
            0.0, 1.0, rtol=1e-8)
    else 
        if v == 0 || v == vMax-1
            B,_ = quadgk(x -> (conj(DgFun(sArr[i,v+1],SQuad(x,j,struc),k₀)) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
                0.0+es, 1.0-es, rtol=1e-8)
    
            B = B + cte
        else
            B1,_ = quadgk(x -> (conj(DgFun(sArr[i,v+1],SQuad(x,j,struc),k₀)) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
                0.0+es,rpos[v+1]-es, rtol=1e-8)
            B2,_ = quadgk(x -> (conj(DgFun(sArr[i,v+1],SQuad(x,j,struc),k₀)) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x), rpos[v+1]+es,1.0-es, rtol=1e-8)
    
            B = B1+B2 + cte
        end
    end

    B
end

function BelementSum(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,Grf::GreenFunctions,
    sArr::Matrix{Vector{Float64}},struc::Structure)

    nᵤ(u) = u == 1 ? n₁ : n₂
    tGr = 0:0.05:1.0
    DgFun = Grf.GrFuncs[u].DgTot
    s = 0.0*im

    for t ∈ tGr
        s += ((conj(DgFun(sArr[i,v+1],SQuad(t,j,struc),k₀))) ⋅ nQuad(t,j,struc))*fpol(m,vˡ,t)
    end

    Δ(j,sArr)*s/length(tGr)

end

function GetMatrixInt(u::Int64,m::Int64,v::Int64,vˡ::Int64,N::Int64,str::Structure,
    k₀::Number,Grf::GreenFunctions,dThr::Number,layer::layerstructure,Opt::String)
    Amat,Bmat = zeros(N,N)*im, zeros(N,N)*im
    rp,sArr = getSvec(m,str)

    for j ∈ 1:N
        for i ∈ 1:N
            dis = norm(sArr[i,v+1]-SQuad(0,j,str))
        
            if dis ≤ dThr
                Amat[i,j] = AelementQuad(i,j,u,m,v,vˡ,k₀,Grf,sArr,rp,str)
                Bmat[i,j] = BelementQuad(i,j,u,m,v,vˡ,k₀,Grf,sArr,rp,str,layer,Opt)
            else
                Amat[i,j] = AelementSum(i,j,u,m,v,vˡ,k₀,Grf,sArr,str)
                Bmat[i,j] = BelementSum(i,j,u,m,v,vˡ,k₀,Grf,sArr,str)
            end
        end
    end

    Amat,Bmat
end

function getFullMatr(u::Int64,m::Int64,str::Structure,k₀::Number,Grf::GreenFunctions,dThr::Number,layer::layerstructure,Opt::String)
    N = numEl(str)
    Dₕ,Dϕ = getD(N)
    AMatrix, BMatrix = zeros(m*N,m*N)*im,zeros(m*N,m*N)*im

    ind = 0

    for v ∈ 1:m
        ind1 = 0
        for vˡ ∈ 1:m
            if vˡ ≠ 1
                 AMatrix[ind+1:ind+N,ind1+1:ind1+N],BMatrix[ind+1:ind+N,ind1+1:ind1+N] = GetMatrixInt(u,m,v-1,vˡ-1,N,str,k₀,Grf,dThr,layer,Opt)
            else
                A0,B0 = GetMatrixInt(u,m,v-1,vˡ-1,N,str,k₀,Grf,dThr,layer,Opt)
                Am,Bm = GetMatrixInt(u,m,v-1,m,N,str,k₀,Grf,dThr,layer,Opt)
            
                AMatrix[ind+1:ind+N,ind1+1:ind1+N] = A0 + Am*Dϕ 
                BMatrix[ind+1:ind+N,ind1+1:ind1+N] = B0 + Bm*Dₕ 
            end
            ind1 += N
        end
        ind += N
    end

    AMatrix, BMatrix
end

function getHϕ(m::Int64,str::Structure,k₀::Number,Grf::GreenFunctions,n₂::Function,n₁::Function,dThr::Number,α::Number,
    layer::layerstructure=layerstructure([material(n₁(k₀)^2,1.0),material(n₁(k₀)^2,1.0)],[0.0,],"up"),Opt::String="Homo")
    N = numEl(str)
    Dₕ,Dϕ = getD(N)


    A1,B1 = getFullMatr(1,m,str,k₀,Grf,dThr,layer,Opt)
    A2,B2 = getFullMatr(2,m,str,k₀,Grf,dThr,layer,Opt)

    _,SArr = getSvec(m,str)
    
    if Opt == "Homo"
        H₀Arr = GetH0Arr(N,SArr[:,1:m],k₀,n₁(k₀),α)
    else
        H₀Arr = GetH0Arr(N,SArr[:,1:m],k₀,n₁(k₀),layer,α)
    end

    M = [(I(m*N) - B1) A1;
            (I(m*N) + B2) -(n₂(k₀)^2/n₁(k₀)^2)*A2]
    b = vcat(reshape(H₀Arr,(:,1)),zeros(m*N))

    solM = M\b

    H = solM[1:m*N]; ϕ = solM[m*N+1:end]
    H = vcat(H,Dₕ*H[1:N]); ϕ = vcat(ϕ,Dϕ*ϕ[1:N]) 
    H = reshape(H,(:,m+1)); ϕ = reshape(ϕ,(:,m+1))

    H,ϕ
end

# =====================================================================================
# =====================================================================================

function Mint(i::Int64,m::Int64,v::Int64,vˡ::Int64,SArr::Matrix{Vector{Float64}})
    Inval,_ = quadgk(x->fpol(m,v,x)*fpol(m,vˡ,x),0,1,rtol=1e-8)
    Δ(i,SArr)*Inval
end

function getHₛ(r::Vector{Float64},R::Number,ϕ::Matrix{ComplexF64},H::Matrix{ComplexF64},k₀::Number,n₁::Number,m::Int64,str::Structure)
    _,SArr= getSvec(m,str)
    N = size(SArr,1)
    rᵤ = r/norm(r)

    sum1=0.0
    for i in 1:N
        for v in  0:m
            valQuad,_ = quadgk(x->exp(-im*k₀*n₁*(rᵤ ⋅ SQuad(x,i,str)))*(ϕ[i,v+1]+H[i,v+1]*im*k₀*n₁*(nQuad(x,i,str) ⋅ rᵤ))*fpol(m,v,x)
                            ,0.0,1.0,rtol=1e-8)
            sum1 += Δ(i,SArr)*valQuad
        end
    end

    (-(1/4)*sqrt(2/(π*k₀*n₁*R))*exp(im*k₀*n₁*R)*exp(im*π/4)*sum1)

end

function getσₐ(m::Int64,k₀::Vector{Float64},matScatter::Vector{MaterialParams},str::Structure,dThr::Float64,α::Float64,Opt::String="Homo";
        layer::layerstructure = layerstructure(matScatter,[0. ,],"up"))
    _,SArr= getSvec(m,str)
    N = size(SArr,1)

    σₐArr = zeros(length(k₀))

    for i ∈ eachindex(k₀)
        k0 = k₀[i]
        n₁ = matScatter[1].n
        n₂ = matScatter[2].n

        if Opt == "Homo"

            Grf = GreenFunctions(matScatter,Opt)
        else 
            SParms = SommerfieldParams(layer,2k0)
            Grf = GreenFunctions(matScatter,Opt;k0=k0,SParms=SParms)
        end


        H,ϕ = getHϕ(m,str,k0,Grf,n₂,n₁,dThr,α,layer,Opt)

        sum1 = 0.0
        for i ∈ 1:N
            for v ∈ 0:m
                for vˡ ∈ 0:m
                    sum1 += Mint(i,m,v,vˡ,SArr)*imag(ϕ[i,v+1]*conj(H[i,vˡ+1]))
                end
            end
        end

        σₐArr[i] = -(1/(k0*n₁(k0)))*sum1
    end
    σₐArr
end

function getσₑ(m::Int64,k₀::Vector{Float64},matScatter::Vector{MaterialParams},str::Structure,dThr::Float64,r::Number,α::Float64,Opt::String="Homo";
    layer::layerstructure = layerstructure(matScatter,[0. ,],"up"))

    σₑArr = zeros(length(k₀))

    Nₛ = 500
    Θₛ = [(2π/Nₛ)*(i-(1/2)) for i ∈ 1:Nₛ]
    R= r*[cos.(Θₛ) sin.(Θₛ)]

    for i ∈ eachindex(k₀)
        k0 = k₀[i]
        n₁ = matScatter[1].n
        n₂ = matScatter[2].n

        if Opt == "Homo"

            Grf = GreenFunctions(matScatter,Opt)
        else 
            SParms = SommerfieldParams(layer,2k0)
            Grf = GreenFunctions(matScatter,Opt;k0=k0,SParms=SParms)
        end

        H,ϕ = getHϕ(m,str,k0,Grf,n₂,n₁,dThr,α,layer,Opt)

        αₛ = (α*π)/180
        Rₛ = r*[cos(αₛ), sin(αₛ)]
        Hₛ = getHₛ(Rₛ,r,ϕ,H,k0,n₁(k0),m,str)

        σₑArr[i] = -2*sqrt(2π/(k0*n₁(k0)))*real(conj(Hₛ)*sqrt(r)*exp(im*k0*n₁(k0)*r)*exp(-im*π/4))
    end
    σₑArr
end

function getσₛ(m::Int64,k₀::Vector{Float64},matScatter::Vector{MaterialParams},str::Structure,dThr::Float64,r::Number,α::Float64,Opt::String="Homo";
    layer::layerstructure = layerstructure(matScatter,[0. ,],"up"))

    σₛArr = zeros(length(k₀))

    Nₛ = 500
    Θₛ = [(2π/Nₛ)*(i-(1/2)) for i ∈ 1:Nₛ]; ΔΘₛ = Θₛ[2] - Θₛ[1]
    R= r*[cos.(Θₛ) sin.(Θₛ)]

    for i ∈ eachindex(k₀)
        k0 = k₀[i]
        n₁ = matScatter[1].n
        n₂ = matScatter[2].n

        if Opt == "Homo"

            Grf = GreenFunctions(matScatter,Opt)
        else 
            SParms = SommerfieldParams(layer,2k0)
            Grf = GreenFunctions(matScatter,Opt;k0=k0,SParms=SParms)
        end

        H,ϕ = getHϕ(m,str,k0,Grf,n₂,n₁,dThr,α,layer,Opt)

        HₛArr = zeros(Nₛ)*im

        for j ∈ eachindex(HₛArr)
            HₛArr[j] = getHₛ(R[j,:],r,ϕ,H,k0,n₁(k0),m,str)
        end

        sInt = 0.0
        for j ∈ 1:Nₛ
            sInt += norm(HₛArr[j])^2*r*ΔΘₛ
        end

        σₛArr[i] = sInt
    end
    σₛArr
end