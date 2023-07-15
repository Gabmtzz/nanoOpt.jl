function Get_H0(k₀::Number,n₁::Number,X::Vector{Float64},Y::Vector{Float64},α::Number)
    return exp.(im*k₀*n₁.*(cos(α)*X + sin(α)*Y))
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


function getD(N::Int64)
    D = Matrix(Bidiagonal(zeros(N),ones(N-1),:U))
    D[N,1] = 1.0
    
    Dₕ = D; Dϕ = D
    Dₕ,Dϕ
end


function Δ(i::Int64,SArr::Matrix{Vector{Float64}})
    S1a = SArr[:,1];  
    if i == length(S1a)
        return norm(S1a[1] - S1a[end])
    else
        return norm(S1a[i+1] - S1a[i])
    end
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

function AelementQuad(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,n₂::Number,n₁::Number,
    sArr::Matrix{Vector{Float64}},rpos::Vector{Float64},struc::Structure)

vMax = size(sArr,2)
nᵤ(u) = u == 1 ? n₁ : n₂

if i ≠ j
A,_ = quadgk(x -> (im/4)*hankelh1(0,k₀*nᵤ(u)* norm(sArr[i,v+1]-SQuad(x,j,struc)))*Δ(j,sArr)*fpol(m,vˡ,x), 
        0.0, 1.0, rtol=1e-8)
else
if v == 0 || v == vMax-1
    A,_ = quadgk(x -> (im/4)*hankelh1(0,k₀*nᵤ(u)* norm(sArr[i,v+1]-SQuad(x,j,struc)))*Δ(j,sArr)*fpol(m,vˡ,x), 
        0.0, 1.0, rtol=1e-8)
else
    A,_ = quadgk(x -> (im/4)*hankelh1(0,k₀*nᵤ(u)* norm(sArr[i,v+1]-SQuad(x,j,struc)))*Δ(j,sArr)*fpol(m,vˡ,x), 
        0.0, rpos[v+1],1.0, rtol=1e-8)
end
end

A
end

function AelementSum(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,n₂::Number,n₁::Number,
    sArr::Matrix{Vector{Float64}},struc::Structure)

nᵤ(u) = u == 1 ? n₁ : n₂
tGr = 0:0.05:1.0

s = 0.0*im

for t ∈ tGr
s += (im/4)*hankelh1(0,k₀*nᵤ(u)* norm(sArr[i,v+1]-SQuad(t,j,struc)))*fpol(m,vˡ,t)
end

Δ(j,sArr)*s/length(tGr)

end

function BelementQuad(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,n₂::Number,n₁::Number,
    sArr::Matrix{Vector{Float64}},rpos::Vector{Float64},struc::Structure)

vMax = size(sArr,2)
nᵤ(u) = u == 1 ? n₁ : n₂
δ(i,j) = i == j ? 1 : 0

if i ≠ j
B,_ = quadgk(x -> -(im/4)*hankelh1(1,k₀*nᵤ(u)*norm(sArr[i,v+1]-SQuad(x,j,struc)))*k₀*nᵤ(u)*(((SQuad(x,j,struc)-sArr[i,v+1])/(norm(SQuad(x,j,struc)-sArr[i,v+1]))) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
    0.0, 1.0, rtol=1e-8)
else
if v == 0 || v == vMax-1
    B,_ = quadgk(x -> -(im/4)*hankelh1(1,k₀*nᵤ(u)*norm(sArr[i,v+1]-SQuad(x,j,struc)))*k₀*nᵤ(u)*(((SQuad(x,j,struc)-sArr[i,v+1])/(norm(SQuad(x,j,struc)-sArr[i,v+1]))) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
    0.0, 1.0, rtol=1e-8)
    
    B = B + δ(v,vˡ)*(1/2)*(δ(u,1) - δ(u,2))
else
    B,_ = quadgk(x -> -(im/4)*hankelh1(1,k₀*nᵤ(u)*norm(sArr[i,v+1]-SQuad(x,j,struc)))*k₀*nᵤ(u)*(((SQuad(x,j,struc)-sArr[i,v+1])/(norm(SQuad(x,j,struc)-sArr[i,v+1]))) ⋅ nQuad(x,j,struc))*Δ(j,sArr)*fpol(m,vˡ,x),
    0.0,rpos[v+1],1.0, rtol=1e-8)
    
    B = B + δ(v,vˡ)*(1/2)*(δ(u,1) - δ(u,2))
end
end

B
end

function BelementSum(i::Int64,j::Int64,u::Int64,m::Int64,v::Int64,vˡ::Int64,k₀::Number,n₂::Number,n₁::Number,
    sArr::Matrix{Vector{Float64}},struc::Structure)

nᵤ(u) = u == 1 ? n₁ : n₂
tGr = 0:0.05:1.0

s = 0.0*im

for t ∈ tGr
s += -(im/4)*hankelh1(1,k₀*nᵤ(u)*norm(sArr[i,v+1]-SQuad(t,j,struc)))*k₀*nᵤ(u)*(((SQuad(t,j,struc)-sArr[i,v+1])/(norm(SQuad(t,j,struc)-sArr[i,v+1]))) ⋅ nQuad(t,j,struc))*fpol(m,vˡ,t)
end

Δ(j,sArr)*s/length(tGr)

end



function GetMatrixInt(u::Int64,m::Int64,v::Int64,vˡ::Int64,N::Int64,str::Structure,
    k₀::Number,n₂::Number,n₁::Number,dThr::Number)
Amat,Bmat = zeros(N,N)*im, zeros(N,N)*im
rp,sArr = getSvec(m,cyl)

for i ∈ 1:N
    for j ∈ 1:N
        dis = norm(sArr[i,v+1]-SQuad(0,j,str))
    
        if dis ≤ dThr
            Amat[i,j] = AelementQuad(i,j,u,m,v,vˡ,k₀,n₂,n₁,sArr,rp,str)
            Bmat[i,j] = BelementQuad(i,j,u,m,v,vˡ,k₀,n₂,n₁,sArr,rp,str)
        else
            Amat[i,j] = AelementSum(i,j,u,m,v,vˡ,k₀,n₂,n₁,sArr,str)
            Bmat[i,j] = BelementSum(i,j,u,m,v,vˡ,k₀,n₂,n₁,sArr,str)
        end
    end
end

Amat,Bmat
end

function getFullMatr(u::Int64,m::Int64,str::Structure,k₀::Number,n₂::Number,n₁::Number,dThr::Number)
N = numEl(str)
Dₕ,Dϕ = getD(N)
AMatrix, BMatrix = zeros(m*N,m*N)*im,zeros(m*N,m*N)*im

ind = 0

for v ∈ 1:m
    ind1 = 0
    for vˡ ∈ 1:m
        if vˡ ≠ 1
             AMatrix[ind+1:ind+N,ind1+1:ind1+N],BMatrix[ind+1:ind+N,ind1+1:ind1+N] = GetMatrixInt(u,m,v-1,vˡ-1,N,str,k₀,n₂,n₁,dThr)
        else
            A0,B0 = GetMatrixInt(u,m,v-1,vˡ-1,N,str,k₀,n₂,n₁,dThr)
            Am,Bm = GetMatrixInt(u,m,v-1,m,N,str,k₀,n₂,n₁,dThr)
        
            AMatrix[ind+1:ind+N,ind1+1:ind1+N] = A0 + Am*Dϕ 
            BMatrix[ind+1:ind+N,ind1+1:ind1+N] = B0 + Bm*Dₕ 
        end
        ind1 += N
    end
    ind += N
end

AMatrix, BMatrix
end

function getHϕ(m::Int64,str::Structure,k₀::Number,n₂::Number,n₁::Number,dThr::Number,α::Number)
N = numEl(str)
Dₕ,Dϕ = getD(N)

A1,B1 = getFullMatr(1,m,str,k₀,n₂,n₁,dThr)
A2,B2 = getFullMatr(2,m,str,k₀,n₂,n₁,dThr)

_,SArr = getSvec(m,str)
H₀Arr = GetH0Arr(N,SArr[:,1:m],k₀,n₁,α)

M = [(I(m*N) - B1) A1;
        (I(m*N) + B2) -(n₂^2/n₁^2)*A2]
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

function getσₐ(m::Int64,k₀::Vector{Float64},mat1::MaterialParams,mat2::MaterialParams,str::Structure,dThr::Float64,α::Float64)
_,SArr= getSvec(m,str)
N = size(SArr,1)

σₐArr = zeros(length(k₀))

for i ∈ eachindex(k₀)
    k0 = k₀[i]
    n₁ = mat1.n(k0)
    n₂ = mat2.n(k0)

    λ = 2π/k0
    H,ϕ = getHϕ(m,str,k0,n₂,n₁,dThr,α)

    sum1 = 0.0
    for i ∈ 1:N
        for v ∈ 0:m
            for vˡ ∈ 0:m
                sum1 += Mint(i,m,v,vˡ,SArr)*imag(ϕ[i,v+1]*conj(H[i,vˡ+1]))
            end
        end
    end

    σₐArr[i] = -(1/(k0*n₁))*sum1
end
σₐArr
end

function getσₑ(m::Int64,k₀::Vector{Float64},mat1::MaterialParams,mat2::MaterialParams,str::Structure,dThr::Float64,r::Number,α::Float64)
_,SArr= getSvec(m,str)
N = size(SArr,1)

σₑArr = zeros(length(k₀))

Nₛ = 500
Θₛ = [(2π/Nₛ)*(i-(1/2)) for i ∈ 1:Nₛ]; ΔΘₛ = Θₛ[2] - Θₛ[1]
R= r*[cos.(Θₛ) sin.(Θₛ)]

for i ∈ eachindex(k₀)
    k0 = k₀[i]
    n₁ = mat1.n(k0)
    n₂ = mat2.n(k0)

    λ = 2π/k0
    H,ϕ = getHϕ(m,str,k0,n₂,n₁,dThr,α)

    αₛ = (α*π)/180
    Rₛ = r*[cos(αₛ), sin(αₛ)]
    Hₛ = getHₛ(Rₛ,r,ϕ,H,k0,n₁,m,str)

    σₑArr[i] = -2*sqrt(2π/(k0*n₁))*real(conj(Hₛ)*sqrt(r)*exp(im*k0*n₁*r)*exp(-im*π/4))
end
σₑArr
end

function getσₛ(m::Int64,k₀::Vector{Float64},mat1::MaterialParams,mat2::MaterialParams,str::Structure,dThr::Float64,r::Number,α::Float64)
_,SArr= getSvec(m,str)
N = size(SArr,1)

σₛArr = zeros(length(k₀))

Nₛ = 500
Θₛ = [(2π/Nₛ)*(i-(1/2)) for i ∈ 1:Nₛ]; ΔΘₛ = Θₛ[2] - Θₛ[1]
R= r*[cos.(Θₛ) sin.(Θₛ)]

for i ∈ eachindex(k₀)
    k0 = k₀[i]
    n₁ = mat1.n(k0)
    n₂ = mat2.n(k0)

    λ = 2π/k0
    H,ϕ = getHϕ(m,str,k0,n₂,n₁,dThr,α)

   HₛArr = zeros(Nₛ)*im

    for j ∈ eachindex(HₛArr)
        HₛArr[j] = getHₛ(R[j,:],r,ϕ,H,k0,n₁,m,str)
    end

    sInt = 0.0
    for j ∈ 1:Nₛ
        sInt += norm(HₛArr[j])^2*r*ΔΘₛ
    end

    σₛArr[i] = sInt
end
σₛArr
end