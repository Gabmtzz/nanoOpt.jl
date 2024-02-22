abstract type greenFunct end

struct greenHomo <: greenFunct
    gTot::Function
    DgTot::Function

    function greenHomo(k::Function)

        gTot = (r::Vector{Float64},rl::Vector{Float64},k0::Number) -> (im/4)*hankelh1(0,k(k0)*norm(r-rl))
        DrTot = (r::Vector{Float64},rl::Vector{Float64},k0::Number) -> (-im/4)*hankelh1(1,k(k0)*norm(r-rl))*((rl-r)/norm(rl-r))*k(k0)

        new(gTot,DrTot)
    end
end

struct greenAnalytic <: greenFunct
    gTot::Function
    DgTot::Function

    function greenAnalytic(k0::Number,layer::layerstructure,opt::String="up")
        
        mVec(r::Vector{Float64}) = r.*[1,-1]
    
        matA = layer.mat
        ε1f,εLf = matA[1].ε,matA[end].ε
        k1 = matA[1].k

        gHomo = greenHomo(k1)

        if opt == "up"
            gA  = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> gHomo.gTot(r,rl,k)  + ((εLf(k)-ε1f(k))/(εLf(k)+ε1f(k)))*(im/4)*hankelh1(0,k1(k)*norm(mVec(r)-rl))
            dgA = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> gHomo.DgTot(r,rl,k) - ((εLf(k)-ε1f(k))/(εLf(k)+ε1f(k)))*(im/4)*hankelh1(1,k1(k)*norm(mVec(r)-rl))*((rl-mVec(r))/norm(rl-mVec(r)))*k1(k)
        elseif opt == "down"
            gA  = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> ((2εLf(k))/(εLf(k)+ε1f(k)))*gHomo.gTot(r,rl,k)
            dgA = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> ((2εLf(k))/(εLf(k)+ε1f(k)))*gHomo.DgTot(r,rl,k)
        end

        new(gA,dgA)
    end
end

struct SommerfieldParams <: greenFunct
    EL::Float64
    EH::Float64
    layer::layerstructure
    
    function SommerfieldParams(layer::layerstructure,EL::Float64,EH::Float64 = 1e-4)
        new(EL,EH,layer)
    end
end

struct greenFunLayer <: greenFunct
    gTot::Function
    DgTot::Function
    
    function greenFunLayer(k0::Number, SParms::SommerfieldParams,xP::Tuple{Float64, Int64},yP::Tuple{Float64,Float64, Int64})
        mVec(r::Vector{Float64}) = r.*[1,-1]

        gfIn = greenFunLayerInd(k0,SParms,xP,yP)
        matA = SParms.layer.mat
        ε1f,εLf = matA[1].ε,matA[end].ε
        k1 = matA[1].k

        gHomo = greenHomo(k1)

        gtot  = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> gHomo.gTot(r,rl,k) + gfIn.gInd(abs(r[1]-rl[1]),r[2]+rl[2]) + ((εLf(k)-ε1f(k))/(εLf(k)+ε1f(k)))*(im/4)*hankelh1(0,k1(k)*norm(mVec(r)-rl))
        dgif = (r::Vector{Float64},rl::Vector{Float64},k::Number=k0) -> gHomo.DgTot(r,rl,k) + [gfIn.DxgInd(abs(r[1]-rl[1]),r[2]+rl[2])*sign(rl[1]-r[1]); gfIn.DygInd(abs(r[1]-rl[1]),r[2]+rl[2])] - ((εLf(k)-ε1f(k))/(εLf(k)+ε1f(k)))*(im/4)*hankelh1(1,k1(k)*norm(mVec(r)-rl))*((rl-mVec(r))/norm(rl-mVec(r)))*k1(k)
        
        new(gtot,dgif)
    end
end

struct greenFunLayerInd <: greenFunct
    gInd::Any
    DxgInd::Any
    DygInd::Any
    
    function greenFunLayerInd(k0::Number, SParms::SommerfieldParams,xP::Tuple{Float64, Int64},yP::Tuple{Float64,Float64, Int64})
        
        EL,EH,layer = SParms.EL,SParms.EH,SParms.layer
        
        xE,nx = xP
        yB,yE,ny = yP
        
        xdA = collect(LinRange(0.,xE,nx)) 
        ysA = collect(LinRange(yB,yE,ny))
        
        yT = 0:1:yE
        xT = 0:1:xE
        
        arrGr = [Evalgreenrefns1(k0,xdA[i],ysA[j],EL,EH,layer) for i in eachindex(xdA), j in eachindex(ysA)]

        method = BSpline(Cubic(Natural(OnGrid())))

        itpGind = scale(interpolate(arrGr, method),LinRange(0.,xE,nx),LinRange(0.,yE,ny))
        
        
        MatGFa = [itpGind(x,y) for x ∈ xT, y ∈ yT];
        
        mtdx = diff(MatGFa,dims=1)
        mtdy = diff(MatGFa,dims=2)

        itpGindDx = scale(interpolate(mtdx, method),xT[1:end-1],yT)
        itpGindDy = scale(interpolate(mtdy, method),xT,yT[1:end-1])
        
        gifn = (x,y)-> itpGind(x,y)
        gifnDx = (x,y)-> itpGindDx(x,y)
        gifnDy = (x,y)-> itpGindDy(x,y)
        
        new(gifn,gifnDx,gifnDy)
    end
end

struct GreenFunctions <: greenFunct
    GrFuncs::Array{greenFunct}
    
    function GreenFunctions(matScatter::Vector{MaterialParams},Opt::String;k0::Number=2π/600,
                SParms::SommerfieldParams = SommerfieldParams(layerstructure(matScatter,[0. ,],"up"),2k0),
                xP::Tuple{Float64, Int64} = (200. , 5),yP::Tuple{Float64,Float64, Int64} = (0.,40. , 5))
        
        mt1,mt2 = matScatter[1],matScatter[2]
        
        if Opt == "Homo"
            gf1 = greenHomo(mt1.k)
        else Opt == "Layer"
            gf1 = greenFunLayer(k0,SParms,xP,yP)
        end
        
        gf2 = greenHomo(mt2.k)
          
        new([gf1; gf2])
    end
end

kyi(β::Number,kᵢ::Number) = zsqrt(Complex(kᵢ^2 - β^2))

function greenRefnsInt(kx::Number,k0::Number,xdif::Float64,ysum::Float64,layer::layerstructure)
    εₗₑ = layer.mat[end].ε(k0)
    ε₁ = layer.mat[1].ε(k0)
    k₁ = layer.mat[1].k(k0)
    rtC = rtcoeffs(layer,k0,[kx,],"up")
    rp = rtC.r.TM[1]

    ky1 = kyi(kx,k₁)
    
    (cos(kx*xdif)*(rp-(εₗₑ - ε₁)/(εₗₑ + ε₁))*exp(im*ky1*ysum))/ky1
end

function kAng(α::Float64,Eₗ::Float64,Eₕ::Float64)
    (1+cos(α))*Eₗ/2 + im*Eₕ*sin(α)
end

function dkAng(α::Float64,Eₗ::Float64,Eₕ::Float64)
    (-sin(α))*Eₗ/2 + im*Eₕ*cos(α)
end

function Evalgreenrefns1(k0::Number,xdif::Float64,ysum::Float64,EL::Float64,EH::Float64,layer::layerstructure)
    fn1(x) = greenRefnsInt(kAng(x,EL,EH),k0,xdif,ysum,layer)*dkAng(x,EL,EH)
    fn2(x) = greenRefnsInt(x,k0,xdif,ysum,layer)
    
    I1,_ = quadgk(x -> fn1(x),-π,0)
    I2,_ = quadgk(x -> fn2(x),EL,Inf)
    
    (im/(2π))*(I1+I2)
end