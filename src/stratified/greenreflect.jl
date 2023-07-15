mutable struct greenreflect <: Field
    layer::layerstructure
    semi::Number
    ratio::Int
    k0::Number
    k1::Number
    kn::Number
    ϕ::Vector{Float64}
    ρ::Vector{Float64}
    z::Vector{Float64}
    zp::Vector{Float64}

    function greenreflect(layer::layerstructure,semi::Number=0.1,ratio::Int=2)
        new(layer,semi,ratio,0.0,0.0,0.0,[0.0,],[0.0,],[0.0,],[0.0,])
    end
end

function GreenEval(green::greenreflect,pos1::Matrix{Float64},pos2::Matrix{Float64},k0::Number)
    green.k0 = k0
    z1,zn = green.layer.z[1],green.layer.z[end]

    k = [green.layer.mat[i].k(k0) for i ∈ eachindex(green.layer.mat)]
    green.k1,green.kn,kmax = k[1],k[end],maximum(real.(k))

    cVal = cart2pol.(pos1[:,1] - pos2[:,1],pos1[:,2] - pos2[:,2])
    ϕ,ρ = first.(cVal),last.(cVal)
    z,zp = pos1[:,3],pos2[:,3]

    ind1 = z .> zn .&& zp .> zn .&& (z + zp .- 2zn) .≥ green.ratio*ρ; Index1 = findall(ind1)
    ind2 = z .> zn .&& zp .> zn .&& (z + zp .- 2zn) .< green.ratio*ρ; Index2 = findall(ind2)
    ind3 = z .< z1 .&& zp .> zn; Index3 = findall(ind3)

    g = zeros(size(pos1,1),9)*im

    if any(ind1)
        green.ρ,green.ϕ,green.z,green.zp = ρ[Index1],ϕ[Index1],z[Index1],zp[Index1]
    
        fInt = (kr,kz,mode) -> ifun1(green,[kr,],mode)
        g[Index1,:] = isommerfield(green.kn,fInt,green.semi,kmax+k0,"bessel")
    end
    
    if any(ind2)
        green.ρ,green.ϕ,green.z,green.zp = ρ[Index2],ϕ[Index2],z[Index2],zp[Index2]
    
        fInt = (kr,kz,mode) -> ifun1(green,[kr,],mode)
        g[Index2,:] = isommerfield(green.kn,fInt,green.semi,kmax+k0,"hankel")
    end
    
    if any(ind3)
        green.ρ,green.ϕ,green.z,green.zp = ρ[Index3],ϕ[Index3],z[Index3],zp[Index3]
    
        fInt = (kr,kz,mode) -> ifun2(green,[kr,],mode)
        g[Index3,:] = isommerfield(green.kn,fInt,green.semi,kmax+k0,"bessel")
    end
    
    g
end

function ifun1(green::greenreflect,kr::Vector{Float64},mode::String)
    knz = zsqrt(green.kn^2 - kr[1]^2)
    kz1 = knz
    kz2 = knz*Int(green.z > green.zp) - knz*Int(green.z ≤ green.zp)

    refld = rtcoeffs(green.layer,green.k0,kr,"down")
    r = refld.r
    te,tm = average(green,kr[1],green.kn,green.kn,kz1,kz2,mode)
    fac = exp.(im*knz*(green.z + green.zp .- 2*green.layer.z[end]))/(4π)
    
    (r.TM.*tm + r.TE.*te ).*fac
end

function ifun2(green::greenreflect,kr::Vector{Float64},mode::String)
    knz = zsqrt(green.kn^2 - kr[1]^2)
    k1z = zsqrt(green.k1^2 - kr[1]^2)

    refld = rtcoeffs(green.layer,green.k0,kr,"down")
    t = refld.t

    te,tm = average(green,kr[1],green.k1,green.kn,-k1z,-knz,mode)
    
    fac = exp.(im*k1z.*(green.layer.z[1] .- green.z) .+ im*knz.*(green.zp.-green.layer.z[end])) / (4π)

    (t.TM.*tm + t.TE.*te ).*fac
end

function ifun1(green::greenreflect,kr::Vector{ComplexF64},mode::String)
    knz = zsqrt(green.kn^2 - kr[1]^2)
    kz1 = knz
    kz2 = knz*Int(green.z > green.zp) - knz*Int(green.z ≤ green.zp)

    refld = rtcoeffs(green.layer,green.k0,kr,"up")
    r = refld.r
    te,tm = average(green,kr[1],green.kn,green.kn,kz1,kz2,mode)
    fac = exp.(im*knz*(green.z + green.zp .- 2*green.layer.z[end]))/(4π)
    
    (r.TM.*tm + r.TE.*te ).*fac
end

function ifun2(green::greenreflect,kr::Vector{ComplexF64},mode::String)
    knz = zsqrt(green.kn^2 - kr[1]^2)
    k1z = zsqrt(green.k1^2 - kr[1]^2)

    refld = rtcoeffs(green.layer,green.k0,kr,"up")
    t = refld.t

    te,tm = average(green,kr[1],green.k1,green.kn,-k1z,-knz,mode)
    
    fac = exp.(im*k1z.*(green.layer.z[1] .- green.z) .+ im*knz.*(green.zp.-green.layer.z[end])) / (4π)

    (t.TM.*tm + t.TE.*te ).*fac
end

function average(green::greenreflect,kr::Number,k1::Number,k2::Number,kz1::Number,kz2::Number,mode::String)
    ϕ,ρ = green.ϕ,green.ρ;

    if mode == "bessel"
        g0 = besselj.(0,kr*ρ)
        g1 = besselj.(1,kr*ρ)
        g2 = besselj.(2,kr*ρ)
    else mode == "hankel"
        g0 = hankelh1.(0,kr*ρ)
        g1 = hankelh1.(1,kr*ρ)
        g2 = hankelh1.(2,kr*ρ)
    end

    te,tm = zeros(length(green.z),9)*im,zeros(length(green.z),9)*im

    te[:,1] = g0 .+ g2.*cos.(2ϕ)
    te[:,2] = g2.*sin.(2ϕ)
    te[:,4] = g2.*sin.(2ϕ)
    te[:,5] = g0 .- g2.*cos.(2ϕ)

    tm[:,1] = @. kz1*kz2*(g0 - g2*cos(2ϕ))
    tm[:,2] = @. -kz1*kz2*g2*sin(2ϕ)
    tm[:,3] = @. -2*im*kr*kz2*g1*cos(ϕ)
    tm[:,4] = @. -kz1*kz2*g2*sin(2ϕ)
    tm[:,5] = @. kz1*kz2*(g0 + g2*cos(2ϕ))
    tm[:,6] = @. -2*im*kr*kz2*g1*sin(ϕ)
    tm[:,7] = @. -2*im*kr*kz1*g1*cos(ϕ)
    tm[:,8] = @. -2*im*kr*kz1*g1*sin(ϕ)
    tm[:,9] = @. 2*kr*kr*g0

    0.5*te,0.5*tm / (k1*k2)
end

# ===============================================================================================

struct pramsSF
    kmax::Number
    semi::Number
    mode::String
end

function ifun1(x,k::Number,fInt::Function,pS::pramsSF)
    kr = @. pS.kmax*(1-cos(x) - im*pS.semi*sin(x))
    kz = zsqrt.(k^2 .-kr.^2)

    pS.kmax.*(sin.(x) - im*pS.semi.*cos.(x)) .*kr./kz .*reshape(fInt.(kr,kz,"bessel"),(:,1))
end

function ifun2(x,k::Number,fInt::Function,pS::pramsSF)
    kr = @. pS.kmax/x
    kz = zsqrt.(k^2 .-kr.^2)

    -2im*pS.kmax ./ x.^2 .*kr./kz .*reshape(fInt.(kr,kz,"bessel"),(:,1))
end

function ifun3(x,k::Number,fInt::Function,pS::pramsSF)
    kr1 = @. 2*pS.kmax*(1 - im + im/x)
    kr2 = @. 2*pS.kmax*(-1 - im + im/x)

    kz1 = zsqrt.(k^2 .-kr1.^2)
    kz2 = zsqrt.(k^2 .-kr2.^2)

    pS.kmax ./ x.^2 .* (kr1./kz1 .*reshape(fInt.(kr1,kz1,"hankel"),(:,1)) - kr2./kz2 .*reshape(fInt.(kr2,kz2,"hankel"),(:,1)))
end

function isommerfield(k::Number,fInt::Function,semi::Number=0.1,kmax::Number=k,mode::String="bessel")
    pS = pramsSF(kmax,semi,mode)

    siz = size(fInt(1e-10,k,"bessel"))
    n = prod(siz)
    y = zeros(n)*im

    function integrand(du,u,p,t)
        u = ifun1(t,k,fInt,pS)
        
        du[1:n] = u
    end
    
    function integrandB(du,u,p,t)
        u = ifun2(t,k,fInt,pS)
        
        du[1:n] = u
    end

    function integrandI(du,u,p,t)
        u = ifun3(t,k,fInt,pS)
        
        du[1:n] = u
    end

    prob = ODEProblem(integrand, y, (0,π))
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
    y = sol.u[end]

    if mode == "bessel"
        prob = ODEProblem(integrandB, y, (1,1e-10))
        sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
        y = sol.u[end]
    elseif mode == "hankel"
        prob = ODEProblem(integrandI, y, (1,1e-10))
        sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
        y = sol.u[end]
    end

    reshape(y,(siz))
end