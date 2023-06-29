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