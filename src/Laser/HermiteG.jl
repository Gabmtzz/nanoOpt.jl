function hermiteGauss(mat::MaterialParams,k0::Number,w0::Number,m::Int64,n::Int64,x::Vector{Float64},y::Vector{Float64},z::Number,pol::Vector{Float64})
    z = z * ones(length(x))

    k =  mat.k(k0)
    e = fun(k,w0,m,x,z).* fun(k,w0,n,y,z)
    im*(e*(pol*im)')
end

function fun(k::Number,w0::Number,m::Int64,x::Vector{Float64},z::Vector{Float64})
    zr = 0.5*w0^2*k
    wz = w0*sqrt.(1.0 .+ (z/zr).^2)

    q0i = -2im/(w0^2*k)
    qzi = z ./ (z .^2 .+ zr^2) - 2im ./(wz.^2*k)

    q0,qz = 1.0/q0i, 1.0 ./qzi

    H = hermite.(m,sqrt(2.0)*x./wz)
    fac = sqrt(sqrt(2/π) / (2^m*factorial(m)*w0))
    
    fac*sqrt.(q0 ./ qz).*(-conj(qz)./q0) .^ (m/2).*H.*exp.(-0.5im*k*x.^2 ./ qz)
end

function hermite(m::Int64,x::Number)
    H = 0.0
    if m == 0
        H = 0*x+1.
    elseif m == 1
        H = 2.0 *x
    elseif m == 2
        H = 4.0 *x^2 - 2.
    elseif m == 3
        H = 8.0*x^3 - 12.0 * x
    elseif m ==4
        H = 16.0*x^4 - 48.0 * x^2 + 12.
    end
    return H
end