function waveselect(kvec::Vector{Matrix{ComplexF64}},NA::Number)
    kvec = [reshape(kvec[1],(:,1)) reshape(kvec[2],(:,1)) reshape(kvec[3],(:,1))]
    k = sqrt(real(sum(kvec[1,:].^2)))
    kz = kvec[:,3]
    ind = findall(imag(kz) .== 0)
    
    dir = kvec[ind,:] /k
    
    return(ind,dir)
end

function waveselect(kvec::Vector{Matrix{ComplexF64}},n::Number,NA::Number)
    kvec = [reshape(kvec[1],(:,1)) reshape(kvec[2],(:,1)) reshape(kvec[3],(:,1))]
    k = sqrt(real(sum(kvec[1,:].^2)))
    kz = kvec[:,3]
    ind = findall(imag(kz) .== 0)

    sint = NA/n
    cost = kz ./ k
    ind  = findall((real.(1 .-cost.^2) .< sint ^ 2).&(imag(kz) .== 0))

    dir = kvec[ind,:] /k
    
    return(ind,dir)
end

function cart2sph(x::Number,y::Number,z::Number)
    ρ = sqrt(x^2+y^2+z^2)
    θ = atan(z,sqrt(x^2+y^2))
    ϕ = atan(y,x)
    
    return (ρ,θ,ϕ)
end

function cart2unit(vec::Matrix{Number})
    siz = size(vec)
    vec = reshape(vec,(:,3))
    vec = real.(vec)
    sphVec = zeros(size(vec,1),2)
    for i in axes(sphVec,1)
        _,sphVec[i,1],sphVec[i,2] = cart2sph(vec[i,1],vec[i,2],vec[i,3])
    end

    θ,ϕ = sphVec[:,1],sphVec[:,2];
    θ = π/2 .- θ
    
    sinθ, cosθ = sin.(θ),cos.(θ)
    sinϕ, cosϕ = sin.(ϕ),cos.(ϕ)
    uϕ = [-sinϕ cosϕ 0*ϕ]
    ur = [cosϕ sinϕ 0*ϕ]
    uθ = [cosϕ.*cosθ sinϕ.*cosθ -sinθ]
    
    return (uϕ,uθ,ur,ϕ,θ)
end

function unitproject(x,u1)
    sum(x.*u1, dims=ndims(x))
end

function unitproject(x,u1,u2)
    y = sum(x.*u1, dims=ndims(x))
    u2.*y
end