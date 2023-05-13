struct Angular <: Field
    mat::MaterialParams
    iefield::Array{ComplexF64, 3}
    k0::Number
    dir::Number
    x::Vector{Number}
    y::Vector{Number}
    kvec::Vector{Matrix{ComplexF64}}
    kx::Vector{Number}
    ky::Vector{Number}
    
    function Angular(x,y,efield, mat, k0,dir,npad)
        k = mat.k(k0)
        Mx,My = npad*length(x),npad*length(y)
        iefield,kx,ky = fourier2(x,y,efield,Mx,My)
        kxg,kyg= ndgrid(kx,ky)
        kz = dir*sqrt.(Complex.(k^2 .- kxg.^2 .- kyg.^2))
        kvec = [Complex.(kxg),Complex.(kyg),kz]
        
        new(mat,iefield,k0,dir,collect(x),collect(y),kvec,kx,ky)
    end
end

struct far <: Field
    k0::Number
    kx::Vector{Number}
    ky::Vector{Number}
    x::Vector{Number}
    y::Vector{Number}
    far::Array{ComplexF64, 3}
    
    function far(angular::Angular)
        kz= angular.kvec[3]
        ind = findall(imag(kz[:]) .!= 0.0)

        farf = copy(reshape(angular.iefield,(:,3)))
        farf[ind,:] .= 0.0*im

        farf = -im*2π*kz.*reshape(farf,(size(angular.iefield)))
        
        new(angular.k0,angular.kx,angular.ky,angular.x,angular.y,farf)
    end
end

function efieldC(angular::Angular,z::Number)
    ie = reshape(angular.iefield,(:,3))
    kz = reshape(angular.kvec[3],(:,1))
    ie = ie.*exp.(im*kz*z)
    
    ifourier2(angular.x,angular.y,reshape(ie,size(angular.iefield)))
end