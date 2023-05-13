struct lensimage2 <: Field
    ang1::Angular
    mat1::MaterialParams
    mat2::MaterialParams
    k0::Number
    NA::Number
    
    function lensimage2(angular::Angular,mat2::MaterialParams,NA::Number)
       new(angular,angular.mat,mat2,angular.k0,NA) 
    end
end

function efieldC(lens::lensimage2,farf::Array{ComplexF64, 3})
    n1,n2 = lens.mat1.n(lens.k0),lens.mat2.n(lens.k0)
    k1,k2 = lens.mat1.k(lens.k0),lens.mat2.k(lens.k0)

    kx,ky = lens.ang1.kx,lens.ang1.ky;

    kx,ky = ndgrid(kx,ky)

    kz = sqrt.(Complex.(k1^2 .- kx.^2 .- ky.^2 ))

    ind,_ = waveselect([Complex.(kx),Complex.(ky),kz],n1,lens.NA)
    
    siz = size(farf)
    farf = reshape(farf,(:,3))
    farf,kx,ky,kz = farf[ind,:],kx[ind],ky[ind],kz[ind]
    
    uϕ,uθ,ur,_,θ = cart2unit(hcat(kx,ky,kz))
    farf = unitproject(farf,uθ,ur) + unitproject(farf,uϕ,uϕ)
    farf = sqrt(n1/n2)*@. farf*sqrt(1/cos(θ))
    
    et = reshape(zeros(siz)*im,(:,3))
    et[ind,:] = farf
    
    M = n1/n2
    (im/(2π*k2)) *(4π/M)^2 .*ifourier2(lens.ang1.x,lens.ang1.y,reshape(et,siz))
end