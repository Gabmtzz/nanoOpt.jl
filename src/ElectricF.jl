using LazyGrids

function eplane(p::dipole,x,y, npad::Int64, z::Number)
    hx = x[2]-x[1]; Mx = npad*length(x)
    hy = y[2]-y[1]; My = npad*length(y);

    kx,ky = wavefourier(Mx) /hx, wavefourier(My) /hy
    kx,ky = ndgrid(kx,ky)
    kx,ky = kx[:],ky[:]

    k = p.mat.k(p.k₀)

    kz = sqrt.(Complex.(k^2 .- kx.^2 .- ky.^2))

    pos,dip = p.pos,p.dip
    x0,y0,z0= pos[1],pos[2],pos[3].-z
    px,py,pz = dip[1],dip[2],dip[3]
    
    g =  @. im* exp(-im*(kx*x0+ky*y0+kz*z0))/(8π^2*kz)
    gp = im*(kx*px+ky*py+kz*pz).*g /k^2
    gx = sum(g.*px, dims=2)
    gy = sum(g.*py, dims=2)
    gz = sum(g.*pz, dims=2)
    gp = sum(gp,dims=2)
    
    ie=[gx+im*kx.*gp;gy+im*ky.*gp;gz+im*kz.*gp]
    
    p.k₀^2 * ifourier2(x,y,reshape(ie,(Mx,My,3)))
end

function efieldC(angular,z)
    ie = reshape(angular.iefield,(:,3))
    kz = reshape(angular.kvec[3],(:,1))
    ie = ie.*exp.(im*kz*z)
    
    ifourier2(angular.x,angular.y,reshape(ie,size(angular.iefield)))
end