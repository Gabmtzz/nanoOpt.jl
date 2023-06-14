struct dipole <: Field
    mat::MaterialParams
    k₀::Number
    pos::Vector{Number}
    dip::Vector{Number}
    
    function dipole(material::MaterialParams,k₀::Number,pos,dip)
        new(material,k₀,pos,dip)
    end
end

function Dipoles(mat::MaterialParams,k0::Number,dip::Matrix{Float64},pos::Matrix{Float64})
    [dipole(mat,k0,pos[i,:],dip[i,:]) for i in axes(dip,1)]
end

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

function eplane(Dipols::Vector{dipole},x,y, npad::Int64, z::Number)
    Mx,My = npad*length(x),npad*length(y)
    ef = zeros(Mx,My,3).*im

    for dip ∈ Dipols
        e = eplane(dip,x,y,npad,z)
        ef = ef + e
    end
    
    ef
end

function farfield(dip::dipole,dir::Matrix{Float64})
    k = dip.mat.k(dip.k₀)
    position,dipol=dip.pos,dip.dip

    dir = dir ./ sum(dir.*dir,dims=2)
    
    f = dip.k₀^2*exp.(-im*k*dir*position)/(4π)
    
    ex = f .* (dipol[1].-dir*dipol.*dir[:,1])
    ey = f .* (dipol[2].-dir*dipol.*dir[:,2])
    ez = f .* (dipol[3].-dir*dipol.*dir[:,3])
    hcat(ex,ey,ez)
end

function farfield(Dipols::Vector{dipole},dir::Matrix{Float64})
    farf = zeros(size(lens.dir))
    for dip ∈ Dipols
        e = farfield(dip,lens.dir)
        farf = farf.+e
    end
    
    farf
end