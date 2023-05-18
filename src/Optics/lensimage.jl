struct lensimage <: Field
    mat1::MaterialParams
    mat2::MaterialParams
    k0::Number
    NA::Number
    dir::Matrix{Float64}
    ϕ::Vector{Float64}
    θ::Vector{Float64}
    
    function lensimage(mat1::MaterialParams,mat2::MaterialParams,k0::Number,NA::Number)
        nϕ,nθ = 51,51
        
        n1 = mat1.n(k0)
        θ = asin(NA/n1)
        ϕa = collect(range(0,2π,length=nϕ))
        θa = collect(range(0,θ,length=nθ))
        
        ϕm,θm = ndgrid(ϕa,θa)

        dir = hcat(cos.(ϕm[:]).*sin.(θm[:]),sin.(ϕm[:]).*sin.(θm[:]),cos.(θm[:]))
        
        new(mat1,mat2,k0,NA,dir,ϕa,θa)
    end
end

function efieldC(lens::lensimage,farf::Matrix{ComplexF64},x::Vector{Float64},y::Vector{Float64})
    n1 = lens.mat1.n(lens.k0); k1 = lens.k0*n1
    n2 = lens.mat2.n(lens.k0); k2 = lens.k0*n2
    
    uϕ1,uθ1,ur1,ϕ1,θ1 = cart2unit(lens.dir)
    
    farf=unitproject(farf,uθ1,ur1) + unitproject(farf,uϕ1,uϕ1)
    ϕ1,θ1 = reshape(ϕ1,(1,:)),reshape(θ1,(1,:))
    
    xx,yy = ndgrid(x,y)
    elm = cart2pol.(xx[:],yy[:])
    ϕ2,r2 = first.(elm),last.(elm)
    
    trans=exp.(im*k1*r2*sin.(θ1).*cos.( ϕ1 .- ϕ2))
    trans = im*k2*sqrt(n2/n1)/(2π).*trans.*sin.(θ1).*sqrt.(cos.(θ1))
    
    dϕ1 = diff(lens.ϕ[[1,2]])[1]
    dθ1 = diff(lens.θ[[1,2]])[1]
    
    e = trans*farf*dϕ1*dθ1
    reshape(e,(length(x),length(y),3))
end