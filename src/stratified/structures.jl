struct layerstructure <: Structure
    mat::Vector{MaterialParams}
    z::Vector{Float64}
    dir::String
    
    function layerstructure(mat::Vector{MaterialParams},z::Vector{Float64},dir::String)
        new(mat,z,dir)
    end
    
end

struct slabstructure <: Field
    mat::Vector{MaterialParams}
    d::Number

    function slabstructure(mat::Vector{MaterialParams},d::Number)
        new(mat,d)
    end
end

struct fieldTransVec <: Field
    TE::Vector{ComplexF64}
    TM::Vector{ComplexF64}
end

struct fresnel <: Field
    r::fieldTransVec
    t::fieldTransVec
    
    function fresnel(layer::layerstructure,k0::Number,kpar,ind::Int64)
        ε₁,μ₁ = layer.mat[ind].ε(k0),layer.mat[ind].μ(k0)
        ε₂,μ₂ = layer.mat[ind+1].ε(k0),layer.mat[ind+1].μ(k0)

        kz1 = zsqrt.(Complex.(μ₁*ε₁*k0^2 .- kpar.^2))
        kz2 = zsqrt.(Complex.(μ₂*ε₂*k0^2 .- kpar.^2))

        Z₁ = zsqrt(Complex(μ₁/ε₁))
        Z₂ = zsqrt(Complex(μ₂/ε₂))
        
        rTE = (μ₂.*kz1 - μ₁.*kz2) ./ (μ₂.*kz1 + μ₁.*kz2) 
        tTE = (2*μ₂.*kz1) ./ (μ₂.*kz1 + μ₁.*kz2) 

        rTM = (ε₂.*kz1 - ε₁.*kz2) ./ (ε₂.*kz1 + ε₁.*kz2)
        tTM = (2*ε₂.*kz1) ./ (ε₂.*kz1 + ε₁.*kz2)*Z₂/Z₁
        
        r = fieldTransVec(rTE,rTM)
        t = fieldTransVec(tTE,tTM)
        
        new(r,t)
    end

    function fresnel(slab::slabstructure,k0::Number,kpar,ind::Int64)
        ε₁,μ₁ = slab.mat[ind].ε(k0),slab.mat[ind].μ(k0)
        ε₂,μ₂ = slab.mat[ind+1].ε(k0),slab.mat[ind+1].μ(k0)

        kz1 = zsqrt.(Complex.(μ₁*ε₁*k0^2 .- kpar.^2))
        kz2 = zsqrt.(Complex.(μ₂*ε₂*k0^2 .- kpar.^2))

        Z₁ = zsqrt(Complex(μ₁/ε₁))
        Z₂ = zsqrt(Complex(μ₂/ε₂))
        
        rTE = (μ₂.*kz1 - μ₁.*kz2) ./ (μ₂.*kz1 + μ₁.*kz2) 
        tTE = (2*μ₂.*kz1) ./ (μ₂.*kz1 + μ₁.*kz2) 

        rTM = (ε₂.*kz1 - ε₁.*kz2) ./ (ε₂.*kz1 + ε₁.*kz2)
        tTM = (2*ε₂.*kz1) ./ (ε₂.*kz1 + ε₁.*kz2)*Z₂/Z₁
        
        r = fieldTransVec(rTE,rTM)
        t = fieldTransVec(tTE,tTM)
        
        new(r,t)
    end
end

struct transfer <: Field
    TE::Array{ComplexF64, 3}
    TM::Array{ComplexF64, 3}
    
    function transfer(layer::layerstructure,k0::Number,kpar,ind::Int64)
        fresnelC = fresnel(layer,k0,kpar,ind);
        r,t = fresnelC.r,fresnelC.t
        
        mTE = reshape([1 ./ t.TE r.TE ./ t.TE; r.TE ./ t.TE 1 ./ t.TE],(:,2,2))
        mTM = reshape([1 ./ t.TM r.TM ./ t.TM; r.TM ./ t.TM 1 ./ t.TM],(:,2,2))
        
        new(mTE,mTM)
    end
end

function propagate(layer::layerstructure,k0::Number,kpar,ind::Int64)
    ε,μ = layer.mat[ind].ε(k0), layer.mat[ind].μ(k0) 

    kz = zsqrt.(Complex.(μ*ε*k0^2 .- kpar.^2))
    d = layer.z[ind] - layer.z[ind-1]

    reshape([exp.(-im*kz*d) 0 .*kz; 0 .*kz exp.(im*kz*d)],(:,2,2))
end

function propagate(slab::slabstructure,k0::Number,kpar)
    ε,μ = slab.mat[2].ε(k0), slab.mat[2].μ(k0) 

    kz = zsqrt.(Complex.(μ*ε*k0^2 .- kpar.^2))
    d = slab.d

    reshape([exp.(-im*kz*d) 0 .*kz; 0 .*kz exp.(im*kz*d)],(:,2,2))
end

function fun(a::Array{ComplexF64, 3},b::Array{ComplexF64, 3})
    c = zeros(Complex,size(a))

    c[:,1,1] = a[:,1,1].* b[:,1,1]+a[:,1,2].* b[:,2,1]
    c[:,1,2] = a[:,1,1].* b[:,1,2]+a[:,1,2].* b[:,2,2]
    c[:,2,1] = a[:,2,1].* b[:,1,1]+a[:,2,2].* b[:,2,1]
    c[:,2,2] = a[:,2,1].* b[:,1,2]+a[:,2,2].* b[:,2,2]
    
    c
end

function fun(a::Array{Complex, 3},b::Array{ComplexF64, 3})
    c = zeros(Complex,size(a))

    c[:,1,1] = a[:,1,1].* b[:,1,1]+a[:,1,2].* b[:,2,1]
    c[:,1,2] = a[:,1,1].* b[:,1,2]+a[:,1,2].* b[:,2,2]
    c[:,2,1] = a[:,2,1].* b[:,1,1]+a[:,2,2].* b[:,2,1]
    c[:,2,2] = a[:,2,1].* b[:,1,2]+a[:,2,2].* b[:,2,2]
    
    c
end

struct transfertot <: Field
    TE::Array{ComplexF64, 3}
    TM::Array{ComplexF64, 3}
    
    function transfertot(layer::layerstructure,k0::Number,kpar)
         n = length(layer.z)
         mTE,mTM = zeros(ComplexF64,length(kpar),2,2),zeros(ComplexF64,length(kpar),2,2)

        
        for ind ∈ 1:n
            m = transfer(layer,k0,kpar,ind)
    
            if ind == 1
                mTE,mTM = m.TE,m.TM
            else
                mTE = fun(mTE,m.TE)
                mTM = fun(mTM,m.TM)
            end
        
            if ind ≠ n
                p = propagate(layer,k0,kpar,ind+1)
                mTE = fun(mTE,p)
                mTM = fun(mTM,p)
            end
        end
            
        new(mTE,mTM)
    end

    function transfertot(slab::slabstructure,k0::Number,kpar)

        mTE,mTM = zeros(ComplexF64,length(kpar),2,2),zeros(ComplexF64,length(kpar),2,2)

        Fr12 = fresnel(slab,k0,kpar,1)
        r12,t12 = Fr12.r,Fr12.t
        Fr23 = fresnel(slab,k0,kpar,2)
        r23,t23 = Fr23.r,Fr23.t

        p = propagate(slab,k0,kpar)
        p1,p2 = p[:,1,1],p[:,2,2]

        mTE[:,1,1] = (p1+p2.*r12.TE.*r23.TE) ./ (t12.TE .* t23.TE)
        mTE[:,1,2] = ( p1 .* r23.TE + p2 .* r12.TE ) ./ (t12.TE .* t23.TE)
        mTE[:,2,1] = ( p1 .* r12.TE + p2 .* r23.TE ) ./ (t12.TE .* t23.TE)
        mTE[:,2,2] = ( p1 .* r12.TE .* r23.TE + p2 ) ./ (t12.TE .* t23.TE)

        mTM[:,1,1] = ( p1 + p2 .* r12.TM .* r23.TM ) ./ ( t12.TM .* t23.TM )
        mTM[:,1,2] = ( p1 .* r23.TM + p2 .* r12.TM ) ./ ( t12.TM .* t23.TM )
        mTM[:,2,1] = ( p1 .* r12.TM + p2 .* r23.TM ) ./ ( t12.TM .* t23.TM )
        mTM[:,2,2] = ( p1 .* r12.TM .* r23.TM + p2 ) ./ ( t12.TM .* t23.TM )

        new(mTE,mTM)
    end
end
    
struct rtcoeffs <: Field
    r::fieldTransVec
    t::fieldTransVec
    
    function rtcoeffs(layer::layerstructure,k0::Number,kpar,dir::String)
        mtot = transfertot(layer,k0,kpar)
        
        if dir == "up"
            rTE = mtot.TE[:,2,1] ./mtot.TE[:,1,1]
            tTE = 1.0 ./ mtot.TE[:,1,1]

            rTM = mtot.TM[:,2,1] ./mtot.TM[:,1,1]
            tTM = 1.0 ./ mtot.TM[:,1,1]
    
        elseif dir == "down"
            rTE = -mtot.TE[:,1,2] ./mtot.TE[:,1,1]
            tTE =  mtot.TE[:,2,2] + mtot.TE[:,2,1] .*rTE

            rTM = -mtot.TM[:,1,2] ./mtot.TM[:,1,1]
            tTM = mtot.TM[:,2,2] + mtot.TM[:,2,1] .*rTM
        end
    
        r,t = fieldTransVec(rTE,rTM),fieldTransVec(tTE,tTM)
            
        new(r,t)
    end

    function rtcoeffs(slab::slabstructure,k0::Number,kpar,dir::String)
        mtot = transfertot(slab,k0,kpar)
        
        if dir == "up"
            rTE = mtot.TE[:,2,1] ./mtot.TE[:,1,1]
            tTE = 1.0 ./ mtot.TE[:,1,1]

            rTM = mtot.TM[:,2,1] ./mtot.TM[:,1,1]
            tTM = 1.0 ./ mtot.TM[:,1,1]
    
        elseif dir == "down"
            rTE = -mtot.TE[:,1,2] ./mtot.TE[:,1,1]
            tTE =  mtot.TE[:,2,2] + mtot.TE[:,2,1] .*rTE

            rTM = -mtot.TM[:,1,2] ./mtot.TM[:,1,1]
            tTM = mtot.TM[:,2,2] + mtot.TM[:,2,1] .*rTM
        end
    
        r,t = fieldTransVec(rTE,rTM),fieldTransVec(tTE,tTM)
            
        new(r,t)
    end
end

function zsqrt(x::Number)
    y = sqrt(Complex(x))
    y*sign(imag(y+im*eps()))
end