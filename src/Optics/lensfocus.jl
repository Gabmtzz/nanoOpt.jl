struct lensfocus <: Field
    mat::MaterialParams
    k0::Number
    NA::Number
    x::Vector{Number}
    y::Vector{Number}
    rad::Number
    ϕ::Vector{Number}
    θ::Vector{Number}
    ρ::Vector{Number}
    
    function lensfocus(mat::MaterialParams,k0::Number,NA::Number)
        nrad,nphi = 51,51
        
        n = mat.n(k₀)
        θ = asin(NA/n)
        rad = sin(θ)
        ϕa = collect(range(0,2π,length=nphi))
        θa = collect(range(0,θ,length=nrad))
        ρa = sin.(θa)
        
        ϕm,ρm = ndgrid(ϕa,ρa)
        ϕm,ρm = ϕm[:],ρm[:]
        xa,ya = zeros(length(ϕm)),zeros(length(ϕm))

        for i ∈ eachindex(xa)
            xa[i],ya[i] = pol2cart(ϕm[i],ρm[i])
        end
        
        new(mat,k0,NA,xa,ya,rad,ϕa,θa,ρa)
    end
end
