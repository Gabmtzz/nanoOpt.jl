struct lensfocus <: Field
    mat::MaterialParams
    k0::Number
    NA::Number
    x::Vector{Float64}
    y::Vector{Float64}
    rad::Number
    ϕ::Vector{Float64}
    θ::Vector{Float64}
    ρ::Vector{Float64}
    
    function lensfocus(mat::MaterialParams,k0::Number,NA::Number)
        nrad,nphi = 51,51
        
        n = mat.n(k0)
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
