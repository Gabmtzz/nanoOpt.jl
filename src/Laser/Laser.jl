function pol2cart(θ::Number,ρ::Number)
    x = ρ*cos(θ)
    y = ρ*sin(θ)
    
    return(x,y)
end

function pol2cart(θ::Number,ρ::Number,z::Number)
    x = ρ*cos(θ)
    y = ρ*sin(θ)
    
    return(x,y,z)
end

function cart2pol(x::Number,y::Number)
    θ = atan(y,x)
    ρ = sqrt(x^2+y^2)
    
    return(θ,ρ)
end

function cart2pol(x::Number,y::Number,z::Number)
    θ = atan(y,x)
    ρ = sqrt(x^2+y^2)
    
    return(θ,ρ,z)
end

function lglnodes(N)
    N1=N+1
    
    x = cos.(π*(0:N)/N)
    P = zeros(N1,N1)
    xold=2;
    ε = 1e-6
    
    while maximum(abs.(x.-xold))> ε
        xold=x
        P[:,1] .= 1; P[:,2] = x
    
        for k in 2:N
            P[:,k+1] = ((2*k-1)*x.*P[:,k]-(k-1)*P[:,k-1])/k
        end
        x = xold .-(x.*P[:,N1].-P[:,N]) ./ (N1*P[:,N1])
    
    end

    w = 2 ./ (N*N1*P[:,N1].^2)
    
    x,w,P
end

include("focusHG00.jl")