struct miesolver <: Field
    mat1::MaterialParams
    mat2::MaterialParams
    diameter::Number
    lmax::Integer
    
    function miesolver(mat1::MaterialParams,mat2::MaterialParams,diameter::Number,lmax::Integer = 40)
        new(mat1,mat2,diameter,lmax)
    end
end

function riccatibessel(z::Number,ltab::Vector{Int64})
    l = 1:maximum(ltab)

    j₀ = sin(z)/z; j =  √(π/(2*z))*besselj.(l[:].+0.5,z)
    y₀ = -cos(z)/z; y =   √(π/(2*z))*bessely.(l[:].+0.5,z)

    h₀ = j₀ + im*y₀; h = j + im*y
    zjp = z*[j₀; j[1:length(l)-1]] - (l[:].+1).*j
    zhp = z*[h₀; h[1:length(l)-1]] - (l[:].+1).*h

    j = j[ltab]; zjp = zjp[ltab]
    h = h[ltab]; zhp = zhp[ltab]
    
    j,h,zjp,zhp
end

function miecoefficients(mie::miesolver,k0::Number)
    k₁,Z₁ = mie.mat1.k(k0),mie.mat1.Z(k0)
    k₂,Z₂ = mie.mat2.k(k0),mie.mat2.Z(k0)
    
    z₁,z₂ = 0.5k₁*mie.diameter,0.5k₂*mie.diameter
    
    j₁,_,zjp₁,_ = riccatibessel(z₁,collect(1:mie.lmax))
    j₂,h₂,zjp₂,zhp₂ = riccatibessel(z₂,collect(1:mie.lmax))
    
    a = (Z₂*z₁*j₁.*zjp₂ - Z₁*z₂*j₂.*zjp₁) ./ (Z₂*z₁*j₁.*zhp₂ - Z₁*z₂*h₂.*zjp₁)
    b = (Z₂*z₂*j₂.*zjp₁ - Z₁*z₁*j₁.*zjp₂) ./ (Z₂*z₂*h₂.*zjp₁ - Z₁*z₁*j₁.*zhp₂)
    
    a,b
end