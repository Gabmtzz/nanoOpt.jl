#kyi(β::Number,kᵢ::Number) = zsqrt(Complex(kᵢ^2 - β^2))

function funSolv(β::Number,k₀::Number,layer::layerstructure)
    
    Mt = transfertot(layer,k₀,[β,]);
    
    Mt.TM[1,1,1]
end


function NRSolver(βᵢ::Number,k₀::Number,layer::layerstructure,thr::Float64 = 1e-6)
    δ = 1e-6*(1+im)

    flag = true
    βₙ = 0 + im*0
    itr = 1
    
    while flag
        βₙ = βᵢ - δ*funSolv(βᵢ,k₀,layer)/(funSolv(βᵢ+δ,k₀,layer)-funSolv(βᵢ,k₀,layer))
        
        if abs(βₙ - βᵢ) ≤ thr || itr ≥ 300
            flag = false
        else
            βᵢ = βₙ
            itr += 1
        end
    end
    
    βₙ
end

function IsDifπ(βreal::Vector{Float64},βImag::Vector{Float64},k₀::Number,layer::layerstructure)
    flag = false
    βℜ,βℑ = ndgrid(βreal,βImag)

    βArr = βℜ+im*βℑ

    βAm = βArr[:]

    FM = reshape(angle.([funSolv(β,k₀,layer) for β in βAm]),(size(βArr)))

    difF = [abs(FM[1,1]-FM[1,2]), abs(FM[1,2]-FM[2,2]),abs(FM[2,2]-FM[2,1]),abs(FM[2,1]-FM[1,1])]

    arrπ = findall(difF .> π)


    if ~isempty(arrπ) 
        flag = length(arrπ) == 1 
    end

    flag
end

function getβc(βRec::Vector{Float64},βImc::Vector{Float64},N::Int64,k₀::Number,layer::layerstructure)
    βRRr = collect(LinRange(βRec[1],βRec[2],N))
    βRRi = collect(LinRange(βImc[1],βImc[2],N))

    βMMr,βMMi = ndgrid(βRRr,βRRi)
    βMM = βMMr + im*βMMi

    nx,ny = size(βMM)
    βc = 0+0*im
    flag = false
    for i in 1:nx-1
        for j in 1:ny-1

            βI= βMM[i:i+1,j:j+1]

            βreal = real(βI[:,1])
            βImag = imag(βI[1,:])

            flag = IsDifπ(βreal,βImag,k₀,layer)
            if flag
                βRc = sum(βreal)/2
                βIc = sum(βImag)/2

                βc = βRc+im*βIc
            end
        end
    end
    
    βc
end

function Λedge(θ::Number,ϕ::Number,λ₀::Number,n::Number)
    θᵣ = (θ*π)/180
    ϕᵣ = (ϕ*π)/180
    
    λ₀/(-sin(θᵣ)*sin(ϕᵣ)+sqrt((sin(θᵣ)^2)*(sin(ϕᵣ)^2)-sin(θᵣ)^2+abs(n)^2))
end

function KsspSlab(k₀::Number,matd::MaterialParams,matM::MaterialParams)
    ε₁ = matd.ε(k₀)
    ε₂ = matM.ε(k₀)

    sqrt((ε₂*ε₁)/(ε₂+ε₁))*k₀
end