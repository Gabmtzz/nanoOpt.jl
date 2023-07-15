struct circle <: Line
    a::Number
    θArr::Vector{Float64}
    xC::Vector{Float64}
    pos::Matrix{Float64}
    nV::Matrix{Float64}
    
    function circle(a::Number,N::Int64,xC::Tuple{Float64, Float64} = (0.0, 0.0) ,θlims::Tuple{Float64, Float64} = (0.0, 2π))
        θs,θe = θlims
        xc,yc = xC
        
        θArr = [((θe-θs)/N)*(i-0.5)+θs for i ∈ 1:N]
        
        xA,yA = a*cos.(θArr).+xc,a*sin.(θArr).+yc
        pos = [xA yA]
        
        nx,ny = cos.(θArr),sin.(θArr)
        nV = [nx ny]
        
        new(a,θArr,[xc,yc],pos,nV)
    end
end

struct lineR <: Line
    pos::Matrix{Float64}
    nV::Matrix{Float64}
    
    function lineR(xP::Tuple{Float64, Float64},yP::Tuple{Float64, Float64},N::Int64,nVp::Vector{Float64} = [1.0,0.0])
        xs,xe = xP
        ys,ye = yP
        
        xPos = collect(LinRange(xs,xe,N))
        yPos = collect(LinRange(ys,ye,N))
        
        pos = [xPos yPos]
        
        new(pos,nVp'.*ones(N,2))
    end
end

# ==========================================================================================
# ==========================================================================================

struct cylinder <: Structure
    circl::circle
    
    function cylinder(a::Number,N::Int64,xC::Tuple{Float64, Float64} = (0.0, 0.0))
        circl = circle(a,N,xC)
        new(circl)
    end
end

function getSurfPoints(cyl::cylinder)
    X,Y = cyl.circl.pos[:,1], cyl.circl.pos[:,2]

    Xn,Yn = zeros(length(X)+1),zeros(length(Y)+1)

    Xn[1:length(X)],Yn[1:length(Y)] = X,Y

    Xn[end],Yn[end] = X[1],Y[1]
    
    Xn,Yn
end

function numEl(cyl::cylinder)
    Pos = cyl.circl.pos
    size(Pos,1)
end

function getSvec(m::Int64,cyl::cylinder)
    rpos = collect(LinRange(0,1,m+1))[:]

    a = cyl.circl.a
    θa = cyl.circl.θArr
    n = length(θa)
    xₚ,yₚ = cyl.circl.xC[1],cyl.circl.xC[2]
    SArr = Array{Array{Float64,1}}(undef,n,length(rpos))

    for i ∈ eachindex(θa)
        if i == n
            Δθ = 2π - (θa[end] - θa[1])
        else
            Δθ = θa[i+1] - θa[i]
        end

        for j ∈ eachindex(rpos)
            r = rpos[j]

            θᵣ = Δθ*r + θa[i]

            SArr[i,j] = [a*cos(θᵣ)+xₚ; a*sin(θᵣ)+yₚ]
        end
    end
    
    rpos,SArr
end

function SQuad(t::Number,i::Int64,cyl::cylinder)
    a = cyl.circl.a
    θa = cyl.circl.θArr
    n = length(θa)
    xₚ,yₚ = cyl.circl.xC[1],cyl.circl.xC[2]
    
    if i == n
        Δθ = 2π - (θa[end] - θa[1])
    else
        Δθ = θa[i+1] - θa[i]
    end
    
    θᵣ = Δθ*t + θa[i]
    
    [a*cos(θᵣ)+xₚ; a*sin(θᵣ)+yₚ]
end

function nQuad(t::Number,i::Int64,cyl::cylinder)
    θa = cyl.circl.θArr
    n = length(θa)
    xₚ,yₚ = cyl.circl.xC[1],cyl.circl.xC[2]
    
    if i == n
        Δθ = 2π - (θa[end] - θa[1])
    else
        Δθ = θa[i+1] - θa[i]
    end
    
    θᵣ = Δθ*t + θa[i]
    
    [cos(θᵣ); sin(θᵣ)]
end

# ==========================================================================================
# ==========================================================================================

struct Rod <: Structure
    RodELems::Vector{Line}
    boundary::Vector{Float64}
    
    function Rod(w::Number,d::Number,rc::Number,Nc::Int64,Nw::Int64,Nd::Int64,xC::Tuple{Float64, Float64} = (0.,0.))
        xₚ,yₚ = xC 
        wᵣ,dᵣ = w-2*rc,d-2*rc

        dSr,dIr = yₚ+dᵣ/2,yₚ-dᵣ/2
        wLr,wRr = xₚ-wᵣ/2,xₚ+wᵣ/2
        wL,wR = xₚ-w/2,xₚ+w/2
        dS,dI = yₚ+d/2,yₚ-d/2

        xC1,yC1,xC2,yC2 = wLr,dSr,wRr,dIr
        
        wUp = lineR((wLr,wRr),(dS,dS),Nw,[0.,1.])
        dRg = lineR((wR,wR),(dSr,dIr),Nd,[1.,0.])
        wDw = lineR((wRr,wLr),(dI,dI),Nw,[0.,-1.])
        dLf = lineR((wL,wL),(dIr,dSr),Nd,[-1.,0.])
        
        circ1 = circle(rc,Nc,(xC2,yC1),(π/2,0.))
        circ2 = circle(rc,Nc,(xC2,yC2),(0.,-π/2))
        circ3 = circle(rc,Nc,(xC1,yC2),(-π/2,-π))
        circ4 = circle(rc,Nc,(xC1,yC1),(1π,π/2))
        
        boundary = [π/2, 0., 0., -π/2, -π/2, -1π, 1π, π/2]
        
        RodELems = [wUp,circ1,dRg,circ2,wDw,circ3,dLf,circ4]
        
        new(RodELems,boundary)
    end
end

function numEl(rod::Rod)
    sis = 0
    rEl = rod.RodELems
    for i ∈ eachindex(rEl)
       sis +=  size(rEl[i].pos,1)
    end

    sis
end

function getSurfPoints(rod::Rod)
    N = numEl(rod1)

    X,Y = zeros(N+1),zeros(N+1)

    pos=1
    rEl = rod.RodELems

    for i ∈ eachindex(rEl)
        xE,yE = rEl[i].pos[:,1],rEl[i].pos[:,2]
        nx,ny = length(xE),length(yE)
    
        X[pos:(pos-1)+nx] = xE
        Y[pos:(pos-1)+ny] = yE
    
        pos += nx
    end

    X[end] = X[1]
    Y[end] = Y[1]
    
    X,Y
end

function getSvec(m::Int64,rod::Rod)
    el =1
    ind = 1

    rpos = collect(LinRange(0,1,m+1))[:]
    n = numEl(rod1)

    ElRod = rod.RodELems 
    bound = rod.boundary

    SArr = Array{Array{Float64,1}}(undef,n,length(rpos)) 
    
    for i ∈ 1:n
        iRod = ElRod[el]
        Pos = iRod.pos
        if typeof(iRod) == lineR
            if ind == size(Pos,1)
                θᵢ = bound[el]
                el += 1
                θf = ElRod[el].θArr[1]
                a =  ElRod[el].a
                xₚ,yₚ = ElRod[el].xC[1],ElRod[el].xC[2]
                Δθ = θf - θᵢ
        
                for j ∈ eachindex(rpos)
                    r = rpos[j]

                    θᵣ = Δθ*r + θᵢ

                    SArr[i,j] = [a*cos(θᵣ)+xₚ; a*sin(θᵣ)+yₚ]
                end
                ind = 1
            else
        
                x,y = Pos[ind,1], Pos[ind,2]
                xn,yn = Pos[ind+1,1], Pos[ind+1,2]
                Δx,Δy = xn-x,yn-y
    
                for j ∈ eachindex(rpos)
                    r = rpos[j]
    
                    SArr[i,j] = [Δx*r+x; Δy*r+y]
                end
                ind += 1
            end
        elseif typeof(iRod) == circle
            θa = ElRod[el].θArr
            a =  ElRod[el].a
            xₚ,yₚ = ElRod[el].xC[1],ElRod[el].xC[2]
    
            θi = θa[ind]
    
    
            if ind == size(Pos,1)
                    θf  = bound[el]
            
                    Δθ = θf-θi
                    for j ∈ eachindex(rpos)
                        r = rpos[j]
                         θᵣ = Δθ*r + θa[ind]
                         SArr[i,j] = [a*cos(θᵣ)+xₚ; a*sin(θᵣ)+yₚ]
                    end
                    ind = 1
                    if el != 8
                        el +=1
                    end
                else
                    θf  = θa[ind+1]
                    Δθ = θf-θi
                    for j ∈ eachindex(rpos)
                        r = rpos[j]
                         θᵣ = Δθ*r + θa[ind]
                         SArr[i,j] = [a*cos(θᵣ)+xₚ; a*sin(θᵣ)+yₚ]
                    end
                ind += 1
            end
        
        end
    end
    rpos,SArr
end