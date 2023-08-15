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

function getnVPoints(cyl::cylinder)
    X,Y = cyl.circl.nV[:,1], cyl.circl.nV[:,2]

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
    Nels::Vector{Int64}
    
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
        
        Nels = [Nw,Nc,Nd,Nc,Nw,Nc,Nd,Nc]
        
        new(RodELems,boundary,Nels)
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
    N = numEl(rod)

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

function getnVPoints(rod::Rod)
    N = numEl(rod)

    X,Y = zeros(N+1),zeros(N+1)

    pos=1
    rEl = rod.RodELems

    for i ∈ eachindex(rEl)
        xE,yE = rEl[i].nV[:,1],rEl[i].nV[:,2]
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
    n = numEl(rod)

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

function GetIndPos(i::Int64,Nels::Vector{Int64})
    sElms = cumsum(Nels)
    
    flag = true
    ind = 1
    beg = 0
    while flag
        if i ≥ beg && i ≤ sElms[ind]
            flag = false
        else
            beg = sElms[ind]
            ind += 1
        end
    end
    
    ind,i-beg
end

function SQuad(t::Number,i::Int64,rod::Rod)
    ElRod = rod.RodELems 
    bound = rod.boundary
    Nels = rod.Nels
    nEl,ind = GetIndPos(i,Nels)
    
    iRod = ElRod[nEl]
    Pos = iRod.pos
    sVecQ = [0.;0.]
    if typeof(iRod) == lineR
        if ind == size(Pos,1)
            θi = bound[nEl]
            θf = ElRod[nEl+1].θArr[1]
            a =  ElRod[nEl+1].a
            xₚ,yₚ = ElRod[nEl+1].xC[1],ElRod[nEl+1].xC[2]
            Δθ = θf - θi
        
            θₜ = Δθ*t + θi
            sVecQ = [a*cos(θₜ)+xₚ; a*sin(θₜ)+yₚ]
        else
            xi,yi = iRod.pos[ind,1],iRod.pos[ind,2]
            xe,ye = iRod.pos[ind+1,1],iRod.pos[ind+1,2]
    
            Δx,Δy = xe-xi,ye-yi
            sVecQ = [Δx*t+xi; Δy*t+yi]
        end
    elseif typeof(iRod) == circle
        θa = ElRod[nEl].θArr
        a =  ElRod[nEl].a
        xₚ,yₚ = ElRod[nEl].xC[1],ElRod[nEl].xC[2]
    
        θi = θa[ind]
    
        if ind == size(Pos,1)
            θf = bound[nEl]
        else
            θf = θa[ind+1]
        end
    
        Δθ = θf-θi
        θₜ = Δθ*t + θi
        sVecQ = [a*cos(θₜ)+xₚ; a*sin(θₜ)+yₚ]
    end  
    
    sVecQ
end

function nQuad(t::Number,i::Int64,rod::Rod)
    ElRod = rod.RodELems 
    bound = rod.boundary
    Nels = rod.Nels
    nEl,ind = GetIndPos(i,Nels)
    
    iRod = ElRod[nEl]
    Pos = iRod.nV
    nVecQ = [0.;0.]
    if typeof(iRod) == lineR
        if ind == size(Pos,1)
            θi = bound[nEl]
            θf = ElRod[nEl+1].θArr[1]
            Δθ = θf - θi
        
            θₜ = Δθ*t + θi
            nVecQ = [cos(θₜ); sin(θₜ)]
        else
            nxi,nyi = iRod.nV[ind,1],iRod.nV[ind,2]
            nxe,nye = iRod.nV[ind+1,1],iRod.nV[ind+1,2]
    
            Δnx,Δny = nxe-nxi,nye-nyi
            nVecQ = [Δnx*t+nxi; Δny*t+nyi]
        end
    elseif typeof(iRod) == circle
        θa = ElRod[nEl].θArr
    
        θi = θa[ind]
    
        if ind == size(Pos,1)
            θf = bound[nEl]
        else
            θf = θa[ind+1]
        end
    
        Δθ = θf-θi
        θₜ = Δθ*t + θi
        nVecQ = [cos(θₜ); sin(θₜ)]
    end  
    
    nVecQ
end

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

struct StructureC <: Structure
    elements::Vector{Structure}
    boundEls::Vector{Int64}
    
    function StructureC(ss)
        bEl = cumsum(numEl.(ss))
        
        new(ss,bEl)
    end
end

function getIndexStr(i::Int64,Structures::StructureC)
    indexB = Structures.boundEls
    
    flag = true
    indB = 1
    i0 = 0
    indEl = 1

    while flag
        elmB = indexB[indB]
    
        if i > i0  && i ≤ elmB
           indEl = i - i0
            flag = false
        else
            i0 = elmB
        
            if indB < length(indexB)
                indB += 1
            end
        end
    end
    indB,indEl
end

function numEl(Structures::StructureC)
    aRReL = Structures.elements
    sn = 0

    for i ∈ eachindex(aRReL)
    
        sn += numEl(aRReL[i])
    end
    
    sn
end

function getSvec(m::Int64,Structures::StructureC)
    aRReL = Structures.elements
    
    rpos,SArr =  getSvec(m,aRReL[1])
    
    for i ∈ 2:length(aRReL)
    
        _,SArrE =  getSvec(m,aRReL[i])
        SArr=vcat(SArr,SArrE)
    end

    rpos,SArr
end    

function SQuad(t::Number,i::Int64,Structures::StructureC)
    indB,indEl = getIndexStr(i,Structures)

    strEl = Structures.elements[indB]

    SQuad(t,indEl,strEl)
end

function nQuad(t::Number,i::Int64,Structures::StructureC)
    indB,indEl = getIndexStr(i,Structures)

    strEl = Structures.elements[indB]

    nQuad(t,indEl,strEl)
end

# ==========================================================================================
# ==========================================================================================
# ==========================================================================================
# ==========================================================================================

function IsInside(x::Number,y::Number,str::Structure)
    N = numEl(str)

    px,py = getSurfPoints(str)
    Apos = [px py]
    nx,ny = getnVPoints(str)
    An = [nx ny]
    lZeroArr = Array{Bool}(undef,N)
    
    for i ∈ 1:N
        S,n = Apos[i,:],An[i,:]
        r = [x; y]
        vec = (r-S)/norm(r-S)
        lZeroArr[i] = n ⋅ vec ≤ 0.0
    end
    
    isempty(findall(lZeroArr .== false))
end

function getLineSep(j::Int64,xx::Matrix{Float64},yy::Matrix{Float64},str::Structure)
    n = size(xx,1)

    Xout,Yout = Array{Float64}(undef,0),Array{Float64}(undef,0)
    Xin,Yin = Array{Float64}(undef,0),Array{Float64}(undef,0)


    for i ∈ 1:n
        x,y = xx[i,j],yy[i,j]
        if IsInside(x,y,str)
            push!(Xin,x); push!(Yin,y)
        else
            push!(Xout,x); push!(Yout,y)
        end
    end
    
    Xin,Yin,Xout,Yout
end

function getMesh(xA::Vector{Float64},yA::Vector{Float64},str::Structure)
    xx,yy = ndgrid(xA,yA)
    xx,yy = collect(xx),collect(yy)
    Xin,Yin,Xout,Yout = Array{Float64}(undef,0),Array{Float64}(undef,0),Array{Float64}(undef,0),Array{Float64}(undef,0)
    
    for j ∈ axes(xx,2)
        Xi,Yi,Xo,Yo = getLineSep(j,xx,yy,str)
        Xin,Xout,Yin,Yout = vcat(Xin,Xi),vcat(Xout,Xo),vcat(Yin,Yi),vcat(Yout,Yo)
    end
    
    Xin,Xout,Yin,Yout
end