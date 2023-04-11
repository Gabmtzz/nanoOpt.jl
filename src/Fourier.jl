using FFTW

function wavefourier(M::Number)
    k = fftshift(2π/M * (0:M-1))
    k[k .> π] = k[k .> π] .- 2π 

    k
end

function ifourier2(x,y,ifun)
    Mx,My = size(ifun,1),size(ifun,2)
    hx = x[2]-x[1]; hy = y[2]-y[1]

    kx,ky = wavefourier(Mx) /hx, wavefourier(My) /hy

    z = x[1].*kx .+ y[1].*ky'

    ifun = reshape(ifun,(Mx*My,:))

    ifun = ifun.*exp.(im*z[:])/(hx*hy)
    ifun = reshape(ifun,(Mx,My,:))
    fun = ifft(ifft(ifftshift(ifftshift(ifun,2),1),2),1)
    fun[1:length(x),1:length(y),:]
end

function fourier2(x,y,fun,Mx=size(fun,1),My=size(fun,2))
    hx,hy = x[2]-x[1],y[2]-y[1]

    if length(size(fun)) ==2
        funfft = zeros(Mx,My)
        funfft[1:size(fun,1),1:size(fun,2)] = fun
    else
        funfft = zeros(Mx,My,3)*im
        funfft[1:size(fun,1),1:size(fun,2),:] = fun
    end
    ifun = fftshift(fftshift(fft(fft(funfft,1),2),1),2)

    kx,ky = wavefourier(Mx) /hx, wavefourier(My) /hy
    z = x[1].*kx .+ y[1].*ky'

    ifun = reshape(ifun,(Mx*My,:))

    ifun = ifun .* (hx*hy*exp.(-im*(z[:])))

    ifun = reshape(ifun,(Mx,My,:))

    if size(ifun,3) ==1
        ifun = ifun[:,:,1]
    end
    return(ifun,kx,ky)
end