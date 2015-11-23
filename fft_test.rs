# This file serves to understand the concept of FFT and the derivatives in it

L = 2.0;

# The number of samples must be odd to avoid singularities from evaluating
# 1/sinc(2*k/N)
N = 100001;
dx = L/N;

p = float(pi)

using PyPlot;

# Single derivative
function difk!(arr)
    for k in 0:length(arr)-1
        arr[k+1] = arr[k+1] * (im * k * (2*p/L) * sinc(2*k/N))
    end
end

# Double derivative
function ddifk!(arr)
    for k in 0:length(arr)-1
        arr[k+1] = -arr[k+1] * (k^2 * (2*p/L)^2 * sinc(k/N)^2)
    end
end

# TODO: Implement correct integration of the first fourier term
# (i.e. the constant term). For that, we'll have to evaluate the integral
# of the fourier series. The constant term, when integrated, will contribute
# an additional term of a_0*x which has to be Fourier-decomposed and added to
# the final series

# Single integral
function intk!(arr)
    # The DC term (it has to be accounted for in a different way, see below)
    a0 = arr[1]
    
    arr[1] = 0
    for k in 1:length(arr)-1
        arr[k+1] = arr[k+1] / (im * k * (2*p/L) * sinc(2*k/N)) + im*a0*L/(p*k)
    end
end

# Double integral (integrates everything except the average term)
function iintk!(arr)
    # The DC term (it has to be accounted for in a different way, see below)
    a0 = arr[1]
    
    arr[1] = 0
    for k in 1:length(arr)-1
        arr[k+1] = -(arr[k+1]) / (k * (2*p/L) * sinc(k/N))^2
        #arr[k+1] = arr[k+1] / (im * k * (2*p/L) * sinc(2*k/N)) + im*a0*L/(p*k)
        #arr[k+1] = arr[k+1] / (im * k * (2*p/L) * sinc(2*k/N)) + im*a0*L/(p*k)
    end
end

z = convert(Array{Float64, 1}, [(x/N) for x in 1:N])
zk = fft(z)
iintk!(zk)
ddifk!(zk)

# Test of a linear function decomposition
#zk = convert(Array{Complex128, 1}, zeros(N))
#zk[2:end] = [im*N/(p*k) for k in 1:N-1]

zd = real(ifft(zk))
plot(zd[1000:end-1000])
