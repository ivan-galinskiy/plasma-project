# This file serves to understand the concept of FFT and the derivatives in it

L = 2.0;

# The number of samples must be odd to avoid singularities from evaluating
# 1/sinc(2*k/N)
N = 10005;
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
        if k > length(arr)-100
            arr[k+1] = 0
        end
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

# Double integral
function iintk!(arr)
    # The DC term (it has to be accounted for in a different way, see below)
    a0 = arr[1]
    
    arr[1] = 0
    for k in 1:length(arr)-1
        arr[k+1] = -arr[k+1] / (k * (2*p/L) * sinc(2*k/N))^2 + (im*a0*L/(pi*k))^2
    end
end

z = convert(Array{Float64, 1}, [(2x - L)^2 for x in linspace(0, L, N)])
zk = fft(z)
#intk!(zk)
iintk!(zk)
#intk!(zk)
#difk!(zk)
#zk = fft(real(ifft(zk)))
#intk!(zk)
#ddifk!(zk)

# Test of a linear function decomposition
#zk = convert(Array{Complex128, 1}, zeros(N))
#zk[2:end] = [im*N/(p*k) for k in 1:N-1]

zd = real(ifft(zk))
#zd = real(z)
#plot(zd[9500:end])
plot(linspace(0, L, N-109), zd[100:end-10])
println(zd[9990] - zd[5000])
