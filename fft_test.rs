# This file serves to understand the concept of FFT and the derivatives in it

L = 1.0;

# The number of samples must be odd to avoid singularities from evaluating
# 1/sinc(2*k/N)
N = 1001;
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

# Single integral (integrates everything except the average term)
function intk!(arr)
    arr[1] = 0
    for k in 1:length(arr)-1
        arr[k+1] = arr[k+1] / (im * k * (2*p/L) * sinc(2*k/N))
    end
end

# Double integral (integrates everything except the average term)
function iintk!(arr)
    arr[1] = 0
    for k in 1:length(arr)-1
        arr[k+1] = -arr[k+1] / (k * (2*p/L) * sinc(k/N))^2
    end
end

z = convert(Array{Float64, 1}, [sin(10*p*(x/N)) for x in 1:N])
zk = fft(z)

iintk!(zk)

zd = real(ifft(zk))
plot(zd[10:end-10])
