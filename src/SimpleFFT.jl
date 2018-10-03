module SimpleFFT

using AbstractFFTs

import AbstractFFTs: plan_fft
import LinearAlgebra: mul!
import Base: *

mutable struct SimplePlan <: AbstractFFTs.Plan{Complex}
    pinv::AbstractFFTs.Plan
    SimplePlan() = new()
end

function plan_fft(x::AbstractVector, region; kws...)
    println(region)
    @assert region == 1:1
    SimplePlan()
end

function fft2(x)
    N = length(x)
    if N == 1
        return x
    else
        xe = fft2(x[1:2:end])
        xo = fft2(x[2:2:end])
        f = [exp(-2pi*im*k/N) for k in 0:N/2-1]
        vcat(xe+f.*xo, xe - f.*xo)
    end
end

function mul!(y, p::SimplePlan, x::AbstractVector)
    xhat = fft2(x)
    copyto!(y, xhat)
end

*(p::SimplePlan, x::AbstractVector) = mul!(similar(x, Complex), p, x)

end # module
