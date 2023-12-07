module ComplexColor

export complex_color

using Images

"""
    complex_color(A)

Convert an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

# Example
```jldoctest
julia> z = [0 im; -1 -im]
2×2 Matrix{Complex{Int64}}:
  0+0im  0+1im
 -1+0im  0-1im

julia> complex_color(z)
2×2 Array{RGB{Float64},2} with eltype ColorTypes.RGB{Float64}:
 RGB{Float64}(0.0,0.0,0.0)  RGB{Float64}(0.0,0.5,1.0)
 RGB{Float64}(1.0,0.0,1.0)  RGB{Float64}(1.0,0.5,0.0)
```
"""
function complex_color(A::AbstractArray{<:Complex{<:Real}})
    r = abs.(A)
    ϕ = angle.(A)
    H = @. rad2deg(mod(ϕ + 2π/3, 2π))
    L = @. 2atan(r)/π
    S = ones(eltype(H), size(H))
    clamp01nan!(map(RGB, HSL.(H, S, L)))
end

end
