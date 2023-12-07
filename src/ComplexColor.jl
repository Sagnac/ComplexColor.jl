module ComplexColor

export complex_color, complex_plot, hsl, latexstring, var"@L_str"

using Images
using GLMakie
using .Makie: Axis, latexstring

const Lims = Union{AbstractRange, Tuple, AbstractVector}
const ComplexArray = AbstractArray{<:Complex{<:Real}}
const Str = Union{AbstractString, Makie.LaTeXString}

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
function complex_color(A::ComplexArray)
    r = abs.(A)
    ϕ = angle.(A)
    H = @. rad2deg(mod(ϕ + 2π/3, 2π))
    L = @. 2atan(r)/π
    S = ones(eltype(H), size(H))
    clamp01nan1!(map(RGB, HSL.(H, S, L)))
end

"""
    complex_plot(x, y, s; [title])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL color space.
"""
function complex_plot(x::Lims, y::Lims, s::ComplexArray; title::Str = L"s")
    color_matrix = complex_color(s)
    xlen, ylen = size(s)
    nticks = 5
    xticklabels = latexstring.(range(x[begin], x[end], nticks))
    yticklabels = latexstring.(range(y[begin], y[end], nticks))
    xticks = (range(1, xlen, nticks), xticklabels)
    yticks = (range(1, ylen, nticks), yticklabels)
    arg_tick_labels = [L"-\pi", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi"]
    arg_ticks = (-π:π/2:π, arg_tick_labels)
    fig = Figure()
    axis = Axis(fig[1,1]; title, titlesize = 21,
                xlabel = L"Re(z)", xlabelsize = 16,
                ylabel = L"Im(z)", ylabelsize = 16,
                xticks, yticks)
    plot = image!(axis, color_matrix)
    Colorbar(fig[1,2]; colormap = hsl, limits = (-π, π), ticks = arg_ticks,
             label = "Arg(s)")
    fig
end

function clamp01nan1!(A)
    for (i, v) ∈ pairs(A)
        A[i] = mapc(v -> isnan(v) ? oneunit(v) : clamp(v, zero(v), oneunit(v)), v)
    end
    A
end

"HSL colormap"
hsl = map(RGB, HSL(i, 1.0, 0.5) for i = range(-60, 300, 2^10))

end
