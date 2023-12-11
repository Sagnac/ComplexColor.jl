module ComplexColor

export complex_color, complex_plot, var"@L_str"

using Colors
using GLMakie
using .Makie: latexstring

const Lims = Union{Tuple, AbstractVector}
const ComplexArray = AbstractArray{<:Complex{<:Real}}
const RealArray = AbstractArray{<:Real}

polar(s::ComplexArray; squared = false) = squared ? abs2.(s) : abs.(s), angle.(s)

"""
    complex_color(s)

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
function complex_color(s::ComplexArray; septaphase = false)
    r2, ϕ = polar(s; squared = true)
    complex_color_(r2, ϕ; septaphase)
end

function complex_color(r::RealArray, ϕ::RealArray; septaphase = false)
    complex_color_(abs2.(r), ϕ; septaphase)
end

# internal function (note the trailing underscore)
function complex_color_(r2::RealArray, ϕ::RealArray; septaphase = false)
    H = @. rad2deg(mod(ϕ + 2π/3, 2π))
    L = @. r2 / (r2 + 1)
    S = ones(eltype(H), size(H))
    septaphase && map!(hue -> 60 * fld(hue, 60), H, H)
    clamp01nan1!(map(RGB, HSL.(H, S, L)))
end

"""
    complex_plot(x, y, s; [title], [contours = true], [septaphase = false])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL color space.
`septaphase = true` will plot the phase using only 6 colors (green, cyan, blue, magenta, red, yellow).
"""
function complex_plot(x::Lims, y::Lims, s::ComplexArray;
                      title::AbstractString = L"s",
                      contours::Bool = true,
                      septaphase::Bool = false)
    xlen, ylen = size(s)
    nticks = 5
    xticklabels = latexstring.(range(x[begin], x[end], nticks))
    yticklabels = latexstring.(range(y[begin], y[end], nticks))
    xticks = (range(1, xlen, nticks), xticklabels)
    yticks = (range(1, ylen, nticks), yticklabels)
    arg_ticks = (-π:π:π, [L"-\pi", L"0", L"\pi"])
    fig = Figure()
    axis = Axis(fig[1,1]; title, titlesize = 21,
                xlabel = L"Re(z)", xlabelsize = 16,
                ylabel = L"Im(z)", ylabelsize = 16,
                xticks, yticks)
    if contours
        r, ϕ = polar(s)
        color_matrix = complex_color(r, ϕ; septaphase)
        image!(axis, color_matrix)
        contour!(axis, r; levels = [2.0^n for n = -3:8], colormap = Reverse(:acton))
        contour!(axis, ϕ; levels = -π:π/2:π, colormap = hsl)
    else
        color_matrix = complex_color(s; septaphase)
        image!(axis, color_matrix)
    end
    Colorbar(fig[1,2]; colormap = hsl, limits = (-π, π), ticks = arg_ticks,
             label = "Arg(s)")
    fig
end

function complex_plot(z::ComplexArray, s::ComplexArray; kw...)
    x = z |> real |> extrema
    y = z |> imag |> extrema
    complex_plot(x, y, s; kw...)
end

function clamp01nan1!(img::AbstractArray{<:Colorant})
    for (i, v) ∈ pairs(img)
        img[i] = mapc(v -> isnan(v) ? oneunit(v) : clamp(v, zero(v), oneunit(v)), v)
    end
    img
end

"HSL colormap"
const hsl = map(RGB, HSL(i, 1.0, 0.5) for i = range(-60, 300, 2^10))

end
