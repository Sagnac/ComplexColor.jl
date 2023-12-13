module ComplexColor

export complex_color, complex_plot, var"@L_str"

using Colors
using GLMakie
using .Makie: latexstring

const ComplexArray = AbstractArray{<:Complex{<:Real}}
const RealArray = AbstractArray{<:Real}

struct Septaphase end

polar2(s::ComplexArray) = abs2.(s), angle.(s)

polar3(s::ComplexArray) = abs.(s), abs2.(s), angle.(s)

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
complex_color(s::ComplexArray) = HSL_to_RGB(polar2HSL(polar2(s)...))

complex_color(r::RealArray, ϕ::RealArray) = HSL_to_RGB(polar2HSL(abs2.(r), ϕ))

function complex_color(r2::RealArray, ϕ::RealArray, ::Septaphase)
    H, S, L = polar2HSL(r2, ϕ)
    HSL_to_RGB(H, S, L), septaphases(H, S, L)...
end

function polar2HSL(r2::RealArray, ϕ::RealArray)
    H = @. rad2deg(mod(ϕ + 2π/3, 2π))
    L = @. r2 / (r2 + 1)
    S = ones(eltype(H), size(H))
    return H, S, L
end

HSL_to_RGB(H, S, L) = map(RGB, HSL.(H, S, L)) |> clamp01nan1!

function septaphases(H, S, L)
        rounded = HSL_to_RGB(map(hue -> 60 * div(hue, 60, RoundNearest), H), S, L)
    thresholded = HSL_to_RGB(map(hue -> 60 * fld(hue, 60), H), S, L)
    return rounded, thresholded
end

"""
    complex_plot(x, y, s; [title], [contours = true], [septaphase = false])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL color space.
`septaphase = true` will plot the phase using only 6 colors (green, cyan, blue, magenta, red, yellow).
"""
function complex_plot(x::AbstractVector, y::AbstractVector, s::ComplexArray;
                      title::AbstractString = L"s",
                      contours::Bool = true,
                      septaphase::Bool = false)
    xlen = length(x)
    ylen = length(y)
    (xlen, ylen) == size(s) || error("Length mismatch.") 
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
    r, r2, ϕ = polar3(s)
    grad_phase, thresh_phase, round_phase = complex_color(r2, ϕ, Septaphase())
    img = image!(axis, grad_phase)
    contour!(axis, r; levels = [2.0^n for n = -3:8], colormap = Reverse(:acton))
    contour!(axis, ϕ; levels = -π:π/2:π, colormap = hsl)
    Colorbar(fig[1,2]; colormap = hsl, limits = (-π, π), ticks = arg_ticks,
             label = "Arg(s)")
    fig
end

function complex_plot(xlims::T, ylims::T, s::ComplexArray;
                      kw...) where {S <: Real, T <: Tuple{S, S}}
    xlen, ylen = size(s)
    x = range(xlims..., xlen)
    y = range(ylims..., ylen)
    complex_plot(x, y, s; kw...)
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
