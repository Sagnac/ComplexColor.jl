module ComplexColor

export complex_color, complex_plot,
       @L_str, @iv_str, Interval, OpenInterval, ±, (..),
       RGB, HSL, Oklch

using Printf
using Colors
using GLMakie
using .Makie: latexstring
using IntervalSets

const ComplexArray = AbstractArray{<:Complex}
const RealArray = AbstractArray{<:Real}

const Spaces = Union{Type{Oklch}, Type{HSL}}

const chroma = 0.35

# base interval range length
const n = 1000

struct Septaphase end

const colors_compat = pkgversion(Colors) <= v"0.12.11"

const default = colors_compat ? HSL : Oklch

__init__() = colors_compat && @warn "Oklch color is not available in this version."

"""
    complex_color(s [, color = `$default`])

Convert an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

# Example
```jldoctest
julia> z = [0 im; -1 -im]
2×2 Matrix{Complex{Int64}}:
  0+0im  0+1im
 -1+0im  0-1im

julia> complex_color(z)
2×2 Array{RGB{Float64},2} with eltype RGB{Float64}:
 RGB{Float64}(0.0,0.0,0.0)  RGB{Float64}(0.0,0.5,1.0)
 RGB{Float64}(1.0,0.0,1.0)  RGB{Float64}(1.0,0.5,0.0)
```
"""
function complex_color(s::ComplexArray, color::Spaces = default)
    to_RGB(to_color(s, color), color)
end

function complex_color(s::ComplexArray, ::Septaphase, color::Spaces)
    t = to_color(s, color)
    color_matrices = Vector{Matrix{RGB}}(undef, 5)
    t1, t2, t3 = revif(t, color)
    rounded, thresholded = septaphase(t3, color)
    _0_ = zeros(eltype(t1), size(t1))
    ϕ01 = normalize_phase(t3, color)
    color_matrices[1] = to_RGB(t, color)
    color_matrices[2] = rev_to_RGB(t1, t2, rounded; color)
    color_matrices[3] = rev_to_RGB(t1, t2, thresholded; color)
    color_matrices[4] = to_RGB((_0_, _0_, t1), HSL)
    color_matrices[5] = to_RGB((_0_, _0_, ϕ01), HSL)
    return color_matrices
end

revif(t::NTuple{3, RealArray}, color::Spaces) = color == HSL ? reverse(t) : t

function Λ(s)
    r2 = abs2.(s)
    _1_ = one(eltype(r2))
    @. r2 / (_1_ + r2)
end

degrees(s) = @. rad2deg(angle(s))

mod360(s::RealArray, offset) = @. mod(s + offset, 360)

mod360(s::ComplexArray, offset) = mod360(degrees(s), offset)

normalize_phase(H, color::Type{Oklch}) = mod360(H, -150) / 360

normalize_phase(H, color::Type{HSL}) = mod360(H, -120) / 360

function to_color(s::ComplexArray, color::Type{Oklch})
    L = Λ(s)
    C = fill(chroma, size(L))
    H = mod360(s, 150)
    return L, C, H
end

function to_color(s::ComplexArray, color::Type{HSL})
    H = mod360(s, 120)
    S = ones(size(H))
    L = Λ(s)
    return H, S, L
end

to_RGB(t, color) = map(RGB, color.(t...)) |> clamp01nan1!

function rev_to_RGB(t...; color)
    t = revif(t, color)
    to_RGB(t, color)
end

function septaphase(H, color::Spaces)
    if color == Oklch
        H = mod360(H, -30)
    end
    φ = 60
        rounded = map(hue -> φ * round(hue / φ), H)
    thresholded = map(hue -> φ * floor(hue / φ), H)
    return rounded, thresholded
end

function draw_modulus_contours(axis, r)
    levels = exp2.(-3:15)
    colormap = Reverse(:acton)
    contour!(axis, r; levels, colormap, inspectable = false)
end

function draw_phase_contours(axis, ϕ, colormap)
    contour!(axis, ϕ; levels = -180:90:180, colormap, inspectable = false)
end

switch_button_color(state) = state ? RGBf(0.0, 1.0, 1.0) : RGBf(0.94, 0.94, 0.94)

"""
    complex_plot(x, y, s [,color = `$default`]; [title])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the `HSL` or `Oklch` color spaces.

The three-valued `septaphase` slider option on the plot will partition the phase using only 6 colors (green, cyan, blue, magenta, red, yellow) either by thresholding or by rounding, depending on the setting; the default off position plots a gradient phase.
"""
function complex_plot(x::AbstractVector, y::AbstractVector, s::ComplexArray,
                      color::Spaces = default; title::AbstractString = L"f(z)")
    xlen = length(x)
    ylen = length(y)
    (xlen, ylen) == size(s) || error("Length mismatch.") 
    nticks = 5
    xtickrange = range(x[begin], x[end], nticks)
    ytickrange = range(y[begin], y[end], nticks)
    xticklabels = map(i -> latexstring(@sprintf("%.4G", i)), xtickrange)
    yticklabels = map(i -> latexstring(@sprintf("%.4G", i)), ytickrange)
    xticks = (range(1, xlen, nticks), xticklabels)
    yticks = (range(1, ylen, nticks), yticklabels)
    arg_ticks = (-π:π:π, [L"-\pi", L"0", L"\pi"])
    fig = Figure(size = (600, 532))
    axis = Axis(fig[1,2]; title, titlesize = 21,
                xlabel = L"Re", xlabelsize = 16,
                ylabel = L"Im", ylabelsize = 16,
                xticks, yticks)
    r = abs.(s)
    ϕ = degrees(s)
    function inspector(_, inds, _)
        i, j = round.(Int, inds)
        str = @sprintf(
            """
            x: %.4G, y: %.4G
            r: %.4G, ϕ: %.2f\u00b0
            """,
            x[i], y[j], r[i,j], ϕ[i,j]
        )
        replace(str, '-' => '\u2212')
    end
    color_matrices = complex_color(s, Septaphase(), color)
    img = image!(axis, color_matrices[1]; inspector_label = inspector)
    prev_img = color_matrices[1]
    local phase_contours
    modulus_contours = draw_modulus_contours(axis, r)
    colormap = colormaps[color]
    Colorbar(fig[1,3]; colormap, limits = (-π, π), ticks = arg_ticks,
             label = "Arg")
    grid1 = GridLayout(fig[1,1]; tellheight = false, valign = :bottom)
    grid2 = GridLayout(fig[2,2])
    modulus_on = Observable(false)
    phase_on = Observable(false)
    modulus_btn_clr = lift(switch_button_color, modulus_on)
    phase_btn_clr = lift(switch_button_color, phase_on)
    modulus_button = Button(fig; buttoncolor = modulus_btn_clr, label = "modulus",
                            halign = :right)
    phase_button = Button(fig; buttoncolor = phase_btn_clr, label = "phase",
                          halign = :right)
    modulus_contours_toggle = Toggle(fig; active = true)
    phase_contours_toggle = Toggle(fig; active = false)
    septaphase_slider = Slider(fig; range = 1:3, width = 48,
                               linewidth = 12, halign = :left)
    grid1[1,1] = modulus_button
    grid1[2,1] = phase_button
    grid1[3,1] = grid!([Label(fig, "septaphase") septaphase_slider])
    grid2[1,1] = grid!([Label(fig, "phase contours") phase_contours_toggle])
    grid2[1,2] = grid!([Label(fig, "modulus contours") modulus_contours_toggle])
    on(septaphase_slider.value) do value
        modulus_on[] = false
        phase_on[] = false
        prev_img = img[3][] = color_matrices[value]
    end
    on(modulus_button.clicks) do _
        phase_on[] = false
        modulus_on[] = !modulus_on[]
        img[3][] = modulus_on[] ? color_matrices[4] : prev_img
    end
    on(phase_button.clicks) do _
        modulus_on[] = false
        phase_on[] = !phase_on[]
        img[3][] = phase_on[] ? color_matrices[5] : prev_img
    end
    on(phase_contours_toggle.active) do active
        if active
            phase_contours = draw_phase_contours(axis, ϕ, colormap)
        else
            delete!(axis, phase_contours)
        end
    end
    on(modulus_contours_toggle.active) do active
        if active
            modulus_contours = draw_modulus_contours(axis, r)
        else
            delete!(axis, modulus_contours)
        end
    end
    ar = true
    function aspect_control()
        axis.aspect = ar ? DataAspect() : nothing
        colsize!(fig.layout, 2, ar ? Aspect(1, xlen / ylen) : Auto(false, 1.0))
        resize_to_layout!(fig)
    end
    onmouserightup(addmouseevents!(fig.scene)) do _
        ispressed(fig, Keyboard.left_control) || return
        ar = !ar
        aspect_control()
    end
    aspect_control()
    DataInspector(fig)
    fig
end

function complex_plot(xlims::T, ylims::T, s::ComplexArray, color::Spaces = default;
                      kw...) where {S <: Real, T <: Tuple{S, S}}
    xlen, ylen = size(s)
    x = range(xlims..., xlen)
    y = range(ylims..., ylen)
    complex_plot(x, y, s, color; kw...)
end

function complex_plot(z::ComplexArray, s::ComplexArray,
                      color::Spaces = default; kw...)
    x = extrema(real, z)
    y = extrema(imag, z)
    complex_plot(x, y, s, color; kw...)
end

# generic methods where f is a complex function
function complex_plot(x::AbstractVector, y::AbstractVector, f,
                      color::Spaces = default; kw...)
    z = complex.(x, y')
    s = f.(z)
    complex_plot(x, y, s, color; kw...)
end

Δ(x::AbstractVector) = x[end] - x[begin]
Δ(x::Interval) = IntervalSets.width(x)

adjust_len(x, y, len) = round(Int, len * Δ(x) / Δ(y))

function complex_plot(x::Interval, y::Interval, f, color::Spaces = default; kw...)
    x = range(x, adjust_len(x, y, n))
    y = range(y, n)
    complex_plot(x, y, f, color; kw...)
end

function complex_plot(x::AbstractVector, y::Interval, f, color::Spaces = default;
                      kw...)
    complex_plot(x, range(y, adjust_len(y, x, length(x))), f, color; kw...)
end

function complex_plot(x::Interval, y::AbstractVector, f, color::Spaces = default;
                      kw...)
    complex_plot(range(x, adjust_len(x, y, length(y))), y, f, color; kw...)
end

function complex_plot(z::ComplexArray, f, color::Spaces = default; kw...)
    complex_plot(z, f.(z))
end

function clamp01nan1!(img::AbstractArray{<:Colorant})
    for (i, v) ∈ pairs(img)
        img[i] = mapc(v -> isnan(v) ? oneunit(v) : clamp(v, zero(v), oneunit(v)), v)
    end
    img
end

const colormaps = Dict{Spaces, Vector{RGB{Float64}}}(
    HSL => map(RGB, HSL(i, 1.0, 0.5) for i = range(-60, 300, 2^10)),
)

if !colors_compat
    colormaps[Oklch] = map(RGB, Oklch(0.5, chroma, i) for i = range(-30, 330, 2^10))
end

end
