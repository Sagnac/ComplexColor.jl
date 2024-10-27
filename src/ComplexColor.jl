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

const VecOrInterval = Union{AbstractVector, Interval}

const Spaces = Union{Type{Oklch}, Type{HSL}}

const chroma = 0.35

# base interval range length
const n = 500

const modulus_levels = exp2.(-7:15)
const modulus_colormap = Reverse(:acton)

const default = Oklch

@kwdef struct Coordinates
    x::Vector{Float64}
    y::Vector{Float64}
    r::Matrix{Float64}
    ϕ::Matrix{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
end

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
    to_RGB(to_color(polar_coords(s)..., color), color)
end

function complex_color(r, ϕ, u, v, color::Spaces)
    t = to_color(r, ϕ, color)
    color_matrices = Vector{Matrix{RGB}}(undef, 7)
    t1, t2, t3 = revif(t, color)
    rounded, thresholded = septaphase(t3, color)
    _0_ = zeros(eltype(t1), size(t1))
    ϕ = mod.(ϕ, 360) / 360
    color_matrices[1] = to_RGB(t, color)
    color_matrices[2] = rev_to_RGB(t1, t2, rounded; color)
    color_matrices[3] = rev_to_RGB(t1, t2, thresholded; color)
    color_matrices[4] = to_RGB((_0_, _0_, t1), HSL)
    color_matrices[5] = to_RGB((_0_, _0_, ϕ), HSL)
    color_matrices[6] = to_RGB((_0_, _0_, ξ(u)), HSL)
    color_matrices[7] = to_RGB((_0_, _0_, ξ(v)), HSL)
    return color_matrices
end

revif(t::NTuple{3, RealArray}, color::Spaces) = color == HSL ? reverse(t) : t

function Λ(r)
    r2 = r .^ 2
    r2 ./ (one(eltype(r2)) .+ r2)
end

ξ(w) = inv.(one(eltype(w)) .+ exp.(-w))

mod360(a::RealArray, offset) = @. mod(a + offset, 360)

function polar_coords(s::ComplexArray)
    r = abs.(s)
    ϕ = @. rad2deg(angle(s))
    return r, ϕ
end

function coords(s::ComplexArray)
    r, ϕ = polar_coords(s)
    u, v = reim(s)
    return r, ϕ, u, v
end

function to_color(r, ϕ, color::Type{Oklch})
    L = Λ(r)
    C = fill(chroma, size(L))
    H = mod360(ϕ, 150)
    return L, C, H
end

function to_color(r, ϕ, color::Type{HSL})
    H = mod360(ϕ, 120)
    S = ones(size(H))
    L = Λ(r)
    return H, S, L
end

to_RGB(t, color) = map(RGB, color.(t...)) |> clamp01nan1!

function rev_to_RGB(t...; color)::Matrix{RGB}
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

function draw_modulus_contours(axis, r, levels)
    colormap = modulus_colormap
    contour!(axis, r; levels, colormap, inspectable = false)
end

function draw_phase_contours(axis, ϕ, colormap)
    contour!(axis, ϕ; levels = -180:90:180, colormap, inspectable = false)
end

switch_button_color(state) = state ? RGBf(0.0, 1.0, 1.0) : RGBf(0.94, 0.94, 0.94)

function inspector(inds, params, uv)
    (; x, y, r, ϕ, u, v) = params
    i, j = round.(Int, inds)
    if uv[]
        select = @sprintf("u: %.4G, v: %.4G", u[i,j], v[i,j])
    else
        select = @sprintf("r: %.4G, ϕ: %.2f\u00b0", r[i,j], ϕ[i,j])
    end
    str = @sprintf("x: %.4G, y: %.4G\n%s", x[i], y[j], select)
    replace(str, '-' => '\u2212')
end

function aspect_control!(ar, axis, fig, xlen, ylen)
    axis.aspect = ar ? DataAspect() : nothing
    colsize!(fig.layout, 2, ar ? Aspect(1, xlen / ylen) : Auto(false, 1.0))
    resize_to_layout!(fig)
end

"""
    complex_plot(x, y, s [,color = `$default`]; [title])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the `HSL` or `Oklch` color spaces.

The three-valued `septaphase` slider option on the plot will partition the phase using only 6 colors (green, cyan, blue, magenta, red, yellow) either by thresholding or by rounding, depending on the setting; the default off position plots a gradient phase.
"""
function complex_plot(x::AbstractVector, y::AbstractVector, s::ComplexArray,
                      color::Spaces = default; title::AbstractString = L"f(z)",
                      xlabel = L"Re", ylabel = L"Im", levels = modulus_levels,
                      nticks = 5, res = (600, 532), titlesize = 21,
                      xlabelsize = 16, ylabelsize = 16)
    xlen = length(x)
    ylen = length(y)
    (xlen, ylen) == size(s) || error("Length mismatch.")
    xtickrange = range(x[begin], x[end], nticks)
    ytickrange = range(y[begin], y[end], nticks)
    xticklabels = map(i -> latexstring(@sprintf("%.4G", i)), xtickrange)
    yticklabels = map(i -> latexstring(@sprintf("%.4G", i)), ytickrange)
    xticks = (range(1, xlen, nticks), xticklabels)
    yticks = (range(1, ylen, nticks), yticklabels)
    arg_ticks = (-π:π:π, [L"-\pi", L"0", L"\pi"])
    fig = Figure(size = res)
    axis::Axis = Axis(fig[1,2]; title, titlesize,
                      xlabel, xlabelsize,
                      ylabel, ylabelsize,
                      xticks, yticks)
    r, ϕ, u, v = coords(s)
    color_matrices = complex_color(r, ϕ, u, v, color)
    img = image!(axis, color_matrices[1];
                 inspector_label = (_, inds, _) -> inspector(inds, params, uv))
    prev_img = color_matrices[1]
    local phase_contours, modulus_contours
    r_const, ϕ_const = [all(isapprox(first(i)), i) for i in (r, ϕ)]
    colormap = colormaps[color]
    Colorbar(fig[1,3]; colormap, limits = (-π, π), ticks = arg_ticks,
             label = "Arg")
    grid1 = GridLayout(fig[1,1]; tellheight = false, valign = :center)
    grid2 = GridLayout(fig[2,2])
    n_btns = 1:4
    btn_active = [Observable(false) for i = n_btns]
    btn_clr = [lift(switch_button_color, btn) for btn in btn_active]
    btn_labels = ["modulus", "phase", "real", "imaginary"]
    buttons = [Button(fig; buttoncolor = btn_clr[i], label = btn_labels[i],
                      halign = :right) for i = n_btns]
    uv = @lift($(btn_active[3]) || $(btn_active[4]))
    modulus_contours_toggle = Toggle(fig; active = true)
    phase_contours_toggle = Toggle(fig; active = false)
    septaphase_slider = Slider(fig; range = 1:3, width = 48,
                               linewidth = 12, halign = :left)
    for i = n_btns
        grid1[i,1] = buttons[i]
    end
    fig[1,1] = grid!([Label(fig, "septaphase") septaphase_slider];
                     tellheight = false, valign = :bottom)
    grid2[1,1] = grid!([Label(fig, "phase contours") phase_contours_toggle])
    grid2[1,2] = grid!([Label(fig, "modulus contours") modulus_contours_toggle])
    on(septaphase_slider.value) do value
        setindex!.(btn_active, false)
        prev_img = img[3][] = color_matrices[value]
    end
    for i = n_btns
        on(buttons[i].clicks) do _
            btn = btn_active[i]
            setindex!.(btn_active[n_btns .!= i], false)
            btn[] = !btn[]
            img[3][] = btn[] ? color_matrices[i+3] : prev_img
        end
    end
    on(phase_contours_toggle.active) do active
        ϕ_const && return
        if active
            phase_contours = draw_phase_contours(axis, ϕ, colormap)
        else
            delete!(axis, phase_contours)
        end
        return
    end
    on(modulus_contours_toggle.active; update = true) do active
        r_const && return
        if active
            modulus_contours = draw_modulus_contours(axis, r, levels)
        else
            delete!(axis, modulus_contours)
        end
        return
    end
    ar::Bool = true
    onmouserightup(addmouseevents!(fig.scene)) do _
        ispressed(fig, Keyboard.left_control) || return
        ar = !ar
        aspect_control!(ar, axis, fig, xlen, ylen)
    end
    params = Coordinates(; x, y, r, ϕ, u, v)
    aspect_control!(ar, axis, fig, xlen, ylen)
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
    s = @. f(complex(x, y'))
    complex_plot(x, y, s, color; kw...)
end

function complex_plot(x::AbstractVector, y::AbstractVector, f, ::Type{Real},
                      color::Spaces = default; kw...)
    s = f.(x, y')
    complex_plot(x, y, s, color; kw...)
end

Δ(x::AbstractVector) = x[end] - x[begin]
Δ(x::Interval) = IntervalSets.width(x)

adjust_len(x, y, len) = round(Int, len * Δ(x) / Δ(y))

ranges(x::Interval, y::Interval) = (range(x, adjust_len(x, y, n)), range(y, n))
ranges(x::AbstractVector, y::Interval) = (x, range(y, adjust_len(y, x, length(x))))
ranges(x::Interval, y::AbstractVector) = (range(x, adjust_len(x, y, length(y))), y)

function complex_plot(x::VecOrInterval, y::VecOrInterval, f,
                      color::Spaces = default; kw...)
    x, y = ranges(x, y)
    complex_plot(x, y, f, color; kw...)
end

function complex_plot(x::VecOrInterval, y::VecOrInterval, f, ::Type{Real},
                      color::Spaces = default; kw...)
    x, y = ranges(x, y)
    complex_plot(x, y, f, Real, color; kw...)
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
    HSL   => map(RGB, HSL(i, 1.0, 0.5) for i = range(-60, 300, 2^10)),
    Oklch => map(RGB, Oklch(0.5, chroma, i) for i = range(-30, 330, 2^10))
)

end
