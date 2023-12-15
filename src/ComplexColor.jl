module ComplexColor

export complex_color, complex_plot, var"@L_str"

using Printf
using Colors
using GLMakie
using .Makie: latexstring

const ComplexArray = AbstractArray{<:Complex{<:Real}}
const RealArray = AbstractArray{<:Real}

struct Septaphase end

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
complex_color(s::ComplexArray) = HSL_to_RGB(complex_to_HSL(s)...)

function complex_color(s::ComplexArray, ::Septaphase)
    gradphase, S, L = complex_to_HSL(s)
    [HSL_to_RGB(H, S, L) for H ∈ (gradphase, septaphase(gradphase)...)]
end

function complex_to_HSL(s::ComplexArray)
    r2 = abs2.(s)
    H = @. mod(rad2deg(angle(s)) + 120, 360)
    L = @. r2 / (r2 + 1)
    S = ones(eltype(H), size(H))
    return H, S, L
end

HSL_to_RGB(H, S, L) = map(RGB, HSL.(H, S, L)) |> clamp01nan1!

function septaphase(H)
        rounded = map(hue -> 60 * div(hue, 60, RoundNearest), H)
    thresholded = map(hue -> 60 * fld(hue, 60), H)
    return rounded, thresholded
end

function draw_modulus_contours(axis, r)
    levels = map(n -> ldexp(1.0, n), -3:8)
    colormap = Reverse(:acton)
    contour!(axis, r; levels, colormap, inspectable = false)
end

function draw_phase_contours(axis, ϕ)
    contour!(axis, ϕ; levels = -180:90:180, colormap = hsl, inspectable = false)
end

"""
    complex_plot(x, y, s; [title])

Plot a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL color space.

The three-valued `septaphase` slider option on the plot will partition the phase using only 6 colors (green, cyan, blue, magenta, red, yellow) either by thresholding or by rounding, depending on the setting; the default off position plots a gradient phase.
"""
function complex_plot(x::AbstractVector, y::AbstractVector, s::ComplexArray;
                      title::AbstractString = L"s")
    xlen = length(x)
    ylen = length(y)
    (xlen, ylen) == size(s) || error("Length mismatch.") 
    nticks = 5
    xticklabels = latexstring.(range(x[begin], x[end], nticks))
    yticklabels = latexstring.(range(y[begin], y[end], nticks))
    xticks = (range(1, xlen, nticks), xticklabels)
    yticks = (range(1, ylen, nticks), yticklabels)
    arg_ticks = (-π:π:π, [L"-\pi", L"0", L"\pi"])
    fig = Figure(size = (600, 532))
    axis = Axis(fig[1,1]; title, titlesize = 21,
                xlabel = L"Re(z)", xlabelsize = 16,
                ylabel = L"Im(z)", ylabelsize = 16,
                xticks, yticks, aspect = DataAspect())
    r = abs.(s)
    ϕ = @. rad2deg(angle(s))
    function inspector(_, inds, _)
        i, j = round.(Int, inds)
        str = @sprintf(
            """
            x: %.2f, y: %.2f
            r: %.2f, ϕ: %.2f\u00b0
            """,
            x[i], y[j], r[i,j], ϕ[i,j]
        )
        replace(str, '-' => '\u2212')
    end
    color_matrices = complex_color(s, Septaphase())
    img = image!(axis, color_matrices[1]; inspector_label = inspector)
    local phase_contours
    modulus_contours = draw_modulus_contours(axis, r)
    Colorbar(fig[1,2]; colormap = hsl, limits = (-π, π), ticks = arg_ticks,
             label = "Arg(s)")
    grid = GridLayout(fig[2,:])
    modulus_toggle = Toggle(fig; active = true)
    phase_toggle = Toggle(fig; active = false)
    septaphase_slider = Slider(fig; range = 1:3, width = phase_toggle.width[] * 1.5)
    grid[1,1] = grid!([Label(fig, "septaphase") septaphase_slider])
    grid[1,2] = grid!([Label(fig, "phase contours") phase_toggle])
    grid[1,3] = grid!([Label(fig, "modulus contours") modulus_toggle])
    on(septaphase_slider.value) do value
        img[3][] = color_matrices[value]
    end
    on(phase_toggle.active) do active
        if active
            phase_contours = draw_phase_contours(axis, ϕ)
        else    
            delete!(axis, phase_contours)
        end
    end
    on(modulus_toggle.active) do active
        if active
            modulus_contours = draw_modulus_contours(axis, r)
        else    
            delete!(axis, modulus_contours)
        end
    end
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    DataInspector(fig)
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
