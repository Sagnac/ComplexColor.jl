# ComplexColor.jl

`complex_color(A)` converts an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

# Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Sagnac/ComplexColor.jl")
```

# Plotting example

```julia
using ComplexColor
using GLMakie # assumes previous installation of this plotting package

# construct the complex array
lims = -2:0.01:2
z = [complex(x, y) for x in lims, y in lims]
s = map(z -> 1/z, z)
# transform the complex array into an image
color_matrix = complex_color(s)
# plot
nticks = 5
ticklabels = Makie.latexstring.(range(lims[begin], lims[end], nticks))
xticks = yticks = (range(1, length(lims), nticks), ticklabels)
arg_ticks = (-π:π/2:π, [L"-\pi", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi"])
fig = Figure()
axis = Axis(fig[1,1]; title = L"s = z^{-1}", titlesize = 21,
            xlabel = L"Re(z)", xlabelsize = 16, ylabel = L"Im(z)", ylabelsize = 16,
            xticks, yticks)
plot = image!(axis, color_matrix)
grid = GridLayout(fig[1,2]; tellheight = false)
Colorbar(grid[1,1]; colormap = :gray1, ticks = (0:1, [L"0", L"\infty"]))
Label(grid[2,1], L"|s|"; fontsize = 15)
Colorbar(grid[1,2]; colormap = hsl, limits = (-π, π), ticks = arg_ticks)
Label(grid[2,2], L"Arg(s)"; fontsize = 15)
fig
```
