# ComplexColor.jl

`complex_color(s)` converts an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

`complex_color(r, ϕ)` accepts modulus and phase arrays as input.

`complex_plot(x, y, s)` plots a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL color space.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Sagnac/ComplexColor.jl")
```

## Plotting example

```julia
using ComplexColor
# construct the complex array
x = y = -2:0.01:2
z = complex.(x, y')
s = map(z -> 1/z, z)
# transform the complex array into an image
fig = complex_plot(x, y, s; title = L"s = z^{-1}")
```
