# ComplexColor.jl

<img src="./images/complex_plot.png" width="556">

`complex_color(s, color = HSL)` converts an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

`complex_plot(x, y, s, color = HSL)` plots a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL or perceptual OKLCH color spaces.

***Note***: The phase contour lines are kinda bugged at the moment at around `±π`; these contours are not displayed by default, but can be toggled on the plot.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Sagnac/ComplexColor.jl")
```

## Plotting examples

```julia
using ComplexColor
# construct the complex array
x = y = -2:0.01:2
z = complex.(x, y')
s = map(z -> 1/z, z)
# transform the complex array into an image
fig = complex_plot(x, y, s; title = L"s = z^{-1}")

# convenient methods if you don't want to explicitly construct an array

# for a complex function f
f(z) = inv(z)
complex_plot(x, y, f)

# with intervals
x = y = -2..2
complex_plot(x, y, f)
```
