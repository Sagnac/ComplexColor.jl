# ComplexColor.jl

<img src="./images/complex_plot.png" width="556">

`complex_color(s, color = Oklch)` converts an array of complex numbers into an image matrix of RGB values using a hue-lightness color mapping for the phase and modulus.

`complex_plot(x, y, s, color = Oklch)` plots a complex number array `s` within the `x` and `y` limits using domain coloring in the HSL or perceptual OKLCH color spaces.

***Note***: The phase contour lines are kinda bugged at the moment at around `±π`; these contours are not displayed by default, but can be toggled on the plot.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Sagnac/ComplexColor.jl")
```

## Plotting examples

```julia
using ComplexColor
# for a complex array
x = y = -2:0.01:2
z = complex.(x, y')
s = map(z -> 1/z, z)
fig = complex_plot(x, y, s; title = L"s = z^{-1}")
```

```julia
# for a complex function over intervals
x = y = -2..2
f(z) = inv(z)
complex_plot(x, y, f)
```

```julia
# for an arbitrary complex function of two real variables (ℝ² → ℂ)
x = y = -3:0.01:3
f(x, y) = ei2pi(x * y)
complex_plot(x, y, f, Real; xlabel = "x", ylabel = "y")
```
