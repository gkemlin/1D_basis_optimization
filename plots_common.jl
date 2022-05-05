using PGFPlots
using LaTeXStrings
using JSON3
using LinearAlgebra
using DataFrames

pushPGFPlotsPreamble("\\usepackage{amsmath}")

# colors for plots
define_color("myred",    [228,26,28])
define_color("myblue",   [55,124,184])
define_color("mygreen",  [77,175,74])
define_color("myorange", [255,127,14])
define_color("myviolet", [148,103,189])
mycolors=["myred", "myblue", "mygreen", "myorange", "myviolet"]
linestyles=["dashed", "dash dot", "dash dot dot"]
