using JSON3
using LinearAlgebra
using Optim

using Polynomials, SpecialPolynomials
import Polynomials: basis
using KrylovKit
using ProgressMeter

dir=joinpath("common/")
include.(joinpath.(Ref(dir),readdir(dir)[contains.(readdir(dir), ".jl")]))

