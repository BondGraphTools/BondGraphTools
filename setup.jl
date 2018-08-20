ENV["PYTHON"] = "python3"

Pkg.update()
Pkg.add("PyCall")
Pkg.add("DifferentialEquations")

using PyCall
using DifferentialEquations