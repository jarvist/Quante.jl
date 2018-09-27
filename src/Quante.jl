# Quante.jl
## A fork of PqQuante in Julia; experimenting with quantum chemistry in the Julia language.
module Quante

# Pull in necessary STDLIB stuff for Julia 0.7+
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
	using Test # so we can have tests distributed in source
	using Distributed
	using SharedArrays
	using LinearAlgebra
	using SpecialFunctions
end

include("physicalconstants.jl")

include("utility.jl")
include("gamma.jl")
include("basis.jl")
include("atomsmolecules.jl")

include("basissets.jl")

include("overlap.jl")
include("nuclear.jl")
include("kinetic.jl")

include("2eints.jl")
include("hrr.jl")

include("rhf.jl")

end # module

