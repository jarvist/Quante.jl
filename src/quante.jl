"""
# Quante.jl
A fork of PqQuante in Julia; experimenting with quantum chemistry in the Julia language.

"""

module quante

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

end

