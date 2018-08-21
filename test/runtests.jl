push!(LOAD_PATH, "../src")

using Quante
if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "Quante.jl" begin

#include("C60.jl")
include("SzaboOstlundReference.jl") # further references from Szabo and Ostlund
#include("vrrspeedtest.jl") # speed test of two-centre integrals

function microtest() # collection of microscopic tests, within scope of functions
    Quante.test_utils()
    Quante.test_pgbf()
    Quante.test_cgbf()
    Quante.test_overlap()
    Quante.test_kinetic()
    Quante.test_a_terms()
    Quante.test_gamma()
    Quante.test_na()
    Quante.test_fgamma()
    Quante.test_one()
    Quante.test_na2()
    Quante.test_two_terms()
    Quante.test_coul1()
    Quante.test_vrr()
    Quante.test_hrr()
    Quante.test_geo_basis()
end

microtest()

# OK! Let's try the full blown restricted Hartree Fock...
function test_h2() # via PyQuante
    @time Energy, E, U = rhf(Quante.h2,verbose=true)
    @test Energy ≈ -1.1170996
end

function test_lih() # via PyQuante
    @time Energy, E, U = rhf(Quante.lih,verbose=true)
    @test Energy ≈ -7.860745582058085 #-7.86073270525799) <-- old PyQuante version
end

function test_h2o() # via PyQuante
    @time Energy,E,U = rhf(Quante.h2o)
    @test Energy ≈ -74.9598566070913 #-74.9597609118851) <-- old PyQuante version
end

println("Testing simple molecules, sto-3g energy...")
test_h2()
test_h2o()
test_lih()

end # testset
println("All tests passed! 8^)")

