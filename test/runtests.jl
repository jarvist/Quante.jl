push!(LOAD_PATH, "../src")

using Quante

include("SzaboOstlundReference.jl") # further references from Szabo and Ostlund
#include("vrrspeedtest.jl") # speed test of two-matrix integrals

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
    @assert isapprox(Energy,-1.1170996)
end

function test_lih() # via PyQuante
    @time Energy, E, U = rhf(Quante.lih,verbose=true)
    @assert isapprox(Energy,-7.86073270525799)
end

function test_h2o() # via PyQuante
    @time Energy,E,U = rhf(Quante.h2o)
    @assert isapprox(Energy,-74.9597609118851)
end

test_h2()
test_h2o()
test_lih()
