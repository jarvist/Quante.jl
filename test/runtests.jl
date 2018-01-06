push!(LOAD_PATH, "../src")

using quante

function microtest() # collection of microscopic tests, within scope of functions
    quante.test_utils()
    quante.test_pgbf()
    quante.test_cgbf()
    quante.test_overlap()
    quante.test_kinetic()
    quante.test_a_terms()
    quante.test_gamma()
    quante.test_na()
    quante.test_fgamma()
    quante.test_one()
    quante.test_na2()
    quante.test_two_terms()
    quante.test_coul1()
    quante.test_vrr()
    quante.test_hrr()
    quante.test_geo_basis()
end

microtest()

# OK! Let's try the full blown restricted Hartree Fock...
function test_h2()
    @time Energy, E, U = rhf(quante.h2)
    @assert isapprox(Energy,-1.1170996)
end

function test_lih()
    @time Energy, E, U = rhf(quante.lih)
    @assert isapprox(Energy,-7.86073270525799)
end

function test_h2o()
    @time Energy,E,U = rhf(quante.h2o)
    @assert isapprox(Energy,-74.9597609118851)
end

test_h2()
test_lih()
test_h2o()

