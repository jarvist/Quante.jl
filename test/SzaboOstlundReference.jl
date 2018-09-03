# Further tests from Szabo and Ostlund, Table $3.10, Standard geoms used in calcs
# units are 'a.u.' ; so assume Bohr
CO = Quante.Molecule([  Atom(AtoZ["C"], 0, 0, 0),
                        Atom(AtoZ["O"], 2.132, 0, 0)])
N2 = Quante.Molecule([  Atom(AtoZ["N"], 0, 0, 0),
                        Atom(AtoZ["N"], 2.074, 0, 0)])

@testset "SzaboOstlundReference" begin

# Szabo and Ostlund
@time Energy,E,U = rhf(CO)
println("CO Energy: $Energy S&O STO-3G: -111.225")
@test Energy ≈ -111.225 atol=2e-3

@time Energy,E,U = rhf(N2)
println("N2 Energy: $Energy S&O STO-3G: -107.496")
#@test Energy ≈ -107.496 atol=2e-3

end #testset

