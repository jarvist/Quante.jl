# ## Atoms and Molecules

type Atom
    atno::Int64
    x::Float64
    y::Float64
    z::Float64
end

type Molecule
    atomlist::Array{Atom,1}
end

function push!(mol::Molecule,at::Atom)
    Base.push!(atomlist,at)
end

tobohr(x::Float64) = x/0.52918
function tobohr!(at::Atom)
    at.x /= 0.52918
    at.y /= 0.52918
    at.z /= 0.52918
end
function tobohr!(mol::Molecule)
    for at in mol.atomlist
        tobohr!(at)
    end
end

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Molecule)
    nr = 0
    for (i,j) in spairs(nat(mol))
        nr += nuclear_repulsion(mol.atomlist[i],mol.atomlist[j])
    end
    return nr
end

nel(mol::Molecule) = sum([at.atno for at in mol.atomlist])
nat(mol::Molecule) = length(mol.atomlist)

# Other molecule methods to implement
# nocc, nclosed, nopen, nup, ndown, stoich, mass,
# center_of_mass, center!

# Array of symbols, masses


# Sample molecules for tests
h2 = Molecule([Atom(1,  0.00000000,     0.00000000,     0.36628549),
               Atom(1,  0.00000000,     0.00000000,    -0.36628549)])

h2o = Molecule([Atom(8,   0.00000000,     0.00000000,     0.04851804),
                Atom(1,   0.75300223,     0.00000000,    -0.51923377),
                Atom(1,  -0.75300223,     0.00000000,    -0.51923377)])

ch4 = Molecule([Atom(6,   0.00000000,     0.00000000,     0.00000000),
                Atom(1,   0.62558332,    -0.62558332,     0.62558332),
                Atom(1,  -0.62558332,     0.62558332,     0.62558332),
                Atom(1,   0.62558332,     0.62558332,    -0.62558332),
                Atom(1,  -0.62558332,    -0.62558332,    -0.62558332)])

c6h6 = Molecule([ Atom(6,  0.98735329,     0.98735329,     0.00000000),
                  Atom(6,  1.34874967,    -0.36139639,     0.00000000),
                  Atom(6,  0.36139639,    -1.34874967,     0.00000000),
                  Atom(6, -0.98735329,    -0.98735329,     0.00000000),
                  Atom(6, -1.34874967,     0.36139639,     0.00000000),
                  Atom(6, -0.36139639,     1.34874967,     0.00000000),
                  Atom(1,  1.75551741,     1.75551741,     0.00000000),
                  Atom(1,  2.39808138,    -0.64256397,     0.00000000),
                  Atom(1,  0.64256397,    -2.39808138,     0.00000000),
                  Atom(1, -1.75551741,    -1.75551741,     0.00000000),
                  Atom(1, -2.39808138,     0.64256397,     0.00000000),
                  Atom(1, -0.64256397,     2.39808138,     0.00000000)])

lih = Molecule([Atom(3,    0.00000000,     0.00000000,    -0.53999756),
                Atom(1,    0.00000000,     0.00000000,     1.08999756)])

# Convert to atomic units (bohr)
tobohr!(h2)
tobohr!(h2o)
tobohr!(ch4)
tobohr!(c6h6)
tobohr!(lih)

