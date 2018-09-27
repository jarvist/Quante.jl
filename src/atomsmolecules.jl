# ## Atoms and Molecules

export AtoBohr!, Atom, Molecule, AtoZ

mutable struct Atom
    atno::Int64
    x::Float64
    y::Float64
    z::Float64
end

mutable struct Molecule
    atomlist::Array{Atom,1}
end

function push!(mol::Molecule,at::Atom)
    Base.push!(atomlist,at)
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
# First we need to dereference atom name to atomic charge (Z)
AtoZ=Dict(j=>i for (i,j) in enumerate(["H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La",
    "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]))
# Magic comprehension ends up like:
# ZfromAtom = Dict("H" => 1, "He" => 2, "Li" => 3)


# Sample molecules for tests - all from PyQuante
# All units currently Angstrom; converted into Bohr below
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
AtoBohr(x::Float64) = x*BohrInA # These assume units start as Angstrom.
function AtoBohr!(at::Atom)
    at.x *= BohrInA
    at.y *= BohrInA 
    at.z *= BohrInA
end
function AtoBohr!(mol::Molecule)
    for at in mol.atomlist
        AtoBohr!(at)
    end
end

AtoBohr!(h2)
AtoBohr!(h2o)
AtoBohr!(ch4)
AtoBohr!(c6h6)
AtoBohr!(lih)

