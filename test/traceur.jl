# traceur.jl - automatic performance tips from traceur.jl

push!(LOAD_PATH, "../src")
using Traceur

using Revise # to catch edits to Quante
using Quante

mol=Quante.h2 

bfs = Quante.build_basis(mol)

S,T,V = Quante.all_1e_ints_singlecore(bfs,mol)

@trace  ints=Quante.all_2e_ints_singlecore(bfs)

