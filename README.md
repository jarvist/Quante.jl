# Quante.jl

An open-source implementation of quantum-chemistry methods in the Julia
programming language. 
These codes are based on Rick Muller's pyquante2 
(https://github.com/rpmuller/pyquante2/). 

Julia is a high-level dynamic programming language that is capable of fast
execution. In particular it can generate performant code targetting SIMD,
parallelisation and GPU execution, from high level constructions.

### ToDo

- [ ] Debug H2O / LiF RHF energy error
- [ ] Add CCSD 'explicit forloops' - see Psi4Numpy
  https://github.com/psi4/psi4numpy/blob/master/Coupled-Cluster/CCSD.py
- [ ] Add CCSD 'Tensor contractions' - via
  https://github.com/Jutho/TensorOperations.jl
- [ ] Play with C_∞ / D_∞ symmetry for diatomic CCSD interaction energies

