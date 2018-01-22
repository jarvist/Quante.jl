# Quante.jl

[![Build Status](https://travis-ci.org/jarvist/Quante.jl.svg?branch=master)](https://travis-ci.org/jarvist/Quante.jl)
[![Coverage
Status](https://coveralls.io/repos/jarvist/Quante.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jarvist/Quante.jl?branch=master)
[![codecov.io](http://codecov.io/github/jarvist/Quante.jl/coverage.svg?branch=master)](http://codecov.io/github/jarvist/Quante.jl?branch=master)

An open-source implementation of quantum-chemistry methods in the Julia
programming language. 
These codes are based on Rick Muller's pyquante2 
(https://github.com/rpmuller/pyquante2/). 

Julia is a high-level dynamic programming language that is capable of fast
execution. In particular it can generate performant code targetting SIMD,
parallelisation and GPU execution, from high level constructions.

### ToDo

- [x] Debug H2O / LiF RHF energy error
- [ ] Add CCSD 'explicit forloops' - see Psi4Numpy
  https://github.com/psi4/psi4numpy/blob/master/Coupled-Cluster/CCSD.py
- [ ] Add CCSD 'Tensor contractions' - via
  https://github.com/Jutho/TensorOperations.jl
- [ ] Play with C_∞ / D_∞ symmetry for diatomic CCSD interaction energies

### References

- Attila Szabo and Neil S. Ostlund - Modern Quantum Chemistry
(https://www.amazon.co.uk/Modern-Quantum-Chemistry-Introduction-Electronic/dp/0486691861)
- David B. Cook - Handbook of Computational Quantum Chemistry
(https://www.amazon.co.uk/Handbook-Computational-Quantum-Chemistry-Dover/dp/0486443078)
- PyQuante2 (https://github.com/rpmuller/pyquante2)

