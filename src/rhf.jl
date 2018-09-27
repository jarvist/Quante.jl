export rhf

"""
1e integrals; single core version. 
"""
function all_1e_ints_singlecore(bfs::BasisSet,mol::Molecule)
    n = length(bfs.bfs)
    S = Array{Float64}(n,n)
    T = Array{Float64}(n,n)
    V = Array{Float64}(n,n)
    for (i,j) in pairs(n)
        a,b = bfs.bfs[i],bfs.bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
    end
    return S,T,V
end

function all_1e_ints(bfs::BasisSet,mol::Molecule)
    n = length(bfs.bfs)
    S = SharedArray{Float64}(n,n)
    T = SharedArray{Float64}(n,n)
    V = SharedArray{Float64}(n,n)
    @sync @distributed for i = 1:n 
        @sync @distributed for j = 1:i 
    #for (i,j) in pairs(n)
        a,b = bfs.bfs[i],bfs.bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
        end
    end
    return S,T,V
end

"""
2e integrals; single core version.
"""
function all_2e_ints_singlecore(bflist,ERI=coulomb_hgp) 
    n = length(bflist.bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = SharedArray{Float64}(totlen)
    for (i,j,k,l) in iiterator(n)
        ints2e[iindex(i,j,k,l)] = ERI(bflist.bfs[i],bflist.bfs[j],bflist.bfs[k],bflist.bfs[l])
    end
    return ints2e
end

function all_2e_ints(bflist,ERI=coulomb_hgp;verbose::Bool=false) #coloumb
    n = length(bflist.bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = SharedArray{Float64}(totlen)
    #for (i,j,k,l) in iiterator(n)
    if verbose println(n^4/4," ints to calculate. Each '.' is ",n^2/2) end
    @sync @distributed for i=1:n # Nb: need sync on both these parallel, to force completion of Ints before SCF
        @sync @distributed for j=1:i 
            for (k,l) in pairs(n) 
                if triangle(i,j) <= triangle(k,l)
                    ints2e[iindex(i,j,k,l)] = ERI(bflist.bfs[i],bflist.bfs[j],bflist.bfs[k],bflist.bfs[l])
                end
            end
            if verbose print(".") end
        end 
    end
    return ints2e
end

"""
 make2JmK

 Coulomb and exchange contribution to the Fock matrix. 
 The N^4 component of Hartree Fock.

 TODO: Optimise me :^)
 Nb: Now accepts SharedArray (parallel constructed) Ints 
"""
function make2JmK(D::Array{Float64,2},Ints::SharedArray{Float64,1})
    n = size(D,1)
    G = Array{Float64}(undef,n,n) # undef is for 1.0 compat. Unsure what it does.
    D1 = reshape(D,n*n)
    temp = Array{Float64}(undef,n*n)
    for (i,j) in pairs(n)
        kl = 1
        for (k,l) in rpairs(n)
            temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
            kl += 1
        end
        G[i,j] = G[j,i] = dot(D1,temp)
    end
    return G
end

densitymatrix(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'

"""
 rhf

Restricted Hartree-Fock function. Will run a self-consistent field (SCF)
calculation on the supplied molecular geometry.

Algorithm (and comments herein) follow Thijssen 2007 2nd Ed pp70-73

"""
function rhf(mol::Molecule,MaxIter::Int64=40; verbose::Bool=false, Econvergence::Float64=1e-6)
    if verbose println("Building basis") end
    bfs = build_basis(mol)
    if verbose println(length(bfs.bfs)," basis functions") end
    # S = Overlap matrix
    # T = one-electron Kinetic Energy
    # V = one-electron potential
    if verbose println("Calculating 1e ints") end
    S,T,V = all_1e_ints(bfs,mol)
    # Ints = two-electron integrals, <pr|g|qs>
    if verbose println("Calculating 2e ints") end
    Ints = all_2e_ints(bfs,verbose=verbose)
    # h = Form one-eletron Hamiltonian
    h = T+V
    # generalised eigenvalue decomposition of one-electron hamiltonian and overlap matrix 
    # used as starting guess for self-consistent procedure?
    E,U = eigen(h,S)
    # Enuke = (classical) nuclear repulision
    Enuke = nuclear_repulsion(mol)
    # Define occupied and virtual orbitals
    nclosed,nopen = divrem(nel(mol),2)
    D=densitymatrix(U,nclosed) # initial density matrix
    Eold = Inf 
    Energy = 0
    println("Nel=$(nel(mol)) Nclosed=$nclosed")
    if verbose
        println("S=\n$S")
        println("h=\n$h")
        println("T=\n$T")
        println("V=\n$V")
        println("E: $E")
        println("U: $U")
        println("2e ints:\n$Ints")
    end
    println("SCF: Iteration TotalEnergy :=: Enuke + Eone + Etwo")
    for iter in 1:MaxIter
        α=0.4 # mixing of current and prior density matrices
        D = α*D + (1-α)*densitymatrix(U,nclosed) # density matrix; from eigenvalues U
        if verbose
            println("D=\n$D")
        end
        # Coulomb and Exchange contributions to Fock Matrix:
        # This is the N^4, and most time consuming, step of Hartree Fock.
        G = make2JmK(D,Ints)
        # Fock matrix H = one-electron hamiltonian + G matrix
        H = h+G
        # Diagonalise the Fock Matrix ; eig=Generalised eigenvalue decomposition
        E,U = eigen(H,S)
        # expanded versions of trace2 functions
        Eone = sum(D.*h) # one electron contribution
        Etwo = sum(D.*H) # should be = sum of occupied Fock orbitals
        Energy = Enuke + Eone + Etwo #c.f. Thijssen, factor of 2 already taken care of?
        println("HF: $iter  $Energy : $Enuke    $Eone    $Etwo")
        if abs(Energy-Eold)<Econvergence
            break
        end
        Eold  = Energy
    end
    return Energy,E,U
end

