export rhf

function all_1e_ints(bfs::BasisSet,mol::Molecule)
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

function all_twoe_ints(bflist,ERI=coulomb)
    n = length(bflist.bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = Array{Float64}(totlen)
    for (i,j,k,l) in iiterator(n)
        ints2e[iindex(i,j,k,l)] = ERI(bflist.bfs[i],bflist.bfs[j],bflist.bfs[k],bflist.bfs[l])
    end
    return ints2e
end

function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
    n = size(D,1)
    G = Array{Float64}(n,n)
    D1 = reshape(D,n*n)
    temp = Array{Float64}(n*n)
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

dmat(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'


function rhf(mol::Molecule,MaxIter::Int64=8,verbose::Bool=false)
    bfs = build_basis(mol)
    S,T,V = all_1e_ints(bfs,mol)
    Ints = all_twoe_ints(bfs)
    h = T+V
    E,U = eig(h,S)
    Enuke = nuclear_repulsion(mol)
    nclosed,nopen = divrem(nel(mol),2)
    Eold = 0
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
    for iter in 1:MaxIter
        D = dmat(U,nclosed)
        if verbose
            println("D=\n$D")
        end
        G = make2JmK(D,Ints)
        H = h+G
        E,U = eig(H,S)
        Eone = trace2(D,h)
        Etwo = trace2(D,H)
        Energy = Enuke + Eone + Etwo
        println("HF: $iter  $Energy : $Enuke    $Eone    $Etwo")
        if isapprox(Energy,Eold)
            break
        end
        Eold  = Energy
    end
    return Energy,E,U
end

