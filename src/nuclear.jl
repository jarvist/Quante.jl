"""
### Nuclear attraction term
"""
function Aterm(i::Int64,r::Int64,u::Int64,l1::Int64,l2::Int64,ax::Float64,bx::Float64,cx::Float64,gamma::Float64)
    term1 = (-1)^(i+u)*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1::Int64,l2::Int64,a::Float64,b::Float64,c::Float64,g::Float64)
    Imax = l1+l2+1
    A = zeros(Float64,Imax)
    for i in 0:(Imax-1)
        for r in 0:div(i,2)
            for u in 0:div(i-2r,2)
                I = i-2r-u+1
                A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
            end
        end
    end
    return A
end

function nuclear_attraction(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                            aI::Int64,aJ::Int64,aK::Int64,
                            bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                            bI::Int64,bJ::Int64,bK::Int64,
                            cx::Float64,cy::Float64,cz::Float64)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    gamma = aexpn+bexpn
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcp2 = dist2(cx-px,cy-py,cz-pz)
    Ax = Aarray(aI,bI,px-ax,px-bx,px-cx,gamma)
    Ay = Aarray(aJ,bJ,py-ay,py-by,py-cy,gamma)
    Az = Aarray(aK,bK,pz-az,pz-bz,pz-cz,gamma)
    total = 0.0
    for I in 0:(aI+bI)
        for J in 0:(aJ+bJ)
            for K in 0:(aK+bK)
                total += Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
            end
        end
    end
    val=-2pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
    #println(val)
    #println((Ax,Ay,Az,rcp2*gamma,Fgamma(0,rcp2*gamma)))
    return val
end

function nuclear_attraction(a::PGBF,b::PGBF,cx::Float64,cy::Float64,cz::Float64)
    return a.norm*b.norm*nuclear_attraction(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                            b.expn,b.x,b.y,b.z,b.I,b.J,b.K,cx,cy,cz)
end
nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.x,c.y,c.z)
nuclear_attraction(a::PGBF,b::PGBF,m::Molecule) = sum([nuclear_attraction(a,b,c) for c in m.atomlist])

# Need a nested scope to squeeze this into the contract function
function nuclear_attraction(a::CGBF,b::CGBF,cx::Float64,cy::Float64,cz::Float64)
    na(a,b) = nuclear_attraction(a,b,cx,cy,cz)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
    na(a,b) = nuclear_attraction(a,b,c)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,m::Molecule)
    na(a,b) = nuclear_attraction(a,b,m)
    contract(na,a,b)
end

function test_a_terms()
    @testset "a_terms" begin
    @test Aterm(0,0,0,0,0,0.,0.,0.,0.) == 1.0
    @test Aarray(0,0,0.,0.,0.,1.) == [1.0]
    @test Aarray(0,1,1.,1.,1.,1.) == [1.0, -1.0]
    @test Aarray(1,1,1.,1.,1.,1.) == [1.5, -2.5, 1.0]
    @test Aterm(0,0,0,0,0,0.,0.,0.,1.) == 1.0
    @test Aterm(0,0,0,0,1,1.,1.,1.,1.) == 1.0
    @test Aterm(1,0,0,0,1,1.,1.,1.,1.) == -1.0
    @test Aterm(0,0,0,1,1,1.,1.,1.,1.) == 1.0
    @test Aterm(1,0,0,1,1,1.,1.,1.,1.) == -2.0
    @test Aterm(2,0,0,1,1,1.,1.,1.,1.) == 1.0
    @test Aterm(2,0,1,1,1,1.,1.,1.,1.) == -0.5
    @test Aterm(2,1,0,1,1,1.,1.,1.,1.) == 0.5
    end #testset
end

function test_na()
    @testset "test_na" begin
    s = pgbf(1.0)
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @test amplitude(c,0,0,0) ≈ 0.71270547
    @test nuclear_attraction(s,s,0.,0.,0.) ≈ -1.59576912
    @test nuclear_attraction(c,c,0.,0.,0.) ≈ -1.59576912
    end #testset
end
#todo make into a test
function test_one()
    s1 = pgbf(1)
    s2 = pgbf(1,0,1,0)
    x=y=0.
    println("S: $(overlap(s1,s2))")
    println("T: $(kinetic(s1,s2))")
    for z in range(0,stop=1,length=5)
        println("V: $z $(nuclear_attraction(s1,s2,x,y,z))")
    end
end


function test_na2()
    @testset "test_na2" begin
    li,h = lih.atomlist
    bfs = build_basis(lih)
    s1,s2,x,y,z,h1s = bfs.bfs
    @test nuclear_attraction(s1,s1,lih) ≈ -8.307532656
    end #testset
end

