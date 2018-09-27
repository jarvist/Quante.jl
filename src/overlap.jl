## One-electron integrals
### Overlap matrix elements

function overlap(a::PGBF,b::PGBF)
    return a.norm*b.norm*overlap(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

overlap(a::CGBF,b::CGBF) = contract(overlap,a,b)

function overlap(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                 aI::Int64,aJ::Int64,aK::Int64, 
                 bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                 bI::Int64,bJ::Int64,bK::Int64)
    gamma = aexpn+bexpn
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    rab2 = dist2(ax-bx,ay-by,az-bz) 
    pre = (pi/gamma)^1.5*exp(-aexpn*bexpn*rab2/gamma)
    wx = overlap1d(aI,bI,px-ax,px-bx,gamma)
    wy = overlap1d(aJ,bJ,py-ay,py-by,gamma)
    wz = overlap1d(aK,bK,pz-az,pz-bz,gamma)
    return pre*wx*wy*wz
end

function gaussian_product_center(a::PGBF,b::PGBF)
    return (a.expn*[a.x,a.y,a.z]+b.expn*[b.x,b.y,b.z])/(a.expn+b.expn)
end

function gaussian_product_center(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                                    bexpn::Float64,bx::Float64,by::Float64,bz::Float64)
    return (aexpn*[ax,ay,az]+bexpn*[bx,by,bz])/(aexpn+bexpn)    
end

function overlap1d(la::Int64,lb::Int64,ax::Float64,bx::Float64,gamma::Float64)
    total = 0
    for i in 0:div(la+lb,2)
        total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
    end
    return total
end

function binomial_prefactor(s::Int64,ia::Int64,ib::Int64,xpa::Float64,xpb::Float64)
    #println("binomial_prefactor($s,$ia,$ib,$xpa,$xpb)")
    total = 0
    for t in 0:s
        if (s-ia) <= t <= ib
            total += binomial(ia,s-t)*binomial(ib,t)*xpa^(ia-s+t)*xpb^(ib-t)
        end
    end
    return total
end


function test_overlap()
    @testset "test_overlap" begin
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @test overlap1d(0,0,0.,0.,1.) == 1
    @test gaussian_product_center(s,s) == [0,0,0]
    @test isapprox(overlap(s,s),1)
    @test isapprox(overlap(px,px),1)
    @test isapprox(overlap(s,px),0)
    @test binomial_prefactor(0,0,0,0.,0.) == 1
    end #testset
end

