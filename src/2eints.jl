## Two electron integrals

function coulomb(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                 aI::Int64,aJ::Int64,aK::Int64,
                 bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                 bI::Int64,bJ::Int64,bK::Int64,
                 cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,
                 cI::Int64,cJ::Int64,cK::Int64,
                 dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,
                 dI::Int64,dJ::Int64,dK::Int64)
    # This is the slow method of computing integrals from Huzinaga et al.
    # Use the HRR/VRR scheme from Head-Gordon & Pople instead

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)

    g1 = aexpn+bexpn
    g2 = cexpn+dexpn
    delta = 0.25*(1/g1+1/g2)
    
    Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
    By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
    Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)
    
    s = 0
    #println("$(aI+bI+cI+dI),$(aJ+bJ+cJ+dJ),$(aK+bK+cK+dK)")
    for I in 0:(aI+bI+cI+dI)
        for J in 0:(aJ+bJ+cJ+dJ)
            for K in 0:(aK+bK+cK+dK)
                #println("coul: $I,$J,$K,$(Bx[I+1]),$(By[J+1]),$(Bz[K+1])")
                s += Bx[I+1]*By[J+1]*Bz[K+1]*Fgamma(I+J+K,0.25*rpq2/delta)
            end
        end
    end
    return 2pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-aexpn*bexpn*rab2/g1)*exp(-cexpn*dexpn*rcd2/g2)*s
end

function coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*coulomb(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
        b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
        c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
        d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end

fB(i::Int64,l1::Int64,l2::Int64,p::Float64,a::Float64,b::Float64,r::Int64,g::Float64) = binomial_prefactor(i,l1,l2,p-a,p-b)*B0(i,r,g)
B0(i::Int64,r::Int64,g::Float64) = fact_ratio2(i,r)*(4g)^(r-i)
fact_ratio2(a::Int64,b::Int64) = factorial(a)factorial(b)/factorial(a-2b)

function Bterm(i1::Int64,i2::Int64,r1::Int64,r2::Int64,u::Int64,
               l1::Int64,l2::Int64,l3::Int64,l4::Int64,
               Px::Float64,Ax::Float64,Bx::Float64,Qx::Float64,Cx::Float64,Dx::Float64,
               gamma1::Float64,gamma2::Float64,delta::Float64)
    # THO eq. 2.22
    #print("Bterm($i1,$i2,$r1,$r2,$u,$l1,$l2,$l3,$l4,$Px,$Ax,$Bx,$Qx,$Cx,$Dx,$gamma1,$gamma2,$delta)=")
    val = (-1)^(i2+u)*fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)*(
            fact_ratio2(i1+i2-2*(r1+r2),u)*(Qx-Px)^(i1+i2-2*(r1+r2)-2*u)/delta^(i1+i2-2*(r1+r2)-u))
    #println("$val")
    return val
end

function Barray(l1::Int64,l2::Int64,l3::Int64,l4::Int64,p::Float64,a::Float64,b::Float64,
                q::Float64,c::Float64,d::Float64,g1::Float64,g2::Float64,delta::Float64)
    Imax = l1+l2+l3+l4+1
    B = zeros(Float64,Imax)
    for i1 in 0:(l1+l2)
        for i2 in 0:(l3+l4)
            for r1 in 0:div(i1,2)
                for r2 in 0:div(i2,2)
                    for u in 0:(div(i1+i2,2)-r1-r2)
                        I = i1+i2-2*(r1+r2)-u
                        B[I+1] += Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta)
                    end
                end
            end
        end
    end
    return B
end

coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb,a,b,c,d)


function test_two_terms()
    @testset "test_two_terms" begin
    @test fB(0,0,0,0.0,0.0,0.0,0,2.0) == 1
    @test fB(0,0,0,1.0,1.0,1.0,0,2.0) == 1
    @test fB(0,0,0,0.0,0.0,0.0,0,2.0 ) == 1
    @test fB(1,0,1,0.0,0.0,0.0,0,2.0 ) == 0.125
    @test B0(0,0,2.0) == 1
    @test fact_ratio2(0,0) == 1
    @test Bterm(0,0,0,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==1
    @test Bterm(0,1,0,0,0,0,0,0,1,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==0
    end #testset
end


function test_coul1()
    @testset "test_coul1" begin
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @test coulomb(s,s,s,px)==0 # 0
    @test coulomb(s,s,px,px) ≈ 0.9403159725793305 
    end #testset
end

function coulomb_hgp(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*hrr(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
        b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
        c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
        d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end
coulomb_hgp(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb_hgp,a,b,c,d)

