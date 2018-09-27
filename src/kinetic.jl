### Kinetic matrix elements

function kinetic(a::PGBF,b::PGBF)
    return a.norm*b.norm*kinetic(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

function kinetic(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
            bexpn::Float64,bx::Float64,by::Float64,bz::Float64,bI::Int64,bJ::Int64,bK::Int64)
    overlap0 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    overlapx1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI+2,bJ,bK)
    overlapy1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ+2,bK)
    overlapz1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK+2)
    overlapx2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-2,bJ,bK)
    overlapy2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-2,bK)
    overlapz2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-2)
    term0 = bexpn*(2*(bI+bJ+bK)+3)*overlap0
    term1 = -2*(bexpn^2)*(overlapx1+overlapy1+overlapz1)
    term2 = -0.5*(bI*(bI-1)*overlapx2+bJ*(bJ-1)*overlapy2+bK*(bK-1)*overlapz2)
    return term0+term1+term2
end

kinetic(a::CGBF,b::CGBF) = contract(kinetic,a,b)


function test_kinetic()
    @testset "test_kinetic" begin
    s = pgbf(1.0)
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @test amplitude(c,0,0,0) ≈ 0.71270547
    @test kinetic(1.,0.,0.,0.,0,0,0,1.,0.,0.,0.,0,0,0) ≈ 2.9530518648229536
    @test kinetic(s,s) ≈ 1.5
    @test kinetic(c,c) ≈ 1.5
    end #testset
end

