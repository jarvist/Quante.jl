function Fgamma(m::Int64,x::Float64,SMALL::Float64=1e-12)
    #println("Fgamma($m,$x)")
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

function gammainc(a::Float64,x::Float64)
    # This is the series version of gamma from pyquante. For reasons I don't get, it 
    # doesn't work around a=1. This works alright, but is only a stopgap solution
    # until Julia gets an incomplete gamma function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    if x < (a+1.0)
        #Use the series representation
        gam,gln = gser(a,x)
    else 
        #Use continued fractions
        gamc,gln = gcf(a,x)
        gam = 1-gamc
    end
    return exp(gln)*gam
end

function gser(a::Float64,x::Float64,ITMAX::Int64=100,EPS::Float64=3e-9)
    # Series representation of Gamma. NumRec sect 6.1.
    gln=lgamma(a)
    if x == 0
        return 0,gln
    end
    ap = a
    delt = s = 1/a
    for i in 1:ITMAX
        ap += 1
        delt *= (x/ap)
        s += delt
        if abs(delt) < abs(s)*EPS
            break
        end
    end
    return s*exp(-x+a*log(x)-gln),gln
end

function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
    #Continued fraction representation of Gamma. NumRec sect 6.1"
    gln=lgamma(a)
    b=x+1.0-a
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in 1:ITMAX
        an=-i*(i-a)
        b=b+2.0
        d=an*d+b
        if abs(d) < FPMIN
            d=FPMIN
        end
        c=b+an/c
        if abs(c) < FPMIN
            c=FPMIN
        end
        d=1.0/d
        delt=d*c
        h=h*delt
        if abs(delt-1.0) < EPS
            break
        end
    end
    gammcf = exp(-x+a*log(x)-gln)*h
    return gammcf,gln
end

function test_gamma()
    # gammainc test functions. Test values taken from Mathematica
    @testset "test_gamma" begin

    @test maximum([gammainc(0.5,float(x)) for x in 0:10]
            -[0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                1.77151, 1.77213, 1.77234, 1.77241, 1.77244]) < 1e-5

    @test maximum([gammainc(1.5,float(x)) for x in 0:10]
            -[0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                1.77151, 1.77213, 1.77234, 1.77241, 1.77244]) < 1e-5
    @test maximum([gammainc(2.5,float(x)) for x in 0:10]
            -[0, 0.200538, 0.59898, 0.922271, 1.12165, 1.22933, 
                1.2831, 1.30859, 1.32024, 1.32542, 1.32768]) < 1e-5
    end #testset
end

function test_fgamma()
    @testset "test_fgamma" begin
    for (x,res) in [(0.,1),
                    (30.,0.161802159),
                    (60.,0.114411404),
                    (90.,0.0934165203),
                    (120.,0.08090108),
                    (300.,0.051166336)]
        @test res ≈ Fgamma(0,x)
    end
    end # testset
end

