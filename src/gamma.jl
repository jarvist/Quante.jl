"""
    Fgamma(m, x, small)

Returns the value of Boys integral
```math
F_m(x) = \\int_0^1 t^{2m} e^{-xt^2} dt = 
\\frac{1}{2 x^{m + \\frac{1}{2}}} \\gamma \\left( m + \\frac{1}{2}, x \\right),
```
where ``\\gamma(n, x)`` is the lower incomplete gamma function.

Taylor series for Boys integral, of use for small values of ``x`` 
```math
F_m(x) = \\sum_{k=0}^{\\infty} \\frac{(-x)^k}{k! (2(m + k) + 1)}
```

External links: [DLMF](https://dlmf.nist.gov/8.2.4), 
[Wikipedia](https://en.wikipedia.org/wiki/Incomplete_gamma_function)

See also [`gamma_inc(a,x,ind=0)`](@ref SpecialFunctions.gamma_inc) 
"""
function Fgamma(m::Int64, x::Float64, SMALL::Float64 = 1e-12)
    if x < SMALL
        # use small x expansion
        bf = 1/(2*m + 1) - x/(2*m + 3)    
    else
        bf = gamma(0.5 + m) * gamma_inc(0.5 + m, x)[1] / (2 * x^(0.5 + m))
    end
    return bf
end

function test_fgamma()
    @testset "test_fgamma" begin
        for (x, res) in [
            (0.0, 1),
            (30.0, 0.161802159),
            (60.0, 0.114411404),
            (90.0, 0.0934165203),
            (120.0, 0.08090108),
            (300.0, 0.051166336),
        ]
            @test res ≈ Fgamma(0, x)
        end
    end # testset
end

