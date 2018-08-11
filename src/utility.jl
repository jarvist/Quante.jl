factorial2(n::Int64) = prod(n:-2:1) # double factorial !!
dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz # Is there something in the standard library that does this?

pairs(n::Int64) = ((i, j) for i = 1:n for j = 1:i)
rpairs(n::Int64) = ((i,j) for i in 1:n for j in 1:n) # rectangular option to old pairs
spairs(n::Int64) = ((i, j) for i = 1:n for j = 1:(i-1)) # subdiagonal option to old pairs
 
triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i<j ? triangle(j-1)+i : triangle(i-1)+j
                        
iiterator(n::Int64) = ((i,j,k,l) for (i,j) in pairs(n) for (k,l) in pairs(n) if triangle(i,j) <= triangle(k,l))

iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l))
trace2(A,B) = sum(A.*B)

function test_utils()
    @testset "test_utils" begin
    @test factorial2(6)==48
    @test collect(pairs(3)) == [(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)]
    @test collect(spairs(3))== [(2,1),(3,1),(3,2)]
    @test collect(rpairs(2)) == [(1,1),(1,2),(2,1),(2,2)]
    @test iindex(1,1,1,1) == 1
    @test iindex(1,1,1,2) == iindex(1,1,2,1) == iindex(1,2,1,1) == iindex(2,1,1,1) == 2
    @test iindex(1,1,2,2) == iindex(2,2,1,1) == 4
    end #testset
end

