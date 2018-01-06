# via https://github.com/rpmuller/pyquante2/blob/4b05475c3616fb48b2536f07023ddc6820dbfccc/othertests/int2e-speed-julia.ipynb

using BenchmarkTools # more accurate @btime macro

println("vrrspeedtest.jl: Compare times between the test cases run using Python/C code (int2e-speed.ipynb) and Julia code.")

function pythonsig_vrr(ax,ay,az,na,aI,aJ,aK,aexpn,bx,by,bz,nb,bexpn,
    cx,cy,cz,nc,cI,cJ,cK,cexpn,dx,dy,dz,nd,dexpn,M)
    # Helper function for vrr, since I for some reason changed the call signature
    # from the python version. This version has the same signature as the python version
    # and calls the Julia version
    na*nb*nc*nd*Quante.vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
end

function pythonsig_vrr_iter(ax,ay,az,na,aI,aJ,aK,aexpn,bx,by,bz,nb,bexpn,
    cx,cy,cz,nc,cI,cJ,cK,cexpn,dx,dy,dz,nd,dexpn,M)
    # Helper function for vrr, since I for some reason changed the call signature
    # from the python version. This version has the same signature as the python version
    # and calls the Julia version
    na*nb*nc*nd*Quante.vrr_iter(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
        cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
end

# Timing comparisons
# Python/C times are from the int2e-speed.ipynb notebook in the same folder.

println("s orbitals on the same center:")
println("RMuller: 9 usec time Julia compared to 0.75 usec in C and 49 usec in Python")
println("Jarv/Asriel commit c2190d5: 1.292 μs (47 allocations: 2.44 KiB) //  15.969 μs (77 allocations: 3.38 KiB)")
#                       xyza   na   lmna   aa   xyzb    nb ab   xyzc      nc  lmnc   ac  xyzd   nd ad  M
@btime pythonsig_vrr(0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0)
@btime pythonsig_vrr_iter(0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0)

println("sp orbitals on the same center:")
println("RMuller: 0.95 usec C, 58 usec Python, 13 usec Julia")
println("Jarv/Asriel commit c2190d5:   3.413 μs (133 allocations: 7.17 KiB) //  16.397 μs (91 allocations: 3.63 KiB)")
#                       xyza   na  lmna   aa  xyzb   nb ab  xyzc   nc  lmnc   ac  xyzd   nd ad  M
@btime pythonsig_vrr(0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0.,0.,0., 1., 1,0,0, 1., 0.,0.,0., 1.,1., 0)
@btime pythonsig_vrr_iter(0.,0.,0., 1., 0,0,0, 1., 0.,0.,0., 1.,1., 0.,0.,0., 1., 1,0,0, 1., 0.,0.,0., 1.,1., 0)

println("sp integrals, 4 different centers")
println("RMuller: 1.2 usec C, 66 usec Python, 23 Julia")
println("Jarv/Asriel commit c2190d5:   3.800 μs (131 allocations: 7.13 KiB) //  16.791 μs (89 allocations: 3.58 KiB)")
@btime pythonsig_vrr(0.,0.,0., 1., 0,0,0, 1., 0.,0.,1., 1.,1., 0.,1.,0., 1., 1,0,0, 1., 1.,1.,0., 1.,1., 0)
@btime pythonsig_vrr_iter(0.,0.,0., 1., 0,0,0, 1., 0.,0.,1., 1.,1., 0.,1.,0., 1., 1,0,0, 1., 1.,1.,0., 1.,1., 0)

println("pd integrals, 4 different centers")
println("RMuller: 3 usec C, 95 usec Python, 32 usec @time Julia")
println("Jarv/Asriel commit c2190d5:   20.108 μs (743 allocations: 40.44 KiB) //  18.694 μs (117 allocations: 4.31 KiB)")
@btime pythonsig_vrr(0.,0.,0., 1., 1,0,0, 1., 0.,0.,1., 1.,1., 0.,1.,0., 1., 1,1,0, 1., 1.,1.,0., 1.,1., 0)
@btime pythonsig_vrr_iter(0.,0.,0., 1., 1,0,0, 1., 0.,0.,1., 1.,1., 0.,1.,0., 1., 1,1,0, 1., 1.,1.,0., 1.,1., 0)

