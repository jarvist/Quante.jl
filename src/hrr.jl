# Horizontal and Vertical Recurrance Relationship of Head-Gordon and Pople
#
# A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations
# Martin Head-Gordon and John A. Pople
# Department of Chemistry, Carnegie-Mellon University, Pittsburgh, Pennsylvania 15213 (Received 26 May 1988; accepted 20 July 1988)

# An efficient method is presented for evaluating two-electron Cartesian
# Gaussian integrals, and their first derivatives with respect to nuclear
# coordinates. It is based on the recurrence relation (RR) of Obara and Saika
# [J. Chem. Phys. 84, 3963 (1986)], and an additional new RR, which are
# combined together in a general algorithm applicable to any angular momenta.
# This algorithm exploits the fact that the new RR can be applied outside
# contraction loops. It is shown, by floating point operation counts and
# comparative timings, to be generally superior to existing methods,
# particularly for basis sets containing d functions.
#  https://doi.org/10.1063/1.455553

function hrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
    bexpn::Float64,bx::Float64,by::Float64,bz::Float64,bI::Int64,bJ::Int64,bK::Int64,
    cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
    dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,dI::Int64,dJ::Int64,dK::Int64)
    if bI > 0
        return hrr(aexpn,ax,ay,az,aI+1,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
        (ax-bx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ+1,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (ay-by)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK+1,bexpn,bx,by,bz,bI,bJ,bK-1,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (az-bz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-1,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif dI > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI+1,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK) +
            (cx-dx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK)
    elseif dJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ+1,cK,dexpn,dx,dy,dz,dI,dJ-1,dK) +
            (cy-dy)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ-1,dK)
    elseif dK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK+1,dexpn,dx,dy,dz,dI,dJ,dK-1) +
            (cz-dz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK-1)
    end
    return vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,0)
end

"""
Vertical Recurrence Relationship (VRR) 

"""
function vrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
        bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
        cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
        dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,m::Int64)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
    #println("P: $px,$py,$pz, Q: $qx,$qy,$qz, W: $wx,$wy,$wz, $zeta,$eta")
    
    val::Float64 = 0.0
    if cK>0
        val = (qz-cz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m) +
            (wz-qz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val1=$val")
        if cK>1
            val += 0.5*(cK-1)/eta*(
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m+1) )
        #println("val2=$val")
        end
        if aK>0
            val += 0.5*aK/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        end
        #println("val3=$val")
        return val
    elseif cJ>0
        val = (qy-cy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m) +
        (wy-qy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val4=$val")
        if cJ>1
            val += 0.5*(cJ-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m+1)
            )
        #println("val5=$val")
        end
        if aJ>0
            val += 0.5*aJ/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val6=$val")
        return val
    elseif cI>0
        val = (qx-cx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-qx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val7=$val")
        if cI>1
            val += 0.5*(cI-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val8=$val")
        if aI>0
            val += 0.5*aI/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val9=$val")
        return val
    elseif aK>0
        val = (pz-az)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wz-pz)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val10=$val")
        if aK>1
            val += 0.5*(aK-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val11=$val")
        return val
    elseif aJ>0
        val = (py-ay)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m)+
        (wy-py)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val12=$val")
        if aJ>1
            val += 0.5*(aJ-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val13=$val")
        return val
    elseif aI>0
        val = (px-ax)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-px)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val14=$val")
        if aI>1
            val += 0.5*(aI-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val15=$val")
        return val
    end

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    #println("rab2=$rab2,rcd2=$rcd2,rpq2=$rpq2,T=$T,Kab=$Kab,Kcd=$Kcd")
    return Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
end


function vrr_iter(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
        bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
        cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
        dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,M::Int64)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    mtot = aI+aJ+aK+cI+cJ+cK+M

    vrr_terms = zeros(Float64,(aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,mtot+1))
    
    for m in 0:mtot
        vrr_terms[1,1,1, 1,1,1, m+1] = Fgamma(m,T)*Kab*Kcd/sqrt(zeta+eta)
    end
    
    for i in 0:(aI-1)
        for m in 0:(mtot-i-1)
            vrr_terms[i+2,1,1, 1,1,1, m+1] = (
                 (px-ax)*vrr_terms[i+1,1,1, 1,1,1, m+1] + 
                 (wx-px)*vrr_terms[i+1,1,1, 1,1,1, m+2])

            if i>0
                vrr_terms[i+2,1,1, 1,1,1, m+1] += i/2/zeta*(
                    vrr_terms[i,1,1, 1,1,1, m+1] -
                    eta/(zeta+eta)*vrr_terms[i,1,1, 1,1,1, m+2])
            end
        end
    end
#=
    for j in 0:(aJ-1)
        for i in 0:aI
            for m in 0:(mtot-i-j-1)
                println(("b",i,j,m))
                vrr_terms[i+1,j+2,1, 1,1,1, m+1] = (
                (py-ay)*vrr_terms[i+1,j+1,1, 1,1,1, m+1] +
                (wy-py)*vrr_terms[i+1,j+1,1, 1,1,1, m+2])
                if j>0
                    vrr_terms[i+1,j+2,1, 1,1,1, m+1] += j/2/zeta*(
                        vrr_terms[i+1,j,1, 1,1,1, m+1] -
                        eta/(zeta+eta)*vrr_terms[i+1,j,1, 1,1,1, m+2])
                end
            end
        end
    end
    for k in 0:(aK-1)
        for j in 0:aJ
            for i in 0:aI
                for m in 0:(mtot-i-j-k-1)
                    println(("c",i,j,k,m))
                    vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] = (
                    (pz-az)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+1] +
                    (wz-pz)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+2])
                    if k>0
                        vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] += k/2/zeta*(
                        vrr_terms[i+1,j+1,k, 1,1,1, m+1] -
                        eta/(zeta+eta)*vrr_terms[i+1,j+1,k, 1,1,1, m+2])
                    end
                end
            end
        end
    end
    for q in 0:(cI-1)
        for k in 0:aK
            for j in 0:aJ
                for i in 0:aI
                    for m in 0:(mtot-i-j-k-q-1)
                        println(("d",i,j,k,q,m))
                        vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] = (
                        (qx-cx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+1] +
                        (wx-qx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+2])
                        if q>0
                            vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += q/2/eta*(
                            vrr_terms[i+1,j+1,k+1, q,1,1, m+1] -
                            eta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q,1,1, m+2])
                        end
                        if i>0
                            vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += (
                            i/2/(zeta+eta)*vrr_terms[i,j+1,j+1, q+1,1,1, m+2])
                        end
                    end
                end
            end
        end
    end
    for r in 0:(cJ-1)
        for q in 0:cI
            for k in 0:aK
                for j in 0:aJ
                    for i in 0:aI
                        for m in 0:(mtot-i-j-k-q-r-1)
                            println(("e",i,j,k,q,r,m))
                            vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] = (
                            (qy-cy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] +
                            (wy-qy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+2])
                            if r>0
                                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += r/2/eta*(
                                vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] -
                                zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q+1,r,1, m+2])
                            end
                            if j>0
                                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += (
                                j/2/(zeta+eta)*vrr_terms[i+1,j,k+1, q+1,r+1,1, m+2])
                            end
                        end
                    end
                end
            end
        end
    end
    for s in 0:(cK-1)
        for r in 0:cJ
            for q in 0:cI
                for k in 0:aK
                    for j in 0:aJ
                        for i in 0:aI
                            for m in 0:(mtot-i-j-k-q-r-s-1)
                                println(("f",i,j,k,q,r,s,m))
                                vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2, m+1] = (
                                (qz-cz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+1]+
                                (wz-qz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+2])
                                if s>0
                                    vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += s/2/eta*(
                                    vrr_terms[i+1,j+1,k+1,q+1,r+1,s, m+1] -
                                    zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1,q+1,r+1,s,m+2])
                                end
                                if k>0
                                    vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += (
                                    k/2/(zeta+eta)*vrr_terms[i+1,j+1,k, q+1,r+1,s+1,m+2])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("before return")
    =#
    vrr_terms[aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,M+1]
end

function test_vrr()
    @testset "test_vrr" begin
    ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
    aexpn=bexpn=cexpn=dexpn=1.0
    aI=aJ=aK=0
    cI=cJ=cK=0
    M=0

    for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
            (0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
            (0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
            (0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
            (0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

            (0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
            (0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

            (1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
            (1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
            (1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

            (1.,2.,3., 2,0,0, 2,0,0, 0.00022503308545040895),
            (1.,2.,3., 0,2,0, 0,2,0, 0.0006102470883881907),
            (1.,2.,3., 0,0,2, 0,0,2, 0.0013427831014563411),

            (0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
            (0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

            (3.,2.,1., 1,1,0, 1,1,0, 5.976771621486971e-5),
            (3.,2.,1., 0,1,1, 0,1,1, 1.5742904443905067e-6),
            (3.,2.,1., 1,0,1, 1,0,1, 4.00292848649699e-6)
        ]

        val1 = vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
        val2 = vrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
        @test val1 ≈ val2
        @test val1 ≈ result
        val3 = vrr_iter(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
        val4 = vrr_iter(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
    end
    end# testset
end

function test_hrr()
    @testset "test_hrr" begin
    ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
    aexpn=bexpn=cexpn=dexpn=1.0
    aI=aJ=aK=0
    bI,bJ,bK = 1,0,1
    cI=cJ=cK=0
    dI,dJ,dK = 1,0,1


    for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in [
            (0.,0.,0., 0,0,0, 0,0,0, 0.0136667330685),
            (0.,0.,0., 1,0,0, 1,0,0, 0.00821630976139),
            (0.,0.,0., 0,1,0, 0,1,0, 0.00122024402397),
            (0.,0.,0., 0,0,1, 0,0,1, 0.00821630976139),

            (0.,0.,0., 2,0,0, 2,0,0,   0.0039759617781),
            (0.,0.,0., 0,2,0, 0,2,0,   0.000599953311785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.0039759617781),

            (1.,2.,3., 1,0,0, 1,0,0, -1.1851316496333975e-6),
            (1.,2.,3., 0,1,0, 0,1,0,  -4.669991667384835e-6),
            (1.,2.,3., 0,0,1, 0,0,1, -3.474373852654044e-5),

            (1.,2.,3., 2,0,0, 2,0,0, 2.81002247462e-6),
            (1.,2.,3., 0,2,0, 0,2,0, 7.09856891538e-6),
            (1.,2.,3., 0,0,2, 0,0,2, 3.62153023224e-5),

            (0.,0.,0., 1,1,0, 1,1,0, 0.000599953311785),
            (0.,0.,0., 0,1,1, 0,1,1, 0.000599953311785),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0116431617287),

            (3.,2.,1., 1,1,0, 1,1,0, 7.37307761485e-6),
            (3.,2.,1., 0,1,1, 0,1,1, 2.5333243119843164e-7),
            (3.,2.,1., 1,0,1, 1,0,1, 2.452115184675799e-6)
        ]
        #println("hrr($aexpn,$ax,$ay,$az,$aI,$aJ,$aK,$bexpn,$bx,$by,$bz,$bI,$bJ,$bK,")
        #println("    $cexpn,$cx,$cy,$cz,$cI,$cJ,$cK,$dexpn,$dx,$dy,$dz,$dI,$dJ,$dK)")
        val1 = hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
        val2 = hrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
        @test val1 ≈ val2
        @test val1 ≈ result
    end
    end #testset
end

