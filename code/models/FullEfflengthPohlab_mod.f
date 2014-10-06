! ----------------------------------------------------------------------
! Sean Ressler's code w/ magnetic damping
!
! Cleaned up for review, before porting B damping to Python
!
! Statement labels right aligned from column 5
! Continuation ampersands on column 6
! Exclamation points for all comments
! ----------------------------------------------------------------------

      program Fullefflengthab
        double precision cm, ecut
        double precision B0, v0, rs, alpha
        double precision nu0, nu, ab, emis1
        double precision a, nugraph(100)
        double precision r, delr(40), rgraph(10000,40)
        double precision dx, integrand, integral
        Double Precision compratio
        double precision distr, Bfield
        double precision intensitygraph(40,10000)
        double precision xex(100), fex(100)
        double precision nuroll,x, rho
        double precision width(40), advection(40)
        double precision rmin, norms(40)
        double precision intensity, fluxplat(40), fluxpeak(40)
        double precision gammaplat(40), gammapeak(40)
        double precision rtab(1000), inttab(1000), titab(1000)
        double precision brtab(1000), bttab(1000), altab(1000)
        double precision bgraph(1000), en, rad, delrad, c1
        double precision disttab(100, 10000), radtab(10000)
        double precision distgraph(10000), Bmin, B, slopeemis
        double precision disttab1(100,10000), emistab(10000)
        double precision spline,BIATX(1000)
        DOUBLE PRECISION BCOEF(1000), GTAU(1000), Q(112640), T(1006)

        common/xrays/xex, fex, ir

        !Physics Parameters
        compratio = 4d0
        v0 = 5*1.d8/compratio
        ab = .5d0       ! Damping length, between [0, 1]
        Bmin =5d-6      ! Minimum magnetic field
        rs = 2.96e19    !tan(943.7'')*2.2 kpc
        alpha = .6d0
        s = 2d0*alpha+1d0

        !Constants
        cm = 1.82d18
        c1 = 6.27e18
        a = 1.57d-3

        !Options
        icut = 1

        !Grid Parameters
        iradmax=400  !Grid size for e-distribution table
        irmax =1000  !Grid size for intensity
        ixmax = 500  !Grid size for line of sight
        inumax = 4   !Grid size for frequency
        rmin = .96d0  !Minimum radius of intensity profile
        nu0 = 2.417989d17*1d0  !Initial frequency

!       ! rtab, altab, inttab, titab, brtab, bttab all appear unused?!
!       open(8, file='intpb.tab', status='old')  ! intpb.tab has br, bt
!       read(8, 99) (rtab(i), altab(i), inttab(i), titab(i),
!    c              brtab(i), bttab(i), i = 1, 1000)
!99     format(6(1x, e13.6))

        open(225, file = 'FLUX.dat')
        open(242, file = 'width.dat')
        open(243, file = 'distint.dat')
        open(119, file = 'disttab.dat')
        open(219, file = 'disttab2.dat')
        open(1001, file = 'emistab.dat')
        open(1002, file = 'emisint.dat')
        open(1004, file = 'emistab2.dat')

        print *, 'Enter B0'
        Read *, Bfield

        ! eta=1, mu=1 assumed
        if(icut.eq.1) then
          ecut = 8.3d0*(Bfield/(100d-6))**(-0.5d0)
     c           *(v0*4d0/1d8)*1.6021773d0
        else
          ecut = 1d40
        endif
        print *, ecut

        ! Read in Pacholczyk table
        open(9, file='fglists.dat', status = 'old')
        do 90 ir = 1, 300
           read(9,*) xex(ir), fex(ir)
           if (xex(ir).eq.0.) go to 91
 90     continue
 91     continue
        ir = ir-1  ! Set ir = number of lines read
        close(9)

! _________________________________________________
! Start doing things -- perform numerical integrals


        nu = nu0  ! Initial frequency
        ! Loop over each observation frequency
        do inu = 1,inumax
            B0 = Bfield
            nugraph(inu) = nu

            ! Tabulate e- distribution in a grid of r
            ! and the appropriate Es
            ! for each r (E = sqrt(nu/c1/B(r)/x)
            do j = 1, ir  !loop in xex
                delrad = (1d0-rmin)/dble(iradmax-1)
                rad = rmin
                do irad = 1, iradmax
                    B = Bmin + (B0-Bmin)*dexp((rad-1d0)/ab)
                    en = dsqrt(nu/c1/B/xex(j))
                    disttab(j,irad) = distr(En,B0, Ecut,ab, rad, Bmin)
                    radtab(irad) = rad
                    ! write to disttab.dat if inu=3 (4 keV) and j=3
                    ! (i.e., fixed e- energy?)
                    ! AH not fixed.  remember en varies w/ B field
                    if (inu.eq.3.and.j.eq.3) then
                        write (119, *) rad,disttab(j, irad),en
                        ! file is disttab.dat
                    endif
                    rad = rad + delrad
                enddo
            enddo
            ! Finished populating disttab, radtab

            ! disttab2.dat, normalized e- distr #s
            do i = 1, iradmax
              if (inu.eq.3) then
                write (219, *) i, disttab(31, i)/disttab(31,iradmax)
              endif
            enddo

            ! Takes in radtab, disttab; nu, B0=Bmax, ab, Bmin, iradmax
            ! Integrates over e- energy to get emissivity
            ! but, B = B(z) varies w/ radial dist
            ! so must change scaling factor sqrt(nu*B) in emissivity 
            call emisx(nu, B0, emistab, radtab, disttab, ab, Bmin,
     c                 iradmax)
            ! Finished populating emistab

            if (inu.eq.3) then
              do i = 1, iradmax
                write (1001,*), radtab(i), emistab(i)/emistab(iradmax)
                ! emistab.dat
              enddo
            endif

            if (inu.eq.1) then
              do i = 1, iradmax
                write (1004,*), radtab(i), emistab(i)/emistab(iradmax)
                ! emistab2.dat
              enddo
            endif

            ! Spline stuff is all for the intensity integral
            ! All B-damping stuff is contained in the disttab/emisx
            ! calculations

            k = 6  ! order of the spline
            n = iradmax !number of points for r

            !Set up knot points for spline interpolation over gamma
            do i = 1, k
              t(i) = radtab(1)
            enddo

            do i = 1, n-k
              t(k+i) = 0d0
              do j= 1, k-1
                t(k+i) = t(k+i)+radtab(i+j)
              enddo
              t(k+i) = t(k+i)/dble(k-1)
            enddo

            do i = n+1, n+k
              t(i) = radtab(n)
            enddo

            do i = 1, n
              gtau(i) = emistab(i)
            enddo

            Call Splint (radtab, gtau, t, n, k, q, bcoef, iflag)
            ! Gets spline coefficients

            if (inu.eq.inumax) close (119)

            ! Iterate over grid of r values, to compute radial profiles
            ! at each r we compute intensity(r)
            delr(inu) = -(1d0-rmin)/real(irmax)
            r = 1d0
            do i = 1, irmax
                integrand = 0d0
                integral = 0d0
                dx = -(dsqrt(1d0-r**2))/real(ixmax)
                icurrent = 1
                do isp = 1, n+k-1
                  if (r.lt.t(isp+1).and.r.ge.t(isp)) then
                    icurrent = isp
                    go to 189
                  endif
                enddo
189             continue

                Call BSPLVB (T, k, 1, r, icurrent, BIATX)  !Get the b-splines

                spline = 0d0

                do isp= icurrent-k+1,icurrent
                  iprime = isp-(icurrent-k)
                  spline = spline + bcoef(isp)*BIATX(iprime)
                enddo

                emis1 = spline
                x = dsqrt(1-r**2)+dx
                integral = emis1*dx/3d0
        
                ! Integration over line-of-sight to find intensity
                do ix = 2,ixmax-1
                    rho = dsqrt(x**2+r**2)

                    icurrent = 1
                    do isp = 1, n+k-1
                        if (rho.lt.t(isp+1).and.rho.ge.t(isp)) then
                            icurrent = isp
                            go to 190
                        endif
                    enddo
190                 continue

                    Call BSPLVB (T, k, 1, rho, icurrent, BIATX)  !Get the b-splines

                    spline = 0d0
                    do isp= icurrent-k+1,icurrent
                        iprime = isp-(icurrent-k)
                        spline = spline + bcoef(isp)*BIATX(iprime)
                    enddo
                    emis1 = spline

                    if (i.eq.100.and.inu.eq.1) then
                        write (1002, *) rho, emis1
                    endif

                    rgraph(i, inu) = r

                    ! Simpson's rule
                    if (mod(ix,2).eq.0) then
                        integrand =2d0*emis1
                    else
                        integrand =4d0*emis1
                    endif

                    integral = integral + integrand*dx
                    x = x+dx
                enddo
                ! Done integrating over line of sight

                intensity = integral
                if (i.eq.1) intensity = 0d0
                if (intensity.lt.0d0) intensity = -1d0*intensity
                intensitygraph(inu,i) = intensity

                Call Distgraphr(r, B0, radtab, disttab, distgraph(i))

                r = r+delr(inu)
            enddo
            ! Finished computing intensity over all r

!           rmin = 1d0-(1d0-rmin)*((2d0**.25)**(-.05d0/4))
            nu = nu*2d0
        enddo

        close (999)

        Call FWHM(intensitygraph,inumax,irmax,delr,nugraph,
     c   width, fluxpeak, fluxplat, gammapeak, gammaplat)

        do n = 1, inumax
            print 124, nugraph(n)/2.417989d17, 'keV', width(n)*900d0,
     c              'arcseconds',
     c        log10(width(n)/width(n-1))/log10(nugraph(n)/nugraph(n-1))
     c        ,'mnu'

            write(225,*) nugraph(n)/2.417989d17, fluxpeak(n),
     c                   fluxplat(n), gammapeak(n), gammaplat(n)

            write (242, *) nugraph(n)/2.417989d17, width(n)*900d0

124         format (F5.2,2x, A3, 2x, F6.2, 2x,A10, 2x, F4.2, 2x, A3)
        enddo



        open (111, file = 'Intensity.dat')
        do i =1, irmax
            norms = maxval(intensitygraph, dim = 2)
            write (111,*) rgraph(i,1),
     c                    (intensitygraph(j,i)/norms(j), j = 1, inumax)
            write (243, *) rgraph(i,1), distgraph(i)
        enddo

        close(111)
        close(225)
        close(242)
        close(243)

        end



      subroutine emisx(nu, Bmax,emistab,
     c  radtab, disttab,ab,Bmin,iradmax)

        double precision en, en0, dist, argexp, oldsum
        double precision r,xold, xint1, x, rho, nu, ef
        double precision v, ecut,a, c1, nuroll,cm
        double precision efinv,rs, Xi, sum, B, ab
        double precision Bfield, Bmin, z
        double precision fex(100), xex(100)
        double precision rtab(1000), brtab(1000)
        double precision bttab(1000), br, bt, slopebr
        double precision slopebt, Bsed, int, integrand
        double precision xprime, dxprime, oldintegrand
        double precision radtab(10000), disttab(100,10000)
        double precision slopedist, Bmax, emistab(10000)
        integer iradmax
    
        Double precision emis1
    
        common/xrays/xex, fex, ir
    
        do irad = 1, iradmax
            z = 1d0-radtab(irad)
            Bfield = Bmin + (Bmax-Bmin)*dexp(-z/ab)

            xint1 = 0.
            x = xex(1)
            dist = disttab(1, irad)
            oldsum1 = (fex(1))*dist/x**1.5d0
            xold = x

            ! Integrate over y-values (e- energies) in Pacholczyk
            do ix = 2, ir
                x = xex(ix)  ! y = xex(ix) in paper
                dist = disttab(ix,irad)  ! Not interpolating here?
                sum1 = fex(ix)*dist/(x)**1.5d0
                ! Increment trapezoidal sum
                xint1 = xint1 + (sum1 + oldsum1)*(x - xold)
                oldsum1 = sum1
                xold = x
            enddo

            emistab(irad) = xint1*dsqrt(nu*Bfield)
        enddo
        return
      end


      ! Assumes Bohm diffusion (so eta=1, mu=1)
      function distr(E,Bmax, Ecut,absc, rsc, Bmin)

        double precision D0, B, a, b0, alpha
        double precision k0, s, v, Ecut, integrand
        double precision x, n, argexp,pi, dn, rs, r
        double precision integral, E, rsc, dr,vs
        double precision oldintegrand, nmin, dy
        double precision y,q,t, dt, dellog, dE
        double precision nold, nnext
        double precision tmax
        double precision tmin, tnext, told, lad, ldiff
        double precision tau, ndelta, a1, b1, c1
        double precision epsmax, epsmin, eps
        double precision norm, scale, dist
        double precision en0, ef, efinv, xi
        double precision ab, Bmax, Bmin, z
        double precision const, distr, absc, eta

        const = Bmax-Bmin
        eta = 1d0

        D0 = eta*2.083d19/Bmax
        a = 1.57d-3
        k0 = 1d0
        b0 = a*Bmax**2d0
        alpha = 1d0
        pi = 4d0*datan(1d0)

        ! SNR parameters
        s = 2.2d0
        vs = 5d8
        v = vs/4d0
        rs = 2.96d19

        ! Damping lengthscale (cm)
        ab = absc*rs

        tau = 1d0/b0/E
        lad = v*tau
        ldiff = dsqrt(D0*E*tau)

        r = rsc*rs
        z = (rs-r)

        x = (z*(Bmin**2d0)/Bmax**2d0+
     c    ab/2d0*(Bmax-Bmin)**2d0/Bmax**2d0*(1d0-dexp(-2d0*z/ab)) +
     c    2d0*ab*(Bmax-Bmin)*Bmin/Bmax**2d0*(1d0-dexp(-z/ab)))/alpha

        if (absc.gt.1d3) x =z  ! Use z(x)=x for weak damping

        tmin = 0d0
        tmax = 1d0

        ! Advective solution if l_ad/l_diff > 30
        if (lad/ldiff.gt.30d0) then
            efinv = b0/v*x
            ef =1d0/efinv
            en0 = E/(1d0-E/ef)
            argexp = en0/ecut
            Xi = 1d0-E/ef

            if (Xi.gt.0d0) then
              integral = 1./en0**s*(en0/E)**2/dexp(argexp)*8d-9
            else
              integral = 0d0
            endif

            go to 444
        endif

        itmax = 5000
        oldintegrand = 0d0
        integral = 0d0
        dt = (tmax-tmin)/dble(itmax-1)
        t = tmin

        do it = 1, itmax

        n = dexp(t**2d0)

        argexp = n*E/Ecut+(lad/alpha*(1d0-1d0/n)-x)**2d0/
     c      (4d0*ldiff**2d0/alpha*dlog(n))
        if (n.eq.1d0) argexp = 1d35
        integrand = 2d0*dexp((1d0-s)*t**2d0)/dexp(argexp)

        if (integrand.lt.1d45) then
        else
c       print *, x, z, bmax, bmin, n
        endif

        if (it.eq.1) go to 2005
        integral = integral + (integrand + oldintegrand)*dt/2d0
2005    continue
        oldintegrand = integrand

        t = t + dt
        enddo

        inmax=100
        dy = 1d0/dble(inmax)

        y = 0d0
        q = 2d0
        nmin = dexp(1d0)

        oldintegrand = 0d0
        do i = 1, inmax

          n = y/(1d0-y**2d0)**q+nmin


          argexp = n*E/Ecut+(lad/alpha*(1d0-1d0/n)-x)**2d0/
     c        (4d0*ldiff**2d0/alpha*dlog(n))
          integrand = n**(-1d0*s)/dsqrt(dlog(n))/dexp(argexp)*
     c     ((2d0*q-1d0)*y**2d0+1d0)/
     c     (1d0-y**2d0)**(q+1d0)

          if (integrand.lt.1d45) then
          else
c           print *, integrand, argexp
          endif

          if (i.eq.1) go to 114
          integral = integral + (integrand+oldintegrand)*dy/2d0
114       continue
          oldintegrand = integrand

          y= y + dy
        enddo

        integral = k0/2d0/dsqrt(pi*alpha*b0*D0)*
     c   E**(-1d0*(1d0+s))*integral

444     continue

        distr = integral

        return
        end



      SUBROUTINE FWHM(intensitygraph,inumax,irmax,delr,nugraph,width,
     c                  fluxpeak, fluxplat, gammapeak, gammaplat)

        double precision xmax xmin, maxintensity, halfint
        double precision onedint(10000), delr(40), nugraph(40)
        double precision width(40), fluxplat(40), fluxpeak(40)
        double precision intensitygraph(40,10000)
        double precision gammaplat(40), gammapeak(40), fac
        integer iplat(100), ixmax(100), ixmin(100)

        imax = 1
        halfint =0d0
        do i = 1, irmax
          onedint(i) =0d0
        enddo

        do inu=1,inumax
          do i = 1, irmax
            onedint(i) = intensitygraph(inu,i)
          enddo

          xmax =0d0
          xmin = 0d0
          maxintensity = 0d0
          maxintensity = maxval(onedint)

        do i = 1, irmax
          if (onedint(i).eq.maxintensity) imax =i
        enddo

        halfint = .5d0*maxintensity

        do j = 1, imax-1
          i = imax-j
          if (onedint(i).lt.halfint) then
            xmax = 1d0+real(i+1)*delr(inu)
            ixmax(inu) = imax-i    !closest to shock
            goto 993
          endif
        enddo
        if (xmax.eq.0d0) xmax = 1d0
993     continue

        do i = imax, irmax
          if (onedint(i).lt.halfint) then
            xmin = 1d0+real((i-1))*delr(inu)
            ixmin(inu) = i     !Farthest from shock
            goto 994
          endif
        enddo
994     continue

        width(inu) = xmax-xmin

        if (xmax.eq.0d0) then
          print *, 'Box Length Error at', nugraph(inu)/2.417989d17
        elseif (xmin.eq.0d0) then
          print *, 'Box Length Error at', nugraph(inu)/2.417989d17
        elseif ((xmax-xmin)/dabs(delr(inu)).lt.20d0) then
          print *, 'Resolution Error at', nugraph(inu)/2.417989d17
        endif

        fluxpeak(inu) = 0d0
        fluxplat(inu) = 0d0

        iplat(inu) = irmax !2*ixmin(inu)!+(ixmin(inu)-ixmax(inu))

        if (iplat(inu).lt.1) then
          print *, 'Make Box Larger to fit Plateau'
        endif

        enddo

        ifix = 8   !Frequency at which rims are observed

        do inu = 1, inumax
          do i = 2, ixmin(ifix)
            fluxpeak(inu) = fluxpeak(inu) +
     c       (intensitygraph(inu,i)+intensitygraph(inu,i-1))/2d0
          enddo

          do i = ixmin(ifix)+1, iplat(ifix)
            fluxplat(inu) = fluxplat(inu) +
     c        (intensitygraph(inu,i)+intensitygraph(inu,i-1))/2d0
          enddo
        enddo

        fac = 2d0

        do inu = 1, inumax-1

!         if (inu.eq.1) then
            gammapeak(inu) =-1d0*
     c        dlog(fluxpeak(inu+1)/fluxpeak(inu)/fac)/dlog(fac)
            gammaplat(inu) = -1d0*
     c        dlog(fluxplat(inu+1)/fluxplat(inu)/fac)/dlog(fac)
!         elseif (inu.eq.inumax) then
!           gammapeak(inu) =-1d0*
!    c       dlog(fluxpeak(inu)/fluxpeak(inu-1)/fac)/dlog(fac)
!           gammaplat(inu) = -1d0*
!    c       dlog(fluxplat(inu)/fluxplat(inu-1)/fac)/dlog(fac)
!         else
!           gammapeak(inu) =-1d0*
!    c       dlog(fluxpeak(inu+1)/fluxpeak(inu-1)/fac**2d0)/dlog(fac**2d0)
!           gammaplat(inu) = -1d0*
!    c       dlog(fluxplat(inu+1)/fluxplat(inu-1)/fac**2d0)/dlog(fac**2d0)
!         endif
        enddo

        return
      end

      subroutine Distgraphr(r, B, radtab, disttab, dist)

        double precision en, en0, dist, argexp, oldsum
        double precision r,xold, xint1, x, rho, nu, ef
        double precision v, ecut,a, c1, nuroll,cm
        double precision efinv,rs, Xi, sum, B, ab
        double precision Bfield, Bmin, const, z
        double precision fex(100), xex(100)
        double precision rtab(1000), brtab(1000)
        double precision bttab(1000), br, bt, slopebr
        double precision slopebt, Bsed, int, integrand
        double precision xprime, dxprime, oldintegrand
        double precision radtab(1000), disttab(100,1000)
        double precision slopedist

        Double precision emis1

        common/xrays/xex, fex, ir

        do 17 j = 2, 1000
          if (r.gt.radtab(j)) go to 17
          slopedist = (disttab(30,j) - disttab(30,j-1)) /
     c                (radtab(j) - radtab(j-1))
          dist = disttab(30,j-1) + (r - radtab(j-1))*slopedist
          go to 18
 17     continue
 18     continue

        return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c The following collection of subroutines is taken directly out of A Practical Guide
c To Splines by de Boor.  They solve a linear system for the coefficients of the linear
c combination of splines (recursively defined polynomials)
        SUBROUTINE SPLINT (TAU, GTAU, T, N, K, Q, BCOEF, IFLAG)


c       k is the order of the spline, n is the number of data points
c       tau is an array of the coordinate
c       gtau is an array of the values of the function at those coordinates
c       Q contains the LU factorizaiton of the spline matrix
c       iflag = 1 for success
c
        Integer iflag, k, n, i, ilp1mx, j, jj, km1, kpkm2, left, lenq
        integer np1
        DOUBLE PRECISION BCOEF(N), GTAU(N), Q(112640), T(N+K)
        DOUBLE PRECISION TAU(N)
        DOUBLE PRECISION TAUI


        NP1 = N+1
        KM1 = K-1
        KPKM2 = 2*KM1
        LEFT = K

        LENQ = N*(K+KM1)


        DO 5 I = 1, LENQ
5       Q(I) = 0

        DO 30 I = 1, N
                TAUI = TAU(I)
                ILP1MX = MIN0(I+K,NP1)

                LEFT = MAX0(LEFT,I)
                IF (TAUI.LT.T(LEFT))    GO TO 998
15              IF (TAUI.LT.T(LEFT+1)) GO TO 16
                LEFT = LEFT +1
                IF (LEFT.LT.ILP1MX) GO TO 15
                LEFT = LEFT-1
                IF (TAUI.GT.T(LEFT+1)) GO TO 998


16              CALL BSPLVB (T,K,1,TAUI,LEFT,BCOEF)

                JJ = I-LEFT+1+(LEFT-K)*(K+KM1)
                DO 30 J=1, K
                        JJ = JJ+KPKM2
30                      Q(JJ) = BCOEF(J)

                CALL BANFAC(Q,K+KM1,N,KM1, KM1, IFLAG)

                IF (IFLAG.EQ.1) THEN
                        GO TO 40
                ELSE
                        GO TO 999
                ENDIF


40              DO 41 I = 1,N
41                      BCOEF(I) = GTAU(I)
                CALL BANSLV(Q,K+KM1,N,KM1,KM1,BCOEF)


                RETURN

998             IFLAG = 2
999             PRINT *, 'LINEAR SYSTEM IN SPLINT NOT INVERTIBLE'

                                                                RETURN

                END




                SUBROUTINE BANFAC(W,NROWW,NROW,NBANDL,NBANDU,IFLAG)

                INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW
                INTEGER I, IPK, J, JMAX
                INTEGER K, KMAX, MIDDLE, MIDMK, NROWM1
                DOUBLE PRECISION W(NROWW,NROW), FACTOR, PIVOT

                IFLAG = 1
                MIDDLE = NBANDU+1

                NROWM1 = NROW-1
                IF (NROWM1)     999, 900, 1
1               IF (NBANDL.GT.0)        GO TO 10

                DO 5 I = 1,NROWM1
                        IF (W(MIDDLE,I).EQ.0D0) GO TO 999
5               CONTINUE

                GO TO 900

10              IF (NBANDU.GT.0) GO TO 20

                DO 15 i = 1, NROWM1
                        PIVOT = W(MIDDLE,I)
                        IF (PIVOT.EQ.0D0) GO TO 999
                        JMAX = MIN0(NBANDL,NROW-I)
                        DO 15 J = 1,JMAX
15                              W(MIDDLE+J, I) = W(MIDDLE+J,I)/PIVOT
                                GO TO 900


20              DO 50 I = 1, NROWM1

                        PIVOT = W(MIDDLE,I)
                        IF (PIVOT.EQ.0D0) GO TO 999

                JMAX = MIN0(NBANDL,NROW-I)

                DO 32 J  = 1,JMAX
32                      W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT

                KMAX = MIN0(NBANDU,NROW-I)

                DO 40 K=1,KMAX
                        IPK = I+K
                        MIDMK = MIDDLE-K
                        FACTOR = W(MIDMK,IPK)
                        DO 40 J= 1, JMAX
40                              W(MIDMK+J,IPK) = W(MIDMK+J,IPK)-
     C                             W(MIDDLE+J,I)*FACTOR
50              CONTINUE

900             IF (W(MIDDLE,NROW).NE.0) RETURN
999             IFLAG = 2


                RETURN
                END



                SUBROUTINE BANSLV(W,NROWW,NROW,NBANDL,NBANDU,B)

                INTEGER NBANDL, NBANDU, NROW, NROWW, I,J, JMAX, MIDDLE
                INTEGER NROWM1
                DOUBLE PRECISION W(NROWW,NROW), B(NROW)
                MIDDLE = NBANDU+1
                IF (NROW.EQ.1) GO TO 49
                NROWM1 = NROW-1
                IF (NBANDL.EQ.0) GO TO 30


                DO 21 I = 1, NROWM1
                        JMAX = MIN0(NBANDL, NROW-I)
                        DO 21 J=1,JMAX
21                              B(I+J)= B(I+J) - B(I)*W(MIDDLE+J,I)

30              IF (NBANDU.GT.0) GO TO 40

                DO 31 I = 1, NROW
31                      B(I) = B(I)/W(1,I)
                RETURN

40              DO 45 I = NROW, 2, -1
                        B(I) = B(I)/W(MIDDLE,I)
                        JMAX = MIN0(NBANDU,I-1)
                        DO 45 J=1, JMAX
45                              B(I-J) = B(I-J)-B(I)*W(MIDDLE-J,I)
49              B(1) = B(1)/W(MIDDLE,1)

                RETURN
                END



                SUBROUTINE BSPLVB (T, JHIGH, INDEX, X, LEFT, BIATX)

                PARAMETER (JMAX = 20)
                INTEGER INDEX, JHIGH, LEFT, I, J, JP1
                DOUBLE PRECISION BIATX(1000), T(1006), X, DELTAL(JMAX)
                DOUBLE PRECISION SAVED, TERM, DELTAR(JMAX)
                DATA J/1/
                SAVE J, DELTAL, DELTAR

                IF (INDEX.EQ.1) THEN
                        GO TO 10
                ELSE
                        GO TO 20
                ENDIF


10              J = 1
                BIATX(1) = 1D0
                IF (J.GE.JHIGH) GO TO 99
20                  JP1 = J+1
                    DELTAR(J) = T(LEFT+J) -X
                    DELTAL(J) = X-T(LEFT+1-J)
                    SAVED = 0D0
                    DO 26 I = 1,J
                        TERM = BIATX(I)/(DELTAR(I)+DELTAL(JP1-I))
                        BIATX(I) = SAVED+DELTAR(I)*TERM
26                      SAVED = DELTAL(JP1-I)*TERM
                    BIATX(JP1) = SAVED
                    J = JP1
                    IF (J.LT.JHIGH) GO TO 20

99              RETURN
                END



