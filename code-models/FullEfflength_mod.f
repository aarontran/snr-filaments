! ----------------------------------------------------------------------
! Sean Ressler's program for computing radial profiles
! of synchrotron thin rims, given parameters B0, eta2, mu
! B0 = magnetic field, eta2 = scaled diffusion coefficient at 2 keV
! mu = spectral index of diffusion coefficient power law dependence
!
! Reference: Ressler et al., 2014, ApJ
!
! UPDATED by Aaron Tran for application to Tycho's SNR
! Summer 2014
!
! Statement labels right aligned from column 5
! Continuation ampersands on column 6
! Exclamation points for all comments (even though not FORTRAN 77)
! ----------------------------------------------------------------------

      program Fullefflength
          ! implicit none
          double precision cm, ecut
          double precision B0, v0, rs, rsarc, alpha
          double precision nu0, nu, nukev, emis1
          double precision a, nugraph(100)
          double precision r, delr(40), rgraph(10000,40)
          double precision dx, integrand, integral
          double precision compratio
          double precision distrpohl, distrmlt1, distrmgt1,Bfield
          double precision intensitygraph(40,10000)
          double precision xex(100), fex(100)
          double precision x, rho
          double precision width(40), advection(40)
          double precision rmin, rminarc, norms(40)
          double precision intensity, fluxplat(40), fluxpeak(40)
          double precision gammaplat(40), gammapeak(40)
          double precision en, rad, delrad, c1
          double precision disttab(100, 1000), radtab(1000)
          double precision distgraph(1000)
          double precision eta, mu, eta2

          ! ADDITION: explicit variable declarations
          ! to avoid type mismatch (and bugs and things)
          ! BUT I think variables E, B are wrong anyways...
          double precision E, B, s
          double precision emisx  ! function
          double precision Distgraphr  ! function
          integer ir

          integer icut, ifit, iplot

          common /xrays/ xex, fex, ir

! Physical Parameters

          compratio = 4d0          ! Rankine-Hugoniot
          v0 = 5d0*1.d8/compratio  ! As measured by Satoru et al. (cm/s)
          rs = 2.96d19             ! Modified from nufit code....distance is now 2.2 kpc,
                                   ! so rs is found from tan(943.7'')*2.2 kpc 
          rsarc = 900d0            ! From Green's SNR Catalog...rs in arcseconds
          alpha = .6d0             ! Spectral index
          s = 2d0*alpha+1d0        

! Constants
          cm = 1.82d18         ! Constant needed for max synchrotron frequency
          c1 = 6.27e18         ! Constant needed for synchrotron emissivity
          a = 1.57d-3          ! Constant needed for synchrotron loss time
          nukev = 2.417989d17  ! 1keV photon frequency

! TODO: allow options/parameters to be read from external text file
! so that source doesn't have to be recompiled to vary parameters
! low priority though

! Options
          icut = 0    ! 1 for cut-off, 0 for no cut-off 
          ifit = 1    ! 1 for the fitting version of the code
          iplot = 0   ! 1 for the plotting version of the code  

! Grid Parameters

          irmax = 400    ! resolution on intensity profile
          iradmax = 100  ! resolution on tabulated electron distribution
          ixmax =  500   ! resolution along line of sight for integration
 
! Start doing things -- read input parameters

          open(242, file = 'widthcut2.dat')
          open(243, file = 'distint.dat')
          ! open(119, file = 'disttab.dat')  ! Doesn't seem to be in use?
                                             ! No other references to 119
          open(225, file='FLUX.dat')  ! Doesn't seem to be in use, as
                                      ! lines are commented out below
                                      ! Superseded by widthcut2.dat?

          print *, 'Enter B0'
          read *, Bfield        ! Gauss
          print *, 'Enter Eta2'
          read *, eta2          ! dimensionless
          print *, 'Enter mu'
          read *, mu            ! dimensionless

          ! eta2 = dimensionless diffusion coefficient scale at 2 keV
          ! eta = eta2 * (2 keV)**(1-mu)
          ! I think eta maps to \eta*(E_h)^(1-\mu), as in Section 4.2
          ! of the paper
          eta =  eta2 * (2d0*nukev/(cm*Bfield))**(-(mu-1d0)/2d0)

          if (icut.eq.1) then
              ! TODO: Equation (19), but disagrees w/paper? check!
              ecut = (8.3d0*(Bfield/(100d-6))**(-1d0/2d0)*(v0*4d0/1d8)*
     &            1.6021773d0)**(2d0/(1d0+mu))*
     &            (1d0/eta)**(1d0/(mu+1d0))
          else
              ecut = 1d40
          endif

! Read frequency-dep. synchrotron emissivity table to common block
! variables xex, fex
          open(9, file='fglists.dat', status='old')
          do 90 ir = 1, 300
              read(9,*) xex(ir), fex(ir)
              if (xex(ir).eq.0.) go to 91  ! Break at end of file
   90     continue
   91     continue
          ir = ir-1  ! Set ir = number of lines read
          close(9)

! Set up from optional parameters (fit/plot)
! TODO: twiddle-able rminarc, should be placed with other params
! (more generally -- port to Python to vary nu more readily)
          if (ifit.eq.1) then
              nu0 = nukev
              inumax = 2
              rminarc = 60d0  ! close to FWHM being fitted?  SET THIS!!
          endif

          if (iplot.eq.1) then
              nu0 = nukev*.1d0
              inumax = 32
              rminarc = 110d0
          endif

          ! rmin = fraction of SNR r (900arcsec) to go in, range [0,1]
          rmin = 1d0 - rminarc/rsarc  ! The only place rminarc is used
          nu = nu0

! Perform all numerical integrals and calculate intensity profiles
! in each observation band (iterate over nu, which is incremented at
! bottom of the loop).
          do inu = 1, inumax
            B0 = Bfield
            nugraph(inu) = nu  ! Store nu in array

            ! Tabulate electron distribution for fixed nu
            ! generating 2-D grid as function of energy, radial position

            ! Step over particle energies from Pacholczyk table
            do j = 1, ir  ! ir = number of entries in Pacholczyk
              ! Step from r=rmin to r=1 (IRADMAX resolution)
              delrad = (1d0-rmin)/dble(iradmax-1)
              rad = rmin
              do irad = 1, iradmax
                en = dsqrt(nu/c1/B0/xex(j))  ! e- energy for xex(j)
                ! TODO: B->B0 or Bfield?
                ! TODO: dive into details of distribution functions
                if (mu.lt.1d0) then
                  disttab(j,irad) = distrmlt1(en,B0, ecut,rad,
     &                                        eta,mu, rs, v0,s)
                elseif (mu.gt.1d0) then
                  disttab(j,irad) = distrmgt1(en,B0, ecut,rad,
     &                                        eta,mu, rs, v0,s)
                else
                  disttab(j,irad) = distrpohl(en,B0, ecut,rad,
     &                                        eta,mu, rs, v0,s)
                endif
                radtab(irad) = rad  ! Store/increment radial coord
                rad = rad + delrad  ! rmin to 1.0 (scaled coord)
              enddo
            enddo


            ! For each radial position in intensity profile
            ! Integrate over line of sight to get I_\nu(r)
            ! At each step of line-of-sight integration,
            ! Compute synchrotron emissivity at radial coord
            ! rho = s**2 + r**2 (r = "viewer's" radial coord)
            ! (integrates over particle energy via emisx)

            ! Step from r=1 to r=rmin (IRMAX resolution)
            delr(inu) = -(1d0-rmin)/real(irmax)  ! Could depend on inu,
            r = 1d0                              ! but currently doesn't
            do i = 1, irmax

              ! TODO: Are we using Simpson's rule?????
              ! if so, need to 1. divide by 3.0 in loop, and 
              ! 2. handle edge case at ix = ixmax-1

              ! Initialize integral at s=sqrt(1 - r**2) to s=0
              ! argument of j_nu ranges (r to rs), flipped wrt paper
              ! emis1 = j_nu, emisx integrates disttab against G(y) dy
              dx = -(dsqrt(1d0-r**2))/real(ixmax)

              ! originally called emisx(r,...), changed to emisx(1, ...)
              ! SHOULD BE
              ! x = dsqrt(1d0-r**2)
              ! rho = dsqrt(x**2 + r**2)  ! Simply 1d0
              ! emis1 = emisx(rho ...)
              ! x = x + dx
              emis1 = emisx(1d0, nu,B0,radtab,disttab) ! j_nu(r)
              integral = emis1*dx/3d0
              x = dsqrt(1d0-r**2)+dx

              ! Integration to find intensity
              do ix = 2,ixmax-1
                rho = dsqrt(x**2+r**2)
                emis1 = emisx(rho, nu,B0,radtab,disttab)

                if (mod(ix,2).eq.0) then
                  integrand =2d0*emis1     
                else
                  integrand =4d0*emis1 
                endif

                integral = integral + integrand*dx
                x = x+dx
              enddo

              intensity = integral
              ! TODO: this if check shouldn't be necessary
              if (intensity.lt.0d0) intensity = -1d0*intensity
              intensitygraph(inu,i) = intensity
  
              ! Plots, remember i indexes "viewer's" radial coord
              rgraph(i, inu) = r
              distgraph(i) = Distgraphr(r, radtab, disttab)
  
              r = r+delr(inu)  ! Increment radial pos for next loop
            enddo
            ! Finished computing integrals for intensity profile


            ! Update value of nu for next loop
            ! TODO: change this for Tycho, or other SNRs
            if (ifit.eq.1) then
              nu = 2d0*nu
            endif
            if (iplot.eq.1) then
              nu = nu*10d0**(.07d0)
              rmin = 1d0-(1d0-rmin)*((10d0)**(-.25d0/16d0))
            endif

            print *, 'rmin is: ', rmin
            print *, 'ecut is: ', ecut
          enddo
          ! DONE computing stuff for all energy bands!

          Call FWHM(intensitygraph,inumax,irmax,delr,nugraph,
     &              width, fluxpeak, fluxplat, gammapeak, gammaplat,imu)

! Now, just writing output to stdout / files

          do n = 1, inumax
            print 124, nugraph(n)/nukev, 'keV',
     &        width(n)*rsarc, 'arcseconds',
     &        log10(width(n)/width(n-1))/log10(nugraph(n)/nugraph(n-1)),
     &        'mnu'  ! TODO: this doesn't work for n=1
  124       format (F5.2,2x, A3, 2x, F6.2, 2x,A10, 2x, F4.2, 2x, A3)

!           write(225,*) nugraph(n)/nukev, fluxpeak(n), fluxplat(n), 
!    &        gammapeak(n), gammaplat(n)
            write (242,*) nugraph(n)/nukev, width(n)*rsarc
          enddo

          open (111, file = 'Intensity.dat') 
          do i =1, irmax
            norms = maxval(intensitygraph, dim = 2)
            write (111,*) rgraph(i,1), (intensitygraph(j,i)/norms(j)
     &                    , j = 1, inumax)
            write (243,*) rgraph(i,1), distgraph(i)  ! distint.tab
          enddo

          close(111)  ! Intensity.dat (opened just above)
          close(225)  ! FLUX.dat (not currently in use)
          close(242)  ! widthcut2.tab
          close(243)  ! distint.tab

          stop
      end


! ----------------------------------------------------------------
! Integrate 1-particle emissivity with particle distribution over
! radiation frequency (change dE to dy) to get emissivity j_\nu(r)
!
! Cleaned 2014 July 16
! ----------------------------------------------------------------

! j_nu(r, nu, B, radtab, disttab) is a function of radial coordinate r
! noting that r != x ... set by e- dist. functions in terms of r, not x

! TODO: what happened to factor $c_3 * 1/sqrt{c_1} * 1/4$ ???
! Maybe overall normalization doesn't matter...
! If overall normalization doesn't affect FWHMs, why include sqrt(nu*B)?
! It could matter for B = B(x) in magnetically damped model.
! nu would not matter for FWHM at all, except to get comparative
! normalization right for plots of multiple filaments.
! Well, that makes sense too.

      function emisx(r, nu, B, radtab, disttab)
          implicit none

          double precision r, nu, B
          double precision radtab(1000), disttab(100,1000)
          double precision fex(100), xex(100)
          integer ir

          ! Variables to perform integration
          double precision dist, slopedist
          double precision x, xold, intgdold, intgd, trpsum, emisx
          integer j, ix

          common /xrays/ xex, fex, ir
 
          do 14 j = 2, 1000  ! Interpolate e- distr., at position r
              if (r.gt.radtab(j)) go to 14  ! Increment radtab(j)
              slopedist = (disttab(1,j) - disttab(1,j-1))/
     &                    (radtab(j) - radtab(j-1))
              dist = disttab(1,j-1) + (r - radtab(j-1))*slopedist
              go to 15  ! Break out of loop
   14     continue
   15     continue

          ! TODO: catch edge case when r = 1.0
          ! Because I think I'm seeing cases of j=1001
          if (j.eq.1001) then
              print *, 'DEBUG emisx: r=',r,
     &                 ', radtab(1000)=',radtab(1000),
     &                 ', radtab(1)=',radtab(1)
          endif

          ! Initialize trapezoidal sum
          x = xex(1)
          intgdold = (fex(1))*dist/(x**1.5d0)
          xold = x
          trpsum = 0.

          ! Integrate over y-values (e- energies) in Pacholczyk
          do ix = 2, ir
            ! Compute G(y) for new y = xex(ix); j is same as before!
            slopedist = (disttab(ix,j) - disttab(ix,j-1))/
     &                  (radtab(j) - radtab(j-1))
            dist = disttab(ix,j-1) + (r - radtab(j-1))*slopedist

            ! Increment trapezoidal sum
            x = xex(ix)
            intgd = fex(ix)*dist/(x**1.5d0)
            trpsum = trpsum + (intgd + intgdold)*(x - xold)
            xold = x
            intgdold = intgd
          enddo

          emisx = trpsum*dsqrt(nu*B)  ! Why this normalization?
          return
      end


! ------------------------------------------------
! Electron distribution function, for mu < 1
! Equation (15), from Lerche & Schlickeiser (1980)
! ------------------------------------------------

      function distrmlt1(E,B, Ecut,r, eta,mu, rs, v,s)
          double precision v, tau, E, Ecut
          double precision distrmlt1, vs
          double precision D0, B, a, b0, alpha
          double precision k0, s, integrand
          double precision x, n, argexp,pi, dn, rs, r
          double precision integral, rsc, dr
          double precision oldintegrand, nmin, dy
          double precision y,q,t, dt, dellog, dE
          double precision tmax, integrandmax
          double precision tmin, lad, ldiff
          double precision epsmax, epsmin, eps
          double precision en0, efinv, ef,  Xi
          double precision br, bt
          double precision  mu, muvector(1:5)
          double precision Bfield, eta, eh
          integer imu

          D0 = eta*2.083d19/B  ! eta*C_d/B0
          a = 1.57d-3
          k0 = 1d0

          b0 = a*B**2d0
          alpha = 1d0
          pi = 4d0*datan(1d0)

          tau = 1d0/b0/E
          lad = v*tau
          ldiff = dsqrt(D0*E**mu*tau)

          x = (rs-r*rs)/alpha

          tmin = 0d0
          tmax = 1d0

          if (lad/ldiff.gt.30d0) then  !If lad<<ldiff the solution is easy
            efinv = a/v*B**2d0*(rs-r*rs)
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

          itmax = 1000
          oldintegrand = 0d0
          integral = 0d0
          dt = (tmax-tmin)/dble(itmax-1)
          t = tmin

          do it = 1, itmax
            n = 1d0-t**2d0

            argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     &               (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     &               (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)

            if (n.eq.1d0) argexp = 1d35

            integrand = 2d0*(1d0-t**2d0)**((s+mu-2d0)/(1d0-mu))
     &                  /dexp(argexp)
            if (it.eq.1) go to 2222
            integral = integral + (integrand + oldintegrand)*dt/2d0
 2222       continue
            oldintegrand = integrand
            t = t + dt
          enddo

          integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(1d0-mu))*
     &               E**(-1d0*(mu/2d0+1d0/2d0+s))*integral
          distrmlt1 = integral

  444     continue

          distrmlt1 = integral

          return 
      end

! ------------------------------------------
! Electron distribution function, for mu = 1
! Equation (14), from Rettig and Pohl (2012)
! ------------------------------------------

      function distrpohl(E,B, Ecut,r, eta,mu, rs, v,s)
        double precision v, tau, E, Ecut
        double precision distrpohl
        double precision D0, B, a, b0, alpha
        double precision k0, s, integrand
        double precision x, n, argexp,pi, dn, rs, r
        double precision integral, rsc, dr
        double precision oldintegrand, nmin, dy
        double precision y,q,t, dt, dellog, dE
        double precision tmax, integrandmax
        double precision tmin, lad, ldiff
        double precision epsmax, epsmin, eps
        double precision en0, efinv, ef,  Xi
        double precision br, bt
        double precision Bfield, eta, mu

        D0 = eta*2.083d19/B
        a = 1.57d-3
        k0 = 1d0
        b0 = a*B**2d0
        alpha = 1d0
        pi = 4d0*datan(1d0)

        tau = 1d0/b0/E
        lad = v*tau
        ldiff = dsqrt(D0*E*tau)

        x = (rs-r*rs)/alpha

        tmin = 0d0
        tmax = 1d0
        
        if (lad/ldiff.gt.30d0) then 
            efinv = a/v*B**2d0*(rs-r*rs)
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

        itmax = 1000
        oldintegrand = 0d0
        integral = 0d0
        dt = (tmax-tmin)/dble(itmax-1)
        t = tmin
        do it = 1, itmax
            n = dexp(t**2d0)
            argexp = n*E/Ecut+(lad/alpha*(1d0-1d0/n)-x)**2d0/
     &               (4d0*ldiff**2d0/alpha*dlog(n))

            if (n.eq.1d0) argexp = 1d35
            integrand = 2d0*dexp((1d0-s)*t**2d0)/dexp(argexp)

            if (integrand.lt.1d45) then
            else
              print *, n,E,Ecut,lad, alpha, x, ldiff, t, s
            endif
        
            if (it.eq.1) go to 2222
            integral = integral + (integrand + oldintegrand)*dt/2d0
 2222       continue
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
     &               (4d0*ldiff**2d0/alpha*dlog(n))
            integrand = n**(-1d0*s)/dsqrt(dlog(n))/dexp(argexp)*
     &                  ((2d0*q-1d0)*y**2d0+1d0)/
     &                  (1d0-y**2d0)**(q+1d0)
         
            ! This structure looks familiar...
            if (integrand.lt.1d45) then
            else
                print *, integrand, argexp
            endif

            if (i.eq.1) go to 1111
            integral = integral + (integrand+oldintegrand)*dy/2d0
 1111       continue
            oldintegrand = integrand

            y = y + dy
        enddo

        integral = k0/2d0/dsqrt(pi*alpha*b0*D0)*E**(-1d0*(1d0+s))
     &             *integral
        distrpohl = integral

  444   continue
        distrpohl = integral

        return 
      end
    
! ------------------------------------------------
! Electron distribution function, for mu > 1
! Equation (16), from Lerche & Schlickeiser (1980)
! ------------------------------------------------

      function distrmgt1(E,B, Ecut,r, eta,mu, rs, v,s)
          double precision v, tau, E, Ecut
          double precision distrmgt1
          double precision D0, B, a, b0, alpha
          double precision k0, s, integrand
          double precision x, n, argexp,pi, dn, rs, r
          double precision integral, rsc, dr
          double precision oldintegrand, nmin, dy
          double precision y,q,t, dt, dellog, dE
          double precision tmax, integrandmax
          double precision tmin, lad, ldiff
          double precision epsmax, epsmin, eps
          double precision en0, efinv, ef,  Xi
          double precision br, bt
          double precision mu
          double precision Bfield, eta

          D0 = eta*2.083d19/B
          a = 1.57d-3
          k0 = 1d0
          b0 = a*B**2d0
          alpha = 1d0
          pi = 4d0*datan(1d0)

          tau = 1d0/b0/E
          lad = v*tau
          ldiff = dsqrt(D0*E**mu*tau)

          x = (rs-r*rs)/alpha

          tmin = 0d0
          tmax = 1d0

          if (lad/ldiff.gt.30d0) then 
              efinv = a/v*B**2d0*(rs-r*rs)
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

          itmax = 1000
          oldintegrand = 0d0
          integral = 0d0
          dt = (tmax-tmin)/dble(itmax-1)
          t = tmin
          do it = 1, itmax
              n = 1d0+t**2d0
              argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     &                 (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     &                 (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)

              if (n.eq.1d0) argexp = 1d35
              integrand = 2d0*(1d0+t**2d0)**((s+mu-2d0)/(1d0-mu))
     &                    /dexp(argexp)

              if (it.eq.1) go to 2222
              integral = integral + (integrand + oldintegrand)*dt/2d0

 2222         continue
              oldintegrand = integrand
              t = t + dt
          enddo

          inmax=100
          dy = 1d0/dble(inmax)

          y = 0d0
          q = 2d0
          nmin = 2d0

          oldintegrand = 0d0
          do i = 1, inmax
              n = y/(1d0-y**2d0)**q+nmin
              argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     &          (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     &          (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)

              integrand = n**((s+mu-2d0)/(1d0-mu))/dsqrt(n-1d0)
     &                    /dexp(argexp)*
     &                    ((2d0*q-1d0)*y**2d0+1d0)/
     &                    (1d0-y**2d0)**(q+1d0)

          if (integrand.lt.1d45) then  !Prints out numbers when bad things happen
          else
              print *, integrand, argexp
          endif

          if (i.eq.1) go to 1111
          integral = integral + (integrand+oldintegrand)*dy/2d0

 1111     continue
          oldintegrand = integrand
          y = y + dy
        enddo

        integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(mu-1d0))*
     &             E**(-1d0*(mu/2d0+1d0/2d0+s))*integral

        distrmgt1 = integral

  444   continue
        distrmgt1 = integral

        return 
      end

    
! ---------------------------
! FWHM calculation subroutine
! ---------------------------

! Calculates the FWHM of the array intensitygraph, which is a function of frequency.
! Note that the radial coordinate spacing (delr) might also depend on frequency.
! This subroutine also contains the machinery to calculate the flux in different 
! regions behind the shock, broken up by whether they are part of the peak intensity 
! region or not

      subroutine FWHM(intensitygraph,inumax,irmax,delr,nugraph,width, 
     &                fluxpeak, fluxplat, gammapeak, gammaplat,imu)

          double precision xmax, xmin, maxintensity, halfint
          double precision onedint(10000), delr(40), nugraph(40)
          double precision width(40), fluxplat(40), fluxpeak(40)
          double precision intensitygraph(40,10000)
          double precision gammaplat(40), gammapeak(40), fac
          integer iplat(100), ixmax(100), ixmin(100), imu

          ! Initialize variables reused for each energy band
          imax = 1
          halfint =0d0
          do i = 1, irmax
            onedint(i) = 0d0
          enddo

          ! Compute FWHM for each energy band
          do inu = 1, inumax

            ! Copy intensitygraph to onedint
            do i = 1, irmax
              onedint(i) = intensitygraph(inu,i)
            enddo

            ! Is this initialization necessary?
            ! YES: if xmax/xmin are not found
            xmax = 0d0
            xmin = 0d0

            ! Find peak intensity, position (index), and halfmax
            maxintensity = 0d0  ! Does this line do anything?
            maxintensity = maxval(onedint)
            do i = 1, irmax
              if (onedint(i).eq.maxintensity) imax=i
            enddo
            halfint = .5d0*maxintensity

            print *, 'DEBUG: imax=',imax,'irmax=',irmax
            ! Something is terribly wrong, it's not getting the right
            ! profiles...

            ! Search points to left of max
            ! Start from i = imax-1, and move left (downhill)
            ! (I think the labeling of xmax/xmin is backwards?)
            do j = 1, imax-1  ! Could decrement i instead of using j
              i = imax-j
              if (onedint(i).lt.halfint) then
                xmax = 1d0+real(i+1)*delr(inu)  ! 1.0 + (imax-j+1)*delr
                ! This computation seems fishy
                ixmax(inu) = imax-i  ! I think this should be just `i`
                goto 993  ! break out of the loop
              endif
            enddo
            ! Edge case, if no xmax is found
            if (xmax.eq.0d0) xmax = 1d0
  993       continue

            ! Search points to right of max
            ! Start from i = imax, and move right (downhill)
            do i = imax, irmax
              if (onedint(i).lt.halfint) then
                xmin = 1d0+real(i-1)*delr(inu)
                print *, 'DEBUG, xmin=', xmin  ! THIS IS NEVER REACHED
                ! This computation also seems fishy
                ixmin(inu) = i
                goto 994  ! break out of the loop
              endif
            enddo
    
  994       continue

            width(inu) = xmax-xmin

            if (xmax.eq.0d0) then
                print *, 'Box Length Error at', nugraph(inu)/2.417989d17
            elseif (xmin.eq.0d0) then
                print *, 'Box Length Error at', nugraph(inu)/2.417989d17
            elseif ((xmax-xmin)/dabs(delr(inu)).lt.20d0) then
                print *, 'Resolution Error at', nugraph(inu)/2.417989d17
            endif

            ! fluxpeak(inu) = 0d0
            ! fluxplat(inu) = 0d0
            ! iplat(inu) = irmax !2*ixmin(inu)!+(ixmin(inu)-ixmax(inu))
            ! if (iplat(inu).lt.1) then
            !   print *, 'Make Box Larger to fit Plateau'
            ! endif
          enddo

!         ifix = 8
!         do inu = 1, inumax
!           do i = 2, ixmin(ifix)
!             fluxpeak(inu) = fluxpeak(inu) + 
!    &          (intensitygraph(inu,i,imu)+intensitygraph(inu,i-1,imu))/2d0
!           enddo

!           do i = ixmin(ifix)+1, iplat(ifix)
!             fluxplat(inu) = fluxplat(inu) +  
!    &          (intensitygraph(inu,i,imu)+intensitygraph(inu,i-1,imu))/2d0
!           enddo
!         enddo

!         fac = 2d0

!         do inu = 1, inumax-1
!           if (inu.eq.1) then
!             gammapeak(inu) =-1d0*
!    &          dlog(fluxpeak(inu+1)/fluxpeak(inu)/fac)/dlog(fac)
!             gammaplat(inu) = -1d0*
!    &          dlog(fluxplat(inu+1)/fluxplat(inu)/fac)/dlog(fac)

!           elseif (inu.eq.inumax) then
!             gammapeak(inu) =-1d0*
!    &          dlog(fluxpeak(inu)/fluxpeak(inu-1)/fac)/dlog(fac)
!             gammaplat(inu) = -1d0*
!    &          dlog(fluxplat(inu)/fluxplat(inu-1)/fac)/dlog(fac)

!           else
!             gammapeak(inu) =-1d0*
!    &          dlog(fluxpeak(inu+1)/fluxpeak(inu-1)/fac**2d0)/dlog(fac**2d0)
!             gammaplat(inu) = -1d0*
!    &          dlog(fluxplat(inu+1)/fluxplat(inu-1)/fac**2d0)/dlog(fac**2d0)
!           endif
!         enddo

          return 
      end


! ------------------------------------------------------
! Return linearly interpolated value of
! electron distribution at position r
! (Diagnostic routine to plot the calculated distribution)
!
! Cleaned 2014 July 16
! ------------------------------------------------------

      function Distgraphr(r, radtab, disttab)
          implicit none

          double precision r
          double precision radtab(1000), disttab(100,1000)
          double precision fex(100), xex(100)
          integer ir

          ! Variables to perform interpolation
          double precision slopedist, Distgraphr
          integer j

          common /xrays/ xex, fex, ir

          ! Originally disttab(30,...) was hardcoded here... why?
          ! First index of disttab is row # in Pacholczyk table
          ! Looks like we're just picking an arbitrary energy
          ! index 30 --> y=5, i.e. 5x crit. frequency
          do 17 j = 2, 1000
              if (r.gt.radtab(j)) go to 17  ! Increment radtab(j)
              slopedist = (disttab(30,j) - disttab(30,j-1))/
     &                    (radtab(j) - radtab(j-1))
              Distgraphr = disttab(30,j-1) + (r - radtab(j-1))*slopedist
              go to 18  ! Break out of loop
   17     continue
   18     continue

          return
      end


