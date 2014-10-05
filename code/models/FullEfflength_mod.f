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
! Fixed format Fortran 77 with mishmash of useful newer things
! Statement labels right aligned from column 5
! Continuation ampersands on column 6
! Exclamation points for all comments
! ----------------------------------------------------------------------


! --------------------------------------
! Just a driver for the main subroutine
! E.g. for manual fitting, checking, etc
! --------------------------------------

      program Fullefflength
          implicit none
          ! Inputs
          double precision kevs(40)
          integer inumax
          double precision B0, eta2, mu, rminarc
          integer snrflag
          double precision widtharc(40)
          double precision rminarcarray(40)
          ! rminarc array

          ! For spitting out results, optional
          integer i

          ! Common block
          double precision xex(100), fex(100)
          integer ir
          common /xrays/ xex, fex, ir

          ! Load common block first
          call readfglists()

          !print *, 'Enter B0'
          !read *, B0            ! Gauss
          !print *, 'Enter Eta2'
          !read *, eta2          ! dimensionless
          !print *, 'Enter mu'
          !read *, mu            ! dimensionless
          !print *, 'Enter rminarc'
          !read *, rminarc       ! arcsec
          !print *, 'Enter 0 for SN1006, 1 for Tycho'
          !read *, snrflag       ! binary...

          ! For timing test
          B0 = 150d-6
          eta2=1d0
          mu=1.5d0
          rminarc=20d0
          snrflag=0

          kevs(1) = 0.7d0
          kevs(2) = 1d0
          kevs(3) = 2d0
          !kevs(4) = 3d0
          inumax = 3

          ! Same as old code, array of rminarcs used in python wrapper
          rminarcarray(1) = rminarc
          rminarcarray(2) = rminarc
          rminarcarray(3) = rminarc
          !rminarcarray(4) = rminarc - 15d0

!      subroutine Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu,
!     &                            vs, v0, rs, rsarc, s, rminarc, icut,
!     &                            irmax, iradmax, ixmax)

          if (snrflag.eq.0) then
            call Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu,
     &                          5d0*1.d8, 5d0*1.d8/4d0, 2.96d19, 900d0,
     &                          2d0*0.6d0+1d0, rminarcarray, 1, 400,
     &                          100,500, 1, 0.05d0, 5d-6)
          elseif (snrflag.eq.1) then
            call Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu,
     &                            4.7d0*1.d8, 4.7d0*1.d8/4d0, 1.077d19,
     &                            240d0, 2d0*0.65d0+1d0, rminarcarray,
     &                            1, 400, 100,500, 0, 0.05d0, 5d-6)
          endif

          do i = 1, inumax
            print *, kevs(i), ' keV: ', widtharc(i)
          enddo

          stop
      end

! ------------------------------------------
! Main subroutine for generating FWHM values
! using subroutine rather than program for
! f2py wrapping, so we can fit output values
! ------------------------------------------

! Yes, I die a little everytime I write a parameter list this long

! Parameters (SN 1006 values):
!     vs = 5d0*1.d8     ! Shock veloc (only used to compute Ecut)
!     v0 = vs/4d0       ! Plasma veloc measured by Satoru et al. (cm/s)
!     rs = 2.96d19      ! Shock radius, tan(943.7 arcsec) * 2.2 kpc
!     rsarc = 900d0     ! rs in arcsec, from Green's SNR Catalog
!     s = 2d0*.6d0+1d0  ! Spectral index

!     rminarc = 60d0    ! How far back to extend model (array, now!!)
!     icut = 0          ! 1 for cut-off, 0 for no cut-off 
!     irmax = 400       ! Resolution on intensity profile
!     iradmax = 100     ! Resolution on tabulated electron distribution
!     ixmax =  500      ! Resolution along line of sight for integration

!     idamp = ...
!     ab = 0.05d0
!     Bmin = 5d-6

      subroutine Fullefflengthsub(kevs, inumax, widtharc, B0, eta2, mu,
     &                            vs, v0, rs, rsarc, s, rminarc, icut,
     &                            irmax, iradmax, ixmax,
     &                            idamp, ab, Bmin)
          implicit none
          ! Input parameters (widtharc is for output)
          integer inumax
          double precision mu, eta2, B0, ab, Bmin
          double precision kevs(inumax), widtharc(inumax)
          ! Constants, or nominally constant inputs (SNR specific)
          double precision vs, v0, rs, rsarc, s
          double precision rminarc(inumax)
          ! Settings
          integer irmax, ixmax, iradmax
          integer icut, idamp

          ! "Fundamental" constants
          double precision cm, c1, a ! Synchrotron params
          double precision nukev ! nu-kev conversion
          double precision ecut, eta  ! Derived constants
          double precision rmin(inumax) ! From rminarc(inumax)
          ! Emissivity data for common block
          double precision xex(100), fex(100)
          integer ir  ! Should label as ifmax or so, to be consistent

          ! Iterate over nu
          double precision nu, nugraph(inumax)
          integer inu
          ! Tabulate e- distribution over particle energy (1 to ir)
          ! and over radial coord ( (rsarc - rminarc) / iradmax )
          double precision disttab(100, 1000), radtab(1000)
          ! Iterate over radial coord (r, irmax pts), at each point
          ! Integrate over line of sight (x, ixmax pts)
          double precision r, delr(inumax), rgraph(10000,40)
          double precision dx, x, rho
          integer i, ix
          double precision emis1, integrand, integral, intensity
          double precision distgraph(10000),intensitygraph(inumax,10000)
          ! FWHM calculation, outputs
          double precision width(inumax)!, norms(inumax)

          ! Functions used
          !double precision distrpohl, distrmlt1, distrmgt1
          double precision emisx, Distgraphr

          common /xrays/ xex, fex, ir

! f2py directives for fitting

cf2py intent(out) widtharc
cf2py depend(inumax) widtharc

! Constants
          cm = 1.82d18         ! Constant needed for max synchrotron frequency
          c1 = 6.27e18         ! Constant needed for synchrotron emissivity
          a = 1.57d-3          ! Constant needed for synchrotron loss time
          nukev = 2.417989d17  ! 1keV photon frequency

! Calculate derived parameters

          do inu = 1, inumax
            rmin(inu) = 1d0 - rminarc(inu) / rsarc ! Scaled r, range [0,1]
          end do
          !rmin = 1d0 - rminarc/rsarc 

          ! eta = \eta*(E_h)^(1-\mu), per section 4.2 of paper
          eta =  eta2 * (2d0*nukev/(cm*B0))**(-(mu-1d0)/2d0)

          if (icut.eq.1) then
              ! Scaling checks out, and 1.602177 erg = 1 TeV
              ecut = (8.3d0*(B0/(100d-6))**(-1d0/2d0)
     &               *(vs/1d8)*1.6021773d0)**(2d0/(1d0+mu))
     &               *(1d0/eta)**(1d0/(mu+1d0))
          else
              ecut = 1d40
          endif

! _________________________________________________
! Start doing things -- perform numerical integrals

! Perform all numerical integrals and calculate intensity profiles
! in each observation band (iterate over nu, which is incremented at
! bottom of the loop).
          do inu = 1, inumax
            nu = kevs(inu) * nukev ! Convert keV to Hz
            nugraph(inu) = nu  ! Store nu in array TODO redundant

            ! Tabulate electron distribution for fixed nu
            ! generating 2-D grid as function of energy, radial position
            call distr(disttab, radtab, iradmax, nu, rmin(inu), B0,Ecut,
     &                 eta, mu, rs, v0, s, c1, idamp, ab, Bmin)

            ! For each radial position in intensity profile
            ! Integrate over line of sight to get I_\nu(r)
            ! At each step of line-of-sight integration,
            ! Compute synchrotron emissivity at radial coord
            ! rho = s**2 + r**2 (r = "viewer's" radial coord)
            ! (integrates over particle energy via emisx)

            ! Step from r=1 to r=rmin (IRMAX resolution)
            delr(inu) = -(1d0-rmin(inu))/real(irmax)
            r = 1d0
            do i = 1, irmax
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
              integral = emis1*dx ! /3d0  ! Removed division by 3
              x = dsqrt(1d0-r**2)+dx

              ! Integration to find intensity
              do ix = 2,ixmax-1
                rho = dsqrt(x**2+r**2)

                ! Integral over e- distribution!
                emis1 = emisx(rho, nu,B0,radtab,disttab)

                if (mod(ix,2).eq.0) then
                  integrand =4d0*emis1  ! Changed to 4 from 2
                else
                  integrand =2d0*emis1  ! Changed to 2 from 4
                endif

                integral = integral + integrand*dx
                x = x+dx
              enddo

              ! Endpoint - middle
              rho = dsqrt(x**2 + r**2)  ! Should be just r
              emis1 = emisx(rho, nu,B0,radtab,disttab)
              integral = integral + emis1*dx / 2d0  ! No double count

              intensity = integral
              ! as dx < 0, this check is necessary...
              if (intensity.lt.0d0) intensity = -1d0*intensity
              intensitygraph(inu,i) = intensity
  
              ! Plots, remember i indexes "viewer's" radial coord
              rgraph(i, inu) = r
              distgraph(i) = Distgraphr(r, radtab, disttab)
              ! distragph not used for anything except text output
  
              r = r+delr(inu)  ! Increment radial pos for next loop
            enddo
            ! Finished computing integrals for intensity profile

          enddo
          ! DONE computing stuff for all energy bands!

! Now, compute FWHM values and print widths

          Call FWHM(intensitygraph,inumax,irmax,delr,nugraph,width)

          do inu = 1, inumax
            widtharc(inu) = width(inu)*rsarc  ! Convert to arcsec output
          enddo

          ! TODO:  We may want intensity profiles, and/or
          ! electron distributions, for debugging.  How to get out?
          ! From best fit parameters, could get them manually...

!          open (111, file = 'Intensity.dat') 
!          open(243, file = 'distint.dat')
!          do i =1, irmax
!            norms = maxval(intensitygraph, dim = 2)
!            write (111,*) rgraph(i,1), (intensitygraph(j,i)/norms(j)
!     &                    , j = 1, inumax)
!            write (243,*) rgraph(i,1), distgraph(i)
!          enddo
!          close(111)
!          close(243)

          return
      end


! --------------------------------------------
! Read frequency-dep. 1-particle emissivity
! table to common block variables xex, fex, ir
! --------------------------------------------

      subroutine readfglists()
          implicit none
          double precision xex(100), fex(100)
          integer ir
          common /xrays/ xex, fex, ir

          open(9, file='fglists.dat', status='old')
          do 90 ir = 1, 300
              read(9,*) xex(ir), fex(ir)
              if (xex(ir).eq.0.) go to 91  ! Break at end of file
   90     continue
   91     continue
          ir = ir-1  ! Set ir = number of lines read
          close(9)

          return
      end

! ----------------------------------------------------------------
! Integrate 1-particle emissivity with particle distribution over
! radiation frequency (change dE to dy) to get emissivity j_\nu(r)
!
! Cleaned 2014 July 16
! ----------------------------------------------------------------

! j_nu(r, nu, B, radtab, disttab) is a function of radial coordinate r
! noting that r != x ... set by e- dist. functions in terms of r, not x

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

          ! If radtab(iradmax) is slightly less than 1, and r=1.0
          ! It will never find the right place to interpolate
          ! (note radtab(1000) extends well past iradmax, all zeros-ish
          ! But in this edge case, dist should be 0 anyways
          if (j.eq.1001) then
              j = maxloc(radtab, 1)
              dist = disttab(1, j)  ! Assume r=1.0 case
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

          emisx = trpsum*dsqrt(nu*B)  ! Constant prefactors skipped
          return
      end

! --------------------------------------
! Electron distribution tabulation
! Compute e- distr and return full table
! --------------------------------------

! NOTE: THIS IS THE ONLY PART OF THE FORTRAN CODE CURRENTLY IN USE

      subroutine distr(disttab, radtab, iradmax, nu, rmin, B0, Ecut,
     &                                           eta, mu, rs, v0, s, c1,
     &                                           idamp, ab, Bmin)
          implicit none
          ! Input arrays to modify
          double precision disttab(100, 1000), radtab(1000)
          ! Inputs
          integer iradmax ! Resolution for radtab
          double precision nu, rmin, B0, mu, rs, v0, s, c1  ! Input/constants
          double precision Ecut, eta ! Derived from program inputs
          ! Magnetic damping inputs
          integer idamp  ! Flag for magnetic damping
          double precision ab, Bmin  ! Magnetic damping params

          ! Extra variable -- compute energy w/ correct field
          double precision B

          ! Iteration variables
          integer irad, j
          double precision en, delrad, rad
          ! Functions to use
          double precision distrmlt1, distrmgt1, distrpohl
          ! Emissivity data for common block
          double precision xex(100), fex(100)
          integer ir ! Number of entries in Pacholczyk

          common /xrays/ xex, fex, ir

! f2py directives

cf2py intent(out) disttab

          rad = rmin
          delrad = (1d0 - rmin) / dble(iradmax-1)
          ! rad will span rmin to 1.0, inclusive

          do irad = 1, iradmax  ! Fortran-order access

            ! Compute B at radial position, if needed
            if (idamp.eq.1) then
              B = Bmin + (B0-Bmin)*dexp((rad-1d0)/ab) ! -x = (r-1)
            else
              B = B0
            endif

            do j = 1, ir
              en = dsqrt(nu/c1/B/xex(j))
              ! Don't pass B into distr*, use B0 and Bmin
              if (mu.lt.1d0) then
                disttab(j,irad) = distrmlt1(en,B0, ecut,rad,
     &                                      eta,mu, rs, v0,s,
     &                                      idamp, ab, Bmin)
              elseif (mu.gt.1d0) then
                disttab(j,irad) = distrmgt1(en,B0, ecut,rad,
     &                                      eta,mu, rs, v0,s,
     &                                      idamp, ab, Bmin)
              else
                disttab(j,irad) = distrpohl(en,B0, ecut,rad,
     &                                      eta,mu, rs, v0,s,
     &                                      idamp, ab, Bmin)
              endif
            enddo
            radtab(irad) = rad
            rad = rad + delrad
          enddo

          return
      end


! --------------------------------------------------
! Electron distribution function arguments, for all!
!
! E = particle energy, B = magnetic field
! ecut = cut-off energy of injected electrons
! s = spectral index of injected electrons
! r = radial position (scaled, 0 to 1)
! eta = \eta*(E_h)^(1-\mu), from the paper
! mu = diffusion coefficient exponent
! rs, v = shock radius, plasma velocity
! idamp = toggle for magnetic damping
! absc = damping length ab (scaled, 0 to 1)
! Bmin = minimum magnetic field
! --------------------------------------------------


! ------------------------------------------------
! Electron distribution function, for mu < 1
! Equation (15), from Lerche & Schlickeiser (1980)
! ------------------------------------------------

      function distrmlt1(E,B, Ecut,r, eta,mu, rs,v,s, idamp,absc,Bmin)
          implicit none

          ! Inputs
          double precision E, r  ! Variables to grid over
          double precision B, mu, rs, v, s  ! Input/constants
          double precision Ecut, eta ! Derived from program inputs
          integer idamp  ! Flag for damping
          double precision absc, Bmin  ! Magnetic damping stuff

          ! Derived constants etc
          double precision D0, a, k0, b0, alpha, pi
          double precision tau, lad, ldiff, x
          double precision ab, z  ! For damping

          ! Integration
          double precision efinv, ef, en0, Xi  ! For lad >> ldiff case
          double precision t, tmin, tmax, dt
          integer it, itmax
          double precision n, argexp
          double precision integral, integrand, oldintegrand
          double precision distrmlt1


          D0 = eta*2.083d19/B  ! eta (E_h)^(1-mu) * C_d / B0
          a = 1.57d-3  ! This is b in the paper
          k0 = 1d0  ! Normalization constant Q_0
          b0 = a*B**2d0  ! b * B^2 in paper
          alpha = 1d0  ! TODO: not sure what this is, scaling factor?
          pi = 4d0*datan(1d0)  ! Just pi

          tau = 1d0/b0/E  ! Synchrotron cooling time
          lad = v*tau  ! Advective lengthscale
          ldiff = dsqrt(D0*E**mu*tau)  ! Diffusive lengthscale

          if (idamp.eq.0) then
              x = (rs-r*rs)/alpha  ! Convert scaled radial coord to x
          elseif (absc.gt.1d3) then
              x = rs - r*rs  ! Let z(x) = x for extreme damping
          else
              ab = absc*rs  ! convert ab to cm
              z = rs - r*rs  ! This would normally be x
              ! variable x is z(x) in paper
              x = (z*(Bmin**2d0)/B**2d0 +
     &             ab/2d0*(B-Bmin)**2d0/B**2d0 *(1d0-dexp(-2d0*z/ab)) +
     &             2d0*ab*(B-Bmin)*Bmin/B**2d0 *(1d0-dexp(-z/ab))
     &            )/alpha
          endif

          ! Integration limits
          tmin = 0d0
          tmax = 1d0

          !If lad >> ldiff, compute easy solution and skip to end
          if (lad/ldiff.gt.30d0) then  !If lad<<ldiff the solution is easy
            if (idamp.eq.0) then
                efinv = a/v * B**2d0 * (rs-r*rs)
            else
                efinv = b0/v*x
            endif
            ef =1d0/efinv
            en0 = E/(1d0-E/ef)
            argexp = en0/ecut
            ! argexp = E/ecut / (1 - a*B**2d0*E * (rs-r*rs) / v)
            Xi = 1d0-E/ef  ! Only used for if-else structure below
        
            if (Xi.gt.0d0) then
              integral = 1./en0**s*(en0/E)**2/dexp(argexp)*8d-9
            else
              integral = 0d0
            endif

            go to 444
          endif

          ! Set up integral over t, where n = 1-t^2
          ! Using trapezoidal sum
          itmax = 1000
          oldintegrand = 0d0
          integral = 0d0
          dt = (tmax-tmin)/dble(itmax-1)
          t = tmin

          do it = 1, itmax
            n = 1d0-t**2d0

            argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     &               ( lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x )**2d0 /
     &               (4d0*ldiff**2d0*(1d0-n)/alpha) * (1d0-mu)

            if (n.eq.1d0) argexp = 1d35  ! Prevent argexp blowup

            integrand = 2d0*(1d0-t**2d0)**((s+mu-2d0)/(1d0-mu))
     &                  /dexp(argexp)
            if (it.eq.1) go to 2222  ! Skip first point
            integral = integral + (integrand + oldintegrand)*dt/2d0
 2222       continue
            oldintegrand = integrand

            t = t + dt
          enddo

          integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(1d0-mu))*
     &               E**(-1d0*(mu/2d0+1d0/2d0+s))*integral
          distrmlt1 = integral

  444     continue  ! Exit point for lad >> ldiff case
          distrmlt1 = integral

          return 
      end

! ------------------------------------------
! Electron distribution function, for mu = 1
! Equation (14), from Rettig and Pohl (2012)
! ------------------------------------------

      function distrpohl(E,B, Ecut,r, eta,mu, rs,v,s, idamp,absc,Bmin)
          !implicit none

          ! Inputs
          double precision E, r  ! Variables to grid over
          double precision B, mu, rs, v, s  ! Input/constants
          double precision Ecut, eta ! Derived from program inputs
          integer idamp  ! Flag for damping
          double precision absc, Bmin  ! Magnetic damping stuff

          ! Derived constants etc
          double precision D0, a, k0, b0, alpha, pi
          double precision tau, lad, ldiff, x
          double precision ab, z  ! for damping

          ! Integration
          double precision efinv, ef, en0, Xi  ! For lad >> ldiff case
          double precision t, tmin, tmax, dt  ! 1st variable change
          double precision y, q, nmin, dy  ! 2nd variable change
          !integer i, inmax, it, itmax
          double precision n, argexp
          double precision integral, integrand, oldintegrand

          double precision distrpohl


        D0 = eta*2.083d19/B  ! C_D = c/(3e)
        a = 1.57d-3
        k0 = 1d0
        b0 = a*B**2d0
        alpha = 1d0
        pi = 4d0*datan(1d0)

        tau = 1d0/b0/E
        lad = v*tau
        ldiff = dsqrt(D0*E*tau)

        if (idamp.eq.0) then
            x = (rs-r*rs)/alpha
        elseif (absc.gt.1d3) then
            x = rs - r*rs  ! Let z(x) = x for extreme damping
        else
            ab = absc*rs  ! convert ab to cm
            z = rs - r*rs  ! This would normally be x
            ! variable x is z(x) in paper
            x = (z*(Bmin**2d0)/B**2d0 +
     &           ab/2d0*(B-Bmin)**2d0/B**2d0 *(1d0-dexp(-2d0*z/ab)) +
     &           2d0*ab*(B-Bmin)*Bmin/B**2d0 *(1d0-dexp(-z/ab))
     &          )/alpha
        endif

        tmin = 0d0
        tmax = 1d0
        
        if (lad/ldiff.gt.30d0) then 
          if (idamp.eq.0) then
              efinv = a/v*B**2d0*(rs-r*rs)
          else
              efinv = b0/v*x
          endif
          ef =1d0/efinv
          en0 = E/(1d0-E/ef)
          argexp = en0/ecut
          ! argexp = E/ecut / (1 - a*B**2d0*E * (rs-r*rs) / v)
          Xi = 1d0-E/ef  ! Only used for if-else structure below

          if (Xi.gt.0d0) then
            integral = 1./en0**s*(en0/E)**2/dexp(argexp)*8d-9
          else
            integral = 0d0
          endif

          go to 444
        endif

        ! This integral only runs from t=0 to t=1 (n=1 to n=e)
        ! Using trapezoidal sum
        if (idamp.eq.0) then
            itmax = 1000
        else
            itmax = 5000  ! TODO TEMPORARY just copied fr Sean's code
        endif
        oldintegrand = 0d0
        integral = 0d0
        dt = (tmax-tmin)/dble(itmax-1)
        t = tmin
        do it = 1, itmax
            n = dexp(t**2d0)
            argexp = n*E/Ecut+(lad/alpha*(1d0-1d0/n)-x)**2d0/
     &               (4d0*ldiff**2d0/alpha*dlog(n))

            if (n.eq.1d0) argexp = 1d35  ! Prevent argexp blowup
            integrand = 2d0*dexp((1d0-s)*t**2d0)/dexp(argexp)

            ! Spew numbers if bad things happen?
            if (integrand.gt.1d45) then
              print *, n,E,Ecut,lad, alpha, x, ldiff, t, s
            endif
        
            if (it.eq.1) go to 2222 ! Skip first point
            integral = integral + (integrand + oldintegrand)*dt/2d0
 2222       continue
            oldintegrand = integrand

            t = t + dt
        enddo

        ! Second part of integral from n=e to n=\infty
        ! Change variables, n = y / (1 - y^2)^q + n_min
        inmax=100
        dy = 1d0/dble(inmax) ! y ranges JUST shy of y=1 (n=infty)
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
         
            ! Again, if integrand blows up, spew numbers
            if (integrand.gt.1d45) then
                print *, integrand, argexp
            endif

            if (i.eq.1) go to 1111 ! Skip first point
            integral = integral + (integrand+oldintegrand)*dy/2d0
 1111       continue
            oldintegrand = integrand

            y = y + dy
        enddo

        integral = k0/2d0/dsqrt(pi*alpha*b0*D0)*E**(-1d0*(1d0+s))
     &             *integral
        distrpohl = integral

  444   continue  ! Exit point for lad >> ldiff
        distrpohl = integral

        return 
      end


! ------------------------------------------------
! Electron distribution function, for mu > 1
! Equation (16), from Lerche & Schlickeiser (1980)
! ------------------------------------------------

      function distrmgt1(E,B, Ecut,r, eta,mu, rs,v,s, idamp,absc,Bmin)
          implicit none
          double precision E, r  ! Variables to grid over
          double precision B, mu, rs, v, s  ! Input/constants
          double precision Ecut, eta  ! Derived from program inputs
          integer idamp  ! Flag for damping
          double precision absc, Bmin  ! Magnetic damping stuff

          ! Derived constants etc
          double precision D0, a, k0, b0, alpha, pi
          double precision tau, lad, ldiff, x
          double precision ab, z  ! for damping

          ! Integration
          double precision efinv, ef, en0, Xi  ! For lad >> ldiff case
          integer i, inmax, it, itmax
          double precision t, tmin, tmax, dt  ! 1st variable change
          double precision y, q, nmin, dy  ! 2nd variable change
          double precision n, argexp
          double precision integral, integrand, oldintegrand

          double precision distrmgt1

          D0 = eta*2.083d19/B
          a = 1.57d-3
          k0 = 1d0
          b0 = a*B**2d0
          alpha = 1d0
          pi = 4d0*datan(1d0)

          tau = 1d0/b0/E
          lad = v*tau
          ldiff = dsqrt(D0*E**mu*tau)

          if (idamp.eq.0) then
              x = (rs-r*rs)/alpha
          elseif (absc.gt.1d3) then
              x = rs - r*rs  ! Let z(x) = x for extreme damping
          else
              ab = absc*rs  ! convert ab to cm
              z = rs - r*rs  ! This would normally be x
              ! variable x is z(x) in paper
              x = (z*(Bmin**2d0)/B**2d0 +
     &             ab/2d0*(B-Bmin)**2d0/B**2d0 *(1d0-dexp(-2d0*z/ab)) +
     &             2d0*ab*(B-Bmin)*Bmin/B**2d0 *(1d0-dexp(-z/ab))
     &            )/alpha
          endif

          tmin = 0d0
          tmax = 1d0

          ! Same, use easy soln for lad >> ldiff
          if (lad/ldiff.gt.30d0) then 
            if (idamp.eq.0) then
                efinv = a/v*B**2d0*(rs-r*rs)
            else
                efinv = b0/v*x
            endif
            ef =1d0/efinv
            en0 = E/(1d0-E/ef)
            argexp = en0/ecut
          ! argexp = E/ecut / (1 - a*B**2d0*E * (rs-r*rs) / v)
            Xi = 1d0-E/ef  ! Only used for if-else structure below

            if (Xi.gt.0d0) then
                integral = 1./en0**s*(en0/E)**2/dexp(argexp)*8d-9
            else
                integral = 0d0
            endif

            go to 444
          endif

          ! Set up integral over t, n = 1+t^2 (not same as before!)
          ! Using trapezoidal sum
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

              if (n.eq.1d0) argexp = 1d35 ! Catch argexp blowup
              integrand = 2d0*(1d0+t**2d0)**((s+mu-2d0)/(1d0-mu))
     &                    /dexp(argexp)

              if (it.eq.1) go to 2222 ! Skip first point
              integral = integral + (integrand + oldintegrand)*dt/2d0
 2222         continue
              oldintegrand = integrand

              t = t + dt
          enddo

          ! Set up second integral, for n=2 to n=\infty
          ! Same change of variables as for mu=1
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

          !Prints out numbers when bad things happen
          if (integrand.gt.1d45) then
              print *, integrand, argexp
          endif

          if (i.eq.1) go to 1111 ! Skip first point
          integral = integral + (integrand+oldintegrand)*dy/2d0
 1111     continue
          oldintegrand = integrand

          y = y + dy
        enddo

        integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(mu-1d0))*
     &             E**(-1d0*(mu/2d0+1d0/2d0+s))*integral

        distrmgt1 = integral

  444   continue ! Exit point for lad >> ldiff case
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

! TODO: turn this into a function?

      subroutine FWHM(intensitygraph,inumax,irmax,delr,nugraph,width)
          ! implicit none
          integer inumax
          double precision intensitygraph(inumax,10000)
          double precision xmax, xmin, maxintensity, halfint
          double precision onedint(10000), delr(inumax), nugraph(inumax)
          double precision width(inumax)
          

          ! Initialize variables reused for each energy band
          imax = 1
          halfint =0d0
          !do i = 1, irmax
          !  onedint(i) = 0d0
          !enddo

          ! Compute FWHM for each energy band
          do inu = 1, inumax

            ! Copy intensitygraph to onedint
            ! Sean notes: i=1 to irmax goes /inwards/ from shock(r=1)
            ! See main program calculation of intensity profiles
            do i = 1, irmax
              onedint(i) = intensitygraph(inu,i)
            enddo

            ! Initialize in case xmax/xmin not found
            xmax = 0d0
            xmin = 0d0

            ! Find peak intensity, position (index), and halfmax
            maxintensity = maxval(onedint)
            do i = 1, irmax
              if (onedint(i).eq.maxintensity) imax=i
            enddo
            halfint = .5d0*maxintensity

            ! Sean notes: delr(inu) is negative!
            ! xmin, xmax may be better named rmin, rmax
            ! they do NOT correspond to x as defined in paper
            ! See main program calculation of intensity profiles

            ! Search points to left of max
            ! Start at i = imax-1 and move left (down to shock front)
            do j = 1, imax-1  ! Could decrement i instead of using j
              i = imax-j  ! i = index position
              if (onedint(i).lt.halfint) then
                xmax = 1d0 + real(i+1)*delr(inu)  ! r = xmax < 1
                goto 993  ! break once xmax is found
              endif
            enddo
            ! Edge case, if no xmax is found
            if (xmax.eq.0d0) then
                xmax = 1d0
                print *, 'Rim falloff towards shock is weird?'
            endif
  993       continue

            ! Search points to right of max
            ! Start at i = imax and move right (downstream of shock)
            do i = imax, irmax
              if (onedint(i).lt.halfint) then
                xmin = 1d0 + real(i-1)*delr(inu)
                goto 994  ! break once xmin is found
              endif
            enddo
  994       continue

            width(inu) = xmax - xmin ! Yes this is positive

            if (xmax.eq.0d0) then
                print *, 'Box Length Error (xmax) at',
     &                   nugraph(inu)/2.417989d17
                ! This won't ever be reached, with edge case fix
            elseif (xmin.eq.0d0) then
                print *, 'Box Length Error (xmin) at',
     &                   nugraph(inu)/2.417989d17
                width(inu) = 1d0  ! Add-on for error detection
            elseif ((xmax-xmin)/dabs(delr(inu)).lt.20d0) then
                ! Warning, but still returns an approx value
                print *, 'Resolution Error at',
     &                   nugraph(inu)/2.417989d17
                width(inu) = 0d0  ! Add-on for error detection
            endif

          enddo

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

          ! If Distgraphr is not found...
          if (r.gt.radtab(j)) then
              print *, 'Distgraphr error (interpolation fail)'
              Distgraphr = 1d45
          endif

          return
      end


