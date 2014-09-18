	program Fullefflength

	
       double precision cm, ecut
       double precision B0, v0, rs, rsarc, alpha
       double precision nu0, nu, nukev, emis1
       double precision a, nugraph(100)
       double precision r, delr(40), rgraph(10000,40)
       double precision dx, integrand, integral
       Double Precision compratio
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
       
       integer icut, ifit, iplot
    
       common/xrays/xex, fex, ir
   


!Physical Parameters

       compratio = 4d0			!Rankine-Hugoniot
       v0 = 5d0*1.d8/compratio  !As measured by Satoru et al.
       rs = 2.96d19   			!Modified from nufit code....distance is now 2.2 kpc,
       							!so rs is found from tan(943.7'')*2.2 kpc 
       rsarc = 900d0  			!From Green's SNR Catalog...rs in arcseconds
       alpha = .6d0   			!Spectral index
       s = 2d0*alpha+1d0		
       
!Constants
		cm = 1.82d18			!Constant needed for max synchrotron frequency
       	c1 = 6.27e18			!Constant needed for synchrotron emissivity
       	a = 1.57d-3				!Constant needed for synchrotron loss time
       	nukev = 2.417989d17     !1keV photon frequency
       	
!Options
			icut = 0   !1 for cut-off, 0 for no cut-off 
			ifit = 0    !1 for the fitting version of the code
			iplot = 1   !1 for the plotting version of the code  
			
!Grid Parameters

	  irmax = 400       !resolution on intensity profile
	  iradmax = 100     !resolution on the tabulated electron distribution
	  ixmax =  500      !resolution along line of sight for integration
       
 
	open(242, file = 'widthcut2.dat')
	open(243, file = 'distint.dat')
	open(119, file = 'disttab.dat')

      

       
       print *, 'Enter B0'
       Read *, Bfield
       PRINT *, 'Enter Eta2'
       READ *, eta2
       PRINT *, 'Enter mu'
       READ *, mu
       
       eta =  eta2 * 
     c    (2d0*nukev/(cm*Bfield))**(-(mu-1d0)/2d0)
       
       if(icut.eq.1) then
       		ecut = (8.3d0*(Bfield/(100d-6))**(-1d0/2d0)*(v0*4d0/1d8)*
     c  		1.6021773d0)**(2d0/(1d0+mu))*
     c   		(1d0/eta)**(1d0/(mu+1d0))
        else
        	ecut = 1d40
        
        endif
       
 
        open(9, file='fglists.dat', status = 'old')
		open(225, file = 'FLUX.dat')
   

        do 90 ir = 1, 300

           read(9,*) xex(ir), fex(ir)
           if (xex(ir).eq.0.) go to 91
 90     continue   
 91     continue
        ir = ir-1
	
	close(9)


	   	if(ifit.eq.1) then
       			nu0 = nukev
       			inumax = 2
       	        rminarc = 30d0
       	endif
       	
       	if(iplot.eq.1) then
       			nu0 = nukev*.1d0
       			inumax = 32
       			rminarc = 110d0
       	endif
       	        rmin =1d0 - rminarc/rsarc 

		nu = nu0
	   
	    do inu = 1,inumax


       B0 = Bfield
       nugraph(inu) = nu
    
       
       
       do j = 1, ir

       
       delrad = (1d0-rmin)/dble(iradmax-1)
       rad = rmin
       
       
       do irad = 1, iradmax  !Tabulate the electron distribution


       		en = dsqrt(nu/c1/B0/xex(j))
       		
      		if (mu.lt.1d0) then
       
       			disttab(j,irad) = distrmlt1(E,B, Ecut,r, eta,mu, rs, v0,s)
       
      		elseif (mu.gt.1d0) then
       			disttab(j,irad) = distrmgt1(E,B, Ecut,r, eta,mu, rs, v0,s)
 		
 	   		else
 	   
 	   			disttab(j,irad) = distrpohl(E,B, Ecut,r, eta,mu, rs, v0,s)
 
 	   		endif
       
       		radtab(irad) = rad
       
       
       		rad = rad + delrad
       enddo
       enddo
     

  
        delr(inu) = -(1d0-rmin)/real(irmax)
	r = 1d0

	do i = 1, irmax





	integrand = 0d0
	integral = 0d0

	

	dx = -(dsqrt(1d0-r**2))/real(ixmax)

	
	
	call emisx(r,nu,B0,emis1, 
     c    radtab, disttab, imu)

	
	x = dsqrt(1-r**2)+dx
	
	integral = emis1*dx/3d0     
    
	

	
	
	
	
	!Integration to find intensity
        do ix = 2,ixmax-1
	
	
	
	rho = dsqrt(x**2+r**2)


        call emisx(rho, nu,B0, emis1, 
     c     radtab, disttab, imu)
	
	
	rgraph(i, inu) = r
	
	
	if (mod(ix,2).eq.0) then
	               
	integrand =2d0*emis1     
	
   	else
	
 	integrand =4d0*emis1 

	endif

	

	integral = integral + integrand*dx

	

	x = x+dx
	enddo
	
	

	
	intensity = integral    
	


	if (intensity.lt.0d0) intensity = -1d0*intensity




	intensitygraph(inu,i) = intensity


	Call Distgraphr(r, B0, radtab, disttab, distgraph(i))
	   
	
	r = r+delr(inu)

	enddo
	
	
		if(ifit.eq.1) then
			nu = 2d0*nu
		endif
		if(iplot.eq.1) then
			nu = nu*10d0**(.07d0)
	        rmin = 1d0-(1d0-rmin)*((10d0)**(-.25d0/16d0))
	        
	    endif

	print *, 'rmin is: ', rmin
	print *, 'ecut is: ', ecut
	enddo


	close (999)



	
	
	Call FWHM(intensitygraph,inumax,irmax,delr,nugraph,
     c   width, fluxpeak, fluxplat, gammapeak, gammaplat,imu)
     
     
     

	do n = 1, inumax
	
	print 124, nugraph(n)/nukev, 'keV', 
     c    width(n)*rsarc, 
     c          'arcseconds', 
     c    log10(width(n)/width(n-1))/
     c     log10(nugraph(n)/nugraph(n-1))
     c    ,'mnu'
     
!	write(225,*) nugraph(n)/nukev, fluxpeak(n), fluxplat(n), 
!     c   gammapeak(n), gammaplat(n)
     
	write (242, *) nugraph(n)/nukev, 
     c    width(n)*rsarc


124	format (F5.2,2x, A3, 2x, F6.2, 2x,A10, 2x, F4.2, 2x, A3)
	
	enddo



	open (111, file = 'Intensity.dat') 
	do i =1, irmax


	norms = maxval(intensitygraph, dim = 2)
	
	write (111,*) rgraph(i,1), (intensitygraph(j,i)/norms(j)
     c    ,j = 1, inumax)
     
	write (243, *) rgraph(i,1), distgraph(i)

	enddo
	
	
	close(111) 
	close(225)
	close(242)
	close(243)
	    
	    
	    

     
        

	end



      subroutine emisx(r, nu, B, emis1,
     c  radtab, disttab, imu)

      double precision en, en0, dist, argexp, oldsum
      double precision r,xold, xint1, x, rho, nu, ef
      double precision v, ecut,a, c1,cm
      double precision efinv,rs, Xi, sum, B, ab
      double precision Bfield, Bmin, const, z
      double precision fex(100), xex(100)
      double precision  br, bt, slopebr
      double precision slopebt, int, integrand
      double precision xprime, dxprime, oldintegrand
      double precision radtab(1000), disttab(100,1000)
      double precision slopedist
      
      Double precision emis1
      integer imu
      
      
      
      common/xrays/xex, fex, ir
      
      
 

C     Linearly interpolate the electron distribution
	
	
	
		do 14 j = 2, 1000

              if (r.gt.radtab(j)) go to 14
              slopedist = (disttab(1,j) - disttab(1,j-1))/
     c        (radtab(j) - 
     c            radtab(j-1))
            
              dist = disttab(1,j-1) + (r - radtab(j-1))*slopedist
             
              go to 15
	      

 14              continue
 
        
 15              continue




	
! Set up trap. rule integration for emissivity

        xint1 = 0.
        x = xex(1)         
        oldsum1 = (fex(1))*dist/x**1.5d0
        xold = x
	
	


C     Integration:
    
        do 10 ix = 2, ir
	
              
           x = xex(ix)
           do 17 j = 2, 1000

              if (r.gt.radtab(j)) go to 17
              slopedist = (disttab(ix,j) - disttab(ix,j-1))/
     c         (radtab(j) - 
     c            radtab(j-1))
            
             dist = disttab(ix,j-1) + (r - radtab(j-1))*slopedist
             
              go to 18
	      

 17             continue
 
        
 18             continue
 
	   
	   
	   
           
           sum1 = fex(ix)*dist/(x)**1.5d0
           xint1 = xint1 + (sum1 + oldsum1)*(x - xold)
	   	   oldsum1 = sum1
         
           xold = x
	   
	   
	   


 10        continue
 

           emis1 = xint1*dsqrt(nu*B)
	   
   

	

           return
           end
	   
	   
	   



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
	double precision xprime, rprime, dxprime, slobebr
	double precision  mu, muvector(1:5)
	double precision Bfield, eta, eh
	integer imu

	
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
	

	
	if (lad/ldiff.gt.30d0) then    !If lad<<ldiff the solution is easy
	
	
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
	do  it = 1, itmax
	
	n = 1d0-t**2d0

	
	
	argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     c    (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     c      (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)
	if (n.eq.1d0) argexp = 1d35
     	

	
	integrand = 2d0*(1d0-t**2d0)**((s+mu-2d0)/(1d0-mu))/dexp(argexp)
	
	
	if (it.eq.1) go to 2222
	integral = integral + (integrand + oldintegrand)*dt/2d0
2222	continue
	oldintegrand = integrand
	
	
	t = t + dt
	enddo
	
	
	integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(1d0-mu))*
     c   E**(-1d0*(mu/2d0+1d0/2d0+s))*integral

	distrmlt1 = integral
		
444	continue


	distrmlt1 = integral
	

	

	return 
	end
	
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
	double precision xprime, rprime, dxprime, slobebr
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
	do  it = 1, itmax
	
	n = dexp(t**2d0)

	
	
	argexp = n*E/Ecut+(lad/alpha*(1d0-1d0/n)-x)**2d0/
     c      (4d0*ldiff**2d0/alpha*dlog(n))
	if (n.eq.1d0) argexp = 1d35
     
	
	
	integrand = 2d0*dexp((1d0-s)*t**2d0)/dexp(argexp)
	
	
	
	if (integrand.lt.1d45) then
	else
	print *, n,E,Ecut,lad, alpha, x, ldiff, t, s
	endif
	
	if (it.eq.1) go to 2222
	integral = integral + (integrand + oldintegrand)*dt/2d0
2222	continue
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
     c      (4d0*ldiff**2d0/alpha*dlog(n))
	integrand = n**(-1d0*s)/dsqrt(dlog(n))/dexp(argexp)*
     c   ((2d0*q-1d0)*y**2d0+1d0)/
     c   (1d0-y**2d0)**(q+1d0)
     
     
	if (integrand.lt.1d45) then
	else
	print *, integrand, argexp
	endif
	
	
	if (i.eq.1) go to 1111
	integral = integral + (integrand+oldintegrand)*dy/2d0
1111	continue
	oldintegrand = integrand
	
	
	y= y + dy
	enddo
	
	
	
	integral = k0/2d0/dsqrt(pi*alpha*b0*D0)*E**(-1d0*(1d0+s))*integral

	distrpohl = integral
		
		
444	continue


	distrpohl = integral

	return 
	end
	
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
	double precision xprime, rprime, dxprime, slobebr
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
	do  it = 1, itmax
	
	n = 1d0+t**2d0

	
	
	argexp = n**(-1d0/(1d0-mu))*E/Ecut+
     c    (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     c      (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)
	if (n.eq.1d0) argexp = 1d35
	
	
     	

	
	integrand = 2d0*(1d0+t**2d0)**((s+mu-2d0)/(1d0-mu))/dexp(argexp)
     
	
	if (it.eq.1) go to 2222
	integral = integral + (integrand + oldintegrand)*dt/2d0
2222	continue
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
     c    (lad/alpha*(1d0-n**(1d0/(1d0-mu)))-x)**2d0/
     c      (4d0*ldiff**2d0*(1d0-n)/alpha)*(1d0-mu)
     
     
	integrand = n**((s+mu-2d0)/(1d0-mu))/dsqrt(n-1d0)
     c   /dexp(argexp)*
     c   ((2d0*q-1d0)*y**2d0+1d0)/
     c   (1d0-y**2d0)**(q+1d0)
     
     
	if (integrand.lt.1d45) then      !Prints out numbers when bad things happen
	else
	print *, integrand, argexp
	endif
	
	
	if (i.eq.1) go to 1111
	integral = integral + (integrand+oldintegrand)*dy/2d0
1111	continue
	oldintegrand = integrand
	
	
	y= y + dy
	enddo
	
	
	
	integral = k0/2d0/dsqrt(pi*alpha*b0*D0*(mu-1d0))*
     c   E**(-1d0*(mu/2d0+1d0/2d0+s))*integral

	distrmgt1 = integral
		


		
444	continue


	distrmgt1 = integral

	return 
	end
	
	
	
	
	



	SUBROUTINE FWHM(intensitygraph,inumax,irmax,delr,nugraph,width, 
     c   fluxpeak, fluxplat, gammapeak, gammaplat,imu)

	double precision xmax xmin, maxintensity, halfint
	double precision onedint(10000), delr(40), nugraph(40)
	double precision width(40), fluxplat(40), fluxpeak(40)
	double precision intensitygraph(40,10000)
	double precision gammaplat(40), gammapeak(40), fac
	integer	iplat(100), ixmax(100), ixmin(100), imu
	
	!Calculates the FWHM of the array intensitygraph, which is a function of frequency.
	!Note that the radial coordinate spacing (delr) might also depend on frequency.
	!This subroutine also contains the machinery to calculate the flux in different 
	!regions behind the shock, broken up by whether they are part of the peak intensity 
	!region or not
  


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

	
	Do i = 1, irmax
	if (onedint(i).eq.maxintensity) imax =i
	enddo
	
        halfint = .5d0*maxintensity
	
	do j = 1,imax-1
	i = imax-j

	if (onedint(i).lt.halfint) then
	xmax = 1d0+real(i+1)*delr(inu)
	ixmax(inu) = imax-i
	goto 993
	endif

	enddo
	if (xmax.eq.0d0) xmax = 1d0

993   Continue

	do i = imax, irmax

	if (onedint(i).lt.halfint) then
	xmin = 1d0+real((i-1))*delr(inu)
	ixmin(inu) = i
	goto 994
        endif

	enddo
	
994   Continue

	width(inu) = xmax-xmin
	
	
	if (xmax.eq.0d0) then
	print *, 'Box Length Error at', nugraph(inu)/2.417989d17
	elseif (xmin.eq.0d0) then
	print *, 'Box Length Error at', nugraph(inu)/2.417989d17
	elseif ((xmax-xmin)/dabs(delr(inu)).lt.20d0) then
	print *, 'Resolution Error at', nugraph(inu)/2.417989d17
	endif
	
c	fluxpeak(inu) = 0d0
c	fluxplat(inu) = 0d0
	
	
c	iplat(inu) = irmax !2*ixmin(inu)!+(ixmin(inu)-ixmax(inu))
	
c	if (iplat(inu).lt.1) then
c	print *, 'Make Box Larger to fit Plateau'
c	endif
	
	
	
	enddo
	
	
c	ifix = 8
	
c	do inu = 1, inumax
	
	
	
c	do i = 2, ixmin(ifix)
c	fluxpeak(inu) = fluxpeak(inu) + 
c     c   (intensitygraph(inu,i,imu)+intensitygraph(inu,i-1,imu))/2d0
c	enddo
	
c	do i = ixmin(ifix)+1, iplat(ifix)
	
c	fluxplat(inu) = fluxplat(inu) +  
c     c    (intensitygraph(inu,i,imu)+intensitygraph(inu,i-1,imu))/2d0
c	enddo
	
	
c	enddo
	

c	fac = 2d0
	
	
c	do inu = 1, inumax-1
	
	
c	if (inu.eq.1) then
	
c	gammapeak(inu) =-1d0*
c     c   dlog(fluxpeak(inu+1)/fluxpeak(inu)/fac)/dlog(fac)
c	 gammaplat(inu) = -1d0*
c     c   dlog(fluxplat(inu+1)/fluxplat(inu)/fac)/dlog(fac)
     
     
c	elseif (inu.eq.inumax) then
c	
c	 gammapeak(inu) =-1d0*
c     c   dlog(fluxpeak(inu)/fluxpeak(inu-1)/fac)/dlog(fac)
c	 gammaplat(inu) = -1d0*
c     c   dlog(fluxplat(inu)/fluxplat(inu-1)/fac)/dlog(fac)
     
     
c	else
c	 gammapeak(inu) =-1d0*
c     c   dlog(fluxpeak(inu+1)/fluxpeak(inu-1)/fac**2d0)/dlog(fac**2d0)
c	 gammaplat(inu) = -1d0*
c     c   dlog(fluxplat(inu+1)/fluxplat(inu-1)/fac**2d0)/dlog(fac**2d0)
     
c	endif
	
	
c	enddo
	
	



	return 
	end
	

	
	subroutine Distgraphr(r, B, radtab, disttab, dist)

      double precision r,dist
      double precision fex(100), xex(100)
      double precision radtab(1000), disttab(100,1000)
      double precision slopedist
      double precision B
      
      !Diagnostic routine to plot the calculated distribution
            
      
      common/xrays/xex, fex, ir

	

		do 17 j = 2, 1000

              if (r.gt.radtab(j)) go to 17
              slopedist = (disttab(30,j) - disttab(30,j-1))/(radtab(j) - 
     c            radtab(j-1))
            
             dist = disttab(30,j-1) + (r - radtab(j-1))*slopedist
             
              go to 18
	      

 17             continue
 
        
 18             continue

	
	
	return
	end



