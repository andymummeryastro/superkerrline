program wrapper
  ! gfortran amodules.f90 superkerrline.f90
    implicit none
    integer ne,i,ifl, jmax, nr, np
    parameter (nr=1000) 
    parameter (np=1000)
    parameter (ne=2**13,jmax=8)
    real Emax,Emin,ear(0:ne),photar(ne),E,dE
    real sk_param(8)
    real c_param(7)
    real pi
    pi = acos(-1.0)
    
! Set parameters
    sk_param(1) = 1.968 !a
    sk_param(2) = 75.0  !inc (deg)
    sk_param(3) = 6.4   !Eline (keV)
    sk_param(4) = 4.0   !index_i
    sk_param(5) = 3.0   !index_o
    sk_param(6) = 30.0  !r_br (if negative then in units of ISCO)
    sk_param(7) = -1.0  !r_i  (if negative then in units of ISCO)
    sk_param(8) = 50.0  !r_o
  
! Set energy grid
  Emax  = 400.0 !10.0
  Emin  = 0.01  !0.5
  do i = 0,ne
     ear(i) = Emin * (Emax/Emin)**(real(i)/real(ne))
  end do
  
! Call line profile model
  call skline(ear,ne,sk_param,ifl,photar)

! Write out
  open(99, file = 'sk_line_output.dat')
  ! write(99,*)"skip on"  
  do i = 1,ne
    E  = 0.5 * ( ear(i) + ear(i-1) )
    dE =         ear(i) - ear(i-1)
    write(99,*)E,  E**2 * photar(i) / dE
  end do
  ! write(99,*)"no no"
  
  ! Now for convolution model
  ! Read in xillver spectrum
  open(90,file='xillver_mo.dat')
  do i = 1,3
    read(90,*)
  end do
  do i = 1,ne
    read(90,*)E,dE,photar(i)
    photar(i) = photar(i) * ( ear(i) - ear(i-1) )
  end do
  close(90)
  
  ! Write out xillver spectrum
  open(98, file = 'input_spectrum.dat')
  do i = 1,ne
    E  = 0.5 * ( ear(i) + ear(i-1) )
    dE =         ear(i) - ear(i-1)
    write(98,*)E,E**2*photar(i)/dE
  end do
  ! write(99,*)"no no"
  
  ! Set parameters
  c_param(1:2) = sk_param(1:2)
  c_param(3:7) = sk_param(4:8)
  
  ! Call convolution model
  call skconv(ear,ne,c_param,ifl,photar)
  
  ! Write out
  open(97, file = 'sk_conv_output.dat')
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE =         ear(i) - ear(i-1)
     write(97,*)E,  E**2 * photar(i) / dE
  end do
  ! write(99,*)"no no"

  end program wrapper
  
! Can't use an include for compiling within XSPEC
! include 'amodules.f90'

!=======================================================================
subroutine skline(ear,ne,param,ifl,photar)
! Calculates observed disk spectrum
  use internal_grids
  implicit none
  integer ne,ifl,i,j,n
  real ear(0:ne),param(8),photar(ne)
  double precision a,inc,pi,rin,rout,mu0,index_i, index_o, r_br, r_i, r_o
  double precision rnmin,rnmax,d,rfunc,disco,mudisk,re
  double precision alpha(nro,nphi),beta(nro,nphi),dOmega(nro,nphi)
  double precision alphan(nro,nphi),betan(nro,nphi),dOmegan(nro,nphi)
  double precision g, dlgfac, dFe, E_line
  real diskline(nec),sum,E,dE
  logical needtrace
  
  pi  = acos(-1.d0)
  ifl = 1

! Parameters
  a       = dble( param(1) )                !Spin parameter
  inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
  E_line  = dble( param(3) )                !Line energy in keV
  index_i = dble( param(4) )                !inner Emissivity index
  index_o = dble( param(5) )                !outer Emissivity index
  r_br    = dble( param(6) )                !break radius (rg)
  r_i     = dble( param(7) )                !inner radius (rg)
  r_o     = dble( param(8) )                !outer radius (rg)

  ! Derived and hardwired quantities
  if (r_i .lt. 0) then 
    rin     = disco(a) * abs(r_i)
  else
    rin = max(r_i, disco(a))
  end if 

  if (r_br .lt. 0) then 
    r_br     = disco(a) * abs(r_br)
  end if 

  rout    = min( r_o, 300.0 )
  mu0     = cos(inc)

  
! Initialize
  if( firstcall )then
     firstcall = .false.
     !Define coarse internal energy grid
     Emax  = 20.0 !2.0 * E_line!! Might need to change me for XSPEC. 
     Emin  = 0.1  !0.01 * E_line!! Might need to change me for XSPEC. 
     dloge = log10( Emax / Emin ) / real(nec)
     do i = 0,nec
       earc(i) = Emin * (Emax/Emin)**(real(i)/real(nec))
     end do
     !Assign impossible initial values to previous parameters
     !Note that the only impossible parameter is now cos(inc). 
     aprev   = 10.d0
     mu0prev = 10.d0
  end if

! Set up full GR grid
  rnmax = min( 300d0, 2 * rout )          !Sets outer boundary of full GR grid
  rnmin = rfunc(a,mu0)                    !Sets inner boundary of full GR grid
  call impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
  d     = max( 1.0d4 , 2.0d2 * rnmax**2 ) !Sensible distance to BH  

! Set up `straight lines' grid
  rnmax = rout                            !Sets outer boundary of Newtonian grid
  rnmin = 300.d0                          !Sets inner boundary of Newtonian grid
  call impactgrid(rnmin,rnmax,mu0,nro,nphi,alphan,betan,dOmegan)
  
! Do the ray tracing in full GR
  needtrace = .false.
  if( abs( a - aprev ) .gt. tiny(a) ) needtrace = .true.
  if( abs(mu0 - mu0prev) .gt. tiny(mu0) ) needtrace = .true.
  mudisk = 0.d0       !razor thin disk
  if( needtrace )then
     call dGRtrace(nro,nphi,alpha,beta,mu0,a,rin,rout,mudisk,d,pem1,re1)
  end if
  aprev = a
  mu0prev = mu0

! Now calculate GR "line profile"
  diskline = 0.0
! Loop through inner relativistic grid
  do j = 1,nphi
     do i = 1,nro
       if( pem1(i,j) .gt. 0.d0 )then
          re = re1(i,j)
           if( re .gt. rin .and. re .le. rout )then
              g = dlgfac( a,mu0,alpha(i,j),re )
              !Add to line profile
              !Calculate contribution to line profile
              dFe = g**3 * dOmega(i,j) * re**(-index_i)
              if (re .gt. r_br) then 
                dFe = dFe * re**(index_i) * r_br**(-index_i) * (re/r_br)**(-index_o)
              end if 

              !Work out what bin this goes into
              n = ceiling( nec * log10(g*E_line/Emin) / log10(Emax/Emin) )
              n = max( 1 , n   )
              n = min( n , nec )
              !Add to line profile
              diskline(n) = diskline(n) + real( dFe )      
           end if
       end if
     end do
  end do

! Convert to specific photon flux to interpolate
do i = 1,nec
  diskline(i) = diskline(i) / ( earc(i) - earc(i-1) )
end do

! Rebin onto input grid
call myinterp(nec,earc,diskline,ne,ear,photar)

! Sort out edge effects and normalise
sum = 0.0
do i = 1,ne
   E  = 0.5 * ( ear(i) + ear(i-1) )
   dE = ear(i) - ear(i-1)
   if( E .lt. Emin ) photar(i) = 0.0
   if( E .gt. Emax ) photar(i) = 0.0
   sum = sum + E * photar(i) * dE
end do
photar = photar / sum

! Convert back to photons per energy bin
do i = 1,ne
  photar(i) = photar(i) * ( ear(i) - ear(i-1) )
end do
  
  return
end subroutine skline
!=======================================================================



!=======================================================================
subroutine skconv(ear,ne,param,ifl,photar)
! Superkerrline convolution model
! Parameters
! param(1) = a
! param(2) = inc (deg)
! param(3) = index_i
! param(4) = index_o
! param(5) = r_br
! param(6) = r_i
! param(7) = r_o
  implicit none
  integer ne,ifl,i,nex
  parameter (nex=2**13)
  real ear(0:ne),param(7),photar(ne)
  real earx(0:nex),Emin,Emax,dloge
  real restframex(nex),sk_param(8),linex(nex)
  real convx(nex),px(nex),p(ne)
  real E,dE
  complex FTlinex(4*nex),FTrestframex(4*nex),FTconvx(4*nex)

! Define fine internal logarithmic energy grid
  Emax  = max( 1.1*ear(ne) , 12.0 )
  Emin  = min( 0.9*ear(1)  , 1.0  )
  dloge = log10( Emax / Emin ) / real(nex)
  do i = 0,nex
     earx(i) = Emin * (Emax/Emin)**(real(i)/real(nex))
  end do
  
! Rebin input onto internal logarithmic energy grid
  do i = 1,ne
     dE = ear(i) - ear(i-1)
     photar(i) = photar(i) / dE
  end do
  call myinterp(ne,ear,photar,nex,earx,restframex)
  do i = 1,nex
     dE = earx(i) - earx(i-1)
     restframex(i) = restframex(i) * dE
  end do
  
! Calculate line profile
  sk_param(1:2) = param(1:2)   !physical parameters the same
  sk_param(3)   = earx(nex/2)  !Eline = centre of internal grid
  sk_param(4:8) = param(3:7)   !physical parameters the same
  call skline(earx,nex,sk_param,ifl,linex)

! Convolve the two
  !Fourier transform
  call pad4FFT(nex,linex,FTlinex)
  call pad4FFT(nex,restframex,FTrestframex)
  !Multiply together
  FTconvx = FTlinex * FTrestframex
  !Inverse Fourier transform
  call pad4invFFT(1e-7,nex,FTconvx,convx)
  convx = 2.0 * convx
  
! Rebin onto input array
  do i = 1,nex
    dE    = earx(i) - earx(i-1)
    px(i) = convx(i) / dE
  end do
  call myinterp(nex,earx,px,ne,ear,p)
  do i = 1,ne
     dE        = ear(i) - ear(i-1)
     photar(i) = p(i) * dE
  end do

  return
end subroutine skconv  
!=======================================================================


!=======================================================================
subroutine raytrace_grid(alpha, beta, param, g_fac, rs)
  ! Calculates observed disk spectrum
    use internal_grids
    implicit none
    integer ifl, i, j
    real param(2)
    double precision a,inc,pi,rin,rout,mu0
    double precision rnmin,rnmax,d,rfunc,disco,mudisk,re
    double precision alpha(nro,nphi),beta(nro,nphi),dOmega(nro,nphi), g_fac(nro, nphi)
    double precision alphan(nro,nphi),betan(nro,nphi),dOmegan(nro,nphi)
    double precision dlgfac, rs(nro, nphi)
    logical needtrace
    
    pi  = acos(-1.d0)
    ifl = 1
  
  ! Parameters
    a       = dble( param(1) )                !Spin parameter
    inc     = dble( param(2) ) * pi / 180.d0  !Inclination (degrees)
  
    ! Derived and hardwired quantities
    rin     = disco(a)
    rout    = 30
    mu0     = cos(inc)
    
  
  ! Set up full GR grid
    rnmax = min( 300d0, rout )              !Sets outer boundary of full GR grid
    rnmin = rfunc(a,mu0)                    !Sets inner boundary of full GR grid
    call impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
    d     = max( 1.0d4 , 2.0d2 * rnmax**2 ) !Sensible distance to BH  
  
  ! Set up `straight lines' grid
    rnmax = rout                            !Sets outer boundary of Newtonian grid
    rnmin = 300.d0                          !Sets inner boundary of Newtonian grid
    call impactgrid(rnmin,rnmax,mu0,nro,nphi,alphan,betan,dOmegan)
    
  ! Do the ray tracing in full GR
    needtrace = .false.
    if( abs( a - aprev ) .gt. tiny(a) ) needtrace = .true.
    if( abs(mu0 - mu0prev) .gt. tiny(mu0) ) needtrace = .true.
    mudisk = 0.0d0       !razor thin disk
    if( needtrace )then
       call dGRtrace(nro,nphi,alpha,beta,mu0,a,rin,rout,mudisk,d,pem1,re1)
    end if

    do j = 1,nphi
      do i = 1,nro
        if( pem1(i,j) .gt. 0.d0 )then
           re = re1(i,j)
            if( re .gt. rin .and. re .le. rout )then
               g_fac(i, j) = dlgfac( a,mu0,alpha(i,j),re )
               rs(i, j) = re
            end if 
        end if
      end do
    end do 


 return 

end subroutine raytrace_grid

!-----------------------------------------------------------------------
subroutine myinterp(nfx,farx,Gfx,nf,far,Gf)
! Interpolates the function Gfx from the grid farx(0:nfx) to the
! function Gf on the grid far(0:nf)
  implicit none
  integer nfx,nf
  real farx(0:nfx),Gfx(nfx),far(0:nf),Gf(nf)
  integer ix,j
  real fx(nfx),f,fxhi,Gxhi,fxlo,Gxlo
! Define grid of central input frequencies
  do ix = 1,nfx
     fx(ix) = 0.5 * ( farx(ix) + farx(ix-1) )
  end do
! Run through grid of central output frequencies
  ix = 1
  do j = 1,nf
     !Find the input grid frequencies either side of the current
     !output grid frequency
     f = 0.5 * ( far(j) + far(j-1) )
     do while( fx(ix) .lt. f .and. ix .lt. nfx )
        ix = ix + 1
     end do
     ix = max( 2 , ix )
     fxhi = fx(ix)
     Gxhi = Gfx(ix)
     ix = ix - 1
     fxlo = fx(ix)
     Gxlo = Gfx(ix)
     !Interpolate
     Gf(j) = Gxlo + ( Gxhi - Gxlo ) * ( f - fxlo ) / ( fxhi - fxlo )
  end do
  return
end subroutine myinterp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dlgfac(a,mu0,alpha,r)
!c Calculates g-factor for a disk in the BH equatorial plane
  implicit none
  double precision dlgfac,a,mu0,alpha,r
  double precision sin0,omega,Delta,Sigma2,gtt,gtp,gpp
  sin0   = sqrt( 1.0 - mu0**2 )
  omega  = 1. / (r**1.5+a)
  Delta  = r**2 - 2*r + a**2
  Sigma2 = (r**2+a**2)**2 - a**2 * Delta
  gtt    = 4*a**2/Sigma2 - r**2*Delta/Sigma2
  gtp    = -2*a/r
  gpp    = Sigma2/r**2
  dlgfac = sqrt( -gtt - 2*omega*gtp - omega**2.*gpp )
  dlgfac = dlgfac / ( 1.+omega*alpha*sin0 )
  return
end function dlgfac
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine dGRtrace(nro,nphi,alpha,beta,mu0,spin,rmin,rout,mudisk,d,pem1,re1)
! Traces rays in the Kerr metric for a camera defined by the impact
! parameters at infinity: alpha(nro,nphi) and beta(nro,nphi).
! Traces back to a disk defined by mudisk = cos(theta_disk), where
! theta_disk is the angle between the vertical and the disk surface.
! i.e. tan( theta_disk ) = 1 / (h/r)
! OUTPUT:
! pem1(nro,nphi)
! pem > 1: there is a solution
! pem = -1 photon goes to infinity without hitting disk surface
! pem = -2 photon falls into horizon without hitting disk surface
! re1(nro,nphi)      radius that the geodesic hits the disc
  use blcoordinate     ! This is a YNOGK module
  implicit none
  integer nro,nphi,i,j
  double precision alpha(nro,nphi),beta(nro,nphi),mu0,spin,rmin,rout,mudisk,d
  double precision pem1(nro,nphi),re1(nro,nphi)
  double precision cos0,sin0,scal,velocity(3),f1234(4),lambda,q
  double precision pem,re,mucros,phie,taudo,sigmacros      
  cos0  = mu0
  sin0  = sqrt(1.0-cos0**2)
  scal     = 1.d0
  velocity = 0.d0
  re1      = 0.0
  do i = 1,nro
    do j = 1,NPHI
      call lambdaq(-alpha(i,j),-beta(i,j),d,sin0,cos0,spin,scal,velocity,f1234,lambda,q)
      pem = Pemdisk(f1234,lambda,q,sin0,cos0,spin,d,scal,mudisk,rout,rmin)  !Can try rin instead of rmin to save an if statement
      pem1(i,j) = pem
      !pem > 1 means there is a solution
      !pem < 1 means there is no solution
      if( pem .gt. 0.0d0 )then
        call YNOGK(pem,f1234,lambda,q,sin0,cos0,spin,d,scal,re,mucros,phie,taudo,sigmacros)
        re1(i,j)    = re! Should also check mu_cross 
      end if
    end do
  end do
  return
end subroutine dGRtrace
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine impactgrid(rnmin,rnmax,mu0,nro,nphi,alpha,beta,dOmega)
! Calculates a grid of impact parameters
! INPUT:
! rnmin        Sets inner edge of impact parameter grid
! rnmax        Sets outer edge of impact parameter grid
! mu0          Sets `eccentricity' of the grid
! nro          Number of steps in radial impact parameter (b)
! nphi         Number of steps in azimuthal impact parameter (phi)
! OUTPUT:
! alpha(nro,nphi)   Horizontal impact parameter
! beta(nro,nphi)    Vertical impact parameter
! dOmega(nro,nphi)  dalpha*dbeta
  implicit none
  integer nro,nphi,i,j
  double precision rnmin,rnmax,mu0,alpha(nro,nphi),beta(nro,nphi)
  double precision dOmega(nro,nphi),mueff,pi,rar(0:nro),dlogr,rn(nro)
  double precision logr,phin
  pi     = acos(-1.d0)

  mueff = max( mu0 , 0.3d0 )
  
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    do j = 1,nphi
       domega(i,j) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
       phin       = (j-0.5) * 2.d0 * pi / dble(nphi) 
       alpha(i,j) = rn(i)  * sin(phin)
       beta(i,j)  = rn(i) * cos(phin) * mueff
    end do
  end do
  
  return
end subroutine impactgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dISCO(a)
  !ISCO in Rg 
  implicit none
  double precision a,dISCO,z1,z2
  if(abs(a).ge.1.0)then
    z1 = -(a**2.0 - 1.0)**(1.0/3.0)
    z1 = z1 * ( (1.0+abs(a))**(1.0/3.0)-(abs(a)-1.0)**(1.0/3.0))+1.0
  else
    z1 = ( 1.0 - a**2.0 )**(1.0/3.0)
    z1 = z1 * ( (1.0+a)**(1.0/3.0)+(1.0-a)**(1.0/3.0))+1.0
  end if 
	  
  z2 = sqrt( 3.0 * a**2.0 + z1**2.0 )
  if(a.ge.0.0)then
    dISCO = 3.0 + z2 - sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  else
    dISCO = 3.0 + z2 + sqrt( (3.0-z1) * (3.0 + z1 + 2.0*z2) )
  end if
  return
end function dISCO
!-----------------------------------------------------------------------

  

!-----------------------------------------------------------------------
function rfunc(a,mu0)
! Sets minimum rn to use for impact parameter grid depending on mu0
! This is just an analytic function based on empirical calculations:
! I simply set a=0.998, went through the full range of mu0, and then
! calculated the lowest rn value for which there was a disk crossing.
! The function used here makes sure the calculated rnmin is always
! slightly lower than the one required.
!
! Lack of event horizon makes r_min much smaller for a > 1.  Hard coded
! Simple lower bound for now. 
!
  implicit none
  double precision rfunc,mu0,a
  if (a .gt. 1.0) then 
    rfunc = 0.1   
  else if( a .gt. 0.8 )then
    rfunc = 1.5d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.1d0 + 5.6d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  else
    rfunc = 3.0d0 + 0.5d0 * mu0**5.5d0
    rfunc = min( rfunc , -0.2d0 + 10.0d0*mu0 )
    rfunc = max( 0.1d0 , rfunc )
  end if
  end function rfunc
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
  implicit none
  integer nro,nphi,i
  double precision rnmin,rnmax,mueff,rn(nro),domega(nro)
  double precision rar(0:nro),dlogr,logr,pi
  pi     = acos(-1.d0)
  rar(0) = rnmin
  dlogr  = log10( rnmax/rnmin ) / dble(nro)
  do i = 1,NRO
    logr = log10(rnmin) + dble(i) * dlogr
    rar(i)    = 10.d0**logr
    rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
    domega(i) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
  end do
  return
end subroutine getrgrid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine drandphithick(alpha,beta,cosi,costheta,r,phi)
!
! A disk with an arbitrary thickness
! The angle between the normal to the midplane and the disk surface is theta
! The inclination angle is i
      implicit none
      double precision alpha,beta,cosi,sini,r,phi
      double precision pi,costheta,sintheta,x,a,b,c,det
      double precision mu,sinphi
!      double precision muplus,muminus,ra,rb,rab,xplus1,xminus1,xplus2,xminus2
      pi = acos(-1.d0)
      sintheta = sqrt( 1.d0 - costheta**2 )
      sini     = sqrt( 1.d0 - cosi**2 )
      x        = alpha / beta
      if( abs(alpha) .lt. abs(tiny(alpha)) .and. abs(beta) .lt. abs(tiny(beta))  )then
        mu = 0.d0
        r  = 0.d0
      else if( abs(beta) .lt. abs(tiny(beta)) )then
        mu     = sini*costheta/(cosi*sintheta)
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      else if( abs(alpha) .lt. abs(tiny(alpha)) )then
        mu     = 1.d0
        sinphi = 0.d0
        r      = beta / ( sini*costheta - cosi*sintheta )
      else
        a      = sintheta**2 + x**2*cosi**2*sintheta**2
        b      = -2*x**2*sini*cosi*sintheta*costheta
        c      = x**2*sini**2*costheta**2-sintheta**2
        det    = b**2 - 4.d0 * a * c
        if( det .lt. 0.d0 ) write(*,*)"determinant <0!!!"
        if( beta .gt. 0.d0 )then
          mu     = ( -b + sqrt( det ) ) / ( 2.d0 * a )
        else
          mu     = ( -b - sqrt( det ) ) / ( 2.d0 * a )
        end if
        sinphi = sign( 1.d0 , alpha ) * sqrt( 1.d0 - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      end if
      phi = atan2( sinphi , mu )
      return
      end subroutine drandphithick
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine pad4FFT(ne,photar,padFT)
! Takes spectrum photar(1:ne), pads out with zeros to make it a length
! of 4*ne, and Fourier transforms to padFT(1:4*ne), which is a function
! of 1/E
  implicit none
  integer ne,i
  real photar(ne),padphot(4*ne)
  complex padFT(4*ne)

! Pad out the array
  padphot = 0.0
  do i = 1,ne
     padphot(i+2*ne) = photar(i)
  end do

! Call the energy Fourier transform code
  call E_FT(4*ne,padphot,padFT)

  return
end subroutine pad4FFT      
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine pad4invFFT(dyn,ne,padFT,conv)
! Takes padFT(1:4*ne), and zero-padded function of 1/E and inverse
! Fourier transforms to conv(ne), which is a non-zero padded spectrum
  implicit none
  integer ne,i
  real dyn,conv(ne),padconv(4*ne),photmax
  complex padFT(4*ne)

! Inverse Fourier transform padded FT
  call E_invFT(4*ne,padFT,padconv)

! Populate output array
  photmax = 0.0
  do i = 1,ne
    conv(i) = padconv(i+5*ne/2)
    photmax = max( photmax , conv(i) )
  end do

! Clean any residual edge effects
  do i = 1,ne
    if( abs(conv(i)) .lt. abs(dyn*photmax) ) conv(i) = 0.0
  end do

  return
end subroutine pad4invFFT
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine E_invFT(nex,cc,conv)
! Takes the complex array cc(1:nex), which is a function of 1/E
! and inverse Fourier transforms to get back a real spectrum as a
! function of E, conv(1:nex)
  implicit none
  integer nex,i
  real conv(nex),cdata(2*nex)
  complex cc(nex)

! Put back into four1 style arrays
  do i = 1,nex
    cdata(2*i-1) =  real( cc(i) )
    cdata(2*i  ) = aimag( cc(i) )
  end do
      
! Then transform back
  call ourfour1(cdata,nex,1)
      
! Move arrays back into original format
  !-ve frequencies
  do i = 1,nex/2-1
    conv(i) = cdata(2*i+nex+1)
  end do
  !DC component
  conv(nex/2) = cdata(1)
  !+ve frequencies
  do i = nex/2,nex
    conv(i) = cdata(2*i-nex+1)
  end do
  return
end subroutine E_invFT
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine E_FT(nex,photarx,bc)
! Takes the real array photarx(1:nex), which is a spectrum as a
! function of photon energy E and Fourier transforms to bc(1:nex),
! which is complex and a function of 1/E.
! Uses FFTs, so nex must be a power of 2.
! Uses the inverse transform of four1.
  implicit none
  integer nex,i
  real photarx(nex)
  real bdata(2*nex)
  complex bc(nex)

! Move arrays into arrays for four1
  bdata = 0.0
  !-ve frequencies
  do i = 1,nex/2-1
    bdata(2*i+nex+1) = photarx(i)
  end do
  !DC component
  bdata(1) = photarx(nex/2)
  !+ve frequencies
  do i = nex/2,nex
    bdata(2*i-nex+1) = photarx(i)
  end do
      
! Now do the inverse FFT
  call ourfour1(bdata,nex,-1)
      
! Now put into complex arrays
  do i = 1,nex
    bc(i) = complex( bdata(2*i-1) , bdata(2*i) ) / sqrt(float(nex))
  end do

  return
end subroutine E_FT
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      SUBROUTINE ourfour1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!-----------------------------------------------------------------------

      
