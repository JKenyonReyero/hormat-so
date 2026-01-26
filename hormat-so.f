c+mu-----------------------------------------------------------------------
      module potential_parameters
      Save
      integer nr2,ploop
      real *8 v0,rr,ar,rcoul,vcoul,Wv,rrv,av,Ws,rrs,as,Vso,rrso,aso
      real*8 beta,cnl,drwf
      complex *16, dimension(:),allocatable :: uloc0,uloc1,d1uloc,d2uloc
     ,           , perey  
      end module potential_parameters
c-----------------------------------------------------------------------
      module temporary
      Save
      integer, dimension(:),allocatable :: ipiv
c     real *8, dimension(:),allocatable :: wk
      complex*16, dimension(:),allocatable :: wk
      end module temporary
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
C     Main  program starts
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      use temporary
      use potential_parameters
      implicit real *8 (a-h,o-z)
      real *8,allocatable, dimension(:)  :: phin0,g,f,r2int,coul_0,angle
     , ,f2,r4int,ho,flog,fc,gc,dfc,dgc,coulphas,f21,x1,x2,hbeta,ho_r
c     real *8,allocatable, dimension(:,:) :: vlsmat1
      real *8,allocatable, dimension(:,:) :: fcr,gcr,pleg,pleg1,vlsmat
      complex *16,allocatable, dimension(:,:,:) :: d2mat
      real *8 lnab,xsection(0:2)
      complex*16,allocatable,dimension(:) :: vector_C
      complex*16,allocatable,dimension(:,:) :: vlint,dmat,vlsmat1
      complex*16 im1,smatrix,delta,S0,Sl,Slm,Slp,dchiext,wf,Vloc(0:2)
     ,          ,hm,hp,dhm,dhp,chi,chip,rmatrix,famp(180,0:2),fcoul
     ,          ,gamp(180,0:2),Vls1
      dimension eta_store(1),plegndr(180),plegndrm(180)
      logical need_wfs
      pi=3.141592653589793d0
      sqrpi=sqrt(pi)
      sqr4pi=sqrt(4*pi)
      im1=(0.d0,1.d0)
      h2c=20.9008
c     open(unit=104,file='rmatrix.n',status='old',position='append')
c     open(unit=105,file='smatrix-nmax',status='old',position='append')

      ! Reading information about colliding nuclei

      print *,' enter target mass, target charge, projectile mass, ',
     &        ' projectile  charge'
      read *,at,zt,ap,zp
      fmu=at*ap/(at+ap)
      if(ap==1) jp=1

      print*,' enter lab energy of the projectile in MeV'
      read*,elab
      ecm=elab*at/(at+ap)
      print *,' elab = ',elab,' e_c.m. = ',ecm
      q=sqrt(fmu*ecm/h2c)
      print *,' q = ',q
      eta=1.43997*zt*zp*fmu/(2*h2c)/q
      print *,' eta = ',eta
      eta_store=eta

      ! Reading information about  interaction potentials

      print*,'enter real optical potential v0,r0,ar,rc,beta'
      read *,v0,r0,ar,rc,beta
      rr=r0*at**(1.d0/3.d0)
      rcoul=rc*at**(1.d0/3.d0)
      vcoul=1.43997*zt*zp
      print*,'enter imaginary potentials Wv,rv,av,Ws,rs,as'
      read *, Wv,r0v,av,Ws,r0s,as
      rrv=r0v*at**(1.d0/3.d0)
      rrs=r0s*at**(1.d0/3.d0)
      print*,'enter spin-orbit potential Vso,rso,aso'
      read *, Vso,rso,aso
      rrso=rso*at**(1.d0/3.d0)
      if(vso==0) jp=0
      print  *,'nonlocality range is ',beta
      an=at/(at+1)*beta**2/2.d0*931.51/197.32858**2
      cnl=beta**2/16.d0
      if(beta==0) npot=0
      if(beta/=0) npot=2

      ! How many partial waves do you need?

      print *,' enter orbital momenta  Lmin and Lmax'
      read *,lmin,lmax
      print*,' using Lmin and Lmax as ',lmin,lmax

      ! Reading information about parameters of the basis

      print *,' Enter channel radius a'
      read*,a
      print  *,' a = ',a
      a2=a**2
      print *,' Enter harmonic oscillator radius b'
      read*,b
      print  *,' b = ',b
      b32=1.d0/b**1.5d0
      h32=1.d0/(1+(beta/(2*b))**2)**1.5d0
      hbb=1-2*beta**2/(beta**2+4*b**2)

      ! Reading information about harmonic oscillator basis size

      print *,' Enter max N for HO-expansion'
      read*,nmax
      print  *,' nmax = ',nmax

      ! Are wave functions needed?

      need_wfs=.true.
c     need_wfs=.false.
      rmax=20.
      drwf=0.1
      nr1=a/drwf
      nr2=rmax/drwf

      ! Writting top comment to TwoFNR wfn files
      write(211,*) '# Channel distorted waves'
      write(212,*) '# Channel distorted waves'
      write(213,*) '# Channel distorted waves'

      ! Now starting the caclulations

      Nm=nmax+1
      ntotmax=2*nmax+lmax+1

      erfab=erf(a/b)
      expab=exp(-a**2/b**2)
      ab=a/b
      ab2=ab**2
      lnab=log(ab)
      b2=b**2

      h2m=h2c/fmu
      h2ma=h2m*a
      hmb=2*h2m/b**2
      qa=q*a
      print *,' q*a = ',qa


      ! calculate and store HO wave functions phi_{n0} at r = 0

      allocate(ho(0:ntotmax),phin0(0:ntotmax),hbeta(0:ntotmax))
      if(need_wfs) allocate(ho_r(0:ntotmax))
      call phi_ho(ntotmax,0,0.d0,b,phin0)
      hbeta=phin0*h32*[(hbb**i,i=0,ntotmax)]
c     do i=0,ntotmax
c     print *,i,hbeta(i),phin0(i),hbb**i
c     enddo

      ! Calcuate Int(x^(2+2m)Exp(-x^2)) on dx from 0 to  a/b and then
      ! multiply by(-2)^m/m!/(m+1/2)! using recurrence relations
      ! Store as a function of m in array f

      imax=ntotmax
      allocate(f(0:imax),g(0:imax),f2(0:imax),flog(0:imax))
      f(0)=0.5*erfab-exp(-ab2+lnab)/sqrpi
      g(0)=ab**3*exp(-ab2)*2/sqrpi
      do i=1,imax
       denom=dfloat(i)*(i+0.5d0)
       f(i)=-((2*i+1)*f(i-1)-g(i-1))/denom
       g(i)=-2*ab2*g(i-1)/denom
       f2(i-1)=-i*(2*i+1.d0)*f(i)/4.d0
      enddo
      denom=dfloat(imax+1)*(imax+1.5d0)
      fi=-((2*imax+3)*f(imax)-g(imax))/denom
      f2(imax)=-fi*(imax+1)*(2*imax+3.d0)/4.d0


      allocate(r2int(0:ntotmax),r4int(0:ntotmax),vlint(0:ntotmax,0:5))
      r2int=f(0)
      r4int=f2(0)

      pn=sqrt(2.d0*gamma(1.5d0))
      r2int=f(0)!*pn
      pnl=log(pn)
      do i=1,imax,2
      f(I)=-f(I)
      enddo
      flog=log(f)

      do n=0,ntotmax
         bc=n
         if(n/=0) bcl=log(dfloat(n))
c        r2int(n)=r2int(n)*pn
         do m=1,n
c        if(n==34) print *,n,m,(-1)**m*(f(m)*bc*pn**0)
c        print *,n,m,f(m),bc,pn,(f(m)*bc*pn**0),r2int(n)
            r2int(n)=r2int(n)+(-1)**m*f(m)*bc!*pn
            r4int(n)=r4int(n)+f2(m)*bc
c           print *,n,m,flog(m),bcl,pnl,exp(flog(m)+bcl+pnl),r2int(n)
c           r2int(n)=r2int(n)+(-1)**m*exp(flog(m)+bcl+pnl)
            bc=bc*(n-m)
            bcl=bcl+log(dfloat(n-m))
         enddo
         r2int(n)=r2int(n)*pn
         r4int(n)=r4int(n)*pn
         write(511,*)n,r2int(n)/b32!,r4int(n)
c        print *,'--',n,r2int(n)!,r4int(n)
c        print *,'--',n,r4int(n)!,r4int(n)
         if(n/=0) pn=pn*sqrt((1.5d0+n)/dfloat(n+1))
         if(n==0) pn=sqrt(2.d0*gamma(2.5d0))
         pnl=log(pn)
      enddo

      r2int=r2int/b32
      r4int=r4int/b32

      print*,'-------------------'

      gn=sqrt(gamma(1.5d0)/2.d0)/b32
      do n=0,ntotmax
      sum_k=0
      bagam=1.d0/gamma(1.5d0+n)
        do k=n,0,-1
        sum_m=0
c       bc=(-2.d0)**n
        bc=1.d0
          do m=n,k,-1
          sum_m=sum_m+bc
c         print*,'bc=',n,k,m,bc
          bc=-bc/2.d0*dfloat(m+0)/dfloat(n-m+1)
          enddo
        sum_k=sum_k+bagam*sum_m
c         print*,'mmmm',n,k,sum_m,bagam,sum_k
        bagam=bagam/ab2*(k+0.5d0)
        enddo
c       print*,n,k,gn,sum_k,bagam
c     print *,    (-1)**n*erfab,exp(-ab2+(2*n+1)*lnab)*sum_k
      r2int(n)=gn*(-1)**n*(erfab-1.d0/sqrt(2.d0)*exp(-ab2+(2*n+1)
     ,        *(lnab+log(sqrt(2.d0))))*sum_k)
c     r2int(n)=gn*((-1)**n*erfab-exp(-ab2+(2*n+1)*lnab)*sum_k)
c     r2int(n)=gn*((-1)**n*erfab-exp(-ab2+(2*n+1)*lnab+sign(1.d0,sum_k)
c    ,        *log(abs(sum_k))))
      gn=gn*sqrt((n+1.5d0)/dfloat(n+1))
      write(512,*) n,r2int(n)
      enddo

      ! Calculate contribution from point-Coulomb potential (if needed)

      allocate(coul_0(0:ntotmax),f21(0:ntotmax))
      coul_0=0
      if(zt*zp/=0) then
         f21(0)=1.d0/gamma(1.5d0)**1
         f21(1)=-1.d0/3.d0/gamma(1.5d0)**1
         do i=1,ntotmax-1
         f21(i+1)=(i*f21(i-1)-0.5d0*f21(i))/(1.5d0+i)
         enddo
              gn=sqrt(gamma(1.5d0)/2.d0)*sqrt(b)
              do n=0,ntotmax
              gamm=0
              cnm=1.d0/gamma(1.5d0)
              sum_m=0
              do m=0,n
              sum_i=0
              bagam=1.d0  
              do  i=0,m
              sum_i=sum_i+bagam
              bagam=bagam/ab2*dfloat(m-i)
              enddo
              cm=0-1*exp(-ab2+2*m*lnab-gamm)*sum_i
              gamm=gamm+log(dfloat(m+1))
              sum_m=sum_m+cm*cnm
              cnm=-2*cnm/(m+1.5d0)*dfloat(n-m)
              enddo
              coul_0(n)=(1*f21(n)+1*sum_m)*gn
              gn=gn*sqrt((n+1.5d0)/dfloat(n+1))
              enddo
      endif

      
      ! Calculate Coulomb phases

      allocate(coulphas(lmax+1))
      call sigma(lmax+1,eta_store(1),sigma0,coulphas)

      print*,'-------------------'

      allocate(dmat(0:nmax,0:nmax),wk(Nm),ipiv(Nm))
      allocate(d2mat(0:nmax,0:nmax,0:2))
      allocate(vlsmat(0:nmax,0:nmax))
      allocate(vlsmat1(0:nmax,0:nmax))
      dmat=0
      d2mat=0
      allocate(vector_C(0:nmax))

      allocate(fc(lmin:lmax),gc(lmin:lmax)
     &                      ,dfc(lmin:lmax),dgc(lmin:lmax))
      call coulfg(qa,eta,dfloat(lmin),dfloat(lmax),
     &                                fc,gc,dfc,dgc,1,0,ifail)
c     print *,' fc = ',fc
c     print *,' gc = ',gc
      if(ifail/=0) then
      print *,' in coulfg  ifail  = ',ifail
      stop
      endif


      if(need_wfs) then
        allocate(fcr(nr1:nr2,lmin:lmax),gcr(nr1:nr2,lmin:lmax))
        allocate(x1(lmin:lmax),x2(lmin:lmax))
        do  i=nr1,nr2
        r=drwf*dfloat(i)
        qr=q*r
        call coulfg(qr,eta,dfloat(lmin),dfloat(lmax),
     &                    fcr(i,:),gcr(i,:),x1,x2,2,0,ifail)
        if(ifail/=0) then
        print *,' in coulfg  ifail  = ',ifail
        stop
        endif
        enddo
      endif

      call local_integrals(ntotmax,a,b,Vlint,ecm,an)

      do i=0,ntotmax
c     test_v=phin0(i)*sqrt(pi)/4.d0*v0*rr**3
c    ,      /(1+(rr/b)**2)**1.5d0*(1-2*(rr)**2/(rr**2+b**2))**i
c     print *,i,Vlint(i,0),Vlint(i,1),test_v
c     vlint(i,1)=test_v
      enddo

      if(2==1)stop


      call facl(1.d0)

c     do i=0,ntotmax
c     print*,'coul_0', i,coul_0(i)
c     enddo

      coul_0=coul_0*vcoul
      Vlint(:,2)=Vlint(:,2)+coul_0(:)

c     do i=0,npot
c     vlint(:,i)=vlint(:,i)+coul_0(:)
c     enddo


c     do i=0,ntotmax
c     print*,'vlint coul_0', i,vlint(i,2) !,coul_0(i),Vlint(i,2)+coul_0(i)
c     print*,'phin0 ', i,phin0(i)
c     enddo

      th0=1
      theta_max=180
      dtheta=1
      nangle=(theta_max-th0)/dtheta+1
      allocate(angle(nangle))
      angle=[(th0+i*dtheta,i=0,nangle-1)]
      allocate(pleg(nangle,lmin:lmax))
      call calculate_Legendre(lmin,lmax,0,angle,nangle,pleg)
      allocate(pleg1(nangle,lmin:lmax))
      call calculate_Legendre(lmin,lmax,1,angle,nangle,pleg1)


      if(2==1) stop

      famp=0
      gamp=0

      do L=lmin,lmax


      ! Construct  NxN matrices for H+L(0) in HO-basis for each L  

      call phi_ho(Nm,l,a,b,ho)

      cl=(-1)**L/sqrt(2*l+1.d0)
      onorm=0
      do N = 0,nmax
      !                      derivative of phi
      dphi=(dfloat(2*N+L)/a-a/b2)*ho(N)
      if(N /=0) dphi=dphi-2.d0/a*(N+L+0.5d0)*ho(N-1)
     &              *sqrt(dfloat(n)/(n+l+0.5d0)) 

      do Np = 0,nmax
c     do Np = 0,N    

      onorm=0
      r2=0
      Vloc=0
      Vc=0
      Vls=0
      Vls1=0
         do n1=0,N+Np+L
         n2=N+Np+L-n1
         talmi=trb8(n1,0,n2,0,0,Np,L,N,L)
c        print *,l,N,Np ,n1,n2,talmi
         onorm=onorm+talmi*phin0(n1)*r2int(n2)
         r2   =r2   +talmi*phin0(n1)*r4int(n2)
         Vloc(0) =Vloc(0) +talmi*phin0(n1)*Vlint(n2,0)
         Vloc(1) =Vloc(1) +talmi*phin0(n1)*Vlint(n2,4)
         Vloc(2) =Vloc(2) +talmi*hbeta(n1)*Vlint(n2,1)
         Vc =Vc+talmi*phin0(n1)*Vlint(n2,2)
         Vls=Vls+talmi*phin0(n1)*Vlint(n2,3)
         Vls1=Vls1+talmi*phin0(n1)*Vlint(n2,5)
c        print *,N,Np ,phin0(n1)*talmi*Vlint(n2)
c    ,     ,hbeta(n1)*talmi*Vlint(n2)
c        print *,Np,N,n1,n2,talmi,Vlint(n2),hbeta(n1)
         enddo
       tkin=(onorm*(2*N+L+1.5d0)-0.5d0*r2)*hmb
       bloch=h2m*(a*ho(N)*ho(Np)+a2*ho(Np)*dphi)
c      d2mat(Np,N,:)=(tkin+Vloc(:)+Vc-Ecm*onorm)*cl+bloch
       d2mat(Np,N,:)=(tkin+Vloc(:)+Vc-Ecm*onorm)*cl
c      print *,N,Np,tkin,Vloc(npot),Vc,onorm !bloch,d2mat(Np,N,npot),bloch+d2mat(Np,N,npot)
       d2mat(Np,N,:)=d2mat(Np,N,:)+bloch
c      d2mat(Np,N,:)=(tkin+Vloc(:)+Vc-Ecm*onorm)*cl
c      d2mat(N,Np,:)=d2mat(Np,N,:)
c       bloch1=h2m*(a*ho(N)*ho(Np)+a2*ho(Np)*dphi)
c      bloch2=h2m*(a*ho(Np)*ho(N)+a2*ho(Np)*dphi)
c      d2mat(Np,N,:)=bloch
c      print *,Np,N,bloch,d2mat(Np,N,npot)
c      d2mat(Np,N,:)=0
       vlsmat(Np,N)=Vls*cl
       vlsmat1(Np,N)=Vls1*cl
c      d2mat(N,Np,:)=d2mat(Np,N,:)
c      vlsmat(N,Np)=vlsmat(Np,N)
c      vlsmat1(N,Np)=vlsmat1(Np,N)
      enddo


c     print*,'d2mat',d2mat(:,:,npot)
c     do Np=0,nmax
c      bloch=h2m*(a*ho(N)*ho(Np)+a2*ho(Np)*dphi)
c      d2mat(Np,N,:)=d2mat(Np,N,:)+bloch
c      print *,Np,N,bloch,d2mat(Np,N,npot)
c     enddo
c     print*,'d2mat',d2mat(:,:,npot)
      enddo

c     print*,'d2mat',d2mat(:,:,npot)

      do ipot=0,npot
c     do ipot=npot,npot
      do j=abs(2*l-jp),2*l+jp,2
       so_factor=(0.25*j*(j+2)-l*(l+1)-jp*(jp+2)*0.25d0)/2.d0
       dmat=d2mat(:,:,ipot)+vlsmat(:,:)*so_factor
       if(an/=0.and.(ipot==1.or.ipot==0))
     , dmat=d2mat(:,:,ipot)+vlsmat1(:,:)*so_factor

c     print*,'dmat',dmat(:,:)
c     do ipot=0,npot

c     dmat=d2mat(:,:,ipot)

      !  Invert matrix dmat. 
      !  Original  dmat is replaced by its inverse in the process.

c     call dgetrf(Nm,Nm,dmat,Nm,ipiv,info)
      call zgetrf(Nm,Nm,dmat,Nm,ipiv,info)
              if(info/=0) then
              print *,'LU factorization failed with info = ',info
              stop
      endif
c     call dgetri(Nm,dmat,Nm,ipiv,wk,Nm,info)
      call zgetri(Nm,dmat,Nm,ipiv,wk,Nm,info)
              if(info/=0) then
              print *,'matrix inversion failed with info = ',info
              stop
      endif

c     print*,'dmat',dmat(:,:)
c     print  *,'ho',ho

      rmatrix=h2m*a*dot_product(ho(0:nmax),matmul(dmat,ho(0:nmax)))

      !  Obtain S-matrix from R-matrix

      hm=(gc(l)-im1*fc(l))
      hp=(gc(l)+im1*fc(l))
      dhm=(dgc(l)-im1*dfc(l))
      dhp=(dgc(l)+im1*dfc(l))
      Slp=(dgc(l)+im1*dfc(l))/(gc(l)+im1*fc(l))
      Slm=(dgc(l)-im1*dfc(l))/(gc(l)-im1*fc(l))

      smatrix=hm/hp*(1.d0-rmatrix*qa*Slm)/(1.d0-rmatrix*qa*Slp)
      print *,l,j,' smatrix = ',smatrix

      ! Obtain phase shift delta from S-matrix

      delta=log(smatrix)/2.d0/im1

      write(101,*)  l,j,rmatrix
      write(102,*)  l,j,smatrix,abs(smatrix)
      write(103,*)  l,j,real(delta)*180.d0/pi,aimag(delta)*180./pi
c     if(ipot==2) write(104,*)  nmax,abs(rmatrix),l,j
      if(ipot==2) write(105,*)  nmax,abs(smatrix),l,j

      ! construct wf in channel L if needed 

      if(need_wfs) then

      cph=exp(im1*coulphas(l+1))
      dchiext=0.5d0*im1*(dhm-smatrix*dhp)
      dchiext=dchiext*cph   ! multiplying by Coulomb phase

      ! Find expansion coefficients C:

      vector_C=h2m*a*matmul(dmat,ho(0:nmax))*dchiext
      ! Testing writting for twoFNR remove do loop later
      do ploop=1,2
      if(l==0) write(211+ipot,*) '# L=',l,' # J=',0.0d0
      if(l.ge.1) write(211+ipot,*) '# L=',l,' # J=',0.5d0*j
      write(211+ipot,960) 0.d0,0.d0,0.d0
        do i = 1,nr1
        r=i*drwf
        call phi_ho(Nm,l,r,b,ho_r)
        wf=dot_product(vector_C,ho_r(0:nmax))
        if(ipot==1) wf=wf*perey(i)
        write(111+ipot,*) r,q*r*real(wf),q*r*aimag(wf)
        write(211+ipot,960) r,q*r*real(wf),q*r*aimag(wf)
c       write(112+npot,*) r,q*r*abs(wf)
        enddo
        write(211+ipot,960) 0.d0,0.d0,0.d0
        do i = nr1+1,nr2
        r=i*drwf
        wf=(fcr(i,l)+0.5d0*im1*(1-smatrix)*(gcr(i,l)+im1*fcr(i,l)))/r/q
        wf=wf*cph
        if(ipot==1) wf=wf*perey(i)
        write(111+ipot,*) r,q*r*real(wf),q*r*aimag(wf)
        write(211+ipot,960) r,q*r*real(wf),q*r*aimag(wf)
c       write(112+npot,*) r,q*r*abs(wf)
        enddo
        write(111+ipot,*) ' '
c       write(112+npot,*) ' '
        enddo !ploop test

      endif
  960 format(f22.4,3d18.8)

      ! construct contribution ofL-th partial wave to xsec

      if(jp==0) then
      famp(:,ipot)=famp(:,ipot)+(2*l+1)*(smatrix-1)/(2*im1*q)
     ,       *pleg(:,l)*exp((0.d0,1.d0)*2.d0*coulphas(l+1))
      endif

      if(jp==1) then
      famp(:,ipot)=famp(:,ipot)+(l+1-(2*l+1-j)/2)*(smatrix-1)/(2*im1*q)
     ,       *pleg(:,l)*exp((0.d0,1.d0)*2.d0*coulphas(l+1))
       if(l/=0) gamp(:,ipot)=gamp(:,ipot)
     ,                      +(-1)**((2*l+1-j)/2)*(smatrix-1)/(2*q)
     ,       *pleg1(:,l)*exp((0.d0,1.d0)*2.d0*coulphas(l+1))
      endif



      enddo   ! ipot
      enddo   ! jp   

      enddo   ! ending loop on L


      do ith=1,180
      th=(dfloat(ith)*pi/180.)
      fcoul=-eta/(2*q)/(sin(th/2.d0))**2 * exp((0.d0,1.d0)
     ,      *(-eta*log((sin(th/2.d0))**2)+2.d0*sigma0))
      xsecruth=(eta/(2*q))**2/(sin(th/2.d0))**4
      xsecruth=xsecruth*10 !  convert to mb
      xsection(0:npot)=abs(fcoul+famp(ith,0:npot))**2
     ,                +abs(gamp(ith,0:npot))**2
      xsection=xsection*10 !  convert to mb
c     write(123,*) ith,xsection(0:npot)
      write(124,*) ith,xsection(0:npot)/xsecruth
c     write(124,*) ith,xsection(2:2)/xsecruth
c     write(124,*) ith,xsection(npot:npot)/xsecruth
      enddo


c     CALL setrteopts ('cpu_time_type=usertime') 
c     CALL stallingloop 
      CALL CPU_TIME(T2) 
      print *, 'The total user time from the start of the program is',
     ,        T2, 'seconds'
      end

      subroutine phi_ho(nmax,l,r,b,howf)
      implicit real *8 (a-h,o-z)
      dimension howf(0:nmax)        
      b32=1.d0/b**1.5d0
      rb2=(r/b)**2
      exprb=exp(-0.5d0*(r/b)**2)        
      howf(0)=sqrt(2.d0/gamma(l+1.5d0))
      howf(1)=(1.5d0+l-rb2)*sqrt(2.d0/gamma(l+2.5d0))
      do n=1,nmax-1
      rc1=sqrt(dfloat(n+1)/(n+l+1.5d0))
      rc2=rc1*sqrt(dfloat(n)/(n+l+0.5d0))
      howf(n+1)=howf(n)*(2*n+l+1.5d0-rb2)*rc1-(n+l+0.5d0)*howf(n-1)*rc2
      howf(n+1)=howf(n+1)/dfloat(n+1)
      enddo
      howf=howf*exprb*b32
      if(l/=0) howf=howf*(r/b)**l
      return
      end
          


      subroutine local_integrals(nmax,a,b,Vlint,ecm,an)
      use potential_parameters
      implicit real *8 (a-h,o-z)
      complex *16 pot,pot0,uloc,Vlint,denom,pt,delta_U
      dimension ho(0:nmax),vlint(0:nmax,0:5),pot0(0:1)
      sqrt2=sqrt(2.d0)
      dr=drwf*0.5

      vlint=0
      npoint=a/dr*3*5
c     npoint=150

      if(an/=0) then
              print*,' Local equivalent potential is calculated'
      allocate(uloc0(npoint),uloc1(npoint),perey(nr2))
      allocate(d1uloc(npoint),d2uloc(npoint))
      perey=1
      do  i=1,npoint
      r=dfloat(I)*dr
      pot0(1)=pot(r)
      coul0=vcoul/r
      if(r.ge.rcoul) coul=coul0
      if(r.lt.rcoul) coul=vcoul*(3.d0-(r/rcoul)**2)/rcoul*0.5d0

      if(abs(pot0(1))<1.d-10) then
          uloc=pot0(1)*exp(-an*(ec-coul))
      else
          call solve_trans_eq(pot0(1),ecm-coul,an,uloc,pot0(1))
      endif
      uloc0(i)=uloc
      if(MOD(i, 2) == 0  .and. i/2 < nr2) perey(i/2)=exp(0.5d0*an*uloc)
!     this  factor is with a  different step
!     could be calculated somewhere else.
      enddo
      call  derivative_1(npoint,dr,uloc0,d1uloc)
      call  derivative_2(npoint,dr,uloc0,d2uloc)
      do  i=1,npoint
      r=dfloat(I)*dr
        denom=1-an*uloc0(i)
c       delta_U=-btau**2/16.d0*d2uloc(i)-btau**2/8.d0*d1uloc(i)/r
c    ,          -btau**2*alpha/16.d0*d1uloc(i)**2/denom
        delta_U=-cnl*d2uloc(i)-2*cnl*d1uloc(i)/r
     ,          -cnl**2/4.d0*d1uloc(i)**2/denom
        uloc1(i)=uloc0(i)+delta_U
c       write(800,*) r,real(uloc0(i)),real(delta_U),real(uloc1(i))
c       write(801,*) r,aimag(uloc0(i)),aimag(delta_U),aimag(uloc1(i))
      enddo
      deallocate(d1uloc,d2uloc)
      endif

      do i=1,nr2
      print *,drwf*i,perey(i)
      enddo



      do  i=1,npoint
      r=dfloat(I)*dr

      coul0=vcoul/r
      if(r.ge.rcoul) coul=coul0
      if(r.lt.rcoul) coul=vcoul*(3.d0-(r/rcoul)**2)/rcoul*0.5d0

      call phi_ho(nmax,0,sqrt2*r,b,ho)
      pt=pot(r)
      write(805,*) r,real(pt),aimag(pt)
      if(an/=0)  pt=uloc0(i)
      vlint(:,0)=vlint(:,0)+ pt*r**2*ho(:)
      if(an/=0) then
      vlint(:,1)=vlint(:,1)+ pot(r)*r**2*ho(:)
      vlint(:,4)=vlint(:,4)+ uloc1(i)*r**2*ho(:)
      vlint(:,5)=vlint(:,5)+ pot_so(r)*r**2*ho(:)/(1-an*uloc0(i))
      endif
      vlint(:,2)=vlint(:,2)+ (coul-coul0)*r**2*ho(:)
      vlint(:,3)=vlint(:,3)+ pot_so(r)*r**2*ho(:)
      pt=pot_so(r)
      if(an/=0) write(802,*) r,real(pt),real(pt/(1-an*uloc0(i)))
     ,,aimag(pt/(1-an*uloc0(i)))
      enddo
      vlint=vlint*dr

      return
      end

      function pot(r)
      use potential_parameters
      implicit none
      real*8 r
      complex *16 pot

      pot=cmplx(v0/(1+exp((r-RR)/ar)),0)
     ,   +cmplx(0,Wv/(1+exp((r-RRv)/av)))
     ,   +cmplx(0,4*Ws*exp((r-RRs)/as)/(1+exp((r-RRs)/as))**2)

c     pot=v0*exp(-(r/rr)**2)
      return
      end

      function pot_so(r)
      use potential_parameters
      implicit none
      real*8 r,pot_so

      ! This s-o potential is for nucleons and/or helions
      pot_so=4*Vso*exp((r-RRso)/aso)/(1+exp((r-RRso)/aso))**2/r/aso

      return
      end



*     ------------------------------------------------------------------
      subroutine sigma(lmax,eta,sigma0,sigmad)
*     ------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension sigmad(90)
      if(eta.ge.10.) go to 20
      eta2=eta*eta
      t1=eta2+16.0
      t2=t1*t1
      t4=t2*t2
      t6=eta/6.0
      sigma0=-(eta/(12.0*t1))*(1.0+(eta2-48.0)/(30.0*t2)+(eta2**2-160.0*
     1eta2+1280.0)/(t4*105.0))-eta+3.0*t6*log(t1)+3.5*atan (1.5*t6)-at
     2an (eta)-atan (3.0*t6)-atan (2.0*t6)
      go to 25
   20 t1=1.0/eta
      t2=t1*t1
      t3=t1*t2
      t5=t3*t2
      t7=t5*t2
      t9=t7*t2
      sigma0=0.7853981634d0+eta*log(eta)-eta
     1      -(0.08333333333d0*t1+0.00277777777d0*t3
     2       +0.00079365079d0*t5+0.00059523810d0*t7
     3       +0.00084175084d0*t9)
   25 mn=sigma0/6.2831853072d0
      sigma0=sigma0-6.2831853072d0*mn
      sigmad(1)=sigma0
      do 30 l=2,lmax
      fl1=l-1
      sigmad(l)=sigmad(l-1)+atan (eta/fl1)
c     print *,' COulomb phase shift',l,sigmad(l)
   30 continue
      return
      end
*     ------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine calculate_Legendre(lmin,lmax,mu,angle,nangle,pleg)
*-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension angle(nangle),pleg(nangle,lmin:lmax),cs(nangle)
      real*8,allocatable,dimension(:,:) ::temp

      pi=3.141592653589793d0
      
      pleg=0
      nangle=size(angle)
      allocate(temp(nangle,0:lmax))
      cs=cos(angle(:)/180.d0*pi)
      temp(:,0)=-(-1)**mu*(1-cs**2)**(0.5d0*mu)
     ,         *product([(2*i+1,i=-1,mu-1)])
      if(lmax-mu>0) then
      temp(:,1)=temp(:,0)*cs*(2*mu+1)
      if(lmax-mu>2) then
         do i=1,lmax-1
         temp(:,i+1)=(2*(mu+i)+1)*cs(:)*temp(:,i)-(mu+i+mu)*temp(:,i-1)
         temp(:,i+1)=temp(:,i+1)/dfloat(mu+i-mu+1)
         enddo
      endif
      endif
      pleg(:,lmin+mu:lmax)=temp(:,lmin:lmax-mu)
      deallocate(temp)

      return
      end
*-----------------------------------------------------------------------
      subroutine solve_trans_eq(un,e,fm,ul,uu)
*-----------------------------------------------------------------------
      implicit none
      real *8 fm,error,tolerance,e
      complex*16  un,ul,xi,xnew,fx,dfx,uu
c     real *8     un,ul,xi,xnew,fx,dfx,uu
      integer max_iter,i
      xi=uu
      i = 1
      error = 1d9
      tolerance = 1d-4
      max_iter  = 1000
10    if ((error > tolerance).and.(i < max_iter)) then

      fx=xi-un* exp(-fm*(e-xi))
      dfx=1-un* exp(-fm*(e-xi))*fm

        if (DFx.ne.0d0) then
             xnew = xi -  Fx/DFx
                 if (xnew.ne.0d0) then
              error = abs((xnew - xi)/xnew)*100
            endif
           xi = xnew
           i  = i + 1
          else if (DFx == 0d0) then
         print*, 'Newton''s method fails....derivative equal to zero.'
        stop
        endif
        goto 10
        endif
*
        if (error <= tolerance) then
        else
       print*, 'The root is not reached within the error limit after t
     $he prescribed number of iteration,', i
       endif
*
      ul = xi
      return
      end

      subroutine  derivative_2(n,h,f,df2)
      implicit none
      integer i,n
      real *8 h
      complex*16 f(n),df2(n)

      do i=2,n-1
      df2(i)=f(i-1)-2*f(i)+f(i+1)
      enddo

      df2(1)=f(3)-2*f(2)+f(1)
      df2(n)=f(n)-2*f(n-1)+f(n-2)
      df2=df2/h**2

      return
      end

      subroutine  derivative_1(n,h,f,df1)
      implicit none
      integer i,n
      real *8 h
      complex*16 f(n),df1(n)

      do i=2,n-1
      df1(i)=f(i+1)-f(i-1)
      enddo

      df1(1)=-3*f(3)+4*f(2)-f(1)
      df1(n)=3*f(n)-4*f(n-1)+f(n-2)
      df1=0.5d0*df1/h

      return
      end
