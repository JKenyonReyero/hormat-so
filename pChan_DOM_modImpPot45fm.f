*-----------------------------------------------------------------
*     inhomogeneous equation solver
*     [e-t-u]h(r)=source*h(r)
*-----------------------------------------------------------------------
      program main
      implicit real*8 (a-h,o-z)
      complex*16 pot1(600),potd1(1200),ucn,pots(600),det,up(600)
      complex*16 potz(1200),potc(1200),potd(1200),uh(600),ui(600)
      complex*16 hp1,hp2,hm1,hm2,phas,uh1,uh2,ui1,ui2,al,bl,cl
      complex*16 wfel(600),wfbu(600),potd1h(600),potd2h(600),du
      complex*16 potd1t(600),potd2t(600),pereyp(600),pereyr(600)
      complex*16 smatel(75),wfn(600),wfv(600),pl,ql,smatin(75)
      complex*16 wxi,wxxi,smatbu(75),potn(1200),pot2(600),potd2(1200)
      complex*16 VNL(151,151,0:30),smaref(75),famp(180),fcoul
      complex*16 deltaup(600)
      dimension potd1o(1200),potd2o(1200),potd1oh(1200),potd2oh(1200) 
      dimension cf3(90),cg3(90),si2(90),cf4(90),cg4(90)
      real*8 xri(200), wri(200), pleg(60,60),plegndr(180),th,pi
      real*8 xsecruth,sigma0,xxx(90),xr,yr,zr
      character(len=20) :: otfile
      integer :: ptc,ith,nth,ptl,ptx,istep,plo,ios
*-----------------------------------------------------------------------
*     common blocks
*-----------------------------------------------------------------------
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      common/bg/v1r,r1,a1,w1i,r1i,a1i,wv1i,v2r,r2,a2,w2i,r2i,a2i,wv2i,rc
      common/bg2/v1so,r1so,a1so,v2so,r2so,a2so,rv1i,av1i,rv2i,av2i,a13
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      common/potsparam/beta,tmas,hb,ptc,iterms
      common/gquad/xri,wri,pleg,igauss
      common/mean/iterm
*-----------------------------------------------------------------------
*     set perey factor to 1 default
*-----------------------------------------------------------------------
      pereyp = (1.0d0, 0.0d0)
      pereyr = (1.0d0, 0.0d0)
*-----------------------------------------------------------------------
*     menu
*-----------------------------------------------------------------------
      write(6,*) ' inhomogeneous equation solver with Nonlocal '
      write(6,*) '       source term for proton channel        '
      write(6,*)
      write(6,*) '  Local from input              [0]          '
      write(6,*) '  NonLocal source               [1]          '
      read(5,*) isc
      if(isc==1) then
       write(6,*) '  Choose nonlocal potential                 '
       write(6,*) '  Be sure to match deuteron channel         '
       write(6,*) '   Giannini-Ricco (N=Z nuclei) [0]          '
       write(6,*) '   Tian-Pang-Ma   (global)     [1]          ' 
       write(6,*) '   Dispersive Optical Model    [2]          ' 
       write(6,*) '   (needs potential file: dataDOM)          ' 
       read(5,*) ptc
      endif
      if(isc==0) then
       write(6,*) '  Import potential?                         '
       write(6,*) '   Use datap potential         [0]          '
       write(6,*) '   Import potential            [1]          ' 
       write(6,*) '   (needs potential files: fort.117)     '
       read(5,*) ptl
      endif
      if(ptl==1) then
       write(6,*) '  Leading order or next to leading order?   '
       write(6,*) '   Leading order               [0]          '
       write(6,*) '   Next to leading order       [1]          ' 
       read(5,*) plo
       if(plo==1) then
        open(unit=205, status='old', action='read', iostat=ios)
        if (ios /= 0) then
         print *, 'Error opening file: fort.205'
         stop
        end if
        read(205,*) n
        do i = 1, n
         read(205,*) xr, yr, zr
C        Example: pack into complex numbers
         pereyr(i) = xr
         pereyp(i) = dcmplx(yr, zr)
        end do
       endif
      endif
      print *,' what step dr do you need?'
      print *,' 0 is  (a-1)/a*0.1 fm'
      print *,' 1 is  (a-2)/a*0.1 fm'
      read *,istep
      if(ptl==1) then
       write(6,*) '  How many points in potential (e.g 240)    '
       read(5,*) ptx
      endif
      write(6,*)  '  Enter outfile data trailer                '
      read(5,*) otfile       
*-----------------------------------------------------------------------
*     files for input/output
*-----------------------------------------------------------------------
      open(27, file = 'smat_proton.'//otfile)
      open(16, file = 'wfns_proton.'//otfile)
      write(16,*) '# Channel distorted waves'
      open(36, file = 'Eff_pots.'//otfile)
*-----------------------------------------------------------------------
*     potz is vanishing coupling array for homogeneous solutions 
*-----------------------------------------------------------------------
      data potz/1200*(0.d0,0.d0)/
*-----------------------------------------------------------------------
*     input data and calculate constants
*-----------------------------------------------------------------------
      call inputs
      nth=180
      pi = 4.d0*datan(1.d0)
*     if non-local potential option chosen
      if(isc == 1) then
       pleg=0.d0
       tmas=a13
       igauss=48
       hb = h
       call gauss(-1.d0,1.d0,igauss,xri,wri)
*      Legendre Polynomials on Gauss quadrature grid 
       do j=0,lmax+1
        do i=1,igauss
         pleg(i,j+1)=plrec(j,xri(i))*wri(i)
        enddo
       enddo 
c      if(ptc==0) beta=sqrt(0.011d0/(rmtarg/(rmtarg+1))*82.940d0)
       if(ptc==0) beta=0.85 
       if(ptc==1) beta=0.88d0
       if(ptc==0.or.ptc==1) then 
        write(6,*) ' nonlocality range =', beta
        iterms=1
       else
        write(6,*) ' DOM: several nonlocal ranges to be read '
        iterms=7
       endif
       write(6,*) ' Now building non-local potential kernel ...'
       call VNLbuild(VNL)
       write(6,*) ' ... Complete'
      endif
*-----------------------------------------------------------------------
*     calculate coulomb functions for matching (energy ecm)
*-----------------------------------------------------------------------
      phas=ucn/2.d0
      rho1=(maxj-16)*rka*h
      rho2=(maxj- 6)*rka*h
      call coulfn(eta,rho1,lmax+15,cf3,cg3,si2)
      call coulfn(eta,rho2,lmax+15,cf4,cg4,si2)
*-----------------------------------------------------------------------
*     deutpot potentials (at maxj2=2*maxj c.m. points)
*-----------------------------------------------------------------------
      if (ptl==0) then
       call deutpot(potd1,potd1o,1)
       call deutpot(potd2,potd2o,2)
      endif
      if (ptl==1) then
*-----------------------------------------------------------------------
*     importpot potentials (at ulocn-pb_MOD.f c.m. points)
*-----------------------------------------------------------------------
       print*, ' maxj = ',maxj, ' step = ',h
       call importpot(potd1)
       call importpot(potd2)
       potd1o(:)=0
       potd2o(:)=0
*       print*,potd1
      endif
*      write(20,'(5e14.7)') real(potd1)
*      write(21,'(5e14.7)') aimag(potd1)
      
*-----------------------------------------------------------------------
*     and at maxj c.m. points
*-----------------------------------------------------------------------
      do i=1,maxj
       potd1h(i) =potd1(2*i-1)
       potd2h(i) =potd2(2*i-1)
       potd1oh(i)=potd1o(2*i-1)
       potd2oh(i)=potd1o(2*i-1)
      enddo
*      write(22,'(5e14.7)') real(potd1h)
*      write(23,'(5e14.7)') aimag(potd1h)
*-----------------------------------------------------------------------
*     potential 2 calculations - not relevant if nonlocal case
*-----------------------------------------------------------------------
      if(isc.eq.1) goto 8042
*-----------------------------------------------------------------------
*     calculate reference (exact) smatrix - potential 2
*     if the potential is local
*-----------------------------------------------------------------------
      do l=0,lmax
      ball=l*(l+1)
*-----------------------------------------------------------------------
*     temporarily removed l.s force and jval loop
*-----------------------------------------------------------------------
*     jlow=max(l-1,0)
*     if(l.eq.0) jlow=1
*     jhigh=min(l+1,lmax)
*     do jval=jlow,jhigh
*     sowt=(jval*(jval+1)-ball-2.d0)/2.d0
*-----------------------------------------------------------------------
      jval=l
      sowt=0.d0
      hp1=phas*(cg3(l+1)+ucn*cf3(l+1))
      hm1=phas*(cg3(l+1)-ucn*cf3(l+1))
      hp2=phas*(cg4(l+1)+ucn*cf4(l+1))
      hm2=phas*(cg4(l+1)-ucn*cf4(l+1))
      det=hm1*hp2-hm2*hp1
*-----------------------------------------------------------------------
*     calculates scattering from potential 2
*-----------------------------------------------------------------------
      do i=1,maxj2
       potd(i)=potd2(i)+sowt*potd2o(i)
      enddo
*-----------------------------------------------------------------------
*     call potadi for diagonal array required for integration
*-----------------------------------------------------------------------
      call potadi(potd)
*-----------------------------------------------------------------------
*     set up for runge-kutta calls (starting values) starting at h
*-----------------------------------------------------------------------
      istart=2
      naxj=maxj-istart
*-----------------------------------------------------------------------
*     homogeneous solution, regular at the origin ihom=1
*-----------------------------------------------------------------------
      ihom=1
      uh(1)=(0.d0,0.d0)
      uh(2)=(1.d-4,1.d-5)
      du=uh(2)*(l+1)/h
      call runkut(naxj,istart,h,uh,du,potd,potz,ihom)
*-----------------------------------------------------------------------
*     assign uh at the two chosen radii for matching
*-----------------------------------------------------------------------
      uh1=uh(maxj-15)
      uh2=uh(maxj-5)
*-----------------------------------------------------------------------
*     coefficients of h+ (bl) and h- (al) in solution uh
*     calculate elastic s-matrix and
*-----------------------------------------------------------------------
      al=(uh1*hp2-uh2*hp1)/det
      bl=(hm1*uh2-hm2*uh1)/det
      smaref(l+1)=-bl/al
*     write(18,970) l,smaref(l+1),abs(smaref(l+1))
  970 format(1h ,i4,3(3x,e16.8))
*-----------------------------------------------------------------------
      enddo
 8042 continue
*-----------------------------------------------------------------------
*     orbital and total angular momentum (partial wave) loops
*-----------------------------------------------------------------------
      famp=0
      call sigma(lmax,eta,sigma0,xxx)                          
      print *,' sigma0 = ',sigma0
      do 999 l=0,lmax
      ball=l*(l+1)
*-----------------------------------------------------------------------
*     temporarily removed l.s force and jval loop
*-----------------------------------------------------------------------
*     jlow=max(l-1,0)
*     if(l.eq.0) jlow=1
*     jhigh=min(l+1,lmax)
*     do 999 jval=jlow,jhigh
*     sowt=(jval*(jval+1)-ball-2.d0)/2.d0
*-----------------------------------------------------------------------
      jval=l
      sowt=0.d0
      hp1=phas*(cg3(l+1)+ucn*cf3(l+1))
      hm1=phas*(cg3(l+1)-ucn*cf3(l+1))
      hp2=phas*(cg4(l+1)+ucn*cf4(l+1))
      hm2=phas*(cg4(l+1)-ucn*cf4(l+1))
      det=hm1*hp2-hm2*hp1
*-----------------------------------------------------------------------
*     calculate scattering from potential 1
*     remove Coulomb from potentials for source term
*-----------------------------------------------------------------------
      do i=1,maxj2
       potd(i)=potd1(i)+sowt*potd1o(i)
      enddo
      do i=1,maxj
       r=(i-1)*h
       pot1(i)=potd1h(i)+sowt*potd1oh(i)-coul(rc*a13,r)
       pot2(i)=potd2h(i)+sowt*potd2oh(i)-coul(rc*a13,r)
      enddo
*-----------------------------------------------------------------------
*     call potadi for diagonal array required for integration
*-----------------------------------------------------------------------
      call potadi(potd)
*-----------------------------------------------------------------------
*     set up for runge-kutta calls (starting values) starting at h
*-----------------------------------------------------------------------
      istart=2
      naxj=maxj-istart
*-----------------------------------------------------------------------
*     homogeneous solution, regular at the origin ihom=1
*-----------------------------------------------------------------------
      ihom=1
      uh(1)=(0.d0,0.d0)
      uh(2)=(1.d-4,1.d-5)
      du=uh(2)*(l+1)/h
      call runkut(naxj,istart,h,uh,du,potd,potz,ihom)
*-----------------------------------------------------------------------
*     assign uh at the two chosen radii for matching
*-----------------------------------------------------------------------
      uh1=uh(maxj-15)
      uh2=uh(maxj-5)
*-----------------------------------------------------------------------
*     coefficients of h+ (bl) and h- (al) in solution uh
*     calculate elastic s-matrix and
*-----------------------------------------------------------------------
      al=(uh1*hp2-uh2*hp1)/det
      bl=(hm1*uh2-hm2*uh1)/det
      smatel(l+1)=-bl/al
*-----------------------------------------------------------------------
*     renormalise to required physical solution
*-----------------------------------------------------------------------
      do i=1,maxj
       uh(i)=uh(i)/al
       wfel(i)=uh(i)
      enddo
  985 format(1h ,a,2i4,3(3x,e16.8))
*-----------------------------------------------------------------------
*     coefficients of h+ and h- from the elastic wave function
*-----------------------------------------------------------------------
*     pl=(wfel(maxj-15)*hp2-wfel(maxj-5)*hp1)/det
*     ql=(hm1*wfel(maxj-5)-hm2*wfel(maxj-15))/det
*     smatin(l+1)=-ql/pl
*     print 985,'other   s-matrix  ',l,jval,smatin(l+1),abs(smatin(l+1))
*-----------------------------------------------------------------------
*     print 985,'starting  s-matrix',l,jval,smatel(l+1),abs(smatel(l+1))
*     write(19,970) l,smatel(l+1),abs(smatel(l+1))
*-----------------------------------------------------------------------
*     start a loop to iterate second solution
*-----------------------------------------------------------------------
      niter=0
 9538 niter=niter+1
      naxj=maxj-1
      naxj2=maxj2-1
*-----------------------------------------------------------------------
      call source(pot1,pot2,wfel,pots,VNL,niter,isc)
*-----------------------------------------------------------------------
*     shift source potential array down one step for interp
*-----------------------------------------------------------------------
      do j=1,naxj
       pots(j)=pots(j+1)
      enddo
      call interp(pots,h,naxj,potc,h2,naxj2)
      do j=naxj2,1,-1
       potc(j+1)=potc(j)
      enddo
      potc(1)=(0.d0,0.d0)
*-----------------------------------------------------------------------
*     set up for runge-kutta calls (starting values) starting at h
*-----------------------------------------------------------------------
      istart=2
      naxj=maxj-istart
*-----------------------------------------------------------------------
*     inhomogeneous solution, regular at the origin 
*-----------------------------------------------------------------------
      ihom=1
      ui(1)=(0.d0,0.d0)
      ui(2)=wfel(2)
      du=ui(2)*(l+1)/h
      call runkut(naxj,istart,h,ui,du,potd,potc,ihom)
*-----------------------------------------------------------------------
*     assign ui at the two chosen radii for matching
*-----------------------------------------------------------------------
      ui1=ui(maxj-15)
      ui2=ui(maxj-5)
*-----------------------------------------------------------------------
   85 format(1h ,f6.2,3(3x,e13.5),2i5)
*-----------------------------------------------------------------------
*     produce physical outgoing waves solution h
*     solve linear equations for coefficients of h- in ui
*     subtract right amount of homogeneous solution to get uphys(r)
*-----------------------------------------------------------------------
      cl=(ui1*hp2-ui2*hp1)/det
      do i=1,maxj
       up(i)=(ui(i)-cl*uh(i))+uh(i)
       wfel(i)=up(i)
      enddo
      pl=(wfel(maxj-15)*hp2-wfel(maxj-5)*hp1)/det
      ql=(hm1*wfel(maxj-5)-hm2*wfel(maxj-15))/det
      smatin(l+1)=-ql/pl
      if(niter.lt.iterm) goto 9538
*-----------------------------------------------------------------------
      print 985,'s-matrix',l,jval,smatin(l+1),abs(smatin(l+1))
      write(27,970) l,smatin(l+1),abs(smatin(l+1))
  
c     if(l==0) smatin(1)=(-0.0206027930, -0.0988632253 )
c     if(l==1) smatin(2)=(-0.1916937355,  0.2312087261)
c     if(l==2) smatin(3)=(-0.1360159483,  0.0964452529)
c     if(l==3) smatin(4)=(   0.3813415538, -0.2670806921 )

      do ith=1,nth
      th=cos(dfloat(ith)*pi/180.)
      plegndr(ith)=plrec(l,th)
      if(ith==180) write(124,*) l,plegndr(ith)
      enddo


      famp(:)=famp(:)+(2*l+1)*(smatin(l+1)-1.d0)/(2*(0.d0,1.d0)*rka)
     ,       *plegndr(:)*exp((0.d0,1.d0)*2.d0*xxx(l+1))


      if(l==0) then
       write(16,*) '# L=',l
       do i=1,maxj
        write(16,960) (i-1)*h,up(i)
c        write(16,960) (i-1)*h,up(i)*pereyp(i)
c        write(6,*) up(i),pereyp(i),up(i)*pereyp(i)
       enddo
      elseif (l.ge.1) then
       do k = 1,2
        write(16,*) '# L=',l
c        write(6,*) maxj
        do i=1,maxj
         write(16,960) (i-1)*h,up(i)
c        write(16,960) (i-1)*h,up(i)*pereyp(i)
c         write(6,*) up(i),pereyp(i),up(i)*pereyp(i)

        enddo
       enddo
      endif
  960 format(f22.4,3d18.8) 
  999 continue

      

      do ith=1,180
      th=(dfloat(ith)*pi/180.)
      fcoul=-eta/(2*rka)/(sin(th/2.d0))**2 * exp((0.d0,1.d0)
     ,      *(-eta*log((sin(th/2.d0))**2)+2.d0*sigma0))
      xsecruth=10*(eta/(2*rka))**2/(sin(th/2.d0))**4                  ! 10* to go to mb
c     write(123,*) ith,abs(famp(ith))**2*10
c     write(123,*) ith,real(famp(ith)),aimag(famp(ith)) !abs(famp(ith))**2
      write(123,*) ith,10*abs(1*fcoul+famp(ith))**2/xsecruth,xsecruth ! 10* to go to mb
      write(126,*) ith,10*abs(1*fcoul+famp(ith))**2                   ! 10* to go to mb

c     write(125,*) ith,abs(fcoul)**2,xsecruth
      enddo


      stop  
      end
*-----------------------------------------------------------------------
*     initialize data
*-----------------------------------------------------------------------
      subroutine inputs
      implicit real*8 (a-h,o-z)
      complex*16 ucn
      character*12 potname,wfname,wfile,enfile
      common/fnames/potname,wfname,wfile,enfile
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      common/bg/v1r,r1,a1,w1i,r1i,a1i,wv1i,v2r,r2,a2,w2i,r2i,a2i,wv2i,rc
      common/bg2/v1so,r1so,a1so,v2so,r2so,a2so,rv1i,av1i,rv2i,av2i,a13
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      common/mean/iterm
*-----------------------------------------------------------------------
      open(12,file='datap',status='unknown')
*-----------------------------------------------------------------------
      rewind 12
*-----------------------------------------------------------------------
*     read basic data from data set 'data_NL'
*-----------------------------------------------------------------------
      read(12,*) rmtarg,ztarg,elab,qval
      read(12,*) hb,maxj,lmax  ! max J is number of points
      print*,' target mass = ',rmtarg,'  target charge = ',ztarg
      print*,' deuteron lab energy = ',elab,'  Q-value = ',qval
      if(qval.ne.0.d0) h = hb*((rmtarg-1.)/rmtarg)
      if(qval==0.d0) h=hb
      print*,' step = ',h,' nstep = ',maxj,' lmax = ',lmax
      read(12,*) iterm
      print*,' number of iterations = ',iterm
      print*
*     Potential 1
      read(12,*) v1r,r1,a1
      read(12,*) w1i,r1i,a1i
      read(12,*) wv1i,rv1i,av1i
      read(12,*) v1so,r1so,a1so
*     Potential 2
      read(12,*) v2r,r2,a2
      read(12,*) w2i,r2i,a2i
      read(12,*) wv2i,rv2i,av2i
      read(12,*) v2so,r2so,a2so
*     Coulomb radius      
      read(12,*) rc
      if(rc.eq.0.d0) rc=1.3d0
*-----------------------------------------------------------------------
*     constants for (A+1)-p centre of mass motion
*-----------------------------------------------------------------------
      a13=(rmtarg)**0.3333333333d0
      rmu=rmproj*rmtarg/(rmproj+rmtarg)
      if(qval.ne.0.d0) then
       ecmd=(rmtarg-1.)/(rmtarg+rmproj)*elab 
       ecm=ecmd+qval
      else
       ecm=(rmtarg)/(rmtarg+rmproj)*elab
      endif
      const=2.d0*con2/sqrt(con1)*ztarg
      eta=con2*ztarg*sqrt(rmu/ecm)
c     eta=0
c     const=0

      rk2=rmu*ecm*con1
      con=rk2/ecm
      rka=sqrt(rk2)
      h2=h/2.d0
      maxj2=2*maxj
      print*,' ecm = ',ecm,' k = ',rka
      print*,' eta = ',eta
      return
      end
*-----------------------------------------------------------------------
*     block data
*-----------------------------------------------------------------------
      block data
      implicit real*8 (a-h,o-z)
      complex*16 ucn
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      data rmproj/1.0d0/,ucn/(0.d0,1.d0)/
      data con1,con2/0.047845d0,0.1574855d0/
      end
*-----------------------------------------------------------------------
*     calculates diagonal potential for runge-kutta (adiabatic)
*-----------------------------------------------------------------------
      subroutine potadi(potd)
      implicit real*8 (a-h,o-z)
      complex*16 potd(1200),ucn
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      potd(1)=(0.d0,0.d0)
      do i=2,maxj2
       r=(i-1)*h2
       potd(i)=ball/(r*r)+con*potd(i)-rk2
      enddo
      return
      end
*-----------------------------------------------------------------------
*     deuteron potential (at step length/2)
*-----------------------------------------------------------------------
      subroutine deutpot(potd,potdo,ipot)
      implicit real*8 (a-h,o-z)
      real*8 rr,potrr,potii
      complex*16 ucn,potd(1200)
      dimension potdo(1200)
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      common/bg/v1r,r1,a1,w1i,r1i,a1i,wv1i,v2r,r2,a2,w2i,r2i,a2i,wv2i,rc
      common/bg2/v1so,r1so,a1so,v2so,r2so,a2so,rv1i,av1i,rv2i,av2i,a13
*-----------------------------------------------------------------------
*     statement function definitions (nucleon potentials)
*-----------------------------------------------------------------------
      ws(r,vv,r0,aa)=-vv/(1.d0+exp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*exp((r-w0)/wa)/(1.+exp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*-----------------------------------------------------------------------
      print*,'-----------------------------------------'
      print 11,ipot
   11 format(1h ,'   Potential ',i4,/,
     +'    vr     r0     a0     wsi     ri     ai     wvi')
      a13=(rmtarg)**0.3333333333d0
      if(ipot.eq.1) then
      print 12,v1r,r1,a1,w1i,r1i,a1i,wv1i
      print 18,'rc = ',rc,' vso = ',v1so,' rso = ',r1so,' aso = ',a1so
      do i=1,maxj2
c       print*,i
c      need to write r, potr, poti to another file and then read it in
c      again, compare read file to calculated potential to double check
c      the read was done properly.
c      can keep the section that does if r=rr to make sure the step is correct.
       r=(i-1)*h2
       potr=ws(r,v1r,r1*a13,a1)
       poti=wsd(r,w1i,r1i*a13,a1i)
       if(wv1i.gt.0.d0) poti=poti+ws(r,wv1i,r1i*a13,a1i)
c       write(24,*)r,potr,poti
c       if (i/=1) then
c       print*,i,r
c       read(17,*)rr,potr,poti
c       if (r/=rr) then
c       print*, "r/=r for ",i
c       print*, "r rr",r,rr
c       stop
c       endif
c       endif
       potd(i)=cmplx(potr,poti)+coul(rc*a13,r)
       if(r.eq.0) then
        potdo(1)=0.d0
       else
        potdo(i)=wso(r,v1so,r1so*a13,a1so)
       endif
      enddo 
c      print*,potd
      endif
      if(ipot.eq.2) then
      print 12,v2r,r2,a2,w2i,r2i,a2i,wv2i
      print 18,'rc = ',rc,' vso = ',v2so,' rso = ',r2so,' aso = ',a2so
      do i=1,maxj2
       r=(i-1)*h2
       potr=ws(r,v2r,r2*a13,a2)
       poti=wsd(r,w2i,r2i*a13,a2i)
       if(wv2i.gt.0.d0) poti=poti+ws(r,wv2i,r2i*a13,a2i)
       potd(i)=cmplx(potr,poti)+coul(rc*a13,r)
       if(r.eq.0) then
        potdo(1)=0.d0
       else
        potdo(i)=wso(r,v1so,r1so*a13,a1so)
       endif
      enddo 
      endif
   12 format(1h ,14f7.3,/)
   18 format(1h ,4(a7,f8.3))
      print*,'-----------------------------------------'
      return
      end
*-----------------------------------------------------------------------
*     source term - including nonlocality
*-----------------------------------------------------------------------
      subroutine source(pot1,pot2,wf,pots,VNL,niter,isc)
      implicit real*8 (a-h,o-z)
      complex*16 ucn,pot1(600),pot2(600),wf(600),pots(600)
      complex*16 VNL(151,151,0:30),NLS(151),loceqv 
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      common/mean/iterm
      if(isc==1) call NLintagain(VNL,wf,NLS,l)
      do i=1,maxj
       pots(i)=(0.d0,0.d0)
       loceqv=(0.d0,0.d0)
       if(isc==0) pots(i)=con*(pot2(i)-pot1(i))*wf(i)
       if(isc==1.and.i.le.151) pots(i)=con*(NLS(i)-pot1(i)*wf(i))
       if(l==2.and.i.le.151.and.niter==iterm) then
        if(isc==0) loceqv=pot2(i)
        if(isc==1.and.(abs(wf(i)).gt.5.d-2)) then
         loceqv=NLS(i)/wf(i)
         write(36,721) (i-1)*h,loceqv
        endif    
       endif
      enddo                             
  721 format(f8.3,2(3x,e17.6))
      return
      end
*-----------------------------------------------------------------------
      subroutine NLintagain(VNL,wf,NLS,l)
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
*     Integrates Nonlocal potential * Wavefunction to
*     produce nonlocal source term for distorted waves
*     VNL - array of NL-potential
*     wf  - distorted wave wavefunction
*     NLS - resultant NL source term
*     R   - radial coordinate
*-----------------------------------------------------------------------
      complex*16 wf(600),NLS(151),VNL(151,151,0:30),intin(151),intres
      real*8, parameter :: fpi = 16.d0*datan(1.d0)
      integer :: ptc
      common/potsparam/beta,tmas,hb,ptc,iterms
      do i = 1, 151
       R = (i-1)*hb
       intres = (0.d0,0.d0)
       intin = (0.d0,0.d0)
       if(i .lt. 41) then  
        lim = 40+i
       elseif(i .gt. 111) then   
        lim = 40+(152-i)
       else
        lim = 81
       endif 
       do j = 1,lim
        m = j+(i-41) 
        if(i.le.40) m = j
        RP = (m-1)*hb
        intin(j) = VNL(i,m,l)*wf(m)*RP
       end do
       call sim(intin,intres,1,lim,hb)
       NLS(i)= R*intres*fpi
      enddo   
      return
      end subroutine NLintagain
*-----------------------------------------------------------------------
      subroutine VNLbuild(VNL)
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
*     Build array of nonlocal array in proton channel from the users 
*     choice of potential
*-----------------------------------------------------------------------
      complex*16 VNL(151,151,0:30)
      real*8 V0,W0,a0,r0,R_0,UFFRE,UFFIM,R,RP,hb,beta,tmas 
      integer i,j,k,lmax,ptc,ipval
      common/jat/h2,maxj2,lmax,rk2,rka,eta,icalr,iwfnr,ball,con
      common/potsparam/beta,tmas,hb,ptc,iterms
      common/dom/V0,W0,r0,a0,R_0,ipval
      VNL = (0.d0,0.d0) 
*     if DOM - input file dataDOM needed
      if(ptc==2) then
       open(22,file='dataDOM',status='unknown')
       read(22,*) iterms
      endif
      do iterm=1,iterms 
       if(ptc==2) then
        read(22,*) beta,ipval,V0,W0,r0,a0
        if(ipval.ne.3) R_0=r0*tmas   
        print*,' DOM term',iterm,'***********'
        print*,' beta =',beta,' ipval = ',ipval
       endif 
       do i = 1,151
        R = (i-1)*hb
        if(i .lt. 41) then  
         lim = 40+i
        elseif(i .gt. 111) then   
         lim = 40+(152-i)
        else
         lim = 81
        endif 
        do j = 1,lim
         m = j+(i-41) 
         if(i.le.40) m = j
         RP = (m-1)*hb
         do k = 0,lmax   
          call UFF(R,RP,k,UFFRE,UFFIM)
          VNL(i,m,k)=VNL(i,m,k)+dcmplx(UFFRE,UFFIM) 
         enddo
        enddo
        write(6,*) 'R = ',R, 'done'
        if(i==76) write(6,*) '   ...half way...'
       enddo
      enddo
      end subroutine VNLbuild
*-----------------------------------------------------------------------
      subroutine UFF(x,R,il,UFFRE,UFFIM)
      implicit none
      common/potsparam/beta,tmas,hb,ptc,iterms
      common/gquad/xri,wri,pleg,igauss
      common/dom/V0,W0,r0,a0,R_0,ipval
      real*8 xri(200),wri(200),pleg(60,60)
      integer ptc, il, j, igauss, iterms, ipval
      real(kind=8) :: R, UNA, tmas, UpS, UpW, Ureal, Uimag, beta
      real(kind=8) :: UFFRE, UFFIM, V0, W0, a0, r0, R_0, mu ,aw,UNAw
      real(kind=8) :: ap, rp, Vp, R_p, x, Y, Z, hb, plrec, c1, c2
      real(kind=8) :: apS, rpS, WpS, R_pS, apW, rpW, WpW, R_pW
      real(kind=16) :: H, Hval
*     potential parameters
      if(ptc==0) then
c      !Giannini-Ricco
c      a0  = 0.57d0!  woods-saxon diff  
c      r0  = 1.25d0!  woods-saxon rad
c      V0  = -88.6d0! potential depth 
c      W0  = -23.3d0! potential depth         
c      R_0 = r0*tmas - 0.285d0
       !Perey-Buck    
       a0  = 0.65d0!  woods-saxon diff  
       aw  = 0.47d0!  woods-saxon diff  
       r0  = 1.22d0!  woods-saxon rad
       V0  = -71.0d0! potential depth 
       W0  = -15.0d0! potential depth         
       R_0 = r0*tmas 
      endif
      if (ptc==1) then
       !Tian,Pang,Ma  
       ap  = 0.58d0!  woods-saxon diff  
       rp  = 1.25d0!  woods-saxon rad
       Vp  = -70.95d0! potential depth 
       R_p = rp*tmas 
       apS  = 0.50d0!  woods-saxon diff  
       rpS  = 1.24d0!  woods-saxon rad
       WpS  = -9.03d0! potential depth 
       R_pS = rpS*tmas   
       apW  = 0.45d0!  woods-saxon diff  
       rpW  = 1.20d0!  woods-saxon rad
       WpW  = -15.74d0! potential depth 
       R_pW = rpW*tmas 
      endif
      !DOM - parameters read in calling routine - in common/dom/
!-----------------------------------------------------------------------
      UFFRE = 0.d0
      UFFIM = 0.d0
      c1=x*x+R*R
      c2=2.d0*x*R
      do j = 1,igauss
       mu = xri(j)
       Y = sqrt(c1+c2*mu)/2.d0
       Z = sqrt(c1-c2*mu)
       Hval=H(Z)*pleg(j,il+1)
       if(ptc==0) then
        UNA = (1.d0 + exp((Y - R_0)/a0))
        Ureal = V0/UNA
        UNAw= (1.d0 + exp((Y - R_0)/aw))
        Uimag = (4.d0*W0/UNAw)*(1.d0-1.d0/UNAw)   
       else if (ptc==1) then
        Ureal = Vp/(1.d0+exp((Y-R_p)/ap))
        UpS=WpS/(1.d0+exp((Y -R_pS)/apS))
        UpW=WpW*4.d0*exp((Y-R_pW)/apW)/(1.+exp((Y-R_pW)/apW))**2
        Uimag = (UpS + UpW)
       else if (ptc==2) then
        if(ipval.eq.1) then
         UNA = (1.d0 + exp((Y - R_0)/a0))
         Ureal = V0/UNA
         Uimag = W0/UNA
        else if(ipval.eq.2) then
         UNA = (1.d0 + exp((Y - R_0)/a0))
         Ureal = (4.d0*V0/UNA)*(1.d0-1.d0/UNA)
         Uimag = (4.d0*W0/UNA)*(1.d0-1.d0/UNA)
        else if(ipval.eq.3) then
         UNA = exp(-(Y*Y)/(R_0*R_0))
         Ureal=V0*UNA
         Uimag=W0*UNA
        endif
       endif
       UFFRE = UFFRE + Hval*Ureal
       UFFIM = UFFIM + Hval*Uimag
      enddo
      UFFRE = UFFRE/2.d0
      UFFIM = UFFIM/2.d0
      end subroutine UFF
!-----------------------------------------------------------------------
      function H(s)
      implicit real*8(a-h,o-z)
!---------------------------------------------------------------------
!     calculates normalised gaussian with an argument of 2s for a 
!     given nonlocality range beta
!---------------------------------------------------------------------
      real(kind=8), parameter :: pi = 4.d0*datan(1.d0)
      real(kind=8) beta, s, tmas, hb
      real(kind=16) :: H1, H2, H3, H4, H
      integer :: ptc
      common/potsparam/beta,tmas,hb,ptc,iterms
      H1 = (sqrt(pi)*beta)**3
      H3 = (s/beta)**2
      H4 = exp(-H3)
      H  = H4/H1
      return
      end function H
!----------------------------------------------------------------------
      subroutine sim(fa,res,m,n,h)
*-----------------------------------------------------------------------
*     subroutine does the integral of fa stored
*     in the arrays of the same name using simpsons rule. the step length
*     is h and the integral is between the elements m and n of the arrays
*     only. resulting integral value is placed in res.
*-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      complex*16 fa(*),dq(905),res,rq1,rq2
      do 90 i=m,n
       dq(i)=fa(i)
   90 continue
      rq1=dq(m+1)
      rq2=dq(m+2)
      i=m+3
   98 continue
      if(i.ge.n) go to 99
      rq1=rq1+dq(i)
      rq2=rq2+dq(i+1)
      i=i+2
      go to 98
   99 continue
      res=0.33333333333d0*h*(dq(m)+4.d0*rq1+2.d0*rq2-dq(n))
      return
      end
*-------------------------------------------------------------------
*     carries out the differentiation
*-------------------------------------------------------------------
      subroutine pderiv (v,vp,maxj,h)
      implicit real*8 (a-h,o-z)
      complex*16 v(maxj),vp(maxj),x(5)
      j2=maxj-2
      v1=1.d0/(h*12.d0)
      do 10 i=1,5
   10 x(i)=v(i)
      vp(1)=v1*(48.d0*x(2)-25.d0*x(1)-36.d0*x(3)+16.d0*x(4)-3.d0*x(5))
      vp(2)=v1*(18.d0*x(3)-3.d0*x(1)-10.d0*x(2)-6.d0*x(4)+x(5))
      vp(3)=v1*(8.d0*(x(4)-x(2))+x(1)-x(5))
      do 30 j=4,j2
      do 20 i=1,4
   20 x(i)=x(i+1)
      x(5)=v(j+2)
   30 vp(j)=v1*(8.d0*(x(4)-x(2))+x(1)-x(5))
      vp(maxj-1)=v1*(6.d0*x(2)-x(1)-18.d0*x(3)+10.d0*x(4)+3.d0*x(5))
      vp(maxj)=(3.d0*x(1)-16.d0*x(2)+36.d0*x(3)-48.d0*x(4)+25.d0*x(5))
      vp(maxj)=vp(maxj)*v1
      return
      end
*-----------------------------------------------------------------------
*     uniform sphere coulomb interaction.
*-----------------------------------------------------------------------
      real*8 function coul(rc,r)
      implicit real*8 (a-h,o-z)
      complex*16 ucn
      common/a/con1,rmu,ecm,ucn,h,l,maxj,rmtarg,rmproj,ztarg,con2,const
      if(r.le.rc) coul=const*(3.d0-((r/rc)**2))/rc/2.d0
      if(r.gt.rc) coul=const/r
      return
      end
*-----------------------------------------------------------------------
*     lagrange interpolation between calculated values
*-----------------------------------------------------------------------
      subroutine interp (v,h,n,vp,hp,np)
      implicit real*8 (a-h,o-z)
      complex*16 v(600),vp(1200)
c      print*,v
      pp=0.d0
      v120=1.d0/120.d0
      dx=hp/h
      x=0.d0
      imax=n-3
c      print*, np
      do 200 j=1,np
      x=x+dx
      i=x
      i=max0(i,3)
      i=min0(i,imax)
      p=x-dfloat(i)
c      print*, p
      if (p) 120,110,120
  110 vp(j)=v(i)
      go to 200
  120 if (p-pp) 130,140,130
  130 a5=p-3.d0
      a4=(p-2.d0)*a5
      a3=(p-1.d0)*a4
      a2=p*a3
      a1=(p+1.d0)*a2
      a6=p+2.d0
      a2=a2*a6
      a6=(p+1.d0)*a6
      a3=a3*a6
      a6=p*a6
      a4=a4*a6
      a6=(p-1.d0)*a6
      a5=a5*a6
      a6=(p-2.d0)*a6
      pp=p
  140 vp(j)=v120*(a6*v(i+3)-a1*v(i-2)+5.d0*(a2*v(i-1)-a5*v(i+2))
     + +10.d0*(a4*v(i+1)-a3*v(i)))
  200 continue
      return
      end
*-----------------------------------------------------------------------
      subroutine runkut (na,ix,h,u,du,vpot,vcoup,ihom)
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
*     routine integrates function at na points given the function u and
*     its derivative du at point ix. the step length h can be +ve or -ve
*     the number of coupled channels is nc. array vpot contains potentia
*     in steps of h/2 (diagonal) and vcoup(i,j) the coupling potentials.
*-----------------------------------------------------------------------
      complex*16 y(2),f(2),g(2),s(2),du,vpot(1200),u(600),vcoup(1200)
*-----------------------------------------------------------------------
*     set odd y's to functions and even y's to derivatives.
*-----------------------------------------------------------------------
      y(1)=u(ix)
      y(2)=du
      h2=0.5d0*h
      h6=h/6.0d0
      i3=ix
      ih=1
      if (h .lt. 0.0d0)  ih=-1
      i2=2*ix-1
*-----------------------------------------------------------------------
*     loop through the na radial points.
*-----------------------------------------------------------------------
      do 900 in=1,na
      do 202 i=1,2
      g(i)=y(i)
  202 s(i)=y(i)
      is=1
      a200=h6
      b200=h2
      go to 100
  203 do 204 i=1,2
      g(i)=g(i)+a200*f(i)
  204 y(i)=s(i)+b200*f(i)
      go to 100
  205 is=2
      i2=i2+ih
      go to 203
  206 is=3
      a200=2.d0*a200
      go to 203
  207 is=4
      b200=h
      i2=i2+ih
      go to 203
  208 do 209 i=1,2
  209 y(i)=g(i)+h6*f(i)
      go to 533
  100 continue
*-----------------------------------------------------------------------
*     start of integration routine
*     the following defines the coupled equations to be solved.
*-----------------------------------------------------------------------
      f(1)=y(2)
      f(2)=vpot(i2)*y(1)+vcoup(i2)
      go to (205,206,207,208),is
*-----------------------------------------------------------------------
*     end of integration routine, renormalize if necessary
*-----------------------------------------------------------------------
  533 if(ihom.ne.1) goto 530
      if(abs(y(1)).gt.1.e5.or.abs(y(2)).gt.1.e5) then
      do 89 i=1,2
   89 y(i)=y(i)*1.d-5
      do 88 i=ix,i3
   88 u(i)=u(i)*1.d-5
      endif
  530 i3=i3+ih
*-----------------------------------------------------------------------
*     store the function value at this radius, r = (i3-1)*dabs(h).
*-----------------------------------------------------------------------
      u(i3)=y(1)
      du=y(2)
  900 continue
      return
      end
*----------------------------------------------------------------------
*     coulomb function routine
*----------------------------------------------------------------------
      subroutine coulfn (eta,rhoz,lmax,fa,ga,sa)
      implicit real*8 (a-h,o-z)
      dimension sa(90),fa(90),fpp(90),ga(90),gpp(90),gd(5)
      data cg0,cg1,cg2,cg3,cg4,cg5/1.223404016d0,4.959570165d-2,
     +8.888888889d-3,2.455199181d-3,9.108958061d-4,2.534684115d-4/
      data cgp,cp1,cp2,cp3,cp4,cp5/-7.078817734d-1,1.728260369d-1,
     +3.174603174d-4,3.581214850d-3,3.117824680d-4,9.073966427d-4/
      data gd/12.d0,-360.d0,1260.d0,-1680.d0,1188.d0/
      nr=1
      drho=0.d0
      eps=1.0d-11
      eta2=eta+eta
      etas=eta*eta
      lp=max0 (lmax,12)+1
      t=lp
      u=t*t+etas
      v=t/u
      w=eta/u
      x=v*v-w*w
      y=2.*v*w
      u=sqrt (u)
      sig0=eta*(log(u)-1.)+(t-0.5)* atan(eta/t)
      do 20 i=1,5
      sig0=sig0-w/gd(i)
      t=v*x-w*y
      w=v*y+w*x
      v=t
   20 continue
   30 if (lp .le. lmax+1) sa(lp)=sig0
      lp=lp-1
      if (lp .le. 0) go to 100
      t=lp
      sig0=sig0- atan (eta/t)
      go to 30
  100 emax=(1.0d-5/eps)**0.16666667
      if (eta .lt. emax) go to 200
      r=eta2
      t=6.0
      t=eta**(1.0/t)
      w=eta*t
      u=t-t*(cg2+cg4/etas)/etas
      v=(cg1+(cg3+cg5/etas)/etas)/w
      g=cg0*(u+v)
      t=1./t
      w=eta*t
      u=t+t*(cp2+cp4/etas)/etas
      v=(cp1+(cp3+cp5/etas)/etas)/w
      gp=cgp*(u-v)
      go to 300
  200 r=max0 (nr,1)-1
      r=rhoz+r*drho
      t=12.+1.4*eta
      if (t .gt. r) r=t
      fk=1.
      f=1.
      gk=0.
      g=0.
      fsk=0.
      fp=0.
      gsk=1.-eta/r
      gp=gsk
      epss=eps*eps
      n=r+r
      do 210 kp=1,n
      t=kp+kp
      u=t*r
      ak=(t-1.)*eta/u
      v=kp*(kp-1)
      bk=(etas-v)/u
      t=ak*fk-bk*gk
      gk=ak*gk+bk*fk
      fk=t
      t=ak*fsk-bk*gsk-fk/r
      gsk=ak*gsk+bk*fsk-gk/r
      fsk=t
      f=f+fk
      g=g+gk
      fp=fp+fsk
      gp=gp+gsk
      test=fk*fk+gk*gk+fsk*fsk + gsk*gsk
      if (test .lt. epss) go to 220
  210 continue
  220 t=r-eta*log(r+r)+sig0
      u=cos (t)
      v=sin (t)
      g=f*u-g*v
      gp=fp*u-gp*v
  300 rs=r
      rho=rhoz
      f=g
      fp=gp
      is=0
      ir=1
      t=r-rho
      if (t) 600,700,310
  310 if (nr .le. 1) go to 320
      is=t/drho
      is=min0 (is+1,nr)
  320 t=is
      rho=rhoz+t*drho
      gg=rho
      is=is+1
      ir=is
  330 rho=rho-drho
      ir=ir-1
      if (ir .gt. 0) go to 600
      ir=max0 (is,1)
      r=rs
      rho=gg
      f=g
      fp=gp
      go to 350
  340 rho=rho+drho
      ir=ir+1
  350 if (ir .gt. nr) return
  600 h=0.5
      w=r-eta2
      if (r-1.0) 601,602,602
  601 h=0.5*r
  602 if (w) 603,605,605
  603 t=sqrt (-r/(w+w))
      if (t-h) 604,605,605
  604 h=t
  605 last=0
      t=rho-r
      if (t) 606,700,607
  606 h=-h
  607 u=t-h
      if (u*h) 608,608,609
  608 h=t
      last=1
  609 u=0.0
      t=1.0
      b1=0.0
      b2=f
      b3=h*fp
      f=f+b3
      v=0.0
  610 it=0
  620 v=-h*(h*b1+w*b2+u*v)/(r*t)
      fp=fp+v
      u=t
      t=t+1.0
      b1=b2
      b2=b3
      b3=h*v/t
      f=f+b3
      test=b3
      testp=v
      if (w) 630,640,640
  630 test=b3/f
      testp=v/fp
  640 if (abs(test)+abs(testp)-eps) 650,610,610
  650 if (it) 660,660,670
  660 it=1
      go to 620
  670 r=r+h
      if (last) 600,600,700
  700 k=lmax+1
      x=f
      y=fp
      do 710 j=1,k
      ga(j)=x
      al=j
      t=j*j
      u=t/rho+eta
      v=sqrt (t+etas)
      w=(u*x-al*y)/v
      y=(v*x-u*w)/al
      x=w
  710 continue
      lp=rho
      lp=max0 (lp+10, lmax+20)
      b3=0.
      b2=1.0d-20
      w=1.0/rho
      al=lp+1
      v=eta/al
      u=0.
      do 840 j=1,lp
      k=lp+1-j
      al=k
      t=k+k+1
      b1=t*(v/al+w)*b2-u*b3
      v=eta/al
      u=sqrt (1.0+v*v)
      b1=b1/u
      b3=b2
      b2=b1
      if (k-lmax-1) 810,810,820
  810 fa(k)=b1
      go to 840
  820 test=b1
      if (abs(test)-1.) 840,840,830
  830 b2=b2*1.0d-20
      b3=b3*1.0d-20
  840 continue
      t=(w+eta)*b2-u*b3
      u=1./(t*f-b1*fp)
      k=lmax+1
      do 850 j=1,k
      fa(j)=u*fa(j)
  850 continue
      do 400 l=1,k
      fl=l
      flsq=fl*fl
      fac1=eta/fl+fl/rhoz
      fac2=sqrt(eta*eta+flsq)/fl
      fpp(l)=fac1*fa(l)-fac2*fa(l+1)
      gpp(l)=fac1*ga(l)-fac2*ga(l+1)
  400 continue
      if (ir-is) 330,340,340
      end
*-----------------------------------------------------------------------
      complex*16 function wxi(xva,ndata,del,wfn)
      implicit real*8 (a-h,o-z)
      complex*16 wfn(ndata)
      complex*16 y1,y2,y3,y4,y5,y6
      xx=0.d0
      do 30 k=2,ndata-1
      nst=0
      xx=xx+del
      if(xx.lt.xva) goto 30
      nst=max0(k-3,1)
      goto 33
   30 continue
   33 if((nst.eq.0).or.(nst.gt.ndata-5)) nst=ndata-5
      x1=(nst-1)*del
      x2=x1+del
      x3=x2+del
      x4=x3+del
      x5=x4+del
      x6=x5+del
      y1=wfn(nst+0)
      y2=wfn(nst+1)
      y3=wfn(nst+2)
      y4=wfn(nst+3)
      y5=wfn(nst+4)
      y6=wfn(nst+5)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)*(x3-x6)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)*(x4-x6)
      pii5=(x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)*(x5-x6)
      pii6=(x6-x1)*(x6-x2)*(x6-x3)*(x6-x4)*(x6-x5)
      xd1=xva-x1
      xd2=xva-x2
      xd3=xva-x3
      xd4=xva-x4
      xd5=xva-x5
      xd6=xva-x6
      pi1=xd2*xd3*xd4*xd5*xd6
      pi2=pi1*xd1/xd2
      pi3=pi2*xd2/xd3
      pi4=pi3*xd3/xd4
      pi5=pi4*xd4/xd5
      pi6=pi5*xd5/xd6
      wxi=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4+y5*pi5/pii5
     *   +y6*pi6/pii6
      return
      end
*-----------------------------------------------------------------------
      recursive real*8 function plrec(il,c) result(plr)
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
*     Legendre functions
*     l<=2 caluclated explicitly
*     l>2  calculated recursively
*-----------------------------------------------------------------------     
      if(il.eq.0) then
       plr=1.d0
      else if(il.eq.1) then
       plr=c
      else if(il.eq.2) then
       plr=(3.d0*c*c-1.d0)/2.d0
      else if(il.gt.2) then
       plr=((2*il-1)*c*plrec(il-1,c) - (il-1)*plrec(il-2,c))/il
      endif
      return
      end      
!-----------------------------------------------------------------------
!****************** Gaussian Quadratures *******************************
      subroutine gauss(a,b,npoint,xri,wri) 
!-----------------------------------------------------------------------
!     a,b - upper and lower limits 
!     npoint - number of points
!     xri - arguments
!     wri - weighting coefficients
!-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8 xg(200),wg(200),xri(200),wri(200)
      call setmgl(npoint,xg,wg)
      do 20 j=1,npoint
      xri(j) = (a+b)/2.d0 + (b-a)/2.d0*xg(j)
      wri(j) = (b-a)/2.d0*wg(j)  
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine setmgl( n, points, weight )
      implicit real*8  ( a-h, o-z )
      real*8 points(200), weight(200)
      real*8 poin16(300)
      pi = 4.d0*atan(1.d0)
      if (n.gt.200)  write (1,50)
50    format(' setmlg call with too many points')
      m = ( n + 1 ) / 2
      e1 = n * ( n + 1 )
      do 1 i = 1, m
      t = ( 4*i - 1 ) * pi / ( 4*n + 2 )
      x0 = ( 1.d0 - ( 1.d0 - 1.d0/n ) / ( 8.d0*n*n ) ) * cos(t)
      pkm1 = 1.d0
      pk = x0
      do 3 k = 2, n
      t1 = x0 * pk
      pkp1 = t1 - pkm1 - ( t1-pkm1 )/k + t1
      pkm1 = pk
      pk = pkp1
3     continue
      den = 1.d0 - x0*x0
      d1 = n * ( pkm1 - x0*pk )
      dpn = d1 / den
      d2pn = ( 2.d0*x0*dpn - e1*pk ) / den
      d3pn = ( 4.d0*x0*d2pn + (2.d0-e1)*dpn ) / den
      d4pn = ( 6.d0*x0*d3pn + (6.d0-e1)*d2pn ) / den
      u = pk / dpn
      v = d2pn / dpn
      h = -u * ( 1.d0 + 0.5d0*u*(v+u*(v*v-u*d3pn/(3.d0*dpn))))
      p = pk + h*(dpn+0.5d0*h*(d2pn+h/3.d0*(d3pn+0.25d0*h*d4pn)))
      dp = dpn + h*(d2pn+0.5d0*h*(d3pn+h*d4pn/3.d0))
      h = h - p / dp
      poin16(i) = x0 + h
      fx = d1 - h*e1*(pk+0.5d0*h*(dpn+h/3.d0* 
     1     (d2pn+0.25d0*h*(d3pn+0.2d0*h*d4pn))))
      weight(i) = 2.d0 * ( 1.d0 - poin16(i)*poin16(i)) / (fx*fx)
1     continue
      if ( m + m .gt. n ) poin16(m) = 0.d0
      do 10 i = n/2 + 1, n
      poin16(i) = poin16( n + 1 - i )
      weight(i) = weight( n + 1 - i )
      poin16( n + 1 - i ) = -poin16( n + 1 - i )
10    continue
      do 30 i=1,n
 30   points(i)=poin16(i)
      return
      end

      subroutine importpot(potd2)
      logical :: exists
      integer :: ptx,ii,nnr
      complex*16 potd(600),potd2(1200)
      real*16    potr(600),poti(600),dr,r(600)
      character(len=20) :: filename
      potd(:)=(0d0,0d0)
      potd2(:)=(0d0,0d0)


      inquire(file="fort.117", exist=exists)
      if (exists) goto 390
      print*,"fort.117 does not exist. Terminating"
      stop
 390  continue
      open(117, action="read")
      read(117,*) nnr, dr
      print*, ' nnr = ',nnr, ' dr = ',dr
      do ii=1, nnr
       read(117,*) r(ii), potr(ii), poti(ii)
       potd(ii) = dcmplx(potr(ii), poti(ii))
      enddo
      close(117)
c     Now need to interp values since code is expecting maxj *2 points
      call interp(potd,0.1d0,600,potd2,0.05d0,1200)

      open(100,file='potd',status='unknown')
      do ii=1, 600
       write(100,*) potd(ii)
      enddo
      close(100)

      open(101,file='potd2',status='unknown')
      do ii=1, 1200
       write(101,*) potd2(ii)
      enddo
      close(101)
      end
*     ------------------------------------------------------------------
      subroutine sigma(lmax,eta,sigma0,sigmad)                          
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
      print *,' Coulomb phase shift',l,sigmad(l)
   30 continue                                                          
      return                                                            
      end                                                               
*     ------------------------------------------------------------------
