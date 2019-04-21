

subroutine odex(n,fcn,x,y,xend,h,                                 &
     rtol,atol,itol,                                 &
     solout,iout,                                    &
     work,lwork,iwork,liwork,lrcont,licont,idid)
  ! ----------------------------------------------------------
  ! read //venus/user/hairer/programme/dgln_nonstiff/odex_quad.for
  !-----------------------------------------------------------------------
  implicit none
   integer*4 n
   integer*4 lwork
   integer*4 liwork
   integer*4 nfcn
   integer*4 nstep
   integer*4 naccpt
   integer*4 nrejct
   integer*4 iwork
   integer*4 nmax
   integer*4 km
   integer*4 nsequ
   integer*4 iout
   integer*4 mstab
   integer*4 jstab
   integer*4 iderr
   integer*4 mudif
   integer*4 nrdens
   integer*4 icomp
   integer*4 i
   real*16 work
   real*16 uround
   real*16 hmax
   real*16 xend
   real*16 x
   real*16 safe3
   real*16 fac1
   real*16 fac2
   real*16 fac3
   real*16 fac4
   real*16 safe1
   real*16 safe2
   integer*4 lfsafe
   integer*4 iedy
   integer*4 ieyh1
   integer*4 ieyh2
   integer*4 iedz
   integer*4 iescal
   integer*4 iet
   integer*4 iefs
   integer*4 ieys
   integer*4 iehh
   integer*4 iew
   integer*4 iea
   integer*4 iefac
   integer*4 istore
   integer*4 ienj
   integer*4 ieip
   integer*4 lrneed
   integer*4 lrcont
   integer*4 licont
   integer*4 idid
   integer*4 nrd
   real*16 y
   real*16 h
   real*16 rtol
   real*16 atol
   integer*4 itol
   integer*4 nnrd
   integer*4 kmit

   dimension y(n),atol(1000),rtol(1000),work(lwork),iwork(liwork)
  logical arret
  external fcn,solout
  common/stat/nfcn,nstep,naccpt,nrejct
  common /conti/nnrd,kmit,icomp(1000)
  ! *** *** *** *** *** *** ***
  !        setting the parameters
  ! *** *** *** *** *** *** ***
  nfcn=0
  nstep=0
  naccpt=0
  nrejct=0
  arret=.false.
  ! -------- nmax , the maximal number of steps -----
  if(iwork(1)==0)then
     nmax=100000
  else
     nmax=iwork(1)
     if(nmax<=0)then
        write(6,*)' wrong input iwork(1)=',iwork(1)
        arret=.true.
     end if
  end if
  ! -------- km     maximum number of columns in the extrapolation
  if(iwork(2)==0)then
     km=16
  else
     km=iwork(2)
     if(km<=2)then
        write(6,*)' curious input iwork(2)=',iwork(2)
        arret=.true.
     end if
  end if
  ! -------- nsequ     choice of step size sequence
  nsequ=iwork(3)
  if(iwork(3)==0.and.iout<=1) nsequ=1
  if(iwork(3)==0.and.iout>=2) nsequ=3
  if(nsequ<=0.or.nsequ>=5)then
     write(6,*)' curious input iwork(3)=',iwork(3)
     arret=.true.
  end if
  if (nsequ<=2.and.iout>=2) then
     write(6,*)' iwork(3) not compatible with iout'
     arret=.true.
  end if
  ! -------- mstab     parameter for stability check
  if(iwork(4)==0)then
     mstab=1
  else
     mstab=iwork(4)
  end if
  ! -------- jstab     parameter for stability check
  if(iwork(5)==0)then
     jstab=2
  else
     jstab=iwork(5)
  end if
  ! -------- iderr  parameter for error estimation in dense output
  if(iwork(6)==0)then
     if(iout<=1) iderr=1
     if(iout>=2) iderr=0
  else
     iderr=iwork(6)
     if(iout<=1)then
        write(6,*)' error estimation in dense output',              &
             ' not possible, wrong iwork(6)=',iwork(6)
        arret=.true.
     end if
  end if
  ! -------- mudif
  if(iwork(7)==0)then
     mudif=4
  else
     mudif=iwork(7)
     if(mudif<=0.or.mudif>=7)then
        write(6,*)' wrong input iwork(7)=',iwork(7)
        arret=.true.
     end if
  end if
  ! -------- nrdens   number of dense output components
  nrdens=iwork(8)
  !      print *,iwork(8),nrdens
  if(nrdens.lt.0.or.nrdens.gt.n)then
     write(6,*)' curious input iwork(8)=',iwork(8)
     arret=.true.
  end if
  icomp(9)=1

  do  i=1,nrdens
     icomp(i)=iwork(9+i)
  enddo
  ! -------- uround   smallest number satisfying 1.q0+uround>1.q0
  if(work(1)==0.q0)then
     uround=1.q-32
  else
     uround=work(1)
     if(uround<=1.q-35.or.uround>=1.q0)then
        write(6,*)' which machine do you have? your uround was:'    &
             ,work(1)
        arret=.true.
     end if
  end if
  ! -------- maximal step size
  if(work(2)==0.q0)then
     hmax=xend-x
  else
     hmax=abs(work(2))
  end if
  ! -------- step size reduction factor
  if(work(3)==0.q0)then
     safe3=0.5q0
  else
     safe3=work(3)
     if(safe3<=uround.or.safe3>=1.q0)then
        write(6,*)' curious input work(3)=',work(3)
        arret=.true.
     end if
  end if
  ! -------  fac1,fac2     parameters for step size selection
  if(work(4)==0.q0)then
     fac1=0.02q0
  else
     fac1=work(4)
  end if
  if(work(5)==0.q0)then
     fac2=4.0q0
  else
     fac2=work(5)
  end if
  ! -------  fac3, fac4   parameters for the order selection
  if(work(6)==0.q0)then
     fac3=0.8q0
  else
     fac3=work(6)
  end if
  if(work(7)==0.q0)then
     fac4=0.9q0
  else
     fac4=work(7)
  end if
  ! ------- safe1, safe2 safety factors for step size prediction
  if(work(8)==0.q0)then
     safe1=0.65q0
  else
     safe1=work(8)
  end if
  if(work(9)==0.q0)then
     safe2=0.94q0
  else
     safe2=work(9)
  end if
  ! ------- prepare the entry-points for the arrays in work -----
  lfsafe=2*km*km+km
  iedy=11
  ieyh1=iedy+n
  ieyh2=ieyh1+n
  iedz=ieyh2+n
  iescal=iedz+n
  iet=iescal+n
  iefs=iet+km*n
  ieys=iefs+lfsafe*nrdens
  iehh=ieys+km*nrdens
  iew=iehh+km
  iea=iew+km
  iefac=iea+km
  ! ------ total storage requirement -----------
  istore=iefac+2*km-1
  if(istore.gt.lwork)then
     write(6,*)' insufficient storage for work, min. lwork=',istore
     arret=.true.
  end if
  ! ------- entry points for integer workspace -----
  ienj=10+nrdens
  ! --------- total requirement ---------------
  ieip=ienj+km
  istore=ieip+km+1-1
  if(istore.gt.liwork)then
     write(6,*)' insuff. storage for iwork, min. liwork=',istore
     arret=.true.
  end if
  ! --------- control of length of common block "contr" -------
  lrneed=(2*km+5)*nrdens+2
  if(lrcont.lt.lrneed)then
     write(6,*)' insuff. storage for contr, min. lrcont=',lrneed
     arret=.true.
  end if
  ! --------- control of length of common block "conti" -------
  lrneed=nrdens+2
  if(licont.lt.lrneed)then
     write(6,*)' insuff. storage for conti, min. licont=',lrneed
     arret=.true.
  end if
  ! ------ when a fail has occured, we return with idid=-1
  if (arret) then
     idid=-1
     return
  end if
  ! -------- call to core integrator ------------
  nrd=max(1,nrdens)
  call odxcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,km,              &
       solout,iout,idid,nmax,uround,work(iedy),work(ieyh1),           &
       work(ieyh2),work(iedz),work(iescal),work(iefs),                &
       work(ieys),work(iet),work(iehh),work(iew),work(iea),           &
       iwork(ienj),iwork(ieip),nsequ,mstab,jstab,lfsafe,              &
       safe1,safe2,safe3,fac1,fac2,fac3,fac4,iderr,work(iefac),       &
       mudif,nrd)
  ! ----------- return -----------
  return
end subroutine odex
!
!
!
!  ----- ... and here is the core integrator  ----------
!
subroutine odxcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,km,        &
     solout,iout,idid,nmax,uround,dy,yh1,yh2,dz,scal,fsafe,         &
     ysafe,t,hh,w,a,nj,ipoint,nsequ,mstab,jstab,lfsafe,             &
     safe1,safe2,safe3,fac1,fac2,fac3,fac4,iderr,errfac,            &
     mudif,nrd)
  ! ----------------------------------------------------------
  !     core integrator for odex
  !     parameters same as in odex with workspace added
  ! ----------------------------------------------------------
  !         declarations
  ! ----------------------------------------------------------
  implicit none
     integer*4 n
   integer*4 km
   integer*4 lfsafe
   integer*4 nrd
   integer*4 nsequ
   integer*4 i
   integer*4 nj
   real*16 a
   integer*4 itol
   real*16 scal
   real*16 atol
   real*16 rtol
   real*16 y
   integer*4 k
   real*16 h
   real*16 hmax
   real*16 xend
   real*16 x
   integer*4 iout
   integer*4 ipoint
   integer*4 njadd
   integer*4 mu
   real*16 errx
   real*16 prod
   integer*4 j
   real*16 errfac
   integer*4 ipt
   integer*4 nnrd
   integer*4 irtrn
   real*16 xold
   integer*4 naccpt
   real*16 err
   real*16 errold
   real*16 hoptde
   real*16 w
   real*16 h1
   real*16 uround
   integer*4 nstep
   real*16 dz
   integer*4 nfcn
   integer*4 kc
   real*16 dy
   real*16 yh1
   real*16 yh2
   real*16 t
   real*16 hh
   real*16 fac
   real*16 safe1
   real*16 fac1
   real*16 fac2
   real*16 safe2
   real*16 safe3
   integer*4 mstab
   integer*4 jstab
   real*16 fsafe
   real*16 ysafe
   integer*4 nmax
   integer*4 kmit
   integer*4 mudif
   real*16 dens
   integer*4 icomp
   real*16 xoldd
   real*16 hhh
   integer*4 kln
   real*16 dblenj
   integer*4 l
   real*16 factor
   integer*4 krn
   integer*4 kmi
   integer*4 kbeg
   integer*4 kk
   real*16 facnj
   integer*4 lbeg
   integer*4 lend
   integer*4 iderr
   real*16 errint
   integer*4 nrejct
   integer*4 kopt
   real*16 fac3
   real*16 fac4
   integer*4 idid

  logical reject,last,atov
  external fcn
  dimension y(n),dy(n),yh1(n),yh2(n),dz(n),scal(n)
  dimension t(km,n),nj(km),hh(km),w(km),a(km),rtol(1000),atol(1000)
  dimension fsafe(lfsafe,nrd),ysafe(km,nrd),ipoint(km+1)
  dimension errfac(2*km)
  common /contr/xoldd,hhh,dens(45000)
  common /conti/nnrd,kmit,icomp(1000)
  common /stat/nfcn,nstep,naccpt,nrejct
  ! --- define the step size sequence
  if (nsequ==1) then
     do  i=1,km
        nj(i)=2*i
     enddo
  end if
  if (nsequ==2) then
     nj(1)=2
     do  i=2,km
        nj(i)=4*i-4
     enddo
  end if
  if (nsequ==3) then
     do  i=1,km
        nj(i)=4*i-2
     enddo
  end if
  if (nsequ==4) then
     do i=1,km
        nj(i)=4*i
     enddo

  end if
  ! --- define the a(i) for order selection
  a(1)=3.q0
  do  i=2,km
     a(i)=a(i-1)+nj(i)
  enddo
  ! --- initial scaling
  do  i=1,n
     if (itol==0) then
        scal(i)=atol(1)+rtol(1)*abs(y(i))
     else
        scal(i)=atol(i)+rtol(i)*abs(y(i))
     end if
  enddo
  ! --- initial preparations
  k=max(3,min(km-1,int(-log10(rtol(1))*0.6q0+1.5q0)))
  h=max(abs(h),1.q-4)
  h=min(h,hmax,(xend-x)/2.q0)
  if (iout>=1) then
     if (iout>=2) then
        ipoint(1)=0
        do  i=1,km
           njadd=4*i-2
           if (nj(i).gt.njadd) njadd=njadd+1
           ipoint(i+1)=ipoint(i)+njadd
        enddo
        do  mu=1,km*2
           errx=sqrt(mu/(mu+4.q0))*0.5q0
           prod=1.q0/(mu+4.q0)**2
           do  j=1,mu
              prod=prod*errx/j
           enddo
           errfac(mu)=prod
        enddo
        ipt=0
     end if
     nnrd=nrd
     irtrn=0
     xold=x
     call solout (naccpt+1,xold,x,y,n,irtrn)
     if (irtrn.lt.0) goto 120
  end if
  err=0.q0
  errold=1.q20
  hoptde=hmax
  w(1)=0.q0
  reject=.false.
  last=.false.
10 atov=.false.
  ! --- is xend reached in the next step?
  h1=xend-x
  if (h1<=uround) go to 110
  h=min(h,h1,hmax,hoptde)
  if (h>=h1-uround) last=.true.
  if (nstep==0.or.iout.ne.2) call fcn(n,x,y,dz)
  nfcn=nfcn+1
  ! --- the first and last step
  if (nstep==0.or.last) then
     ipt=0
     nstep=nstep+1
     do  j=1,k
        kc=j
        call midex(j,x,y,h,hmax,n,fcn,dy,yh1,yh2,dz,t,nj,hh,w,        &
             err,fac,a,safe1,uround,fac1,fac2,safe2,scal,atov,safe3,     &
             reject,km,rtol,atol,itol,mstab,jstab,errold,fsafe,lfsafe,   &
             iout,ipt,ysafe,nrd)
        if (atov) go to 10
        if (j.gt.1.and.err<=1.q0) go to 60
     enddo
     go to 55
  end if
  ! --- basic integration step
30 continue
  ipt=0
  nstep=nstep+1
  if (nstep>=nmax) go to 120
  kc=k-1
  do  j=1,kc
     call midex(j,x,y,h,hmax,n,fcn,dy,yh1,yh2,dz,t,nj,hh,w,           &
          err,fac,a,safe1,uround,fac1,fac2,safe2,scal,atov,safe3,        &
          reject,km,rtol,atol,itol,mstab,jstab,errold,fsafe,lfsafe,      &
          iout,ipt,ysafe,nrd)
     if (atov) go to 10
40   continue
  enddo 
  ! --- convergence monitor
  if (k==2.or.reject) go to 50
  if (err<=1.q0) go to 60
  if (err.gt.((nj(k+1)*nj(k))/4.q0)**2) go to 100
50 continue
  call midex(k,x,y,h,hmax,n,fcn,dy,yh1,yh2,dz,t,nj,hh,w,           &
       err,fac,a,safe1,uround,fac1,fac2,safe2,scal,atov,safe3,        &
       reject,km,rtol,atol,itol,mstab,jstab,errold,fsafe,lfsafe,      &
       iout,ipt,ysafe,nrd)
  if (atov) go to 10
  kc=k
  if (err<=1.q0) go to 60
  ! --- hope for convergence in line k+1
55 continue
  if (err.gt.(nj(k+1)/2.q0)**2) go to 100
  kc=k+1
  call midex(kc,x,y,h,hmax,n,fcn,dy,yh1,yh2,dz,t,nj,hh,w,          &
       err,fac,a,safe1,uround,fac1,fac2,safe2,scal,atov,safe3,        &
       reject,km,rtol,atol,itol,mstab,jstab,errold,fsafe,lfsafe,      &
       iout,ipt,ysafe,nrd)
  if (atov) go to 10
  if (err.gt.1.q0) go to 100
  ! --- step is accepted
60 xold=x
  x=x+h
  if (iout>=2) then
     ! ---  kmit = mu of the paper
     kmit=2*kc-mudif+1
     do  i=1,nrd
        dens(i)=y(icomp(i))
     enddo
     xoldd=xold
     hhh=h
     do  i=1,nrd
        dens(nrd+i)=h*dz(icomp(i))
     enddo
     kln=2*nrd
     do  i=1,nrd
1       dens(kln+i)=t(1,icomp(i))
     enddo
     ! --- compute solution at mid-point ----
     do  j=2,kc
        dblenj=nj(j)
        do  l=j,2,-1
           factor=(dblenj/nj(l-1))**2-1.q0
           do  i=1,nrd
              ysafe(l-1,i)=ysafe(l,i)+(ysafe(l,i)-ysafe(l-1,i))/factor
           enddo
        enddo
     enddo

     krn=4*nrd
     do  i=1,nrd
        dens(krn+i)=ysafe(1,i)
     enddo
     ! --- compute first derivative at right end ----
     do  i=1,n
        yh1(i)=t(1,i)
     enddo
     call fcn(n,x,yh1,yh2)
     krn=3*nrd
     do  i=1,nrd
        dens(krn+i)=yh2(icomp(i))*h
     enddo
     ! --- the loop ---
     do  kmi=1,kmit
        ! --- compute kmi-th derivative at mid-point ----
        kbeg=(kmi+1)/2
        do  kk=kbeg,kc
           facnj=(nj(kk)/2.q0)**(kmi-1)
           ipt=ipoint(kk+1)-2*kk+kmi
           do  i=1,nrd
              ysafe(kk,i)=fsafe(ipt,i)*facnj
           enddo
        enddo
        do  j=kbeg+1,kc
           dblenj=nj(j)
           do  l=j,kbeg+1,-1
              factor=(dblenj/nj(l-1))**2-1.q0
              do  i=1,nrd
                 ysafe(l-1,i)=ysafe(l,i)+(ysafe(l,i)-ysafe(l-1,i))/factor
              enddo
           enddo
        enddo
        krn=(kmi+4)*nrd
        do  i=1,nrd
           dens(krn+i)=ysafe(kbeg,i)*h
        enddo
        if (kmi==kmit) goto 180
        ! --- compute differences
        do  kk=(kmi+2)/2,kc
           lbeg=ipoint(kk+1)
           lend=ipoint(kk)+kmi+1
           if (kmi==1.and.nsequ==3) lend=lend+2
           do  l=lbeg,lend,-2
              do  i=1,nrd
                 fsafe(l,i)=fsafe(l,i)-fsafe(l-2,i)
              enddo
           enddo

           if (kmi==1.and.nsequ==3) then
              l=lend-2
              do i=1,nrd
                 fsafe(l,i)=fsafe(l,i)-dz(icomp(i))
              enddo
           end if
        enddo
        ! --- compute differences
        do  kk=(kmi+2)/2,kc
           lbeg=ipoint(kk+1)-1
           lend=ipoint(kk)+kmi+2
           do  l=lbeg,lend,-2
              do  i=1,nrd
                 fsafe(l,i)=fsafe(l,i)-fsafe(l-2,i)
              enddo
           enddo
        enddo
180     continue
     enddo
     call interp(nrd,dens,kmit)
     ! --- estimation of interpolation error
     if (iderr==0.and.kmit>=1) then
        errint=0.q0
        do  i=1,nrd
           errint=errint+(dens((kmit+4)*nrd+i)/scal(icomp(i)))**2
        enddo
        errint=sqrt(errint/nrd)*errfac(kmit)
        hoptde=h/max((errint)**(1.q0/(kmit+4)),0.01q0)
        if (errint.gt.10.q0) then
           h=hoptde
           x=xold
           nrejct=nrejct+1
           reject=.true.
           goto 10
        end if
     end if
     do  i=1,n
        dz(i)=yh2(i)
     enddo
  end if
  do  i=1,n
     y(i)=t(1,i)
  enddo
  naccpt=naccpt+1
  if (iout>=1) then
     call solout (naccpt+1,xold,x,y,n,irtrn)
     if (irtrn.lt.0) goto 120
  end if
  ! --- compute optimal order
  if (kc==2) then
     kopt=3
     if (reject) kopt=2
     go to 80
  end if
  if (kc<=k) then
     kopt=kc
     if (w(kc-1).lt.w(kc)*fac3) kopt=kc-1
     if (w(kc).lt.w(kc-1)*fac4) kopt=min(kc+1,km-1)
  else
     kopt=kc-1
     if (kc.gt.3.and.w(kc-2).lt.w(kc-1)*fac3) kopt=kc-2
     if (w(kc).lt.w(kopt)*fac4) kopt=min(kc,km-1)
  end if
  ! --- after a rejected step
80 if (reject) then
     k=min(kopt,kc)
     h=min(h,hh(k))
     reject=.false.
     go to 10
  end if
  ! --- compute stepsize for next step
  if (kopt<=kc) then
     h=hh(kopt)
  else
     if (kc.lt.k.and.w(kc).lt.w(kc-1)*fac4) then
        h=hh(kc)*a(kopt+1)/a(kc)
     else
        h=hh(kc)*a(kopt)/a(kc)
     end if
  end if
  k=kopt
  go to 10
  ! --- step is rejected
100 continue
  k=min(k,kc)
  if (k.gt.2.and.w(k-1).lt.w(k)*fac3) k=k-1
  nrejct=nrejct+1
  h=hh(k)
  reject=.true.
  go to 30
  ! --- solution exit
110 continue
  idid=1
  return
  ! --- fail exit
120 continue
  !write (6,979) x,h
979 format(' exit of odex at x=',d14.7,'   h=',d14.7)
  idid=-1
  return
end subroutine odxcor
!

subroutine midex(j,x,y,h,hmax,n,fcn,dy,yh1,yh2,dz,t,nj,hh,w,     &
     err,fac,a,safe1,uround,fac1,fac2,safe2,scal,atov,safe3,        &
     reject,km,rtol,atol,itol,mstab,jstab,errold,fsafe,lfsafe,      &
     iout,ipt,ysafe,nrd)
  ! --- this subroutine computes the j-th line of the
  ! --- extrapolation table and provides an estimation
  ! --- of the optimal stepsize
  implicit none
   integer*4 n
   integer*4 km
   integer*4 lfsafe
   integer*4 nrd
   real*16 hj
   real*16 h
   integer*4 nj
   integer*4 j
   integer*4 i
   real*16 yh1
   real*16 y
   real*16 yh2
   real*16 dz
   integer*4 m
   integer*4 njmid
   integer*4 mm
   integer*4 iout
   real*16 ysafe
   integer*4 icomp
   real*16 x
   real*16 dy
   integer*4 ipt
   real*16 fsafe
   real*16 ys
   integer*4 mstab
   integer*4 jstab
   real*16 del1
   real*16 scal
   real*16 del2
   real*16 quot
   real*16 uround
   integer*4 nfcn
   real*16 t
   real*16 dblenj
   integer*4 l
   real*16 fac
   real*16 err
   real*16 t1i
   integer*4 itol
   real*16 atol
   real*16 rtol
   real*16 errold
   real*16 expo
   real*16 facmin
   real*16 fac1
   real*16 fac2
   real*16 safe1
   real*16 safe2
   real*16 hh
   real*16 hmax
   real*16 w
   real*16 a
   real*16 safe3
   integer*4 nstep
   integer*4 naccpt
   integer*4 nrejct
   integer*4 nrdd
   integer*4 imit

  external fcn
  logical reject,atov
  dimension y(n),dy(n),yh1(n),yh2(n),dz(n),scal(n)
  dimension t(km,n),nj(km),hh(km),w(km),a(km),rtol(1000),atol(1000)
  dimension fsafe(lfsafe,nrd),ysafe(km,nrd)
  common /stat/nfcn,nstep,naccpt,nrejct
  common /conti/nrdd,imit,icomp(1000)
  hj=h/nj(j)
  ! --- euler starting step
  do  i=1,n
     yh1(i)=y(i)
     yh2(i)=y(i)+hj*dz(i)
  enddo
  ! --- explicit midpoint rule
  m=nj(j)-1
  njmid=nj(j)/2
  do  mm=1,m
     if (iout>=2.and.mm==njmid) then
        do  i=1,nrd
           ysafe(j,i)=yh2(icomp(i))
        enddo
     end if
     call fcn(n,x+hj*mm,yh2,dy)
     if (iout>=2.and.abs(mm-njmid)<=2*j-1) then
        ipt=ipt+1
        do  i=1,nrd
           fsafe(ipt,i)=dy(icomp(i))
        enddo
     end if
     do  i=1,n
        ys=yh1(i)
        yh1(i)=yh2(i)
        yh2(i)=ys+2.q0*hj*dy(i)
     enddo
     if (mm<=mstab.and.j<=jstab) then
        ! --- stability check
        del1=0.q0
        do  i=1,n
           del1=del1+(dz(i)/scal(i))**2
        enddo
        del2=0.q0
        do  i=1,n
           del2=del2+((dy(i)-dz(i))/scal(i))**2
        enddo
        quot=del2/max(uround,del1)
        if (quot.gt.4.q0) then
           nfcn=nfcn+1
           goto 79
        end if
     end if
  enddo
  ! --- final smoothing step
  call fcn(n,x+h,yh2,dy)
  if (iout>=2.and.njmid<=2*j-1) then
     ipt=ipt+1
     do  i=1,nrd
        fsafe(ipt,i)=dy(icomp(i))
     enddo

  end if
  do i=1,n
     t(j,i)=(yh1(i)+yh2(i)+hj*dy(i))/2.q0
  enddo
  nfcn=nfcn+nj(j)
  ! --- polynomial extrapolation
  if (j==1) return
  dblenj=nj(j)
  do  l=j,2,-1
     fac=(dblenj/nj(l-1))**2-1.q0
     do  i=1,n
        t(l-1,i)=t(l,i)+(t(l,i)-t(l-1,i))/fac
     end do
  end do

  err=0.q0
  ! --- scaling
  do i=1,n
     t1i=max(abs(y(i)),abs(t(1,i)))
     if (itol==0) then
        scal(i)=atol(1)+rtol(1)*t1i
     else
        scal(i)=atol(i)+rtol(i)*t1i
     end if
     err=err+((t(1,i)-t(2,i))/scal(i))**2
  end do
  err=sqrt(err/n)
  if (err*uround>=1.q0) goto 79
  if (j.gt.2.and.err>=errold) goto 79
  errold=max(4*err,1.q0)
  ! --- compute optimal stepsizes
  expo=1.q0/(2*j-1)
  facmin=fac1**expo
  fac=min(fac2/facmin,max(facmin,(err/safe1)**expo/safe2))
  fac=1.q0/fac
  hh(j)=min(h*fac,hmax)
  w(j)=a(j)/hh(j)
  return
79 atov=.true.
  h=h*safe3
  reject=.true.
  return
end subroutine midex
!
subroutine interp(n,y,imit)
  ! --- computes the coefficients of the interpolation formula
  implicit none
 integer*4 n
 integer*4 imit
 integer*4 i
 real*16 y0
 real*16 y
 real*16 y1
 real*16 yp0
 real*16 yp1
 real*16 ydiff
 real*16 aspl
 real*16 bspl
 real*16 ph0
 real*16 ph1
 real*16 ph2
  real*16 ph3
  real*16 a
  integer*4 im
  real*16 fac1
  real*16 fac2

  dimension y(n*(imit+5)),a(0:30)
  ! --- begin with hermite interpolation
  do  i=1,n
     y0=y(i)
     y1=y(2*n+i)
     yp0=y(n+i)
     yp1=y(3*n+i)
     ydiff=y1-y0
     aspl=-yp1+ydiff
     bspl=yp0-ydiff
     y(n+i)=ydiff
     y(2*n+i)=aspl
     y(3*n+i)=bspl
     if (imit.lt.0) goto 100
     ! --- compute the derivatives of hermite at midpoint
     ph0=(y0+y1)*0.5q0+0.125q0*(aspl+bspl)
     ph1=ydiff+(aspl-bspl)*0.25q0
     ph2=-(yp0-yp1)
     ph3=6.q0*(bspl-aspl)
     ! --- compute the further coefficients
     if (imit.lt.1) goto 20
     a(1)=16.q0*(y(5*n+i)-ph1)
     if (imit.lt.3) goto 20
     a(3)=16.q0*(y(7*n+i)-ph3+3*a(1))
     if (imit.lt.5) goto 20
     do  im=5,imit,2
        fac1=im*(im-1)/2.q0
        fac2=fac1*(im-2)*(im-3)*2.q0
        a(im)=16.q0*(y((im+4)*n+i)+fac1*a(im-2)-fac2*a(im-4))
     enddo
20   continue
     a(0)=(y(4*n+i)-ph0)*16.q0
     if (imit.lt.2) goto 60
     a(2)=(y(n*6+i)-ph2+a(0))*16.q0
     if (imit.lt.4) goto 60
     do  im=4,imit,2
        fac1=im*(im-1)/2.q0
        fac2=im*(im-1)*(im-2)*(im-3)
        a(im)=(y(n*(im+4)+i)+a(im-2)*fac1-a(im-4)*fac2)*16.q0
     enddo
60   continue
     do  im=0,imit
        y(n*(im+4)+i)=a(im)
     end do
100  continue
  end do
  return
end subroutine interp
!
real*16 function contex(ii,x)
  ! ----------------------------------------------------------
  !     this function can be used for coninuous output in conection
  !     with the output-subroutine for odex. it provides an
  !     approximation to the ii-th component of the solution at x.
  ! ----------------------------------------------------------
  implicit none
  integer*4 ii
  real*16 x
  integer*4 i
  real*16 theta
  real*16 xold
  real*16 h
  real*16 theta1
  real*16 phthet
  real*16 y
  integer*4 n
  integer*4 imit
  real*16 thetah
  integer*4 im
  integer*4 icomp

  common /contr/xold,h,y(45000)
  common /conti/n,imit,icomp(1000)
  ! ----- compute place of ii-th component
  i=ii
  !      i=0
  !c       print *,'i',i,ii,(x-xold)/h,imit,icomp(1),icomp(2)
  !      do 5 j=1,n
  !      if (icomp(j)==ii) i=j
  !   5  continue
  !      if (i==0) then
  !         write (6,*) ' no dense output available for comp.',ii
  !         return
  !      end if
  !

  ! ----- compute the interpolated value
  theta=(x-xold)/h
  theta1=1.q0-theta
  phthet=y(i)+theta*(y(n+i)+theta1*(y(2*n+i)*theta+y(3*n+i)*theta1))
  if (imit.lt.0) then
     contex=phthet

     return
  end if
  thetah=theta-0.5q0
  contex=y(n*(imit+4)+i)
  do  im=imit,1,-1
     contex=y(n*(im+3)+i)+contex*thetah/im
  end do
  contex=phthet+(theta*theta1)**2*contex

  return
end function contex
