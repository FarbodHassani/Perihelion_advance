  
  
  include 'sample.h'
  integer* 4 nfcn,nstep,naccpt,nrejct
  common/stat/nfcn,nstep,naccpt,nrejct


  integer*4 n,ndgl,nrdens,km,lwork,liwork,lrcont,licont
  parameter (ndgl=1000,nrdens=1000,km=20)
  parameter (lwork= ndgl*(km+5)+5*km+10+2*km*(km+1)*nrdens)
  parameter (liwork= 2*km+10+nrdens)
  parameter (lrcont= ( 2 * km + 5 ) * nrdens + 2 )
  parameter (licont= nrdens + 2 )
  real*16 y(ndgl),work(lwork)
  integer*4 iwork(liwork)
  real*16 rcont
  integer*4 icont
  common /contr/rcont(lrcont)
  common /conti/icont(licont)
  real*16 atol(ndgl),rtol(ndgl)
  real*16 tol
  real*16 x,xend,xbeg
  real*16 theta,pi
  common theta
  !     vf is vector field
  integer iangle

  external vf,solout
  real*16 h
  integer*4 i,iout,itol,idid
  parameter(iout=2)
  parameter(n=4)
  !FILE OPEN
  ! * * * * * * * * * * * * * * * * * * * * * * * * *
  open (unit = 14, file = "schw_perihelion.0005.normal.ifort.8_12.fast.dat")


  !  open (unit = 12, file = "schw_txy.dat")

  !  open (unit = 16, file = "Energy.dat")





    do iangle = 0,360,2
  !do iangle = 45,360,90
     theta=3.1415926535897932384626433832795028841971693993751/180.0q0*iangle
     ! could use parameters in function
     !! print *,iangle
     !adapt the parameters lrcont and licont in the subprograms, where they must
     !be excplicit 45002  -> x,y,z(45000)
     !  print *,lrcont,licont
     !  stop
     !## x is an independent variable

!     n=4

     ! --- output routine is used during integration
!     iout=2
     ! --- endpoint of integration
     xbeg=0q0
     xend= period * Norbit * beta **(3./2.) !the rescaling for
     ! --- initial values
     x=xbeg


     ! ------ obligatoire au debut -------
     ! --- initial values
     !  q followed by p
     do i=1,n
        y(i)=0.0q0
     end do


     ! * * * * * * * * * * * * * * * * * * * * * * * * *
     !Initial values FOR SCHWARZCHILD,
     ! we start from perihelion values
     ! * * * * * * * * * * * * * * * * * * * * * * * * *
     !#### y(1) = x
     !#### y(2) = y
     !#### y(3) = x' =v_x =
     !#### y(4) = y' =v_y =
     !#### we start from y(1)=x= r_mercury and let other variables to be 0
     !#### alpha rescales the problem, masses and ...
     !real*16 alpha , vmax

     !Aphelion: speed= 38.86 distance= 69.817079
     ! * * * * * * * * * * * * * * * * * * * * * * * * *
     !IC1=
     ! !y(1) = -69.817079 * beta!(Giga meter)
     !y(4)= +38.86q0 / SQRT(beta) !(y'= vmax = 58.98 for mercury;//(*(Giga meter)/(Mega !second)*)(*vmax=vmax(1/10^9)(*(Gm/m)*)(10^5/1)*)(*(s/d)*)(*Gm/d*))

     !IC2:
     !y(1)=46.001272 * beta!(Giga meter)
     !y(4)= 58.98q0 / SQRT(beta) !(y'= vmax = 58.98 for mercury;//(*(Giga meter)/(Mega !second)*)(*vmax=vmax(1/10^9)(*(Gm/m)*)(10^5/1)*)(*(s/d)*)(*Gm/d*))

     !Rotated IC


     y(1) = COS(theta)* 46.001272 * beta!(Giga meter)
     y(2)=  SIN(theta)* 46.001272 * beta
     y(3)= -SIN(theta)* 58.98q0 / SQRT(beta)  !(x'=0)
     y(4)=  COS(theta)* 58.98q0 / SQRT(beta) !(y'= vmax = 58.98 for mercury;//(*(Giga meter)/(Mega !second)*)(*vmax=vmax(1/10^9)(*(Gm/m)*)(10^5/1)*)(*(s/d)*)(*Gm/d*))
!     y(3)=-y(3)
!     y(4)=-y(4)
     !  print *, iangle
     !  print *,y(1),y(2),y(3),y(4)

     !y(4)= 0.q0 !(y'= vmax = 58.98 for mercury;//(*(Giga meter)/(Mega !second)*)(*vmax=vmax(1/10^9)(*(Gm/m)*)(10^5/1)*)(*(s/d)*)(*Gm/d*))


     ! --- required tolerance
     !       rtol=1.0q-20
     !       atol=1.0q-7*rtol
     ! --- initial step size
     h=1.0q-4

     ! --- required (relative) tolerance
     tol=tolerance
     itol=0
     rtol=tol
     atol=tol
     ! --- default values for parameters
     do 10 i=1,20
        iwork(i)=0
10      work(i)=0.d0
        iwork(1)=10000000
        iwork(8)=n
        iwork(9)=4
        do i=1,n
           iwork(9+i)=i
        enddo
        iwork(6)=0
        ! --- call of the subroutine odex
        call odex(n,vf,x,y,xend,h,                                   &
             &                  rtol,atol,itol,                                 &
             &                  solout,iout,                                    &
             &                  work,lwork,iwork,liwork,lrcont,licont,idid)

     enddo
   end program

   !
   real*16 function derivative(t)
     implicit none
     real*16 t
     real*16 distance
     real*16 dt
     dt=1.q-5 ! is this a good value
     derivative=(distance(t+dt)-distance(t-dt))/(2*dt)
   end function derivative

   real*16 function distance(t)
     implicit none
     real*16 t
     real*16 contex
     external contex
     ! JPE no
     ! you want the derivative of sqrt(x^2+y^2)
     ! so that is
     real*16 temp(10)
     integer*4 i
     do i=1,2 ! since you are in R^2
        !contex gives the i component at time t and other variables should be there!
        temp(i)=contex(i,t)
     enddo
     distance=temp(1)*temp(1) + temp(2)*temp(2)

   end function distance

   !################
   !###Function Phi (i,j)= -GM/sqrt(i^2 dx^2 + j^2 dx^2)
   !################
   real*16 function Phi(ix,iy,dl,GM)
     integer*4 ix,iy
     real*16 distance,dl,GM
     distance = (ix*dl)**2 + (iy*dl)**2
     !     print *,distance,'dd'
     Phi= -GM/sqrt(distance)
   end function Phi



   subroutine solout (nr,xold,x,y,n,irtrn)
     ! --- prints solution at equidistant output-points
     ! --- by using "contex", the continuous collocation solution
     include 'sample.h'
     integer*4 nr,n,irtrn
     real*16 y(n)
     real*16  theta
     common theta
!!!!!! jpadded
     real*16 oldderivative
     external derivative
     real*16 derivative
     save lastx,lasti
     real*16 zeroin
     external zeroin
     !! end jpadded
     real*16 dt,t
     integer*4 nrpts,i,j,increasing
     real*16 x,xold,xoldold
     real*16 contex
     real*16 lastx,zerotime
     save increasing,xoldold

     real*16 rr,angle
     real*16 twopi
     real*16 distance
     external distance
     real*16 x1,y1
     real*16 difference
     integer*4 currenti,lasti,jj
     twopi=2.* 3.1415926535897932384626433832795028841971693993751q0
     !  print *,nr
     if(nr==1)then
        xold=0
        lasti=0
        increasing=1
        return
     endif


     dt=sqrt(2.1q0)*1q0/2048q0
     

     currenti=x/dt
     !!    print *,currenti,lasti,xold,x,'current'
     !   !############PRINTING#####################
     !   !#########################################
!     do i=lasti+1,currenti
!        x1=contex(1,i*dt)
!       y1=contex(2,i*dt)
 !       write(12,'(4e17.8)')x1,y1,theta,i*dt
!        lasti=i
!     enddo
     !         print *,distance(xold),distance(x),contex(1,x),increasing,'dist'

 

!     print *,'****',COS(theta)*contex(1,x)+SIN(theta)*contex(2,x)
        if(increasing==0 .and.(COS(theta)*contex(1,x)+SIN(theta)*contex(2,x))>0.and.&
             distance(xold)<distance(x))then
           !! print *,'min reached'
!           !! print *,distance(xoldold),distance(xold),distance(x)
           zerotime=zeroin(xoldold,x,derivative,0Q00)
           if(zerotime==0)then
              !! print *,distance(xold),distance(xoldold),distance(x)
!              stop
           endif
!           !! print *,'der',derivative(zerotime)
           !! print *,'zerotime',zerotime,xoldold,x
           !! print *,contex(1,zerotime),contex(2,zerotime)
           angle=ATAN2(contex(2,zerotime),contex(1,zerotime))-theta
           
           !! print *,'angle   ',angle*360/twopi
           !           !! print *,'angtheta',mod(angle,theta)
           do while(angle<-1.0q0)
              !! print *,'correcting',angle
              angle=angle+twopi
           enddo
           !! print *,'mod     ', mod(angle,twopi)
           write(14,'(2e25.16,i4,2e17.8)')angle,zerotime,int(theta/3.14155926*180),zerotime/period,dx
        endif

 
     if( distance(xold)<distance(x))increasing=1
     if( distance(xold)>distance(x))increasing=0
     lasti=currenti
     xoldold=xold
     xold=x


     !   print *,'lastx=',lastx
     return
     !close(12)
   end subroutine solout



   !####################################################
   !######Differential Equation Adding
   !####################################################
   subroutine vf(n,x,y,f)
     ! --- right-hand side of equ
     include 'sample.h'
     integer*4 n
     real*16 x
     real*16 y(n),f(n)

     !     y(1) = x, y(2)=y
     !     time is the independent
     !(x'=v_x , v_x'=  -GMx/SQRT(x^2+y^2)^3, y'=v_y , v_y'=  -GMy/sqrt(x^2+y^2)^3 )

     !####################VARIABLES###############
     !###### alpha = rescale the problem, mass of sun and maximum velocity.
     !#######################################
     real*16 G, c, Mearth, Lm, MSun, R_perihelion, Rsch, GML2, vmax, TotalEnergy ,KineticEnergy, PotentialEnergy
     real*16 xpos,ypos,r
     integer*4 ix,iy,i,j
     integer*4 fx,fy
     external Phi
     real*16 Phi,u,v,forcex,forcey
     real*16 Gradphix(0:1,0:1), Gradphiy(0:1,0:1)

     !#####################################################
     Mearth = 5.972*1.d24 !(*kg*)
     MSun = 1.99 * 1.d30 /Mearth  !(*kg*) (*Normalized to earth's mass*)
     R_perihelion=46.0  * beta!(*Giga meter*)
     G = 6.67*1.d-26 * Mearth
     c = 2.998 * 1.d5
     vmax = 58.98/ SQRT(beta)   !(*(Giga meter)/(Mega second)*)
     Lm = R_perihelion * vmax;
     GML2=G * MSun/(Lm*Lm)
     Rsch=2 * G * MSun/(c*c)
     !#####################################################
     !print * ,f(1), f(2)
     !    print *,y(1),x
     !write(16,'(3e17.8)')x1,y1,i*dt

     f(1)= y(3) !(x' = v_x)
     f(2)= y(4) !(y' = v_y)
     xpos=y(1)
     ypos=y(2)
     !print *,fx, xpos, dx, xpos/dx
     fx= floor(xpos/dx) !The number of cell particle belongs to
     fy=floor(ypos/dx) !The number of cell particle belongs to

     u=xpos-fx*dx
     v=ypos-fy*dx


     !DF/Dx  = (1-y) * phi(1,0)-phi(0,0)/dx
     ! sprengel
     !! first compute derivative in the 4 corners of
     !! the square in which the particle lies.

     do i=0,1
        do j=0,1

           gradphix(i,j)=(Phi(fx+i+1,fy+j,dx,G*Msun) - Phi(fx+i-1,fy+j,dx,G*Msun))/(2*dx)
           gradphiy(i,j)=(Phi(fx+i,fy+j+1,dx,G*Msun) - Phi(fx+i,fy+j-1,dx,G*Msun))/(2*dx)
        enddo
     enddo

     !! now interpolate between the 4 values
     !! the x-component
     !! 

     !!

     forcex= (gradphix(0,0)*(1d0-u)+gradphix(1,0)*u)*(1-v) &
            +(gradphix(0,1)*(1d0-u)+gradphix(1,1)*u)*v
     forcey= (gradphiy(0,0)*(1d0-v)+gradphiy(0,1)*v)*(1-u) &
            +(gradphiy(1,0)*(1d0-v)+gradphiy(1,1)*v)*u
            
          
     




     f(3)=  -forcex
     f(4)=  -forcey

     !     print *,'y=',y
     !     print *,'f=',f
   end subroutine vf
