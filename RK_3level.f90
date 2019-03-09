subroutine RK_3level(xini,xfin,xNbin,yini,zini,wini)

use laserParameter
use atomicSystem

implicit none

double precision,intent(in) :: xini,xfin
double complex,intent(in) ::yini,zini,wini
integer,intent(in) :: xNbin				!øÃ§ﬂ ˝

integer::n, Ndisp
double precision::x,dx,Ex,Ey,result_mod
double complex::y,z,w
double complex::k1,k2,k3,k4,q1,q2,q3,q4,r1,r2,r3,r4
open(20,file='odeResults.dat')

dx=(xfin-xini)/dble(xNbin)		!øÃ§‡∑˘§Ú”ãÀ„
print*, "dx=",dx

x=xini 							! x --> t
y=yini							! up
z=zini                          ! um
w=wini                          ! u2s


!•Î•Û•≤•Ø•√•ø∑®§Œ”ãÀ„§Ú§π§Î
print*, "The initial condition is:", x,y,z,w

do n=1,xNbin
    call fEt(x, Ex, Ey)

	k1=dx*funcdy(x,y,z,w,Ex,Ey)
	q1=dx*funcdz(x,y,z,w,Ex,Ey)
	r1=dx*funcdw(x,y,z,w,Ex,Ey)

	k2=dx*funcdy(x+0.50d0*dx,y+0.50d0*k1,z+0.50d0*q1,w+0.50d0*r1,Ex,Ey)
	q2=dx*funcdz(x+0.50d0*dx,y+0.50d0*k1,z+0.50d0*q1,w+0.50d0*r1,Ex,Ey)
	r2=dx*funcdw(x+0.50d0*dx,y+0.50d0*k1,z+0.50d0*q1,w+0.50d0*r1,Ex,Ey)

	k3=dx*funcdy(x+0.50d0*dx,y+0.50d0*k2,z+0.50d0*q2,w+0.50d0*r2,Ex,Ey)
	q3=dx*funcdz(x+0.50d0*dx,y+0.50d0*k2,z+0.50d0*q2,w+0.50d0*r2,Ex,Ey)
	r3=dx*funcdw(x+0.50d0*dx,y+0.50d0*k2,z+0.50d0*q2,w+0.50d0*r2,Ex,Ey)

	k4=dx*funcdy(x+dx,y+k3,z+q3,w+r3,Ex,Ey)
	q4=dx*funcdz(x+dx,y+k3,z+q3,w+r3,Ex,Ey)
	r4=dx*funcdw(x+dx,y+k3,z+q3,w+r3,Ex,Ey)

	y=y+(k1+2d0*k2+2d0*k3+k4)/6d0
	z=z+(q1+2d0*q2+2d0*q3+q4)/6d0
	w=w+(r1+2d0*r2+2d0*r3+r4)/6d0
    x=xini+dx*dble(n)

    if(mod(n,1000)==0) then
        write(20,*) x/Toc,real(y),aimag(y),real(z),aimag(z),real(w),aimag(w)
    end if

end do

close(20)

!-----------------------------------------

stop
contains


!dy/dx
function funcdy(X,Y,Z,W,Ex,Ey) result(dydx)
    double precision,intent(in)::X,Ex,Ey
	double complex,intent(in)::Y,Z,W
	double complex :: dydx,H13

    H13 = (0.0d0,1.0d0)*(Ex*d1x+Ey*d1y)*exp( (0.0d0,-1.0d0)*Delta_Ek1*x )

	dydx=H13*W

end function funcdy

!dz/dx
function funcdz(X,Y,Z,W,Ex,Ey) result(dzdx)
    double precision,intent(in)::X,Ex,Ey
	double complex,intent(in)::Y,Z,W
	double complex dzdx,H23

    H23 = (0.0d0,1.0d0)*(Ex*d2x+Ey*d2y)*exp( (0.0d0,-1.0d0)*Delta_Ek2*x )

	dzdx=H23*W

end function funcdz

!dw/dx
function funcdw(X,Y,Z,W,Ex,Ey) result(dwdx)
    double precision,intent(in)::X,Ex,Ey
	double complex,intent(in)::Y,Z,W
	double complex dwdx,H31,H32

    H31 = (0.0d0,1.0d0)*(Ex*conjg(d1x)+Ey*conjg(d1y))*exp( (0.0d0,1.0d0)*Delta_Ek1*x )
    H32 = (0.0d0,1.0d0)*(Ex*conjg(d2x)+Ey*conjg(d2y))*exp( (0.0d0,1.0d0)*Delta_Ek2*x )

	dwdx = H31*Y + H32*Z

end function funcdw


end subroutine RK_3level
