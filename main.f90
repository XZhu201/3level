program main
use laserParameter
use atomicSystem
implicit none

integer, parameter :: Lt = 100000000
doubleprecision :: dt=tmax/Lt
integer :: nn
doubleprecision :: tt, Ex, Ey

doubleprecision, dimension(Lt) :: Ext, Eyt
doubleprecision, dimension(Lt) :: art

! calls
call show_laserParameter()
call show_atom()

! show electric field
open(10,file='Et.dat')

do nn = 0, Lt-1, 10000
    tt = dt*nn
    call fEt(tt, Ex, Ey)
    write(10,*) tt/Toc, Ex, Ey
end do
close(10)

! ode45
call RK_3level( 0,tmax,Lt,(0.0d0,0.0d0),(1.0d0,0.0d0),(0.0d0,0.0d0) )


end program main



