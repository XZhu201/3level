module atomicSystem
implicit none

   ! d1x d1y d2x d2y Delta_Ek1 Delta_Ek2

   doubleprecision, parameter :: Energy1 = -0.793180324283556
   doubleprecision, parameter :: Energy2 = -0.793180324283556
   doubleprecision, parameter :: Energys = -0.216647750753502

   doubleprecision, parameter :: Delta_Ek1 = Energys-Energy1
   doubleprecision, parameter :: Delta_Ek2 = Energys-Energy2

   complex :: d1x = (0.462975910042748, 0.0)     !  d_ps_x
   complex :: d1y = (0.0, -0.462975910042436)    !  d_ps_y
   complex :: d2x = (0.462975910042748, 0.0)     !  d_ms_x
   complex :: d2y = (0.0, 0.462975910042436)     !  d_ms_y

contains
   subroutine show_atom()

      print*, " ----- show atomic parameters -----"
      print*, "Energy1=", Energy1, "Energy2=", Energy2, "Energys=", Energys
      print*, "Delta_Ek1=", Delta_Ek1, "Delta_Ek2=", Delta_Ek2
      print*, "d1x=", d1x, "d1y", d1y
      print*, "d2x=", d2x, "d2y", d2y
      print*

   end subroutine show_atom

end module atomicSystem
