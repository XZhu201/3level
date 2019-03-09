module laserParameter
implicit none

    ! intensity and wavelength
   doubleprecision, parameter :: E0 = 0.037677686771342
   doubleprecision, parameter :: E1 = 0.075355373542683
   doubleprecision, parameter :: w0 = 0.001627269045612 
   doubleprecision, parameter :: w1 = 0.003254538091224

    ! pulse shape
   doubleprecision, parameter :: pi=3.1415926536
   doubleprecision, parameter :: Ncycle = 100
   doubleprecision, parameter :: Toc=2*pi/w0
   doubleprecision, parameter :: tmax = Ncycle*Toc

   doubleprecision, parameter :: T1 = 40*Toc
   doubleprecision, parameter :: T2 = (Ncycle-40)*Toc


contains
   subroutine show_laserParameter()

      print*, " ----- show laser parameters -----"
      print*, "E0=", E0, "E1=", E1
      print*, "w0=", w0, "w1=", w1
      print*, "T=", Toc
      print*, "T1=", T1, "T2=", T2
      print*

   end subroutine show_laserParameter


   subroutine fEt(t,Ext,Eyt)
        doubleprecision, intent(in) :: t
        doubleprecision, intent(out) :: Ext,Eyt
        doubleprecision :: ft

        if(t<=T1) then
            ft = sin(pi/2*t/T1)**2
        elseif(t<=T2) then
            ft = 1
        elseif(t<=T1+T2) then
            ft = cos( pi/2*(t-T2)/T1 )**2
        else
            ft = 0
        end if

        ! print*, "t=",t,"ft=",ft

        Ext = E0*cos(w0*t)*ft + E1*cos(w1*t)*ft
        Eyt = E0*sin(w0*t)*ft - E1*sin(w1*t)*ft

   end subroutine fEt

end module laserParameter
