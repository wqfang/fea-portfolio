! Material-point driver: calls UMAT the way Abaqus/Standard does for one
! integration point. Reads C10, D1, then any number of deformation gradients
! (9 numbers, row by row) from stdin; prints STRESS(6) and DDSDDE(6,6) row-wise.
program driver
  implicit none
  real(8) :: stress(6), statev(1), ddsdde(6,6), sse, spd, scd, rpl, ddsddt(6), &
       drplde(6), drpldt, stran(6), dstran(6), time(2), dtime, temp, dtemp, &
       predef(1), dpred(1), props(2), coords(3), drot(3,3), pnewdt, celent, &
       dfgrd0(3,3), dfgrd1(3,3), f(9)
  integer :: jstep(4), ios, i
  character(len=80) :: cmname = 'NEOHOOKE'
  read(*,*) props(1), props(2)
  dfgrd0 = 0d0; do i = 1, 3; dfgrd0(i,i) = 1d0; end do
  do
    read(*,*,iostat=ios) f
    if (ios /= 0) exit
    dfgrd1 = transpose(reshape(f, [3,3]))
    stress = 0d0; ddsdde = 0d0; pnewdt = 1d0
    call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
         stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, 3, 3, 6, 0, &
         props, 2, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 1, 1, 1, 1, jstep, 1)
    write(*,'(6es24.16)') stress
    do i = 1, 6
      write(*,'(6es24.16)') ddsdde(i,:)
    end do
  end do
end program driver
