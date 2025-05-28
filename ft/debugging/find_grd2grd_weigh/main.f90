
program main
  use module_remap
  implicit none
  real, allocatable :: src_lon(:,:), src_lat(:,:), dst_lon(:,:), dst_lat(:,:)
  real :: tmp_dst_lon, tmp_dst_lat
  real :: w11, w21, w12, w22
  integer, parameter :: xi_dst = 61, eta_dst = 41, xi_src = 121, eta_src = 101
  integer :: i, j
  integer :: base_i, base_j
  integer :: unit
  logical :: found

  allocate( src_lon(eta_src, xi_src), src_lat(eta_src, xi_src) ) ! 120-135E, 20-30N, 0.25x0.25
  allocate( dst_lon(eta_dst, xi_dst), dst_lat(eta_dst, xi_dst) ) ! 100-160E, 5-55N, 0.5x0.5

  ! Define model grid sample
  do i = 1, eta_dst
    do j = 1, xi_dst
      dst_lon(i, j) = 120 + real( j - 1 ) * 0.25
      dst_lat(i, j) =  20 + real( i - 1 ) * 0.25
    end do
  end do

! Define OGCM grid sample
  do i = 1, eta_src
    do j = 1, xi_src
      src_lon(i, j) = 100 + real( j - 1 ) * 0.5
      src_lat(i, j) =   5 + real( i - 1 ) * 0.5
    end do
  end do

  ! open ASCII file for save weights
  unit=10
  open(unit=unit, file='weight.txt', status='unknown', action='write',form='formatted')
  !Define header
  write(unit, '(A)') 'dst_lon   dst_lat   w11   w21   w21   w22'  

  ! Find nearest cell
  do i = 1, eta_dst
    do j = 1, xi_dst

      tmp_dst_lon = dst_lon(i, j)
      tmp_dst_lat = dst_lat(i, j)  

      call find_enclosing_cell(tmp_dst_lon, tmp_dst_lat, src_lon, src_lat, base_i, base_j, found)  
      ! Debugging -> stop 1
      if (.not. found) then
        print *, '!!! Shell not found !!!'
        print *, i, j
        print *, src_lon(base_i, base_j), src_lat(base_i, base_j)
        print *, tmp_dst_lon, tmp_dst_lat
        stop 1
      end if
      
      ! Calulates weights
      call compute_bilinear_weights( &
          src_lon(base_i, base_j),    src_lat(base_i, base_j),    &
          src_lon(base_i, base_j+1),  src_lat(base_i, base_j+1),  &
          src_lon(base_i+1, base_j),  src_lat(base_i+1, base_j),  &
          src_lon(base_i+1, base_j+1),src_lat(base_i+1, base_j+1),&
          tmp_dst_lon, tmp_dst_lat, &          
          w11, w21, w12, w22 )

      write(unit, '(6F12.6)') tmp_dst_lon, tmp_dst_lat, w11, w21, w12, w22      

    end do
  end do
  
  close(unit)

end program main

