! remap_pipeline.f90 : 완전한 remap 예제 (셀 자동 탐색 + weight 계산)
program remap_pipeline
  use module_remap
  implicit none

  real, allocatable :: src_lon(:,:), src_lat(:,:)
  real :: dst_lon, dst_lat
  integer :: base_i, base_j
  logical :: found
  real :: w11, w21, w12, w22

  integer :: i, j

  ! 1. 그리드 할당 및 생성
  allocate(src_lon(31,31), src_lat(31,31))
  do j = 1, 31
    do i = 1, 31
      src_lon(i,j) = 120.0 + real(j-1)  ! 경도 증가
      src_lat(i,j) = 10.0  + real(i-1)  ! 위도 증가
    end do
  end do

  ! 2. 테스트 지점 설정
  dst_lon = 135.25
  dst_lat = 23.75
  print *, 'start debug 1'
  print *, 'src_lon(1,1) =', src_lon(1,1)
  print *, 'src_lat(1,1) =', src_lat(1,1)
  ! 3. 포함 셀 탐색
  call find_enclosing_cell(dst_lon, dst_lat, src_lon, src_lat, base_i, base_j, found)
  print *, 'found =', found
  if (found) then
    print *, 'base_i, base_j =', base_i, base_j
  else
    print *, '셀을 찾지 못했습니다.'
    stop
  end if

  if (base_i >= 31 .or. base_j >= 31) then
    print *, 'base_i or base_j too large!'
    stop
  end if

  print *, '셀 위치: (i,j) =', base_i, base_j

  ! 4. weight 계산
  call compute_bilinear_weights( & 
      src_lon(base_i,base_j),     src_lat(base_i,base_j),      & ! 좌하단
      src_lon(base_i,base_j+1),   src_lat(base_i,base_j+1),    & ! 우하단
      src_lon(base_i+1,base_j),   src_lat(base_i+1,base_j),    & ! 좌상단
      src_lon(base_i+1,base_j+1), src_lat(base_i+1,base_j+1),  & ! 우상단
      dst_lon, dst_lat, &
      w11, w21, w12, w22)

  ! 5. 출력
  print *, 'Weights:'
  print *, 'w11 (좌하단) =', w11
  print *, 'w21 (우하단) =', w21
  print *, 'w12 (좌상단) =', w12
  print *, 'w22 (우상단) =', w22

  deallocate(src_lon, src_lat)

end program remap_pipeline


!---------------------------

