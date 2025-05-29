program test_remap
  use module_remap
  implicit none

  real, allocatable :: src_lon(:,:), src_lat(:,:)
  real :: dst_lon, dst_lat
  integer :: nearest_i, nearest_j
  integer :: i, j
  real :: x1, x2, y1, y2
  real :: w11, w21, w12, w22

  ! source grid 크기 설정 (31x31)
  allocate(src_lon(31,31), src_lat(31,31))

  ! source grid 초기화
  do j = 1, 31
    do i = 1, 31
      src_lon(i,j) = 120.0 + real(j-1)  ! j방향: 경도 증가
      src_lat(i,j) = 10.0  + real(i-1)  ! i방향: 위도 증가
    end do
  end do

  ! destination point 설정
  dst_lon = 135.5
  dst_lat = 23.5

  ! 가장 가까운 점 찾기
  call find_nearest_point(dst_lon, dst_lat, src_lon, src_lat, nearest_i, nearest_j)

  print*, 'Nearest (i,j):', nearest_i, nearest_j
  print*, 'Nearest (lon,lat):', src_lon(nearest_i,nearest_j), src_lat(nearest_i,nearest_j)

  ! 주변 4개 점 잡기
  ! 기준은 (i,j): 좌하단
  ! (i,j+1): 우하단
  ! (i+1,j): 좌상단
  ! (i+1,j+1): 우상단

  x1 = src_lon(nearest_i, nearest_j)     ! 좌하단
  x2 = src_lon(nearest_i, nearest_j+1)   ! 우하단
  y1 = src_lat(nearest_i, nearest_j)     ! 좌하단
  y2 = src_lat(nearest_i+1, nearest_j)   ! 좌상단

  ! weight 계산
  call compute_bilinear_weights(x1, x2, y1, y2, dst_lon, dst_lat, w11, w21, w12, w22)

  ! 결과 출력
  print*, 'Weights:'
  print*, 'w11 (좌하단) =', w11
  print*, 'w21 (우하단) =', w21
  print*, 'w12 (좌상단) =', w12
  print*, 'w22 (우상단) =', w22

  deallocate(src_lon, src_lat)

end program test_remap

