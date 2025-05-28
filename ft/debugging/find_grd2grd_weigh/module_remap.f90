module module_remap
  implicit none


contains

  function great_circle_distance(lon1,lat1,lon2,lat2) result(dist)
    implicit none
    real, intent(in) :: lon1, lat1, lon2, lat2
    real :: dist
    real :: dlon, dlat, a, c
    real, parameter :: pi = 3.141592
    real, parameter :: R = 6371.0 

    dlon = (lon2 - lon1) * pi / 180.0
    dlat = (lat2 - lat1) * pi / 180.0

    a = sin(dlat/2.0)**2 + cos(lat1*pi/180.0) * cos(lat2*pi/180.0) * sin(dlon/2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0-a))
    dist = R * c

  end function great_circle_distance

  subroutine find_nearest_point(dst_lon, dst_lat, src_lon, src_lat, nearest_i, nearest_j)
    implicit none
    real, intent(in) :: dst_lon, dst_lat
    real, intent(in) :: src_lon(:,:), src_lat(:,:)
    integer, intent(out) :: nearest_i, nearest_j

    integer :: i, j
    integer :: ilo, ihi, jlo, jhi
    real :: min_dist, dist

    ilo = lbound(src_lat,1)
    ihi = ubound(src_lat,1)
    jlo = lbound(src_lon,2)
    jhi = ubound(src_lon,2)

    min_dist = 1.0e30

    do j = jlo, jhi 
      do i = ilo, ihi
        dist = great_circle_distance(dst_lon, dst_lat, src_lon(i,j), src_lat(i,j) )
        if (dist < min_dist) then
          min_dist = dist
          nearest_i = i
          nearest_j = j
        end if
      end do
    end do

  end subroutine find_nearest_point

  subroutine compute_bilinear_weights(x1, y1, x2, y2, x3, y3, x4, y4, &
                                      dst_lon, dst_lat, &
                                      w11, w21, w12, w22)

    implicit none
    real, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
    real, intent(in) :: dst_lon, dst_lat
    real, intent(out) :: w11, w21, w12, w22
    real :: dx, dy, s, t

    dx = x2 - x1
    dy = y3 - y1

    s = (dst_lon - x1) / dx
    t = (dst_lat - y1) / dy

    w11 = (1.0 - s) * (1.0 - t)
    w21 = s * (1.0 - t)
    w12 = (1.0 - s) * t
    w22 = s * t
  end subroutine compute_bilinear_weights


  subroutine find_enclosing_cell(dst_lon, dst_lat, src_lon, src_lat, base_i, base_j, found)
    implicit none
    real, intent(in) :: dst_lon, dst_lat
    real, intent(in) :: src_lon(:,:), src_lat(:,:)
    integer, intent(out) :: base_i, base_j
    logical, intent(out) :: found
    integer :: i, j, imax, jmax

    imax = size(src_lon,1)
    jmax = size(src_lon,2)

    found = .false.

    do j = 1, jmax - 1
      do i = 1, imax - 1
        if (dst_lon >= minval([src_lon(i,j), src_lon(i+1,j), src_lon(i,j+1), src_lon(i+1,j+1)]) .and. &
            dst_lon <= maxval([src_lon(i,j), src_lon(i+1,j), src_lon(i,j+1), src_lon(i+1,j+1)]) .and. &
            dst_lat >= minval([src_lat(i,j), src_lat(i+1,j), src_lat(i,j+1), src_lat(i+1,j+1)]) .and. &
            dst_lat <= maxval([src_lat(i,j), src_lat(i+1,j), src_lat(i,j+1), src_lat(i+1,j+1)])) then
          base_i = i
          base_j = j
          found = .true.
          return
        end if
      end do
    end do

  end subroutine find_enclosing_cell

end module module_remap































