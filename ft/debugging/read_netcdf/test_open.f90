program test_open

  use netcdf
  implicit none

  integer :: ncid, status
  character(len=100) :: filename

  filename = '/mnt/c/Users/shjo9/Bridge/DATA/roms_grd_fennel_5km_smooth_v2.nc'

  status = nf90_open(trim(filename), NF90_NOWRITE, ncid )
  if (status /= NF90_NOERR) then
    print *, "Error opening file:", trim(nf90_strerror(status))
  else
    print *, "File opened successfully ! ncid = ", ncid
  end if

  status = nf90_close(ncid)

end program test_open
