program test_open

  use netcdf
  implicit none

  integer :: ncid, status
  integer :: ndims, nvars, ngatts, unlimdimid
  integer :: i, varid, xtype, ndims_var, natts
  integer, dimension(:), allocatable :: dimids
  character(len=256) :: varname
  character(len=100) :: filename

  filename = '/mnt/c/Users/shjo9/Bridge/DATA/test_clean.nc'

  status = nf90_open(trim(filename), NF90_NOWRITE, ncid )
  if (status /= nf90_noerr) stop 'ERROR opening file'

  status = nf90_inquire(ncid, ndims, nvars, ngatts, unlimdimid)
  if (status /= nf90_noerr) stop 'ERROR in nf90_inquire'
  
  print *, "Number of dimensions: ", ndims
  print *, "Number of variables: ", nvars
  print *, "Number of global attributes: ", ngatts

  status = nf90_inq_varid(ncid, "alpha", varid)
  print *, 'varid for alpha =', varid
! 이 방법으로 varid를 구해서 써야 함
  status = nf90_inq_varid(ncid, "alpha", varid)
  if (status /= nf90_noerr) stop 'ERROR in inquire_variable (step 1)'
  status = nf90_inquire_variable(ncid, varid, varname, xtype, ndims_var, natts=natts)
  if (status /= nf90_noerr) stop 'ERROR in inquire_variable (step 2)'
  
  do i=1, nvars
    varid = i - 1 ! 0 starts netcdf index
    
    status = nf90_inquire_variable(ncid, varid, varname, xtype, ndims_var, natts=natts)
    if (status /= nf90_noerr) then
      print *, 'varid=', varid, '  status=', status
      print *, 'Error: ', trim(nf90_strerror(status))
      stop 'ERROR in inquire_variable (step 3)'
    end if

    allocate(dimids(ndims_var))

    status = nf90_inquire_variable(ncid, varid, varname, xtype, ndims_var, dimids, natts)
    if (status /= nf90_noerr) stop 'ERROR in inquire_variable (step 4)'

    print *, "Variables ", i, ":", trim(varname)
    print *, " - Type: ", xtype
    print *, " - #Dims: ", ndims_var
    print *, " - #Atts: ", natts

    deallocate(dimids)
  end do


  status = nf90_close(ncid)

end program test_open
