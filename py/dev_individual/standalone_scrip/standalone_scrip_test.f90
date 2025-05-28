! standalone_scrip_test.f90
! 간단한 SCRIP 실행 예제: compute_remap_weights_in 파일을 읽고 처리 시도

program standalone_scrip_test
  implicit none

  ! Fortran SCRIP 루틴 호출
  interface
    subroutine compute_remap_weights(config_file)
      character(len=*), intent(in) :: config_file
    end subroutine compute_remap_weights
  end interface

  character(len=256) :: config

  config = 'compute_remap_weights_in'

  print *, '>> Fortran SCRIP 테스트 시작'
  call compute_remap_weights(config)
  print *, '>> 종료'

end program standalone_scrip_test
