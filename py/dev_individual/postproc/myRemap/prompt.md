# Python prompt 

## 목적 
- 15km ROMS model 결과를 5km격자로 인터폴레이션하여 파일 생성

## 사용함수 
- netCDF4, numpy, xesmf, xarray 기반 + 알파

## 코딩 스타일
- 명시적 (표준 협업 스타일로)
- 확장가능
- 유지보수 쉽게
- 포트란 스타일
- 속도 우선
- 디버깅 요소 추가 (조건문 형태로)
- 객체지향
- 주석 필수
- 변수이름은 이해하기 쉽게

## 알고리즘 로직
- Inputs: 
  - grid01.nc, grid02.nc -> 그리드 이름 (type: string), 
  - 15km model average file directory (type: string)
  - output directory name for remapped grid file(type: string)
  - output directory name for weight file (type: string)

- Outpus:
  - interpolated grid file (netcdf file)
  - weight file (netcdf file)

- 함수
  - func read_grid(grid_name: str ) --> tuple: (lon_rho, lat_rho), (lon_u, lat_u), (lon_v, lat_v)  
    - open grid file and read grid with xarray
  - func remap(netcdf4.Dataset, numpy.array) --> numpy.array
    - netcdf4.Dataset -> weight객체 , netcdf4.Dataset -> remap 적용할 변수 객체 
    - 기능 변수 이름 및 차원 확인
    - 2차원->보간, 3차원 -> for문을 사용하여 수심에 따라 보간
    - Output -> remapped numpy array
  - func duplicate_ncStruct(input_file_name:str, output_file_name:str, grd_name:str) -> netcdf file 
    - input_file을 netcdf dataset으로 읽고 변수들과 global, variable attribute를 복사
    - grd_name을 읽고 xi_rho, eta_rho, xi_u,eta_u, xi_v, eta_v를 grd_name파일 기준으로 수정
    - 위 변수가 포함 된 데이터는 0으로 초기화

- 로직
  0. 작성자 및 코드 목적, 날짜, 이메일 정보
  1. 패키지 import 
  2. 변수 설명 및 선언 (그리드 이름, 디렉토리 이름 등)
  3. 함수 선언
  4. 그리드 읽기(xarray 사용)
  5. 그리드간 weight 계산
  6. weight 저장 (filename -> weight_"grid01"_to_"grid02".nc )
  7. Read avgerage netcdf file list (type: list)
  8. Load weight file 
  9. Subroutine01(nc_list,output_directory)
    - for nm in len(filelist)
      - func duplicate_ncStruct()로 output 파일 생성
      - netcdf4 Dataset으로 open file
      - 3D, 4D (혹은 2d,3d -> 확인필요)에 해당하는 변수 리스트 생성
        - 혹은 xi, eta를 차원으로 갖는 변수 확인
      - 변수별로 for문 돌려 func remap 적용
      - duplicate_ncStruc로 생성된 netcdf파일에 remap된 변수 저장

