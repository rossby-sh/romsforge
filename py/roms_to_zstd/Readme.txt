프로그램 환경 설치 방법
------------------------------------------------
pip install --no-index --find-links=pkgs -r requirements.txt
------------------------------------------------
pip install --no-index -r requirements.txt


------------------------------------------------
프로그램 사용방법

    1. First step  : config.yml 파일에서 설정
    2. Second step : 프로그램 실행 

------------------------------------------------

****** config.yml 파일 설정 ****** 
경로 설정 
 > config.yml 파일을 열어서 설정을 변경

 > RDomcfgPath : 수행 모델의 domain (ROMS) 파일 경로를 작성
 > RomsRstHead : 수행 모델의 결과 폴더 경로를 작성
 > RomsRstTail :
    RomsRstHead 경로 안에 결과들이 생성되어 있을때, 파일의 tail 정보
    가 일치하는 파일들을 읽어들임. 
    ex ) RomsRstHead 밑에 
    avg_0001.nc 
    avg_0002.nc
    avg_0003.nc
    dif_0001.nc
    dif_0002.nc
    dif_0003.nc

    이런식으로 존재하는 경우에는 "avg_" 라고 입력하면, 
    avg_0001.nc, avg_0002.nc, avg_0003.nc 만 읽어들임.


  > RomsNamlist 에는 ROMS모델 수행시 사용한 OOOO.in 파일을 경로 작성.
    이 파일에서 아래 정보들을 캐치해서 값을 인식하게 됨.
    Vtransform  
    Vstretching 
    theta_s     
    theta_b     
    tcline      
    nlayer      

표준 영역 설정(Standard)
  > latmin : 표준 영역의 최소 위도
  > latmax : 표준 영역의 최대 위도
  > lonmin : 표준 영역의 최소 경도
  > lonmax : 표준 영역의 최대 경도
  > hreosl : 수평해상도(degree)
  > vdepth : 변환 하고자 하는 연직 레이어 깊이


추출대상 변수 선언(ExtraVar)
  > 여기에 추출 하고자 하는 변수를 추가하면 됨. 
  > 규칙은 
    ExtraVar: : 
      - 변수명 : [포인트정보, 차원정보]
      - 변수명 : [포인트정보, 차원정보]
      - 변수명 : [포인트정보, 차원정보]
    이런식으로 작성하면 됨.

  > ex) 만약 수온이면 아래처럼 작성
    - temp : [rho, 3]

  > rho 는 arakawa C grid 의 중심
  > u   는 arakawa C grid 의 u point
  > v   는 arakawa C grid 의 v point
  > zeta의 경우 2차원이므로 차원정보에 2

* 주의사항으로 모델의 결과는 4차원(Time, Depth, Lat, Lon)이 아닌
  (Depth, Lat, Lon)의 결과파일을 그대로 z로 변환하는 코드임.