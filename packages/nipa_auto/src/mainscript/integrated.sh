#!/usr/bin/bash
통합 실행 코드 

###
# pipe01 -> 변수 및 함수 설정, template들 변환 
# pipe02 -> 데이터 다운 및 전처리 & 진행사항 로그로 표시
# pipe03 -> pipe02 작업 실패시 알림 + 재실행 규칙
#
#
#
cron -

source 파이프 1

./


[pipe01 로직]
- 환경변수 및 함수 설정
- 전처리 탬플릿 및 모델  namelist 탬플릿 변환
- export 환경변수

[pipe02 로직] 
- HYCOM 다운
- CMEMS 생물변수 다운
- GFS 다운 (미구현)
- HYCOM으로 initial 제작 및 cmems 생물변수 추가 (만약 pipe01에서 USE_RST=TRUE이면 initial 제작은 스킵)
- HYCOM으로 boundray 제작 및 cmems 생물변수 추가
- GFS로 기상장 제작 (미구현)


[pipe03 로직]
- pipe log에 failed있으면 3시간 후 다시 pipe02 실행  
- pipe log에 failed있으면 다시 3시간 후 pipe02 실행  
- 여전히 failed 있으면 알람과 다음 로직 실행.
  - Initial은 만약 USE_RST=FALSE일 경우 하루전 HYCOM, CMEMS 파일을 사용 (예측장 받는거라 데이터는 있음) 
  - bry도 마찬가지로  하루전 파일을 사용
  - forcing파일은 ${INPUT_ROOT}/fixed 에 있는 장기 평균장을 사용
  - 메일(가능하다면) 로 해당 내용 발송


[pipe04 로직]
- pbs파일을 EOF로 제작  
- 모델 실행 (run pbs)

[pipe05 로직]
- 모델 success/blow-up여부 조사
- 모델 결과 후처리
- 해당 내용 발송



[추가예정]
- HYCOM 데이터에 이상있는지 확인



















