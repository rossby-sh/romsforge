      MODULE positive_mod
!
!svn $Id: positive.F 523 2011-01-05 03:21:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Although positive definite property is applied in the NPZD model,   !
!  negative values appear due to the advection scheme.                 !
!  Thus, let's apply positive definite property again after step3d_t.  !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: positive
      CONTAINS
      SUBROUTINE positive (ng,tile)
!
!svn $Id: positive.h 523 2011-01-05 03:21:38Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2011 The ROMS/TOMS Group       Jerome Fiechter   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Nutrient-Phytoplankton-Zooplankton-Detritus Model,                  !
!  including Iron Limitation on Phytoplankton Growth.                  !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fiechter, J., A.M. Moore, C.A. Edwards, K.W. Bruland,             !
!      E. Di Lorenzo, C.V.W. Lewis, T.M. Powell, E. Curchitser         !
!      and K. Hedstrom, 2009: Modeling iron limitation of primary      !
!      production in the coastal Gulf of Alaska, Deep Sea Res. II      !
!                                                                      ! 
!  Hajoon : subroutine that applies positive definite only             !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
!  Set header file name.
!
      IF (Lbiofile(iNLM)) THEN
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)="ROMS/Nonlinear/Biology/positive_npzd_iron.h"
      END IF
!
      CALL wclock_on (ng, iNLM, 15)
      CALL positive_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % t)
      CALL wclock_off (ng, iNLM, 15)
      RETURN
      END SUBROUTINE positive 
!
!-----------------------------------------------------------------------
      SUBROUTINE positive_tile (ng, tile,                               &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         rmask,                                   &
     &                         h,                                       &
     &                         Hz, z_r, z_w,                            &
     &                         srflx,                                   &
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology 
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: Iter, i, ibio, isink, itime, itrc, iTrcMax, j, k, ks
      real(r8), parameter :: MinVal = 1.0e-10_r8
      real(r8) :: Att, ExpAtt, Itop, PAR
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: fac
      real(r8) :: FeNudgCoef
      real(r8), dimension(NT(ng),2) :: BioTrc
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set nudging coefficient for dissolved iron over the shelf.
!
      FeNudgCoef=dt(ng)/(FeNudgTime(ng)*86400.0_r8)
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Restrict biological tracer to be positive definite. If a negative
!  concentration is detected, nitrogen is drawn from the most abundant
!  pool to supplement the negative pools to a lower limit of MinVal
!  which is set to 1E-6 above.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
              BioTrc(ibio,nstp)=t(i,j,k,nstp,ibio)
              BioTrc(ibio,nnew)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)
            END DO
!
!  Impose positive definite concentrations.
!
            cff2=0.0_r8
            DO itime=1,2
              cff1=0.0_r8
              iTrcMax=idbio(1)
              DO itrc=1,NBT-2
                ibio=idbio(itrc)
                cff1=cff1+MAX(0.0_r8,MinVal-BioTrc(ibio,itime))
                IF (BioTrc(ibio,itime).gt.BioTrc(iTrcMax,itime)) THEN
                  iTrcMax=ibio
                END IF
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
              IF (BioTrc(iTrcMax,itime).gt.cff1) THEN
                BioTrc(iTrcMax,itime)=BioTrc(iTrcMax,itime)-cff1
              END IF
              DO itrc=NBT-1,NBT
                ibio=idbio(itrc)
                BioTrc(ibio,itime)=MAX(MinVal,BioTrc(ibio,itime))
              END DO
            END DO
!
!  Load biological tracers into local arrays.
!
            DO itrc=1,NBT
              ibio=idbio(itrc)
              Bio(i,k,ibio)=BioTrc(ibio,nnew)
            END DO
!
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
            IF (h(i,j).le.FeHmin(ng)) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,ibio)=Bio(i,k,ibio)*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP
      RETURN
      END SUBROUTINE positive_tile
      END MODULE positive_mod
