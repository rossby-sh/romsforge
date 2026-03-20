      MODULE biology_mod
!
!svn $Id: biology.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!  biology model.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: biology
      CONTAINS
      SUBROUTINE biology (ng,tile)
! !!! Modified by Dr.jhlee 
!svn $Id: npzd_iron.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Jerome Fiechter   !
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
        BIONAME(iNLM)="ROMS/Nonlinear/Biology/npzd_iron_vjhlee01.h"
      END IF
!
      CALL wclock_on (ng, iNLM, 15)
      CALL biology_tile (ng, tile,                                      &
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
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
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
      integer, parameter :: Nsink = 2
      integer :: Iter, i, ibio, isink, itime, itrc, iTrcMax, j, k, ks
      integer, dimension(Nsink) :: idsink
      real(r8), parameter :: MinVal = 1.0e-6_r8
      real(r8) :: Att, ExpAtt, Itop, PAR
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, dtdays
      real(r8) :: cffL, cffR, cu, dltL, dltR
      real(r8) :: fac
!  New: column-wise effective params based on depth h(i,j)
      real(r8) :: AttSW_col, PhyIS_col
      real(r8) :: wPhy_col, KFeC_col
      real(r8) :: FeMax_col, FeHmin_col
      real(r8) :: cff_col
      real(r8) :: Nlimit, FNlim
      real(r8) :: FNratio, FCratio, FCratioE, Flimit
      real(r8) :: FeC2FeN, FeN2FeC
      real(r8) :: cffFe
      real(r8) :: FeNudgCoef
      real(r8), dimension(Nsink) :: Wbio
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(NT(ng),2) :: BioTrc
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Light
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
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
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping size (days) according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set nudging coefficient for dissolved iron over the shelf.
!
      FeNudgCoef=dt(ng)/(FeNudgTime(ng)*86400.0_r8)
!
!  Set Fe:N and Fe:C conversion ratio and its inverse.
!
          FeN2FeC=(16.0_r8/106.0_r8)*1.0E3_r8
          FeC2FeN=(106.0_r8/16.0_r8)*1.0E-3_r8
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iphyt                 ! Phytoplankton
      idsink(2)=iSdet                 ! Small detritus
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! Phytoplankton
      Wbio(2)=wDet(ng)                ! Small detritus
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
              Bio_old(i,k,ibio)=BioTrc(ibio,nstp)
              Bio(i,k,ibio)=BioTrc(ibio,nstp)
            END DO
!  Relax dissolved iron at coast (h <= FeHim) to a constant value
!  (FeMax) over a time scale (FeNudgTime; days) to simulate sources
!  at the shelf.
!
!            IF (h(i,j).le.FeHmin(ng)) THEN
!              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                           &
!     &                       FeNudgCoef*(FeMax(ng)-Bio(i,k,iFdis))
!            END IF
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in a fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!  In the implicit algorithm, we have for example (N: nutrient,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!  Compute light attenuation as function of depth.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
            IF (h(i,j).le.50.0_r8) THEN
            AttSW_col=0.09_r8
            ELSE IF (h(i,j).gt.50.0_r8.and.h(i,j).le.100.0_r8) THEN
            AttSW_col=0.085_r8
            ELSE IF (h(i,j).gt.100.0_r8.and.h(i,j).le.150.0_r8) THEN
            AttSW_col=0.08_r8
            ELSE IF (h(i,j).gt.150.0_r8.and.h(i,j).le.2000.0_r8) THEN
            AttSW_col=0.07_r8
            ELSE IF (h(i,j).gt.2000.0_r8.and.h(i,j).le.3000.0_r8) THEN
            AttSW_col=0.068_r8
            ELSE IF (h(i,j).gt.3000.0_r8) THEN
            AttSW_col=0.065_r8
            ELSE 
            AttSW_col=0.065_r8
            END IF
                Att=(AttSW_col+AttPhy(ng)*Bio(i,k,iPhyt))*              &
     &              (z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                Light(i,k)=0.0_r8
              END DO
            END IF
          END DO
!
!  Phytoplankton photosynthetic growth and nitrate uptake (Vm_NO3 rate).
!  The Michaelis-Menten curve is used to describe the change in uptake
!  rate as a function of nitrate concentration. Here, PhyIS is the
!  initial slope of the P-I curve and K_NO3 is the half saturation of
!  phytoplankton nitrate uptake.
!
!  Growth reduction factors due to iron limitation:
!
!    FNratio     current Fe:N ratio [umol-Fe/mmol-N]
!    FCratio     current Fe:C ratio [umol-Fe/mol-C]
!                  (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1E3 mmol-C/mol-C)
!    FCratioE    empirical  Fe:C ratio
!    Flimit      Phytoplankton growth reduction factor due to Fe
!                  limitation based on Fe:C ratio
!
!
          DO k=1,N(ng)
            DO i=Istr,Iend
            IF (h(i,j).le.50.0_r8) THEN
            PhyIS_col=0.024_r8
            KFeC_col =15.0_r8
!            FeMax_col=0.002_r8
            FeMax_col=2.0_r8
            FeHmin_col=50.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE IF (h(i,j).gt.50.0_r8.and.h(i,j).le.100.0_r8) THEN
            PhyIS_col=0.024_r8
            KFeC_col =18.0_r8
!            FeMax_col=0.0016_r8
            FeMax_col=1.60_r8
            FeHmin_col=100.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE IF (h(i,j).gt.100.0_r8.and.h(i,j).le.150.0_r8) THEN
            PhyIS_col=0.023_r8
            KFeC_col =16.0_r8
!            FeMax_col=0.0014_r8
            FeMax_col=1.4_r8
            FeHmin_col=150.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE IF (h(i,j).gt.150.0_r8.and.h(i,j).le.2000.0_r8) THEN
            KFeC_col =20.0_r8
            PhyIS_col=0.022_r8
!            FeMax_col=0.0010_r8
            FeMax_col=1.0_r8
            FeHmin_col=200.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE IF (h(i,j).gt.2000.0_r8.and.h(i,j).le.3000.0_r8) THEN
            KFeC_col =16.9_r8
            PhyIS_col=0.021_r8
!            FeMax_col=0.0006_r8
            FeMax_col=0.60_r8
            FeHmin_col=200.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE IF (h(i,j).gt.3000.0_r8) THEN
            KFeC_col =14.9_r8
            PhyIS_col=0.020_r8
!            FeMax_col=0.0004_r8
            FeMax_col=0.4_r8
            FeHmin_col=300.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            ELSE 
            KFeC_col =14.9_r8
            PhyIS_col=0.020_r8
!            FeMax_col=0.0004_r8
            FeMax_col=0.40_r8
            FeHmin_col=300.0_r8
            IF (h(i,j).le.FeHmin_col) THEN
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       FeNudgCoef*(FeMax_col-Bio(i,k,iFdis))
            END IF
            END IF
          cff1=dtdays*Vm_NO3(ng)*PhyIS_col
          cff2=Vm_NO3(ng)*Vm_NO3(ng)
          cff3=PhyIS_col*PhyIS_col
!
!  Calculate growth reduction factor due to iron limitation.
!
              FNratio=Bio(i,k,iFphy)/MAX(MinVal,Bio(i,k,iPhyt))
              FCratio=FNratio*FeN2FeC
              FCratioE=B_Fe(ng)*Bio(i,k,iFdis)**A_Fe(ng)
              Flimit=FCratio*FCratio/                                   &
     &               (FCratio*FCratio+KFeC_col*KFeC_col)
!
              Nlimit=1.0_r8/(K_NO3(ng)+Bio(i,k,iNO3_))
              FNlim=MIN(1.0_r8,Flimit/(Bio(i,k,iNO3_)*Nlimit))
              cff4=1.0_r8/SQRT(cff2+cff3*Light(i,k)*Light(i,k))
              cff=Bio(i,k,iPhyt)*                                       &
     &            cff1*cff4*Light(i,k)*FNlim*Nlimit
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                            &
     &                       Bio(i,k,iNO3_)*cff
!
!  Iron uptake proportional to growth.
!
              fac=cff*Bio(i,k,iNO3_)*FNratio/MAX(MinVal,Bio(i,k,iFdis))
              Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+fac)
              Bio(i,k,iFphy)=Bio(i,k,iFphy)+                            &
     &                       Bio(i,k,iFdis)*fac
!
!  Iron uptake to reach appropriate Fe:C ratio.
!
              cff5=dtdays*(FCratioE-FCratio)/T_Fe(ng)
              cff6=Bio(i,k,iPhyt)*cff5*FeC2FeN
              IF (cff6.ge.0.0_r8) THEN
                cff=cff6/MAX(MinVal,Bio(i,k,iFdis))
                Bio(i,k,iFdis)=Bio(i,k,iFdis)/(1.0_r8+cff)
                Bio(i,k,iFphy)=Bio(i,k,iFphy)+                          &
     &                         Bio(i,k,iFdis)*cff
              ELSE
                cff=-cff6/MAX(MinVal,Bio(i,k,iFphy))
                Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
                Bio(i,k,iFdis)=Bio(i,k,iFdis)+                          &
     &                         Bio(i,k,iFphy)*cff
              END IF
            END DO
          END DO
!
!  Grazing on phytoplankton by zooplankton (ZooGR rate) using the Ivlev
!  formulation (Ivlev, 1955) and lost of phytoplankton to the nitrate
!  pool as function of "sloppy feeding" and metabolic processes
!  (ZooEEN and ZooEED fractions).
!  The lost of phytoplankton to the dissolve iron pool is scale by the
!  remineralization rate (FeRR).
!
          cff1=dtdays*ZooGR(ng)
          cff2=1.0_r8-ZooEEN(ng)-ZooEED(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,iZoop)*                                       &
     &            cff1*(1.0_r8-EXP(-Ivlev(ng)*Bio(i,k,iPhyt)))/         &
     &            Bio(i,k,iPhyt)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)/(1.0_r8+cff)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+                            &
     &                       Bio(i,k,iPhyt)*cff2*cff
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iPhyt)*ZooEEN(ng)*cff
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iPhyt)*ZooEED(ng)*cff
              Bio(i,k,iFphy)=Bio(i,k,iFphy)/(1.0_r8+cff)
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       Bio(i,k,iFphy)*cff*FeRR(ng)
            END DO
          END DO
!
!  Phytoplankton mortality to nutrients (PhyMRNro rate), detritus
!  (PhyMRD rate), and if applicable dissolved iron (FeRR rate).
!
          cff3=dtdays*PhyMRD(ng)
          cff2=dtdays*PhyMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iPhyt)*cff2
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iPhyt)*cff3
              Bio(i,k,iFphy)=Bio(i,k,iFphy)*cff1
              Bio(i,k,iFdis)=Bio(i,k,iFdis)+                            &
     &                       Bio(i,k,iFphy)*(cff2+cff3)*FeRR(ng)
            END DO
          END DO
!
!  Zooplankton mortality to nutrients (ZooMRN rate) and Detritus
!  (ZooMRD rate).
!
          cff3=dtdays*ZooMRD(ng)
          cff2=dtdays*ZooMRN(ng)
          cff1=1.0_r8/(1.0_r8+cff2+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iZoop)=Bio(i,k,iZoop)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iZoop)*cff2
              Bio(i,k,iSDet)=Bio(i,k,iSDet)+                            &
     &                       Bio(i,k,iZoop)*cff3
            END DO
          END DO
!
!  Detritus breakdown to nutrients: remineralization (DetRR rate).
!
          cff2=dtdays*DetRR(ng)
          cff1=1.0_r8/(1.0_r8+cff2)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,iSDet)=Bio(i,k,iSDet)*cff1
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+                            &
     &                       Bio(i,k,iSDet)*cff2
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms: Phytoplankton and Detritus
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
           DO i=Istr,Iend
              IF (ibio.eq.iphyt) THEN
                cff_col=dtdays*ABS(wPhy(ng))
              ELSE
                cff_col=dtdays*ABS(wDet(ng))
              END IF
            END DO
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
               IF (t(i,j,k,nstp,itemp).lt.10.0_r8) THEN
               Wbio(1)=wPhy(ng)/3.0_r8
               Wbio(2)=wDet(ng)/3.0_r8
               ELSE IF                                                  &
     &  (t(i,j,k,nstp,itemp).ge.10.0_r8.and.                            &
     &   t(i,j,k,nstp,itemp).le.25.0_r8) THEN
            Wbio(1)=wPhy(ng)/                                           &
     &  (3.0_r8-2.0_r8*((t(i,j,k,nstp,itemp)-10.0_r8)/15.0_r8))
            Wbio(2)=wDet(ng)/                                           &
     &  (3.0_r8-2.0_r8*((t(i,j,k,nstp,itemp)-10.0_r8)/15.0_r8))
               ELSE
                Wbio(1)=wPhy(ng)
                Wbio(2)=wDet(ng)
               END IF
               cff=dtdays*ABS(Wbio(isink))
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff_col
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP
      RETURN
      END SUBROUTINE biology_tile
      END MODULE biology_mod
