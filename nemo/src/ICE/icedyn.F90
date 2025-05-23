MODULE icedyn
   !!======================================================================
   !!                     ***  MODULE  icedyn  ***
   !!   Sea-Ice dynamics : master routine for sea ice dynamics
   !!======================================================================
   !! history :  4.0  ! 2018  (C. Rousset)  original code SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn       : dynamics of sea ice
   !!   ice_dyn_init  : initialization and namelist read
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_icedyn     ! SI3 dynamics parameters
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE icedyn_rhg     ! sea-ice: rheology
   USE icedyn_adv     ! sea-ice: advection
   USE icedyn_rdgrft  ! sea-ice: ridging/rafting
   USE icecor         ! sea-ice: corrections
   USE icevar  , ONLY : ice_var_zapsmall
   USE icectl         ! sea-ice: control prints
   USE icefsd  , ONLY : a_ifsd   ! sea-ice: floe size distribution
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! to use sign with key_nosignedzero
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing
   USE fldread        ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn        ! called by icestp.F90
   PUBLIC   ice_dyn_init   ! called by icestp.F90

   INTEGER ::              nice_dyn   ! choice of the type of dynamics
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_dynALL     = 1   ! full ice dynamics               (rheology + advection + ridging/rafting + correction)
   INTEGER, PARAMETER ::   np_dynRHGADV  = 2   ! pure dynamics                   (rheology + advection)
   INTEGER, PARAMETER ::   np_dynADV1D   = 3   ! only advection 1D - test case from Schar & Smolarkiewicz 1996
   INTEGER, PARAMETER ::   np_dynADV2D   = 4   ! only advection 2D w prescribed vel.(rn_uvice + advection)
   !
   ! ** namelist (namdyn) **
   LOGICAL  ::   ln_dynALL        ! full ice dynamics                      (rheology + advection + ridging/rafting + correction)
   LOGICAL  ::   ln_dynRHGADV     ! no ridge/raft & no corrections         (rheology + advection)
   LOGICAL  ::   ln_dynADV1D      ! only advection in 1D w ice convergence (test case from Schar & Smolarkiewicz 1996)
   LOGICAL  ::   ln_dynADV2D      ! only advection in 2D w prescribed vel. (rn_uvice + advection)
   REAL(wp) ::   rn_uice          !    prescribed u-vel (case np_dynADV1D & np_dynADV2D)
   REAL(wp) ::   rn_vice          !    prescribed v-vel (case np_dynADV1D & np_dynADV2D)

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_icbmsk  ! structure of input grounded icebergs mask (file informations, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_fastmsk ! structure of input landfast ice mask      (file informations, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn( kt, Kmm )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn  ***
      !!
      !! ** Purpose :   this routine manages sea ice dynamics
      !!
      !! ** Action : - calculation of friction in case of landfast ice
      !!             - call ice_dyn_rhg    = rheology
      !!             - call ice_dyn_adv    = advection
      !!             - call ice_dyn_rdgrft = ridging/rafting
      !!             - call ice_cor        = corrections if fields are out of bounds
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      INTEGER, INTENT(in) ::   Kmm    ! ocean time level index
      !!
      INTEGER  ::   ji, jj        ! dummy loop indices
      REAL(wp) ::   zcoefu, zcoefv
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zdivu_i
      !!--------------------------------------------------------------------
      !
      ! controls
      IF( ln_timing )   CALL timing_start('ice_dyn')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn: sea-ice dynamics'
         WRITE(numout,*)'~~~~~~~'
      ENDIF
      !
      ! retrieve thickness from volume for landfast param. and UMx advection scheme
      WHERE( a_i(:,:,:) >= epsi20 )
         h_i(:,:,:) = v_i(:,:,:) / a_i(:,:,:)
         h_s(:,:,:) = v_s(:,:,:) / a_i(:,:,:)
      ELSEWHERE
         h_i(:,:,:) = 0._wp
         h_s(:,:,:) = 0._wp
      END WHERE
      !
      WHERE( a_ip(:,:,:) >= epsi20 )
         h_ip(:,:,:) = v_ip(:,:,:) / a_ip(:,:,:)
         h_il(:,:,:) = v_il(:,:,:) / a_ip(:,:,:)
      ELSEWHERE
         h_ip(:,:,:) = 0._wp
         h_il(:,:,:) = 0._wp
      END WHERE
      !
      ! read grounded icebergs and landfast masks
      ! ----------------------------------------
      CALL fld_read( kt, 1, sf_icbmsk  )
      CALL fld_read( kt, 1, sf_fastmsk )
      icb_tmask (:,:) = sf_icbmsk (1)%fnow(:,:,1)
      fast_tmask(:,:) = sf_fastmsk(1)%fnow(:,:,1)
      !
      ! mask at u-, v-points (computed from tmask)
      DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
         fast_umask(ji,jj) = MAX( fast_tmask(ji,jj), fast_tmask(ji+1,jj  ) )
         icb_umask (ji,jj) = MAX( icb_tmask (ji,jj), icb_tmask (ji+1,jj  ) )
         fast_vmask(ji,jj) = MAX( fast_tmask(ji,jj), fast_tmask(ji  ,jj+1) )
         icb_vmask (ji,jj) = MAX( icb_tmask (ji,jj), icb_tmask (ji  ,jj+1) )
      END_2D

      SELECT CASE( nice_dyn )          !-- Set which dynamics is running

      CASE ( np_dynALL )           !==  all dynamical processes  ==!
         !
         CALL ice_dyn_rhg   ( kt, Kmm )                                     ! -- rheology
         CALL ice_dyn_adv   ( kt )                                          ! -- advection of ice
         CALL ice_dyn_rdgrft( kt )                                          ! -- ridging/rafting
         CALL ice_cor       ( kt , 1 )                                      ! -- Corrections
         !
      CASE ( np_dynRHGADV  )       !==  no ridge/raft & no corrections ==!
         !
         CALL ice_dyn_rhg   ( kt, Kmm )                                     ! -- rheology
         CALL ice_dyn_adv   ( kt )                                          ! -- advection of ice
         CALL Hpiling                                                       ! -- simple pile-up (replaces ridging/rafting)
         CALL ice_var_zapsmall                                              ! -- zap small areas
         !
      CASE ( np_dynADV1D )         !==  pure advection ==!   (1D)
         !
         ! --- monotonicity test from Schar & Smolarkiewicz 1996 --- !
         ! CFL = 0.5 at a distance from the bound of 1/6 of the basin length
         ! Then for dx = 2m and dt = 1s => rn_uice = u (1/6th) = 1m/s
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zcoefu = ( REAL(jpiglo+1)*0.5_wp - REAL(ji+nimpp-1) ) / ( REAL(jpiglo+1)*0.5_wp - 1._wp )
            zcoefv = ( REAL(jpjglo+1)*0.5_wp - REAL(jj+njmpp-1) ) / ( REAL(jpjglo+1)*0.5_wp - 1._wp )
            u_ice(ji,jj) = rn_uice * 1.5_wp * SIGN( 1.0_wp, zcoefu ) * ABS( zcoefu ) * umask(ji,jj,1)
            v_ice(ji,jj) = rn_vice * 1.5_wp * SIGN( 1.0_wp, zcoefv ) * ABS( zcoefv ) * vmask(ji,jj,1)
         END_2D
         ! ---
         CALL ice_dyn_adv   ( kt )                                          ! -- advection of ice
         CALL ice_var_zapsmall                                              ! -- zap small areas
         !
      CASE ( np_dynADV2D )         !==  pure advection ==!   (2D w prescribed velocities)
         !
         u_ice(:,:) = rn_uice * umask(:,:,1)
         v_ice(:,:) = rn_vice * vmask(:,:,1)
         !CALL RANDOM_NUMBER(u_ice(:,:)) ; u_ice(:,:) = u_ice(:,:) * 0.1 + rn_uice * 0.9 * umask(:,:,1)
         !CALL RANDOM_NUMBER(v_ice(:,:)) ; v_ice(:,:) = v_ice(:,:) * 0.1 + rn_vice * 0.9 * vmask(:,:,1)
         ! ---
         CALL ice_dyn_adv   ( kt )                                          ! -- advection of ice
         CALL ice_var_zapsmall                                              ! -- zap small areas

      END SELECT
      !
      !
      ! diagnostics: divergence at T points
      IF( iom_use('icediv') ) THEN
         !
         SELECT CASE( nice_dyn )

         CASE ( np_dynADV1D , np_dynADV2D )

            ALLOCATE( zdivu_i(A2D(0)) )
            DO_2D( 0, 0, 0, 0 )
               zdivu_i(ji,jj) = ( ( e2u(ji,jj) * u_ice(ji,jj) - e2u(ji-1,jj) * u_ice(ji-1,jj) )   &   ! add () for NP repro
                  &             + ( e1v(ji,jj) * v_ice(ji,jj) - e1v(ji,jj-1) * v_ice(ji,jj-1) ) ) * r1_e1e2t(ji,jj)
            END_2D
            CALL iom_put( 'icediv' , zdivu_i )
            DEALLOCATE( zdivu_i )

         END SELECT
         !
      ENDIF
      !
      ! --- Lateral boundary conditions --- !
      !     caution: t_su update needed from itd_reb
      !              plus, one needs ldfull=T to deal with the NorthFold in case of Prather advection
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         CALL lbc_lnk( 'icedyn', a_i , 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp, oa_i, 'T', 1._wp, &
            &                    t_su, 'T', 1._wp, a_ip, 'T', 1._wp, v_ip, 'T', 1._wp, v_il, 'T', 1._wp, ldfull = .TRUE. )
      ELSE
         CALL lbc_lnk( 'icedyn', a_i , 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp, oa_i, 'T', 1._wp, &
            &                    t_su, 'T', 1._wp, ldfull = .TRUE. )
      ENDIF
      CALL lbc_lnk( 'icedyn', e_i, 'T', 1._wp, e_s, 'T', 1._wp, szv_i, 'T', 1._wp, ldfull = .TRUE. )
      IF( ln_fsd ) THEN
         CALL lbc_lnk( 'icedyn', a_ifsd, 'T', 1._wp, ldfull = .TRUE. )
      ENDIF
      !
      ! controls
      IF( ln_timing )   CALL timing_stop ('ice_dyn')
      !
   END SUBROUTINE ice_dyn


   SUBROUTINE Hpiling
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hpiling  ***
      !!
      !! ** Purpose : Simple conservative piling comparable with 1-cat models
      !!
      !! ** Method  : pile-up ice when no ridging/rafting
      !!
      !! ** input   : a_i
      !!-------------------------------------------------------------------
      INTEGER ::   jl         ! dummy loop indices
      !!-------------------------------------------------------------------
      ! controls
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'Hpiling', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
      at_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )
      DO jl = 1, jpl
         WHERE( at_i(A2D(0)) > epsi20 )
            a_i(A2D(0),jl) = a_i(A2D(0),jl) * (  1._wp + MIN( rn_amax_2d(A2D(0)) - at_i(A2D(0)) , 0._wp ) / at_i(A2D(0))  )
         END WHERE
      END DO
      !
      ! controls
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'Hpiling', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
   END SUBROUTINE Hpiling


   SUBROUTINE ice_dyn_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio, ierror   ! Local integer output status for namelist read
      !
      CHARACTER(len=256) ::   cn_dir     ! Root directory for location of ice files
      TYPE(FLD_N)        ::   sn_icbmsk  ! informations about the grounded icebergs field to be read
      TYPE(FLD_N)        ::   sn_fastmsk ! informations about the landfast ice field to be read
      !!
      NAMELIST/namdyn/ ln_dynALL, ln_dynRHGADV, ln_dynADV1D, ln_dynADV2D, rn_uice, rn_vice,  &
         &             rn_ishlat ,                                                           &
         &             ln_landfast_L16, rn_lf_depfra, rn_lf_bfr, rn_lf_relax, rn_lf_tensile, &
         &             sn_icbmsk, sn_fastmsk, cn_dir
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn)
      READ_NML_CFG(numnam_ice,namdyn)
      IF(lwm) WRITE( numoni, namdyn )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn:'
         WRITE(numout,*) '      Full ice dynamics      (rhg + adv + ridge/raft + corr) ln_dynALL       = ', ln_dynALL
         WRITE(numout,*) '      No ridge/raft & No cor (rhg + adv)                     ln_dynRHGADV    = ', ln_dynRHGADV
         WRITE(numout,*) '      Advection 1D only      (Schar & Smolarkiewicz 1996)    ln_dynADV1D     = ', ln_dynADV1D
         WRITE(numout,*) '      Advection 2D only      (rn_uvice + adv)                ln_dynADV2D     = ', ln_dynADV2D
         WRITE(numout,*) '         with prescribed velocity given by   (u,v)_ice = (rn_uice,rn_vice)   = (', rn_uice,',',rn_vice,')'
         WRITE(numout,*) '      lateral boundary condition for sea ice dynamics        rn_ishlat       = ', rn_ishlat
         WRITE(numout,*) '      Landfast: param from Lemieux 2016                      ln_landfast_L16 = ', ln_landfast_L16
         WRITE(numout,*) '         fraction of ocean depth that ice must reach         rn_lf_depfra    = ', rn_lf_depfra
         WRITE(numout,*) '         maximum bottom stress per unit area of contact      rn_lf_bfr       = ', rn_lf_bfr
         WRITE(numout,*) '         relax time scale (s-1) to reach static friction     rn_lf_relax     = ', rn_lf_relax
         WRITE(numout,*) '         isotropic tensile strength                          rn_lf_tensile   = ', rn_lf_tensile
         WRITE(numout,*)
      ENDIF
      !                             !== set the choice of ice dynamics ==!
      ioptio = 0
      !      !--- full dynamics                               (rheology + advection + ridging/rafting + correction)
      IF( ln_dynALL    ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynALL       ;   ENDIF
      !      !--- dynamics without ridging/rafting and corr   (rheology + advection)
      IF( ln_dynRHGADV ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynRHGADV    ;   ENDIF
      !      !--- advection 1D only - test case from Schar & Smolarkiewicz 1996
      IF( ln_dynADV1D  ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynADV1D     ;   ENDIF
      !      !--- advection 2D only with prescribed ice velocities (from namelist)
      IF( ln_dynADV2D  ) THEN   ;   ioptio = ioptio + 1   ;   nice_dyn = np_dynADV2D     ;   ENDIF
      !
      IF( ioptio /= 1 )    CALL ctl_stop( 'ice_dyn_init: one and only one ice dynamics option has to be defined ' )
      !
      !                                      !--- Lateral boundary conditions
      IF     (      rn_ishlat == 0.                ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  free-slip'
      ELSEIF (      rn_ishlat == 2.                ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  no-slip'
      ELSEIF ( 0. < rn_ishlat .AND. rn_ishlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  partial-slip'
      ELSEIF ( 2. < rn_ishlat                      ) THEN   ;   IF(lwp) WRITE(numout,*) '   ===>>>   ice lateral  strong-slip'
      ENDIF
      !                                      !--- Landfast ice
      IF( .NOT.ln_landfast_L16 )   tau_icebfr(:,:) = 0._wp
      !
      !                                      !--- allocate and fill structure for grounded icebergs and landfast masks
      ALLOCATE( sf_icbmsk(1), sf_fastmsk(1), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'ice_dyn_init: unable to allocate sf_icbmsk and sf_fastmsk structures' ) ; RETURN
      ENDIF
      !
      CALL fld_fill( sf_icbmsk , (/ sn_icbmsk /) , cn_dir, 'ice_dyn_init',   &
         &                                                 'landfast ice is a function of read grounded icebergs', 'icedyn' )
      CALL fld_fill( sf_fastmsk, (/ sn_fastmsk /), cn_dir, 'ice_dyn_init',   &
         &                                                 'landfast ice is read in a file', 'icedyn' )
      !
      ALLOCATE( sf_icbmsk(1)%fnow(jpi,jpj,1), sf_fastmsk(1)%fnow(jpi,jpj,1) )
      IF( sf_icbmsk (1)%ln_tint )   ALLOCATE( sf_icbmsk (1)%fdta(jpi,jpj,1,2) )
      IF( sf_fastmsk(1)%ln_tint )   ALLOCATE( sf_fastmsk(1)%fdta(jpi,jpj,1,2) )
      IF( TRIM(sf_icbmsk (1)%clrootname) == 'NOT USED' )   sf_icbmsk (1)%fnow(:,:,1) = 0._wp   ! not used field  (set to 0)
      IF( TRIM(sf_fastmsk(1)%clrootname) == 'NOT USED' )   sf_fastmsk(1)%fnow(:,:,1) = 0._wp   ! not used field  (set to 0)
      !
      !                                      !--- other init
      CALL ice_dyn_rdgrft_init          ! set ice ridging/rafting parameters
      CALL ice_dyn_rhg_init             ! set ice rheology parameters
      CALL ice_dyn_adv_init             ! set ice advection parameters
      !
   END SUBROUTINE ice_dyn_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn
