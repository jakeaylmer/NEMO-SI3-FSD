MODULE sbc_oce
   !!======================================================================
   !!                       ***  MODULE  sbc_oce  ***
   !! Surface module :   variables defined in core memory
   !!======================================================================
   !! History :  3.0  ! 2006-06  (G. Madec)  Original code
   !!             -   ! 2008-08  (G. Madec)  namsbc moved from sbcmod
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!            3.3  ! 2010-10  (J. Chanut, C. Bricaud)  add the surface pressure forcing
   !!            4.0  ! 2012-05  (C. Rousset) add attenuation coef for use in ice model
   !!            4.0  ! 2016-06  (L. Brodeau) new unified bulk routine (based on AeroBulk)
   !!            4.0  ! 2019-03  (F. Lemarié, G. Samson) add compatibility with ABL mode
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) modified wave parameters in namelist
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_oce_alloc : allocation of sbc arrays
   !!   sbc_tau2wnd   : wind speed estimated from wind stress
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_oce_alloc   ! routine called in sbcmod.F90
   PUBLIC   sbc_oce_dealloc ! routine called in nemogcm.F90
   PUBLIC   sbc_tau2wnd     ! routine called in several sbc modules

   !!----------------------------------------------------------------------
   !!           Namelist for the Ocean Surface Boundary Condition
   !!----------------------------------------------------------------------
   !                                   !!* namsbc namelist *
   LOGICAL , PUBLIC ::   ln_usr         !: user defined formulation
   LOGICAL , PUBLIC ::   ln_flx         !: flux      formulation
   LOGICAL , PUBLIC ::   ln_blk         !: bulk formulation
   LOGICAL , PUBLIC ::   ln_abl         !: Atmospheric boundary layer model
   LOGICAL , PUBLIC ::   ln_wave        !: wave in the system (forced or coupled)
   LOGICAL , PUBLIC ::   ln_cpl         !: ocean-atmosphere coupled formulation
   LOGICAL , PUBLIC ::   ln_mixcpl      !: ocean-atmosphere forced-coupled mixed formulation
   LOGICAL , PUBLIC ::   ln_dm2dc       !: Daily mean to Diurnal Cycle short wave (qsr)
   LOGICAL , PUBLIC ::   ln_rnf         !: runoffs / runoff mouths
   LOGICAL , PUBLIC ::   ln_ssr         !: Sea Surface restoring on SST and/or SSS
   LOGICAL , PUBLIC ::   ln_apr_dyn     !: Atmospheric pressure forcing used on dynamics (ocean & ice)
   INTEGER , PUBLIC ::   nn_ice         !: flag for ice in the surface boundary condition (=0/1/2/3)
   LOGICAL , PUBLIC ::   ln_ice_embd    !: flag for levitating/embedding sea-ice in the ocean
   !                                             !: =F levitating ice (no presure effect) with mass and salt exchanges
   !                                             !: =T embedded sea-ice (pressure effect + mass and salt exchanges)
   INTEGER , PUBLIC ::   nn_components  !: flag for sbc module (including sea-ice) coupling mode (see component definition below)
   INTEGER , PUBLIC ::   nn_fwb         !: FreshWater Budget:
   !                                             !:  = 0 unchecked
   !                                             !:  = 1 global mean of e-p-r set to zero at each nn_fsbc time step
   !                                             !:  = 2 annual global mean of e-p-r set to zero
   LOGICAL , PUBLIC ::   ln_icebergs    !: Icebergs
   !
   INTEGER , PUBLIC ::   nn_lsm         !: Number of iteration if seaoverland is applied
   !
   !                                   !!* namsbc_cpl namelist *
   INTEGER , PUBLIC ::   nn_cats_cpl    !: Number of sea ice categories over which the coupling is carried out
   !
   !                                   !!* namsbc_wave namelist *
   LOGICAL , PUBLIC ::   ln_sdw         !: =T 3d stokes drift from wave model
   LOGICAL , PUBLIC ::   ln_stcor       !: =T if Stokes-Coriolis and tracer advection terms are used
   LOGICAL , PUBLIC ::   ln_cdgw        !: =T neutral drag coefficient from wave model
   LOGICAL , PUBLIC ::   ln_tauoc       !: =T if normalized stress from wave is used
   LOGICAL , PUBLIC ::   ln_wave_test   !: =T wave test case (constant Stokes drift)
   LOGICAL , PUBLIC ::   ln_charn       !: =T Chranock coefficient from wave model
   LOGICAL , PUBLIC ::   ln_taw         !: =T wind stress corrected by wave intake
   LOGICAL , PUBLIC ::   ln_phioc       !: =T TKE surface BC from wave model
   LOGICAL , PUBLIC ::   ln_bern_srfc   !: Bernoulli head, waves' inuced pressure
   LOGICAL , PUBLIC ::   ln_breivikFV_2016 !: Breivik 2016 profile
   LOGICAL , PUBLIC ::   ln_vortex_force !: vortex force activation
   LOGICAL , PUBLIC ::   ln_stshear     !: Stoked Drift shear contribution in zdftke
   !
   !                                                     !!* namsbc_wave namelist (continued) *
   LOGICAL , PUBLIC                 ::   ln_wave_spec     !: =T to read full wave spectrum
   INTEGER , PUBLIC                 ::   nn_nwfreq        !: number of spectral (frequency) components (discretised bins)
   LOGICAL , PUBLIC                 ::   ln_wfreq_usr     !: =T if frequency bin bounds are defined by rn_wfreq_usr
   REAL    , PUBLIC, DIMENSION(100) ::   rn_wfreq_usr     !: user-defined frequency bin bounds; up to 100 (arbitrary limit) can be specified
   LOGICAL , PUBLIC                 ::   ln_wfreq_usr_exp !: user-defined frequency bin bounds are exponentially spaced
   REAL    , PUBLIC                 ::   rn_wfreq_0       !: lowest frequency (when calculated internally, ln_wfreq_usr=F)
   REAL    , PUBLIC                 ::   rn_wfreq_k       !: ratio between adjacent frequencies (when calculated internally, ln_wfreq_usr=F)
   !
   !!----------------------------------------------------------------------
   !!           switch definition (improve readability)
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_usr     = 1        !: user defined                  formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_flx     = 2        !: flux                          formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_blk     = 3        !: bulk                          formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_abl     = 4        !: Atmospheric boundary layer    formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_purecpl = 5        !: Pure ocean-atmosphere Coupled formulation
   INTEGER , PUBLIC, PARAMETER ::   jp_none    = 6        !: for OCE when doing coupling via SAS module
   !
   !!----------------------------------------------------------------------
   !!           component definition
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_nemo = 0      !: Initial single executable configuration
   !  (no internal OASIS coupling)
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_oce  = 1      !: Multi executable configuration - OCE component
   !  (internal OASIS coupling)
   INTEGER , PUBLIC, PARAMETER ::   jp_iam_sas  = 2      !: Multi executable configuration - SAS component
   !  (internal OASIS coupling)
   !!----------------------------------------------------------------------
   !!              Ocean Surface Boundary Condition fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC ::  ncpl_qsr_freq = 0        !: qsr coupling frequency per days from atmosphere (used by top)
   !
   !!                                   !!   now    ! before   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau               !: sea surface i-stress (ocean referential) T-pt [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vtau               !: sea surface j-stress (ocean referential) T-pt [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utauU  , utau_b    !: sea surface i-stress (ocean referential) U-pt [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   vtauV  , vtau_b    !: sea surface j-stress (ocean referential) V-pt [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau_icb, vtau_icb !: sea surface (i,j)-stress used by icebergs     [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   taum               !: module of sea surface stress (at T-point)     [N/m2]
   !! wndm is used compute surface gases exchanges in ice-free ocean or leads
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wndm               !: wind speed module at T-point (=|U10m-Uoce|)   [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rhoa               !: air density at "rn_zu" m above the sea        [kg/m3]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr                !: sea heat flux:     solar                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns    , qns_b     !: sea heat flux: non solar                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qsr_tot            !: total     solar heat flux (over sea and ice)  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qns_tot            !: total non solar heat flux (over sea and ice)  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp    , emp_b     !: freshwater budget: volume flux                [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sfx    , sfx_b     !: salt flux                                     [PSS.kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   emp_tot            !: total E-P over ocean and ice                  [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fwfice             !: ice-ocean freshwater budget (>0 to the ocean) [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rnf    , rnf_b     !: river runoff                                  [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fwficb             !: iceberg melting                               [Kg/m2/s]
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  sbc_tsc, sbc_tsc_b  !: sbc content trend                      [K.m/s] jpi,jpj,jpts
   !!
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tprecip            !: total precipitation                           [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sprecip            !: solid precipitation                           [Kg/m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   fr_i               !: ice fraction = 1 - lead fraction       (between 0 to 1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   atm_co2            !: atmospheric pCO2                              [ppm]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: xcplmask         !: coupling mask for ln_mixcpl (warning: allocated in sbccpl)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   cloud_fra          !: cloud cover (fraction of cloud in a gridcell) [-]

   !!---------------------------------------------------------------------
   !! ABL Vertical Domain size
   !!---------------------------------------------------------------------
   INTEGER , PUBLIC            ::   jpka   = 2     !: ABL number of vertical levels (default definition)
   INTEGER , PUBLIC            ::   jpkam1 = 1     !: jpka-1
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   ght_abl, ghw_abl          !: ABL geopotential height (needed for iom)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   ::   e3t_abl, e3w_abl          !: ABL vertical scale factors (needed for iom)

   !!----------------------------------------------------------------------
   !!                     Sea Surface Mean fields
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC                     ::   nn_fsbc   !: frequency of sbc computation (as well as sea-ice model)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssu_m     !: mean (nn_fsbc time-step) surface sea i-current (U-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssv_m     !: mean (nn_fsbc time-step) surface sea j-current (V-point) [m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sst_m     !: mean (nn_fsbc time-step) surface sea temperature     [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   sss_m     !: mean (nn_fsbc time-step) surface sea salinity            [psu]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssh_m     !: mean (nn_fsbc time-step) sea surface height                [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tsk_m     !: mean (nn_fsbc time-step) SKIN surface sea temp.      [Celsius]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   e3t_m     !: mean (nn_fsbc time-step) sea surface layer thickness       [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   frq_m     !: mean (nn_fsbc time-step) fraction of solar net radiation absorbed in the 1st T level [-]

   !!----------------------------------------------------------------------
   !!                     Surface atmospheric fields
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: q_air_zt       !: specific humidity of air at z=zt [kg/kg]ww
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: theta_air_zt   !: potential temperature of air at z=zt [K]


   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_oce_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION sbc_oce_alloc  ***
      !!---------------------------------------------------------------------
      INTEGER :: ierr(8)
      !!---------------------------------------------------------------------
      ierr(:) = 0
      !
      ! ----------------- !
      ! == FULL ARRAYS == !
      ! ----------------- !
      ALLOCATE( utau(jpi,jpj) , utau_b(jpi,jpj) , utauU(jpi,jpj) , &
         &      vtau(jpi,jpj) , vtau_b(jpi,jpj) , vtauV(jpi,jpj) , STAT=ierr(1) )
      !
      ALLOCATE( emp(jpi,jpj) , emp_b(jpi,jpj) ,  &
         &      STAT=ierr(2) )
      !
      IF(ln_rnf)   ALLOCATE( rnf(jpi,jpj), rnf_b(jpi,jpj), STAT=ierr(3) )
      !
      ALLOCATE( fr_i(jpi,jpj) ,     &
         &      ssu_m  (jpi,jpj) , sst_m(jpi,jpj) , frq_m(jpi,jpj) ,      &
         &      ssv_m  (jpi,jpj) , sss_m(jpi,jpj) , ssh_m(jpi,jpj) , e3t_m(jpi,jpj) , STAT=ierr(4) )
      !
      ! -------------------- !
      ! == REDUCED ARRAYS == !
      ! -------------------- !
      ALLOCATE( qns    (A2D(0))    , qns_b  (A2D(0))     , qsr   (A2D(0)),    &
         &      qns_tot(A2D(0))    , qsr_tot(A2D(0))     ,                    &
         &      STAT=ierr(5) )
      !
      ALLOCATE( sbc_tsc(A2D(0),jpts) , sbc_tsc_b(A2D(0),jpts) ,  &
         &      sfx (A2D(0)) , sfx_b(A2D(0)) , emp_tot(A2D(0)), fwfice(A2D(0)), fwficb(A2D(0)), &
         &      wndm(A2D(0)) , taum (A2D(0)) , STAT=ierr(6) )
      !
      ALLOCATE( tprecip(A2D(0)) , sprecip(A2D(0)) ,    &
         &      atm_co2(A2D(0)) , tsk_m  (A2D(0)) , cloud_fra(A2D(0)), STAT=ierr(7) )

      ALLOCATE( rhoa(A2D(0)) , q_air_zt(A2D(0)) , theta_air_zt(A2D(0)) , STAT=ierr(8) )
      !
      !
      sbc_oce_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbc_oce', sbc_oce_alloc )
      IF( sbc_oce_alloc > 0 )   CALL ctl_warn('sbc_oce_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_oce_alloc


   SUBROUTINE sbc_oce_dealloc()
      IF( ALLOCATED(utau) ) THEN
         DEALLOCATE( utau , utau_b, utauU , vtau , vtau_b, vtauV )
         IF(ln_rnf)   DEALLOCATE( emp  , emp_b , rnf   , rnf_b )
         DEALLOCATE( fr_i , ssu_m , sst_m , frq_m, ssv_m , sss_m, ssh_m, e3t_m )
         DEALLOCATE( qns, qns_b, qsr,  qns_tot , qsr_tot )
         DEALLOCATE( sbc_tsc, sbc_tsc_b   ,  sfx , sfx_b , emp_tot, fwfice, fwficb, wndm , taum )
         DEALLOCATE( tprecip , sprecip , atm_co2 , tsk_m, cloud_fra )
         DEALLOCATE( rhoa , q_air_zt , theta_air_zt )
      ENDIF
   END SUBROUTINE sbc_oce_dealloc
   
   !!clem => this subroutine is never used in nemo
   SUBROUTINE sbc_tau2wnd
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE   ***
      !!
      !! ** Purpose : Estimation of wind speed as a function of wind stress
      !!
      !! ** Method  : |tau|=rhoa*Cd*|U|^2
      !!---------------------------------------------------------------------
      USE dom_oce         ! ocean space and time domain
      USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztau, zcoef           ! temporary variables
      INTEGER  ::   ji, jj                ! dummy indices
      !!---------------------------------------------------------------------
      zcoef = 0.5 / ( zrhoa * zcdrag )
      DO_2D( 0, 0, 0, 0 )
         ztau = SQRT( utau(ji,jj)*utau(ji,jj) + vtau(ji,jj)*vtau(ji,jj) )
         wndm(ji,jj) = SQRT ( ztau * zcoef ) * smask0(ji,jj)
      END_2D
      !
   END SUBROUTINE sbc_tau2wnd

   !!======================================================================
END MODULE sbc_oce
