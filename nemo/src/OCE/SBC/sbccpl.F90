MODULE sbccpl
   !!======================================================================
   !!                       ***  MODULE  sbccpl  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in coupled mode
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!            3.4  ! 2011_11  (C. Harris) more flexibility + multi-category fields
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi)  wave coupling updates
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namsbc_cpl      : coupled formulation namlist
   !!   sbc_cpl_init    : initialisation of the coupled exchanges
   !!   sbc_cpl_rcv     : receive fields from the atmosphere over the ocean (ocean only)
   !!                     receive stress from the atmosphere over the ocean (ocean-ice case)
   !!   sbc_cpl_ice_tau : receive stress from the atmosphere over ice
   !!   sbc_cpl_ice_flx : receive fluxes from the atmosphere over ice
   !!   sbc_cpl_snd     : send     fields to the atmosphere
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE trc_oce         ! share SMS/Ocean variables
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE sbcapr          ! Stochastic param. : ???
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE sbcwave         ! surface boundary condition: waves
   USE phycst          ! physical constants
   USE isf_oce , ONLY : l_isfoasis, fwfisf_oasis ! ice shelf boundary condition
#if defined key_si3
   USE par_ice        ! SI3 parameters
   USE ice     , ONLY : a_i, h_i, h_s, h_ip, a_ip_frac, a_ip_eff, cnd_ice, t1_ice, u_ice, v_ice
#endif
   USE cpl_oasis3     ! OASIS3 coupling
   USE geo2ocean      !
   USE oce     , ONLY : ts, uu, vv, ssh, fraqsr_1lev
   USE ocealb         !
   USE eosbn2         !
   USE sbcrnf  , ONLY : l_rnfcpl
#if defined key_cice
   USE ice_domain_size, only: ncat
#endif
#if defined key_si3
   USE icevar         ! for CALL ice_var_snwblow
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! NetCDF library
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

#if defined key_oasis3
   USE mod_oasis, ONLY : OASIS_Sent, OASIS_ToRest, OASIS_SentOut, OASIS_ToRestOut
#endif

   USE sbc_phy, ONLY : pp_cldf, rpref

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_cpl_init      ! routine called by sbcmod.F90
   PUBLIC   sbc_cpl_rcv       ! routine called by icestp.F90
   PUBLIC   sbc_cpl_snd       ! routine called by stprk3.F90 and step.F90 (SAS)
   PUBLIC   sbc_cpl_ice_tau   ! routine called by icestp.F90
   PUBLIC   sbc_cpl_ice_flx   ! routine called by icestp.F90
   PUBLIC   sbc_cpl_alloc     ! routine called in sbcice_cice.F90

   !! received fields are only in the interior (without halos)
   INTEGER, PARAMETER ::   jpr_otx1   =  1   ! 3 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jpr_oty1   =  2   !
   INTEGER, PARAMETER ::   jpr_otz1   =  3   !
!!$   INTEGER, PARAMETER ::   jpr_otx2   =  4   ! 3 atmosphere-ocean stress components on grid 2
!!$   INTEGER, PARAMETER ::   jpr_oty2   =  5   !
!!$   INTEGER, PARAMETER ::   jpr_otz2   =  6   !
   INTEGER, PARAMETER ::   jpr_itx1   =  7   ! 3 atmosphere-ice   stress components on grid 1
   INTEGER, PARAMETER ::   jpr_ity1   =  8   !
   INTEGER, PARAMETER ::   jpr_itz1   =  9   !
!!$   INTEGER, PARAMETER ::   jpr_itx2   = 10   ! 3 atmosphere-ice   stress components on grid 2
!!$   INTEGER, PARAMETER ::   jpr_ity2   = 11   !
!!$   INTEGER, PARAMETER ::   jpr_itz2   = 12   !
   INTEGER, PARAMETER ::   jpr_qsroce = 13   ! Qsr above the ocean
   INTEGER, PARAMETER ::   jpr_qsrice = 14   ! Qsr above the ice
   INTEGER, PARAMETER ::   jpr_qsrmix = 15
   INTEGER, PARAMETER ::   jpr_qnsoce = 16   ! Qns above the ocean
   INTEGER, PARAMETER ::   jpr_qnsice = 17   ! Qns above the ice
   INTEGER, PARAMETER ::   jpr_qnsmix = 18
   INTEGER, PARAMETER ::   jpr_rain   = 19   ! total liquid precipitation (rain)
   INTEGER, PARAMETER ::   jpr_snow   = 20   ! solid precipitation over the ocean (snow)
   INTEGER, PARAMETER ::   jpr_tevp   = 21   ! total evaporation
   INTEGER, PARAMETER ::   jpr_ievp   = 22   ! solid evaporation (sublimation)
   INTEGER, PARAMETER ::   jpr_sbpr   = 23   ! sublimation - liquid precipitation - solid precipitation
   INTEGER, PARAMETER ::   jpr_semp   = 24   ! solid freshwater budget (sublimation - snow)
   INTEGER, PARAMETER ::   jpr_oemp   = 25   ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jpr_w10m   = 26   ! 10m wind
   INTEGER, PARAMETER ::   jpr_dqnsdt = 27   ! d(Q non solar)/d(temperature)
   INTEGER, PARAMETER ::   jpr_rnf    = 28   ! runoffs
   INTEGER, PARAMETER ::   jpr_cal    = 29   ! calving
   INTEGER, PARAMETER ::   jpr_taum   = 30   ! wind stress module
   INTEGER, PARAMETER ::   jpr_co2    = 31
   INTEGER, PARAMETER ::   jpr_topm   = 32   ! topmeltn
   INTEGER, PARAMETER ::   jpr_botm   = 33   ! botmeltn
   INTEGER, PARAMETER ::   jpr_sflx   = 34   ! salt flux
   INTEGER, PARAMETER ::   jpr_toce   = 35   ! ocean temperature
   INTEGER, PARAMETER ::   jpr_soce   = 36   ! ocean salinity
   INTEGER, PARAMETER ::   jpr_ocx1   = 37   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jpr_ocy1   = 38   !
   INTEGER, PARAMETER ::   jpr_ssh    = 39   ! sea surface height
   INTEGER, PARAMETER ::   jpr_fice   = 40   ! ice fraction
   INTEGER, PARAMETER ::   jpr_e3t1st = 41   ! first T level thickness
   INTEGER, PARAMETER ::   jpr_fraqsr = 42   ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jpr_mslp   = 43   ! mean sea level pressure
   !**  surface wave coupling  **
   INTEGER, PARAMETER ::   jpr_hsig   = 44   ! Hsig
   INTEGER, PARAMETER ::   jpr_phioc  = 45   ! Wave=>ocean energy flux
   INTEGER, PARAMETER ::   jpr_sdrftx = 46   ! Stokes drift on grid 1
   INTEGER, PARAMETER ::   jpr_sdrfty = 47   ! Stokes drift on grid 2
   INTEGER, PARAMETER ::   jpr_wper   = 48   ! Mean wave period
   INTEGER, PARAMETER ::   jpr_wnum   = 49   ! Mean wavenumber
   INTEGER, PARAMETER ::   jpr_wstrf  = 50   ! Stress fraction adsorbed by waves
   INTEGER, PARAMETER ::   jpr_wdrag  = 51   ! Neutral surface drag coefficient
   INTEGER, PARAMETER ::   jpr_charn  = 52   ! Chranock coefficient
   INTEGER, PARAMETER ::   jpr_twox   = 53   ! wave to ocean momentum flux
   INTEGER, PARAMETER ::   jpr_twoy   = 54   ! wave to ocean momentum flux
   INTEGER, PARAMETER ::   jpr_tawx   = 55   ! net wave-supported stress
   INTEGER, PARAMETER ::   jpr_tawy   = 56   ! net wave-supported stress
   INTEGER, PARAMETER ::   jpr_bhd    = 57   ! Bernoulli head. waves' induced surface pressure
   INTEGER, PARAMETER ::   jpr_tusd   = 58   ! zonal stokes transport
   INTEGER, PARAMETER ::   jpr_tvsd   = 59   ! meridional stokes tranmport
   INTEGER, PARAMETER ::   jpr_isf    = 60
   INTEGER, PARAMETER ::   jpr_icb    = 61
   INTEGER, PARAMETER ::   jpr_ts_ice = 62   ! Sea ice surface temp
   INTEGER, PARAMETER ::   jpr_qtrice = 63   ! Transmitted solar thru sea-ice
   INTEGER, PARAMETER ::   jpr_wpf    = 64   ! Wave peak frequency

   INTEGER, PARAMETER ::   jprcv      = 64   ! total number of fields received
   
   !! sent fields are only in the interior (without halos)
   INTEGER, PARAMETER ::   jps_fice   =  1   ! ice fraction sent to the atmosphere
   INTEGER, PARAMETER ::   jps_toce   =  2   ! ocean temperature
   INTEGER, PARAMETER ::   jps_tice   =  3   ! ice   temperature
   INTEGER, PARAMETER ::   jps_tmix   =  4   ! mixed temperature (ocean+ice)
   INTEGER, PARAMETER ::   jps_albice =  5   ! ice   albedo
   INTEGER, PARAMETER ::   jps_albmix =  6   ! mixed albedo
   INTEGER, PARAMETER ::   jps_hice   =  7   ! ice  thickness
   INTEGER, PARAMETER ::   jps_hsnw   =  8   ! snow thickness
   INTEGER, PARAMETER ::   jps_ocx1   =  9   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jps_ocy1   = 10   !
   INTEGER, PARAMETER ::   jps_ocz1   = 11   !
   INTEGER, PARAMETER ::   jps_ivx1   = 12   ! ice   current on grid 1
   INTEGER, PARAMETER ::   jps_ivy1   = 13   !
   INTEGER, PARAMETER ::   jps_ivz1   = 14   !
   INTEGER, PARAMETER ::   jps_co2    = 15
   INTEGER, PARAMETER ::   jps_soce   = 16   ! ocean salinity
   INTEGER, PARAMETER ::   jps_ssh    = 17   ! sea surface height
   INTEGER, PARAMETER ::   jps_qsroce = 18   ! Qsr above the ocean
   INTEGER, PARAMETER ::   jps_qnsoce = 19   ! Qns above the ocean
   INTEGER, PARAMETER ::   jps_oemp   = 20   ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jps_sflx   = 21   ! salt flux
   INTEGER, PARAMETER ::   jps_otx1   = 22   ! 2 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jps_oty1   = 23   !
   INTEGER, PARAMETER ::   jps_rnf    = 24   ! runoffs
   INTEGER, PARAMETER ::   jps_taum   = 25   ! wind stress module
   INTEGER, PARAMETER ::   jps_fice2  = 26   ! ice fraction sent to OCE (by SAS when doing SAS-OCE coupling)
   INTEGER, PARAMETER ::   jps_e3t1st = 27   ! first level depth (vvl)
   INTEGER, PARAMETER ::   jps_fraqsr = 28   ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jps_ficet  = 29   ! total ice fraction
   INTEGER, PARAMETER ::   jps_ocxw   = 30   ! currents on grid 1
   INTEGER, PARAMETER ::   jps_ocyw   = 31   ! currents on grid 2
   INTEGER, PARAMETER ::   jps_wlev   = 32   ! water level
   INTEGER, PARAMETER ::   jps_fice1  = 33   ! first-order ice concentration (for semi-implicit coupling of atmos-ice fluxes)
   INTEGER, PARAMETER ::   jps_a_p    = 34   ! meltpond area fraction
   INTEGER, PARAMETER ::   jps_ht_p   = 35   ! meltpond thickness
   INTEGER, PARAMETER ::   jps_kice   = 36   ! sea ice effective conductivity
   INTEGER, PARAMETER ::   jps_sstfrz = 37   ! sea surface freezing temperature
   INTEGER, PARAMETER ::   jps_ttilyr = 38   ! sea ice top layer temp

   INTEGER, PARAMETER ::   jpsnd      = 38   ! total number of fields sent

#if ! defined key_oasis3
   ! Dummy variables to enable compilation when oasis3 is not being used
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1
#endif

   !                                  !!** namelist namsbc_cpl **
   TYPE ::   FLD_C                     !
      CHARACTER(len = 32) ::   cldes      ! desciption of the coupling strategy
      CHARACTER(len = 32) ::   clcat      ! multiple ice categories strategy
      CHARACTER(len = 32) ::   clvref     ! reference of vector ('spherical' or 'cartesian')
      CHARACTER(len = 32) ::   clvor      ! orientation of vector fields ('eastward-northward' or 'local grid')
      CHARACTER(len = 32) ::   clvgrd     ! grids on which is located the vector fields
   END TYPE FLD_C
   !                                   ! Send to the atmosphere
   TYPE(FLD_C) ::   sn_snd_temp  , sn_snd_alb , sn_snd_thick, sn_snd_crt   , sn_snd_co2,  &
      &             sn_snd_thick1, sn_snd_cond, sn_snd_mpnd , sn_snd_sstfrz, sn_snd_ttilyr
   !                                   ! Received from the atmosphere
   TYPE(FLD_C) ::   sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_dqnsdt, sn_rcv_qsr,  &
      &             sn_rcv_qns , sn_rcv_emp   , sn_rcv_rnf, sn_rcv_ts_ice, sn_rcv_qtrice
   TYPE(FLD_C) ::   sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2, sn_rcv_mslp, sn_rcv_icb, sn_rcv_isf
   !                                   ! Send to waves
   TYPE(FLD_C) ::   sn_snd_ifrac, sn_snd_crtw, sn_snd_wlev
   !                                   ! Received from waves
   TYPE(FLD_C) ::   sn_rcv_hsig, sn_rcv_phioc, sn_rcv_sdrfx, sn_rcv_sdrfy, sn_rcv_wper, sn_rcv_wnum, sn_rcv_wpf, &
      &             sn_rcv_wstrf, sn_rcv_wdrag, sn_rcv_charn, sn_rcv_taw, sn_rcv_bhd, sn_rcv_tusd, sn_rcv_tvsd
   !                                   ! Other namelist parameters
   INTEGER     ::   nn_cplmodel           ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
   LOGICAL     ::   ln_usecplmask         !  use a coupling mask file to merge data received from several models
                                          !   -> file cplmask.nc with the float variable called cplmask (jpi,jpj,nn_cplmodel)
   LOGICAL     ::   ln_scale_ice_flux     !  use ice fluxes that are already "ice weighted" ( i.e. multiplied ice concentration)

   TYPE( DYNARR ), SAVE, DIMENSION(jprcv) ::   frcv                ! all fields recieved from the atmosphere

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   alb_oce_mix    ! ocean albedo sent to atmosphere (mix clear/overcast sky)
#if defined key_si3 || defined key_cice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i_last_couple !: Ice fractional area at last coupling time
#endif

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   nrcvinfo           ! OASIS info argument

   !! Substitution
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION sbc_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( alb_oce_mix(A2D(0)), nrcvinfo(jprcv),  STAT=ierr(1) )

#if ! defined key_si3 && ! defined key_cice
      ALLOCATE( a_i(jpi,jpj,1) , STAT=ierr(2) )  ! used in sbcice_if.F90 (done here as there is no sbc_ice_if_init)
#endif
      ALLOCATE( xcplmask(A2D(0),1,0:nn_cplmodel) , STAT=ierr(3) )
      !
      IF( .NOT. ln_apr_dyn ) ALLOCATE( ssh_ib(jpi,jpj), ssh_ibb(jpi,jpj), apr(jpi, jpj), STAT=ierr(4) )

      sbc_cpl_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbccpl', sbc_cpl_alloc )
      IF( sbc_cpl_alloc > 0 )   CALL ctl_warn('sbc_cpl_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_cpl_alloc


   SUBROUTINE sbc_cpl_init( k_ice )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_init  ***
      !!
      !! ** Purpose :   Initialisation of send and received information from
      !!                the atmospheric component
      !!
      !! ** Method  : * Read namsbc_cpl namelist
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_ice   ! ice management in the sbc (=0/1/2/3)
      !
      INTEGER ::   jn          ! dummy loop index
      INTEGER ::   ios, inum   ! Local integer
      REAL(wp), DIMENSION(A2D(0)) ::   zacs, zaos
      !!
      NAMELIST/namsbc_cpl/  nn_cplmodel  , ln_usecplmask, nn_cats_cpl , ln_scale_ice_flux,             &
         &                  sn_snd_temp  , sn_snd_alb   , sn_snd_thick, sn_snd_crt   , sn_snd_co2   ,  &
         &                  sn_snd_ttilyr, sn_snd_cond  , sn_snd_mpnd , sn_snd_sstfrz, sn_snd_thick1,  &
         &                  sn_snd_ifrac , sn_snd_crtw  , sn_snd_wlev , sn_rcv_hsig  , sn_rcv_phioc ,  &
         &                  sn_rcv_w10m  , sn_rcv_taumod, sn_rcv_tau  , sn_rcv_dqnsdt, sn_rcv_qsr   ,  &
         &                  sn_rcv_sdrfx , sn_rcv_sdrfy , sn_rcv_wper , sn_rcv_wnum  , sn_rcv_wstrf ,  &
         &                  sn_rcv_charn , sn_rcv_taw   , sn_rcv_bhd  , sn_rcv_tusd  , sn_rcv_tvsd,    &
         &                  sn_rcv_wdrag , sn_rcv_qns   , sn_rcv_emp  , sn_rcv_rnf   , sn_rcv_cal  ,   &
         &                  sn_rcv_iceflx, sn_rcv_co2   , sn_rcv_icb  , sn_rcv_isf   , sn_rcv_ts_ice, sn_rcv_qtrice, &
         &                  sn_rcv_mslp  , sn_rcv_wpf

      !!---------------------------------------------------------------------
      !
      ! ================================ !
      !      Namelist informations       !
      ! ================================ !
      !
      READ_NML_REF(numnam,namsbc_cpl)
      READ_NML_CFG(numnam,namsbc_cpl)
      IF(lwm) WRITE ( numond, namsbc_cpl )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'sbc_cpl_init : namsbc_cpl namelist '
         WRITE(numout,*)'~~~~~~~~~~~~'
      ENDIF
      IF( lwp .AND. ln_cpl ) THEN                        ! control print
         WRITE(numout,*)'  nn_cplmodel                         = ', nn_cplmodel
         WRITE(numout,*)'  ln_usecplmask                       = ', ln_usecplmask
         WRITE(numout,*)'  ln_scale_ice_flux                   = ', ln_scale_ice_flux
         WRITE(numout,*)'  nn_cats_cpl                         = ', nn_cats_cpl
         WRITE(numout,*)'  received fields (mutiple ice categogies)'
         WRITE(numout,*)'      10m wind module                 = ', TRIM(sn_rcv_w10m%cldes  ), ' (', TRIM(sn_rcv_w10m%clcat  ), ')'
         WRITE(numout,*)'      stress module                   = ', TRIM(sn_rcv_taumod%cldes), ' (', TRIM(sn_rcv_taumod%clcat), ')'
         WRITE(numout,*)'      surface stress                  = ', TRIM(sn_rcv_tau%cldes   ), ' (', TRIM(sn_rcv_tau%clcat   ), ')'
         WRITE(numout,*)'                     - referential    = ', sn_rcv_tau%clvref
         WRITE(numout,*)'                     - orientation    = ', sn_rcv_tau%clvor
         WRITE(numout,*)'      non-solar heat flux sensitivity = ', TRIM(sn_rcv_dqnsdt%cldes), ' (', TRIM(sn_rcv_dqnsdt%clcat), ')'
         WRITE(numout,*)'      solar heat flux                 = ', TRIM(sn_rcv_qsr%cldes   ), ' (', TRIM(sn_rcv_qsr%clcat   ), ')'
         WRITE(numout,*)'      non-solar heat flux             = ', TRIM(sn_rcv_qns%cldes   ), ' (', TRIM(sn_rcv_qns%clcat   ), ')'
         WRITE(numout,*)'      freshwater budget               = ', TRIM(sn_rcv_emp%cldes   ), ' (', TRIM(sn_rcv_emp%clcat   ), ')'
         WRITE(numout,*)'      runoffs                         = ', TRIM(sn_rcv_rnf%cldes   ), ' (', TRIM(sn_rcv_rnf%clcat   ), ')'
         WRITE(numout,*)'      calving                         = ', TRIM(sn_rcv_cal%cldes   ), ' (', TRIM(sn_rcv_cal%clcat   ), ')'
         WRITE(numout,*)'      iceberg                         = ', TRIM(sn_rcv_icb%cldes   ), ' (', TRIM(sn_rcv_icb%clcat   ), ')'
         WRITE(numout,*)'      ice shelf                       = ', TRIM(sn_rcv_isf%cldes   ), ' (', TRIM(sn_rcv_isf%clcat   ), ')'
         WRITE(numout,*)'      sea ice heat fluxes             = ', TRIM(sn_rcv_iceflx%cldes), ' (', TRIM(sn_rcv_iceflx%clcat), ')'
         WRITE(numout,*)'      transmitted solar thru sea-ice  = ', TRIM(sn_rcv_qtrice%cldes), ' (', TRIM(sn_rcv_qtrice%clcat), ')'
         WRITE(numout,*)'      atm co2                         = ', TRIM(sn_rcv_co2%cldes   ), ' (', TRIM(sn_rcv_co2%clcat   ), ')'
         WRITE(numout,*)'      Sea ice surface skin temperature= ', TRIM(sn_rcv_ts_ice%cldes), ' (', TRIM(sn_rcv_ts_ice%clcat), ')'
         WRITE(numout,*)'      surface waves:'
         WRITE(numout,*)'      significant wave heigth         = ', TRIM(sn_rcv_hsig%cldes  ), ' (', TRIM(sn_rcv_hsig%clcat  ), ')'
         WRITE(numout,*)'      wave to oce energy flux         = ', TRIM(sn_rcv_phioc%cldes ), ' (', TRIM(sn_rcv_phioc%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid u     = ', TRIM(sn_rcv_sdrfx%cldes ), ' (', TRIM(sn_rcv_sdrfx%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid v     = ', TRIM(sn_rcv_sdrfy%cldes ), ' (', TRIM(sn_rcv_sdrfy%clcat ), ')'
         WRITE(numout,*)'      Mean wave period                = ', TRIM(sn_rcv_wper%cldes  ), ' (', TRIM(sn_rcv_wper%clcat  ), ')'
         WRITE(numout,*)'      Mean wave number                = ', TRIM(sn_rcv_wnum%cldes  ), ' (', TRIM(sn_rcv_wnum%clcat  ), ')'
         WRITE(numout,*)'      Wave peak frequency             = ', TRIM(sn_rcv_wpf%cldes   ), ' (', TRIM(sn_rcv_wpf%clcat   ), ')'
         WRITE(numout,*)'      Stress frac adsorbed by waves   = ', TRIM(sn_rcv_wstrf%cldes ), ' (', TRIM(sn_rcv_wstrf%clcat ), ')'
         WRITE(numout,*)'      Neutral surf drag coefficient   = ', TRIM(sn_rcv_wdrag%cldes ), ' (', TRIM(sn_rcv_wdrag%clcat ), ')'
         WRITE(numout,*)'      Charnock coefficient            = ', TRIM(sn_rcv_charn%cldes ), ' (', TRIM(sn_rcv_charn%clcat ), ')'
         WRITE(numout,*)'  sent fields (multiple ice categories)'
         WRITE(numout,*)'      surface temperature             = ', TRIM(sn_snd_temp%cldes  ), ' (', TRIM(sn_snd_temp%clcat  ), ')'
         WRITE(numout,*)'      top ice layer temperature       = ', TRIM(sn_snd_ttilyr%cldes), ' (', TRIM(sn_snd_ttilyr%clcat), ')'
         WRITE(numout,*)'      albedo                          = ', TRIM(sn_snd_alb%cldes   ), ' (', TRIM(sn_snd_alb%clcat   ), ')'
         WRITE(numout,*)'      ice/snow thickness              = ', TRIM(sn_snd_thick%cldes ), ' (', TRIM(sn_snd_thick%clcat ), ')'
         WRITE(numout,*)'      total ice fraction              = ', TRIM(sn_snd_ifrac%cldes ), ' (', TRIM(sn_snd_ifrac%clcat ), ')'
         WRITE(numout,*)'      surface current                 = ', TRIM(sn_snd_crt%cldes   ), ' (', TRIM(sn_snd_crt%clcat   ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crt%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crt%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crt%clvgrd
         WRITE(numout,*)'      oce co2 flux                    = ', TRIM(sn_snd_co2%cldes   ), ' (', TRIM(sn_snd_co2%clcat   ), ')'
         WRITE(numout,*)'      ice effective conductivity      = ', TRIM(sn_snd_cond%cldes  ), ' (', TRIM(sn_snd_cond%clcat  ), ')'
         WRITE(numout,*)'      meltponds fraction and depth    = ', TRIM(sn_snd_mpnd%cldes  ), ' (', TRIM(sn_snd_mpnd%clcat  ), ')'
         WRITE(numout,*)'      sea surface freezing temp       = ', TRIM(sn_snd_sstfrz%cldes), ' (', TRIM(sn_snd_sstfrz%clcat), ')'
         WRITE(numout,*)'      water level                     = ', TRIM(sn_snd_wlev%cldes  ), ' (', TRIM(sn_snd_wlev%clcat  ), ')'
         WRITE(numout,*)'      mean sea level pressure         = ', TRIM(sn_rcv_mslp%cldes  ), ' (', TRIM(sn_rcv_mslp%clcat  ), ')'
         WRITE(numout,*)'      surface current to waves        = ', TRIM(sn_snd_crtw%cldes  ), ' (', TRIM(sn_snd_crtw%clcat  ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crtw%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crtw%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crtw%clvgrd
      ENDIF
      IF( lwp .AND. ln_wave) THEN                        ! control print
         WRITE(numout,*)'      surface waves:'
         WRITE(numout,*)'      Significant wave heigth         = ', TRIM(sn_rcv_hsig%cldes  ), ' (', TRIM(sn_rcv_hsig%clcat  ), ')'
         WRITE(numout,*)'      Wave to oce energy flux         = ', TRIM(sn_rcv_phioc%cldes ), ' (', TRIM(sn_rcv_phioc%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid u     = ', TRIM(sn_rcv_sdrfx%cldes ), ' (', TRIM(sn_rcv_sdrfx%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid v     = ', TRIM(sn_rcv_sdrfy%cldes ), ' (', TRIM(sn_rcv_sdrfy%clcat ), ')'
         WRITE(numout,*)'      Mean wave period                = ', TRIM(sn_rcv_wper%cldes  ), ' (', TRIM(sn_rcv_wper%clcat  ), ')'
         WRITE(numout,*)'      Mean wave number                = ', TRIM(sn_rcv_wnum%cldes  ), ' (', TRIM(sn_rcv_wnum%clcat  ), ')'
         WRITE(numout,*)'      Wave peak frequency             = ', TRIM(sn_rcv_wpf%cldes   ), ' (', TRIM(sn_rcv_wpf%clcat   ), ')'
         WRITE(numout,*)'      Stress frac adsorbed by waves   = ', TRIM(sn_rcv_wstrf%cldes ), ' (', TRIM(sn_rcv_wstrf%clcat ), ')'
         WRITE(numout,*)'      Neutral surf drag coefficient   = ', TRIM(sn_rcv_wdrag%cldes ), ' (', TRIM(sn_rcv_wdrag%clcat ), ')'
         WRITE(numout,*)'      Charnock coefficient            = ', TRIM(sn_rcv_charn%cldes ), ' (', TRIM(sn_rcv_charn%clcat ), ')'
         WRITE(numout,*)' Transport associated to Stokes drift u = ', TRIM(sn_rcv_tusd%cldes ), ' (', TRIM(sn_rcv_tusd%clcat ), ')'
         WRITE(numout,*)' Transport associated to Stokes drift v = ', TRIM(sn_rcv_tvsd%cldes ), ' (', TRIM(sn_rcv_tvsd%clcat ), ')'
         WRITE(numout,*)'      Bernouilli pressure head        = ', TRIM(sn_rcv_bhd%cldes   ), ' (', TRIM(sn_rcv_bhd%clcat  ), ')'
         WRITE(numout,*)'Wave to ocean momentum flux and Net wave-supported stress = ', TRIM(sn_rcv_taw%cldes ), ' (', TRIM(sn_rcv_taw%clcat ), ')'
         WRITE(numout,*)'      Surface current to waves        = ', TRIM(sn_snd_crtw%cldes  ), ' (', TRIM(sn_snd_crtw%clcat  ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crtw%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crtw%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crtw%clvgrd
      ENDIF
      !                                   ! allocate sbccpl arrays
      IF( sbc_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_cpl_alloc : unable to allocate arrays' )

      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress

      ! for each field: define the OASIS name                              (srcv(nmodsbc)%fld(:)%clname)
      !                 define receive or not from the namelist parameters (srcv(nmodsbc)%fld(:)%laction)
      !                 define the north fold type of lbc                  (srcv(nmodsbc)%fld(:)%nsgn)

      ! default definitions of srcv
      ALLOCATE( srcv(nmodsbc)%fld(jprcv) )
      srcv(nmodsbc)%fld(:)%laction = .FALSE.   ;   srcv(nmodsbc)%fld(:)%clgrid = 'T'   ;   srcv(nmodsbc)%fld(:)%nsgn = 1.
      srcv(nmodsbc)%fld(:)%nct = 1   ;   srcv(nmodsbc)%fld(:)%nlvl = 1

      !                                                      ! ------------------------- !
      !                                                      ! ice and ocean wind stress !
      !                                                      ! ------------------------- !
      !                                                           ! Name
      srcv(nmodsbc)%fld(jpr_otx1)%clname = 'O_OTaux1'      ! 1st ocean component on grid ONE (T or U)
      srcv(nmodsbc)%fld(jpr_oty1)%clname = 'O_OTauy1'      ! 2nd   -      -         -     -
      srcv(nmodsbc)%fld(jpr_otz1)%clname = 'O_OTauz1'      ! 3rd   -      -         -     -
      !
      srcv(nmodsbc)%fld(jpr_itx1)%clname = 'O_ITaux1'      ! 1st  ice  component on grid ONE (T, F, I or U)
      srcv(nmodsbc)%fld(jpr_ity1)%clname = 'O_ITauy1'      ! 2nd   -      -         -     -
      srcv(nmodsbc)%fld(jpr_itz1)%clname = 'O_ITauz1'      ! 3rd   -      -         -     -
      !
      ! Vectors: change of sign at north fold ONLY if on the local grid
      IF(       TRIM( sn_rcv_tau%cldes ) == 'oce only' .OR. TRIM( sn_rcv_tau%cldes ) == 'oce and ice'  &
           .OR. TRIM( sn_rcv_tau%cldes ) == 'mixed oce-ice' ) THEN ! avoid working with the atmospheric fields if they are not coupled
      !
      IF( TRIM( sn_rcv_tau%clvor ) == 'local grid' )   srcv(nmodsbc)%fld(jpr_otx1:jpr_itz1)%nsgn = -1.

      !                                                            ! Set grid and action
      srcv(nmodsbc)%fld(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1
      srcv(nmodsbc)%fld(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1
      !
      IF( TRIM( sn_rcv_tau%clvref ) == 'spherical' )   &           ! spherical: 3rd component not received
         &     srcv(nmodsbc)%fld( (/jpr_otz1, jpr_itz1/) )%laction = .FALSE.
      !
      IF( TRIM( sn_rcv_tau%cldes ) /= 'oce and ice' ) THEN         ! 'oce and ice' case ocean stress on ocean mesh used
         srcv(nmodsbc)%fld(jpr_itx1:jpr_itz1)%laction = .FALSE.    ! ice components not received
      ENDIF
      ENDIF

      !                                                      ! ------------------------- !
      !                                                      !    freshwater budget      !   E-P
      !                                                      ! ------------------------- !
      ! we suppose that atmosphere modele do not make the difference between precipiration (liquide or solid)
      ! over ice of free ocean within the same atmospheric cell.cd
      srcv(nmodsbc)%fld(jpr_rain)%clname = 'OTotRain'      ! Rain = liquid precipitation
      srcv(nmodsbc)%fld(jpr_snow)%clname = 'OTotSnow'      ! Snow = solid precipitation
      srcv(nmodsbc)%fld(jpr_tevp)%clname = 'OTotEvap'      ! total evaporation (over oce + ice sublimation)
      srcv(nmodsbc)%fld(jpr_ievp)%clname = 'OIceEvap'      ! evaporation over ice = sublimation
      srcv(nmodsbc)%fld(jpr_sbpr)%clname = 'OSubMPre'      ! sublimation - liquid precipitation - solid precipitation
      srcv(nmodsbc)%fld(jpr_semp)%clname = 'OISubMSn'      ! ice solid water budget = sublimation - solid precipitation
      srcv(nmodsbc)%fld(jpr_oemp)%clname = 'OOEvaMPr'      ! ocean water budget = ocean Evap - ocean precip
      SELECT CASE( TRIM( sn_rcv_emp%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(nmodsbc)%fld(jpr_oemp)%laction = .TRUE.
      CASE( 'conservative'  )
         srcv(nmodsbc)%fld( (/jpr_rain, jpr_snow, jpr_ievp, jpr_tevp/) )%laction = .TRUE.
         IF( k_ice <= 1 )  srcv(nmodsbc)%fld(jpr_ievp)%laction = .FALSE.
      CASE( 'oce and ice'   )   ;   srcv(nmodsbc)%fld( (/jpr_ievp, jpr_sbpr, jpr_semp, jpr_oemp/) )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_emp%cldes' )
      END SELECT
      !
      !                                                      ! ------------------------- !
      !                                                      !     Runoffs & Calving     !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_rnf   )%clname = 'O_Runoff'
      IF( TRIM( sn_rcv_rnf%cldes ) == 'coupled' ) THEN
         IF( .NOT. ln_rnf ) CALL ctl_stop( 'STOP', 'sbccpl : coupling runoffs requires ln_rnf = .true.' )
         srcv(nmodsbc)%fld(jpr_rnf)%laction = .TRUE.
         l_rnfcpl              = .TRUE.                      ! -> no need to read runoffs in sbcrnf
         ln_rnf                = nn_components /= jp_iam_sas ! -> force to go through sbcrnf if not sas
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   runoffs received from oasis -> force ln_rnf = ', ln_rnf
      ENDIF
      !
      srcv(nmodsbc)%fld(jpr_cal)%clname = 'OCalving'   ;  IF( TRIM( sn_rcv_cal%cldes) == 'coupled' )   srcv(nmodsbc)%fld(jpr_cal)%laction = .TRUE.
      srcv(nmodsbc)%fld(jpr_isf)%clname = 'OIcshelf'   ;  IF( TRIM( sn_rcv_isf%cldes) == 'coupled' )   srcv(nmodsbc)%fld(jpr_isf)%laction = .TRUE.
      srcv(nmodsbc)%fld(jpr_icb)%clname = 'OIceberg'   ;  IF( TRIM( sn_rcv_icb%cldes) == 'coupled' )   srcv(nmodsbc)%fld(jpr_icb)%laction = .TRUE.

      IF( srcv(nmodsbc)%fld(jpr_isf)%laction ) THEN
         l_isfoasis = .TRUE.  ! -> isf fwf comes from oasis
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   iceshelf received from oasis '
      ENDIF
      !
      !
      !                                                      ! ------------------------- !
      !                                                      !    non solar radiation    !   Qns
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_qnsoce)%clname = 'O_QnsOce'
      srcv(nmodsbc)%fld(jpr_qnsice)%clname = 'O_QnsIce'
      srcv(nmodsbc)%fld(jpr_qnsmix)%clname = 'O_QnsMix'
      SELECT CASE( TRIM( sn_rcv_qns%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(nmodsbc)%fld(               jpr_qnsoce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv(nmodsbc)%fld( (/jpr_qnsice, jpr_qnsmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv(nmodsbc)%fld( (/jpr_qnsice, jpr_qnsoce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(nmodsbc)%fld(               jpr_qnsmix   )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qns%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' .AND. nn_cats_cpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qns%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !
      !                                                      ! ------------------------- !
      !                                                      !    solar radiation        !   Qsr
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_qsroce)%clname = 'O_QsrOce'
      srcv(nmodsbc)%fld(jpr_qsrice)%clname = 'O_QsrIce'
      srcv(nmodsbc)%fld(jpr_qsrmix)%clname = 'O_QsrMix'
      SELECT CASE( TRIM( sn_rcv_qsr%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(nmodsbc)%fld(               jpr_qsroce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv(nmodsbc)%fld( (/jpr_qsrice, jpr_qsrmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv(nmodsbc)%fld( (/jpr_qsrice, jpr_qsroce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(nmodsbc)%fld(               jpr_qsrmix   )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qsr%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qsr%cldes ) == 'mixed oce-ice' .AND. nn_cats_cpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qsr%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !
      !                                                      ! ------------------------- !
      !                                                      !   non solar sensitivity   !   d(Qns)/d(T)
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_dqnsdt)%clname = 'O_dQnsdT'
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'coupled' )   srcv(nmodsbc)%fld(jpr_dqnsdt)%laction = .TRUE.
      !
      ! non solar sensitivity mandatory for mixed oce-ice solar radiation coupling technique
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'none' .AND. TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' )  &
         &   CALL ctl_stop( 'sbc_cpl_init: namsbc_cpl namelist mismatch between sn_rcv_qns%cldes and sn_rcv_dqnsdt%cldes' )
      !
      !                                                      ! ------------------------- !
      !                                                      !      10m wind module      !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_w10m)%clname = 'O_Wind10'   ;   IF( TRIM(sn_rcv_w10m%cldes  ) == 'coupled' )   srcv(nmodsbc)%fld(jpr_w10m)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !   wind stress module      !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_taum)%clname = 'O_TauMod'   ;   IF( TRIM(sn_rcv_taumod%cldes) == 'coupled' )   srcv(nmodsbc)%fld(jpr_taum)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !      Atmospheric CO2      !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_co2 )%clname = 'O_AtmCO2'
      IF( TRIM(sn_rcv_co2%cldes   ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_co2 )%laction = .TRUE.
         l_co2cpl = .TRUE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   Atmospheric pco2 received from oasis '
         IF(lwp) WRITE(numout,*)
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      ! Mean Sea Level Pressure   !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_mslp)%clname = 'O_MSLP'     ;   IF( TRIM(sn_rcv_mslp%cldes  ) == 'coupled' )    srcv(nmodsbc)%fld(jpr_mslp)%laction = .TRUE.
      !
      !                                                      ! --------------------------------- !
      !                                                      !  ice topmelt and conduction flux  !   
      !                                                      ! --------------------------------- !
      srcv(nmodsbc)%fld(jpr_topm )%clname = 'OTopMlt'
      srcv(nmodsbc)%fld(jpr_botm )%clname = 'OBotMlt'
      IF( TRIM(sn_rcv_iceflx%cldes) == 'coupled' ) THEN
         IF( TRIM( sn_rcv_iceflx%clcat ) == 'yes' ) THEN
            srcv(nmodsbc)%fld(jpr_topm:jpr_botm)%nct = nn_cats_cpl
         ELSE
            CALL ctl_stop( 'sbc_cpl_init: sn_rcv_iceflx%clcat should always be set to yes currently' )
         ENDIF
         srcv(nmodsbc)%fld(jpr_topm:jpr_botm)%laction = .TRUE.
      ENDIF
      !                                                      ! --------------------------- !
      !                                                      ! transmitted solar thru ice  !   
      !                                                      ! --------------------------- !
      srcv(nmodsbc)%fld(jpr_qtrice)%clname = 'OQtr'
      IF( TRIM(sn_rcv_qtrice%cldes) == 'coupled' ) THEN
         IF ( TRIM( sn_rcv_qtrice%clcat ) == 'yes' ) THEN
            srcv(nmodsbc)%fld(jpr_qtrice)%nct = nn_cats_cpl
         ELSE
           CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qtrice%clcat should always be set to yes currently' )
         ENDIF
         srcv(nmodsbc)%fld(jpr_qtrice)%laction = .TRUE.
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !    ice skin temperature   !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_ts_ice)%clname = 'OTsfIce'    ! needed by Met Office
      IF( TRIM( sn_rcv_ts_ice%cldes ) == 'ice' )   srcv(nmodsbc)%fld(jpr_ts_ice)%laction = .TRUE.
      IF( TRIM( sn_rcv_ts_ice%clcat ) == 'yes' )   srcv(nmodsbc)%fld(jpr_ts_ice)%nct     = nn_cats_cpl
      IF( TRIM( sn_rcv_emp%clcat    ) == 'yes' )   srcv(nmodsbc)%fld(jpr_ievp)%nct       = nn_cats_cpl

#if defined key_si3
      IF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN
         IF( .NOT.srcv(nmodsbc)%fld(jpr_ts_ice)%laction )  &
            &   CALL ctl_stop( 'sbc_cpl_init: srcv(nmodsbc)%fld(jpr_ts_ice)%laction should be set to true when ln_cndflx=T' )
      ENDIF
#endif
      !                                                      ! ------------------------- !
      !                                                      !      Wave breaking        !
      !                                                      ! ------------------------- !
      srcv(nmodsbc)%fld(jpr_hsig)%clname  = 'O_Hsigwa'    ! significant wave height
      IF( TRIM(sn_rcv_hsig%cldes  ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_hsig)%laction = .TRUE.
         cpl_hsig = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_phioc)%clname = 'O_PhiOce'    ! wave to ocean energy
      IF( TRIM(sn_rcv_phioc%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_phioc)%laction = .TRUE.
         cpl_phioc = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_sdrftx)%clname = 'O_Sdrfx'    ! Stokes drift in the u direction
      IF( TRIM(sn_rcv_sdrfx%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_sdrftx)%laction = .TRUE.
         cpl_sdrftx = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_sdrfty)%clname = 'O_Sdrfy'    ! Stokes drift in the v direction
      IF( TRIM(sn_rcv_sdrfy%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_sdrfty)%laction = .TRUE.
         cpl_sdrfty = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_wper)%clname = 'O_WPer'       ! mean wave period
      IF( TRIM(sn_rcv_wper%cldes  ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_wper)%laction = .TRUE.
         cpl_wper = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_wnum)%clname = 'O_WNum'       ! mean wave number
      IF( TRIM(sn_rcv_wnum%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_wnum)%laction = .TRUE.
         cpl_wnum = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_wpf)%clname = 'O_Wpfreq'       ! wave peak frequency
      IF( TRIM(sn_rcv_wpf%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_wpf)%laction = .TRUE.
         cpl_wpf = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_wstrf)%clname = 'O_WStrf'     ! stress fraction adsorbed by the wave
      IF( TRIM(sn_rcv_wstrf%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_wstrf)%laction = .TRUE.
         cpl_wstrf = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_wdrag)%clname = 'O_WDrag'     ! neutral surface drag coefficient
      IF( TRIM(sn_rcv_wdrag%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_wdrag)%laction = .TRUE.
         cpl_wdrag = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_charn)%clname = 'O_Charn'     ! Chranock coefficient
      IF( TRIM(sn_rcv_charn%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_charn)%laction = .TRUE.
         cpl_charn = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_bhd)%clname = 'O_Bhd'     ! Bernoulli head. waves' induced surface pressure
      IF( TRIM(sn_rcv_bhd%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_bhd)%laction = .TRUE.
         cpl_bhd = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_tusd)%clname = 'O_Tusd'     ! zonal stokes transport
      IF( TRIM(sn_rcv_tusd%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_tusd)%laction = .TRUE.
         cpl_tusd = .TRUE.
      ENDIF
      srcv(nmodsbc)%fld(jpr_tvsd)%clname = 'O_Tvsd'     ! meridional stokes tranmport
      IF( TRIM(sn_rcv_tvsd%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_tvsd)%laction = .TRUE.
         cpl_tvsd = .TRUE.
      ENDIF

      srcv(nmodsbc)%fld(jpr_twox)%clname = 'O_Twox'     ! wave to ocean momentum flux in the u direction
      srcv(nmodsbc)%fld(jpr_twoy)%clname = 'O_Twoy'     ! wave to ocean momentum flux in the v direction
      srcv(nmodsbc)%fld(jpr_tawx)%clname = 'O_Tawx'     ! Net wave-supported stress in the u direction
      srcv(nmodsbc)%fld(jpr_tawy)%clname = 'O_Tawy'     ! Net wave-supported stress in the v direction
      IF( TRIM(sn_rcv_taw%cldes ) == 'coupled' )  THEN
         srcv(nmodsbc)%fld(jpr_twox)%laction = .TRUE.
         srcv(nmodsbc)%fld(jpr_twoy)%laction = .TRUE.
         srcv(nmodsbc)%fld(jpr_tawx)%laction = .TRUE.
         srcv(nmodsbc)%fld(jpr_tawy)%laction = .TRUE.
         cpl_taw = .TRUE.
      ENDIF
      !
      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - rcv by opa !
      !                                                      ! ------------------------------- !
      srcv(nmodsbc)%fld(jpr_sflx)%clname = 'O_SFLX'
      srcv(nmodsbc)%fld(jpr_fice)%clname = 'RIceFrc'
      !
      IF( nn_components == jp_iam_oce ) THEN    ! OCE coupled to SAS via OASIS: force received field by OCE (sent by SAS)
         srcv(nmodsbc)%fld(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         srcv(nmodsbc)%fld(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         srcv(nmodsbc)%fld(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv(nmodsbc)%fld( (/jpr_qsroce, jpr_qnsoce, jpr_oemp, jpr_sflx, jpr_fice, jpr_otx1, jpr_oty1, jpr_taum/) )%laction = .TRUE.
         srcv(nmodsbc)%fld(jpr_otx1)%clgrid = 'T'        ! oce components given at T-point
         srcv(nmodsbc)%fld(jpr_oty1)%clgrid = 'T'        
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv(nmodsbc)%fld( (/jpr_otx1,jpr_oty1/) )%nsgn = -1.
         sn_rcv_tau%clvor = 'local grid'
         sn_rcv_tau%clvref = 'spherical'
         sn_rcv_emp%cldes = 'oce only'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OCE coupling  '
            WRITE(numout,*)'               OCE component  '
            WRITE(numout,*)
            WRITE(numout,*)'  received fields from SAS component '
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V on local grid and sperical coordinates '
            WRITE(numout,*)'                  wind stress module'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! -------------------------------- !
      !                                                      !   OCE-SAS coupling - rcv by sas  !
      !                                                      ! -------------------------------- !
      srcv(nmodsbc)%fld(jpr_toce  )%clname = 'I_SSTSST'
      srcv(nmodsbc)%fld(jpr_soce  )%clname = 'I_SSSal'
      srcv(nmodsbc)%fld(jpr_ocx1  )%clname = 'I_OCurx1'
      srcv(nmodsbc)%fld(jpr_ocy1  )%clname = 'I_OCury1'
      srcv(nmodsbc)%fld(jpr_ssh   )%clname = 'I_SSHght'
      srcv(nmodsbc)%fld(jpr_e3t1st)%clname = 'I_E3T1st'
      srcv(nmodsbc)%fld(jpr_fraqsr)%clname = 'I_FraQsr'
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) srcv(nmodsbc)%fld(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(nmodsbc)%fld(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(nmodsbc)%fld(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv(nmodsbc)%fld( (/jpr_toce, jpr_soce, jpr_ssh, jpr_fraqsr, jpr_ocx1, jpr_ocy1/) )%laction = .TRUE.
         srcv(nmodsbc)%fld( jpr_e3t1st )%laction = .NOT.lk_linssh
         srcv(nmodsbc)%fld(jpr_ocx1)%clgrid = 'U'        ! oce components given at U-point
         srcv(nmodsbc)%fld(jpr_ocy1)%clgrid = 'V'        !           and           V-point
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv(nmodsbc)%fld(jpr_ocx1:jpr_ocy1)%nsgn = -1.
         ! Change first letter to couple with atmosphere if already coupled OCE
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_Runoff received by OCE from SAS and therefore S_Runoff received by SAS from the Atmosphere
         DO jn = 1, jprcv
            IF( srcv(nmodsbc)%fld(jn)%clname(1:1) == "O" ) srcv(nmodsbc)%fld(jn)%clname = "S"//srcv(nmodsbc)%fld(jn)%clname(2:LEN(srcv(nmodsbc)%fld(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OCE coupling  '
            WRITE(numout,*)'               SAS component  '
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  received fields from OCE component '
            ELSE
               WRITE(numout,*)'  Additional received fields from OCE component : '
            ENDIF
            WRITE(numout,*)'               sea surface temperature (Celsius) '
            WRITE(numout,*)'               sea surface salinity '
            WRITE(numout,*)'               surface currents '
            WRITE(numout,*)'               sea surface height '
            WRITE(numout,*)'               thickness of first ocean T level '
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF

      ! =================================================== !
      ! Allocate all parts of frcv used for received fields !
      ! =================================================== !
      DO jn = 1, jprcv
         IF( srcv(nmodsbc)%fld(jn)%laction ) ALLOCATE( frcv(jn)%z3(A2D(0),srcv(nmodsbc)%fld(jn)%nct) )
      END DO
      ! Allocate taum part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(nmodsbc)%fld(jpr_taum)%laction ) ALLOCATE( frcv(jpr_taum)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_taum)%nct) )
      ! Allocate w10m part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(nmodsbc)%fld(jpr_w10m)%laction ) ALLOCATE( frcv(jpr_w10m)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_w10m)%nct) )
      ! Allocate jpr_otx1 part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(nmodsbc)%fld(jpr_otx1)%laction ) ALLOCATE( frcv(jpr_otx1)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_otx1)%nct) )
      IF( .NOT. srcv(nmodsbc)%fld(jpr_oty1)%laction ) ALLOCATE( frcv(jpr_oty1)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_oty1)%nct) )
      ! Allocate itx1 and ity1 as they are used in sbc_cpl_ice_tau even if srcv(nmodsbc)%fld(jpr_itx1)%laction = .FALSE.
      IF( k_ice /= 0 ) THEN
         IF( .NOT. srcv(nmodsbc)%fld(jpr_itx1)%laction ) ALLOCATE( frcv(jpr_itx1)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_itx1)%nct) )
         IF( .NOT. srcv(nmodsbc)%fld(jpr_ity1)%laction ) ALLOCATE( frcv(jpr_ity1)%z3(A2D(0),srcv(nmodsbc)%fld(jpr_ity1)%nct) )
      ENDIF

      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (ssnd(nmodsbc)%fld(:)%clname)
      !                 define send or not from the namelist parameters (ssnd(nmodsbc)%fld(:)%laction)
      !                 define the north fold type of lbc               (ssnd(nmodsbc)%fld(:)%nsgn)
      
      ! default definitions of nsnd
      ALLOCATE( ssnd(nmodsbc)%fld(jpsnd) )
      ssnd(nmodsbc)%fld(:)%laction = .FALSE.   ;   ssnd(nmodsbc)%fld(:)%clgrid = 'T'   ;   ssnd(nmodsbc)%fld(:)%nsgn = 1.
      ssnd(nmodsbc)%fld(:)%nct = 1   ;   ssnd(nmodsbc)%fld(:)%nlvl = 1

      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !
      !                                                      ! ------------------------- !
      ssnd(nmodsbc)%fld(jps_toce)%clname   = 'O_SSTSST'
      ssnd(nmodsbc)%fld(jps_tice)%clname   = 'O_TepIce'
      ssnd(nmodsbc)%fld(jps_ttilyr)%clname = 'O_TtiLyr'
      ssnd(nmodsbc)%fld(jps_tmix)%clname   = 'O_TepMix'
      SELECT CASE( TRIM( sn_snd_temp%cldes ) )
      CASE( 'none'                                 )       ! nothing to do
      CASE( 'oce only'                             )   ;   ssnd(nmodsbc)%fld( jps_toce )%laction = .TRUE.
      CASE( 'oce and ice' , 'weighted oce and ice' , 'oce and weighted ice' )
         ssnd(nmodsbc)%fld( (/jps_toce, jps_tice/) )%laction = .TRUE.
         IF( TRIM( sn_snd_temp%clcat ) == 'yes' )  ssnd(nmodsbc)%fld(jps_tice)%nct = nn_cats_cpl
      CASE( 'mixed oce-ice'                        )   ;   ssnd(nmodsbc)%fld( jps_tmix )%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_temp%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !          Albedo           !
      !                                                      ! ------------------------- !
      ssnd(nmodsbc)%fld(jps_albice)%clname = 'O_AlbIce'
      ssnd(nmodsbc)%fld(jps_albmix)%clname = 'O_AlbMix'
      SELECT CASE( TRIM( sn_snd_alb%cldes ) )
      CASE( 'none'                 )     ! nothing to do
      CASE( 'ice' , 'weighted ice' )   ; ssnd(nmodsbc)%fld(jps_albice)%laction = .TRUE.
      CASE( 'mixed oce-ice'        )   ; ssnd(nmodsbc)%fld(jps_albmix)%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_alb%cldes' )
      END SELECT
      !
      ! Need to calculate oceanic albedo if
      !     1. sending mixed oce-ice albedo or
      !     2. receiving mixed oce-ice solar radiation
      IF( TRIM ( sn_snd_alb%cldes ) == 'mixed oce-ice' .OR. TRIM ( sn_rcv_qsr%cldes ) == 'mixed oce-ice' ) THEN
         CALL oce_alb( zaos, zacs )
         ! Due to lack of information on nebulosity : mean clear/overcast sky
         alb_oce_mix(A2D(0)) = ( zacs(A2D(0)) + zaos(A2D(0)) ) * 0.5
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness !
      !                                                      ! ------------------------- !
      ssnd(nmodsbc)%fld(jps_fice)%clname  = 'OIceFrc'
      ssnd(nmodsbc)%fld(jps_ficet)%clname = 'OIceFrcT'
      ssnd(nmodsbc)%fld(jps_hice)%clname  = 'OIceTck'
      ssnd(nmodsbc)%fld(jps_a_p)%clname   = 'OPndFrc'
      ssnd(nmodsbc)%fld(jps_ht_p)%clname  = 'OPndTck'
      ssnd(nmodsbc)%fld(jps_hsnw)%clname  = 'OSnwTck'
      ssnd(nmodsbc)%fld(jps_fice1)%clname = 'OIceFrd'
      IF( k_ice /= 0 ) THEN
         ssnd(nmodsbc)%fld(jps_fice)%laction  = .TRUE.                 ! if ice treated in the ocean (even in climato case)
         ssnd(nmodsbc)%fld(jps_fice1)%laction = .TRUE.                 ! First-order regridded ice concentration, to be used producing atmos-to-ice fluxes (Met Office requirement)
! Currently no namelist entry to determine sending of multi-category ice fraction so use the thickness entry for now
         IF( TRIM( sn_snd_thick%clcat  ) == 'yes' ) ssnd(nmodsbc)%fld(jps_fice)%nct  = nn_cats_cpl
         IF( TRIM( sn_snd_thick1%clcat ) == 'yes' ) ssnd(nmodsbc)%fld(jps_fice1)%nct = nn_cats_cpl
      ENDIF

      IF(TRIM( sn_snd_ifrac%cldes )  == 'coupled') ssnd(nmodsbc)%fld(jps_ficet)%laction = .TRUE.

      SELECT CASE ( TRIM( sn_snd_thick%cldes ) )
      CASE( 'none'         )       ! nothing to do
      CASE( 'ice and snow' )
         ssnd(nmodsbc)%fld(jps_hice:jps_hsnw)%laction = .TRUE.
         IF( TRIM( sn_snd_thick%clcat ) == 'yes' ) THEN
            ssnd(nmodsbc)%fld(jps_hice:jps_hsnw)%nct = nn_cats_cpl
         ENDIF
      CASE ( 'weighted ice and snow' )
         ssnd(nmodsbc)%fld(jps_hice:jps_hsnw)%laction = .TRUE.
         IF( TRIM( sn_snd_thick%clcat ) == 'yes' ) ssnd(nmodsbc)%fld(jps_hice:jps_hsnw)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_thick%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !      Ice Meltponds        !
      !                                                      ! ------------------------- !
      ! Needed by Met Office
      ssnd(nmodsbc)%fld(jps_a_p)%clname  = 'OPndFrc'
      ssnd(nmodsbc)%fld(jps_ht_p)%clname = 'OPndTck'
      SELECT CASE ( TRIM( sn_snd_mpnd%cldes ) )
      CASE ( 'none' )
         ssnd(nmodsbc)%fld(jps_a_p)%laction  = .FALSE.
         ssnd(nmodsbc)%fld(jps_ht_p)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(nmodsbc)%fld(jps_a_p)%laction  = .TRUE.
         ssnd(nmodsbc)%fld(jps_ht_p)%laction = .TRUE.
         IF( TRIM( sn_snd_mpnd%clcat ) == 'yes' ) THEN
            ssnd(nmodsbc)%fld(jps_a_p)%nct  = nn_cats_cpl
            ssnd(nmodsbc)%fld(jps_ht_p)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_mpnd%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(nmodsbc)%fld(jps_a_p)%laction  = .TRUE.
         ssnd(nmodsbc)%fld(jps_ht_p)%laction = .TRUE.
         IF( TRIM( sn_snd_mpnd%clcat ) == 'yes' ) THEN
            ssnd(nmodsbc)%fld(jps_a_p)%nct  = nn_cats_cpl
            ssnd(nmodsbc)%fld(jps_ht_p)%nct = nn_cats_cpl
         ENDIF
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_mpnd%cldes; '//sn_snd_mpnd%cldes )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !      Surface current      !
      !                                                      ! ------------------------- !
      !        ocean currents              !            ice velocities
      ssnd(nmodsbc)%fld(jps_ocx1)%clname = 'O_OCurx1'   ;   ssnd(nmodsbc)%fld(jps_ivx1)%clname = 'O_IVelx1'
      ssnd(nmodsbc)%fld(jps_ocy1)%clname = 'O_OCury1'   ;   ssnd(nmodsbc)%fld(jps_ivy1)%clname = 'O_IVely1'
      ssnd(nmodsbc)%fld(jps_ocz1)%clname = 'O_OCurz1'   ;   ssnd(nmodsbc)%fld(jps_ivz1)%clname = 'O_IVelz1'
      ssnd(nmodsbc)%fld(jps_ocxw)%clname = 'O_OCurxw'
      ssnd(nmodsbc)%fld(jps_ocyw)%clname = 'O_OCuryw'
      !
      ssnd(nmodsbc)%fld(jps_ocx1:jps_ivz1)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( sn_snd_crt%clvgrd == 'U,V' ) THEN
         ssnd(nmodsbc)%fld(jps_ocx1)%clgrid = 'U' ; ssnd(nmodsbc)%fld(jps_ocy1)%clgrid = 'V'
      ELSE IF( sn_snd_crt%clvgrd /= 'T' ) THEN
         CALL ctl_stop( 'sn_snd_crt%clvgrd must be equal to T' )
      ENDIF
      ssnd(nmodsbc)%fld(jps_ocx1:jps_ivz1)%laction = .TRUE.   ! default: all (ocean and ice, x,y,z components) are send
      IF( TRIM( sn_snd_crt%clvref ) == 'spherical' )   ssnd(nmodsbc)%fld( (/jps_ocz1, jps_ivz1/) )%laction = .FALSE.
      IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) ssnd(nmodsbc)%fld(jps_ocx1:jps_ivz1)%nsgn = 1.
      SELECT CASE( TRIM( sn_snd_crt%cldes ) )
      CASE( 'none'                 )   ;   ssnd(nmodsbc)%fld(jps_ocx1:jps_ivz1)%laction = .FALSE.
      CASE( 'oce only'             )   ;   ssnd(nmodsbc)%fld(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE( 'oce and ice'          )   !   nothing to do
      CASE( 'weighted oce and ice' )   !   nothing to do
      CASE( 'mixed oce-ice'        )   ;   ssnd(nmodsbc)%fld(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_crt%cldes' )
      END SELECT

      ssnd(nmodsbc)%fld(jps_ocxw:jps_ocyw)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( sn_snd_crtw%clvgrd == 'U,V' ) THEN
         ssnd(nmodsbc)%fld(jps_ocxw)%clgrid = 'U' ; ssnd(nmodsbc)%fld(jps_ocyw)%clgrid = 'V'
      ELSE IF( sn_snd_crtw%clvgrd /= 'T' ) THEN
         CALL ctl_stop( 'sn_snd_crtw%clvgrd must be equal to T' )
      ENDIF
      IF( TRIM( sn_snd_crtw%clvor ) == 'eastward-northward' ) ssnd(nmodsbc)%fld(jps_ocxw:jps_ocyw)%nsgn = 1.
      SELECT CASE( TRIM( sn_snd_crtw%cldes ) )
         CASE( 'none'                 )   ; ssnd(nmodsbc)%fld(jps_ocxw:jps_ocyw)%laction = .FALSE.
         CASE( 'oce only'             )   ; ssnd(nmodsbc)%fld(jps_ocxw:jps_ocyw)%laction = .TRUE.
         CASE( 'weighted oce and ice' )   !   nothing to do
         CASE( 'mixed oce-ice'        )   ; ssnd(nmodsbc)%fld(jps_ivx1:jps_ivz1)%laction = .FALSE.
         CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_crtw%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !          CO2 flux         !
      !                                                      ! ------------------------- !
      ssnd(nmodsbc)%fld(jps_co2)%clname = 'O_CO2FLX' ;  IF( TRIM(sn_snd_co2%cldes) == 'coupled' )    ssnd(nmodsbc)%fld(jps_co2 )%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      ! Sea surface freezing temp !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      ssnd(nmodsbc)%fld(jps_sstfrz)%clname = 'O_SSTFrz' ; IF( TRIM(sn_snd_sstfrz%cldes) == 'coupled' )  ssnd(nmodsbc)%fld(jps_sstfrz)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !    Ice conductivity       !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      ! Note that ultimately we will move to passing an ocean effective conductivity as well so there
      ! will be some changes to the parts of the code which currently relate only to ice conductivity
      ssnd(nmodsbc)%fld(jps_ttilyr )%clname = 'O_TtiLyr'
      SELECT CASE ( TRIM( sn_snd_ttilyr%cldes ) )
      CASE ( 'none' )
         ssnd(nmodsbc)%fld(jps_ttilyr)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(nmodsbc)%fld(jps_ttilyr)%laction = .TRUE.
         IF( TRIM( sn_snd_ttilyr%clcat ) == 'yes' ) THEN
            ssnd(nmodsbc)%fld(jps_ttilyr)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_ttilyr%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(nmodsbc)%fld(jps_ttilyr)%laction = .TRUE.
         IF( TRIM( sn_snd_ttilyr%clcat ) == 'yes' ) ssnd(nmodsbc)%fld(jps_ttilyr)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_ttilyr%cldes;'//sn_snd_ttilyr%cldes )
      END SELECT

      ssnd(nmodsbc)%fld(jps_kice )%clname = 'OIceKn'
      SELECT CASE ( TRIM( sn_snd_cond%cldes ) )
      CASE ( 'none' )
         ssnd(nmodsbc)%fld(jps_kice)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(nmodsbc)%fld(jps_kice)%laction = .TRUE.
         IF( TRIM( sn_snd_cond%clcat ) == 'yes' ) THEN
            ssnd(nmodsbc)%fld(jps_kice)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_cond%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(nmodsbc)%fld(jps_kice)%laction = .TRUE.
         IF( TRIM( sn_snd_cond%clcat ) == 'yes' ) ssnd(nmodsbc)%fld(jps_kice)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_cond%cldes;'//sn_snd_cond%cldes )
      END SELECT
      !
      !                                                      ! ------------------------- !
      !                                                      !     Sea surface height    !
      !                                                      ! ------------------------- !
      ssnd(nmodsbc)%fld(jps_wlev)%clname = 'O_Wlevel' ;  IF( TRIM(sn_snd_wlev%cldes) == 'coupled' )   ssnd(nmodsbc)%fld(jps_wlev)%laction = .TRUE.

      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - snd by opa !
      !                                                      ! ------------------------------- !
      ssnd(nmodsbc)%fld(jps_ssh   )%clname = 'O_SSHght'
      ssnd(nmodsbc)%fld(jps_soce  )%clname = 'O_SSSal'
      ssnd(nmodsbc)%fld(jps_e3t1st)%clname = 'O_E3T1st'
      ssnd(nmodsbc)%fld(jps_fraqsr)%clname = 'O_FraQsr'
      !
      IF( nn_components == jp_iam_oce ) THEN
         ssnd(nmodsbc)%fld(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd(nmodsbc)%fld( (/jps_toce, jps_soce, jps_ssh, jps_fraqsr, jps_ocx1, jps_ocy1/) )%laction = .TRUE.
         ssnd(nmodsbc)%fld( jps_e3t1st )%laction = .NOT.lk_linssh
         ! vector definition: not used but cleaner...
         ssnd(nmodsbc)%fld(jps_ocx1)%clgrid  = 'U'        ! oce components given at U-point
         ssnd(nmodsbc)%fld(jps_ocy1)%clgrid  = 'V'        !           and           V-point
         sn_snd_crt%clvgrd = 'U,V'
         sn_snd_crt%clvor = 'local grid'
         sn_snd_crt%clvref = 'spherical'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'  sent fields to SAS component '
            WRITE(numout,*)'               sea surface temperature (T before, Celsius) '
            WRITE(numout,*)'               sea surface salinity '
            WRITE(numout,*)'               surface currents U,V on local grid and spherical coordinates'
            WRITE(numout,*)'               sea surface height '
            WRITE(numout,*)'               thickness of first ocean T level '
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - snd by sas !
      !                                                      ! ------------------------------- !
      ssnd(nmodsbc)%fld(jps_sflx  )%clname = 'I_SFLX'
      ssnd(nmodsbc)%fld(jps_fice2 )%clname = 'IIceFrc'
      ssnd(nmodsbc)%fld(jps_qsroce)%clname = 'I_QsrOce'
      ssnd(nmodsbc)%fld(jps_qnsoce)%clname = 'I_QnsOce'
      ssnd(nmodsbc)%fld(jps_oemp  )%clname = 'IOEvaMPr'
      ssnd(nmodsbc)%fld(jps_otx1  )%clname = 'I_OTaux1'
      ssnd(nmodsbc)%fld(jps_oty1  )%clname = 'I_OTauy1'
      ssnd(nmodsbc)%fld(jps_rnf   )%clname = 'I_Runoff'
      ssnd(nmodsbc)%fld(jps_taum  )%clname = 'I_TauMod'
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) ssnd(nmodsbc)%fld(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd(nmodsbc)%fld( (/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/) )%laction = .TRUE.
         !
         ! Change first letter to couple with atmosphere if already coupled with sea_ice
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_SSTSST sent by OCE to SAS and therefore S_SSTSST sent by SAS to the Atmosphere
         DO jn = 1, jpsnd
            IF( ssnd(nmodsbc)%fld(jn)%clname(1:1) == "O" ) ssnd(nmodsbc)%fld(jn)%clname = "S"//ssnd(nmodsbc)%fld(jn)%clname(2:LEN(ssnd(nmodsbc)%fld(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  sent fields to OCE component '
            ELSE
               WRITE(numout,*)'  Additional sent fields to OCE component : '
            ENDIF
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V components'
            WRITE(numout,*)'                  wind stress module'
         ENDIF
      ENDIF

      !
      ! =================================== !
      !   define variables for the coupler  !
      ! =================================== !
      CALL cpl_var(jprcv, jpsnd, nn_cplmodel, nmodsbc)

      IF(ln_usecplmask) THEN
         CALL iom_open( 'cplmask', inum )
         CALL iom_get( inum, jpdom_unknown, 'cplmask', xcplmask(A2D(0),1,1:nn_cplmodel),   &
            &          kstart = (/ mig(Nis0,0),mjg(Njs0,0),1 /), kcount = (/ Ni_0,Nj_0,nn_cplmodel /) )
         CALL iom_close( inum )
         xcplmask(A2D(0),1,0) = 1. - SUM( xcplmask(A2D(0),1,1:nn_cplmodel), dim = 3 )
      ELSE
         xcplmask(A2D(0),1,:) = 1.
      ENDIF
      !
   END SUBROUTINE sbc_cpl_init


   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in frcv array). called at each time step.
      !!                OASIS controls if there is something do receive or not. nrcvinfo contains the info
      !!                to know if the field was really received or not
      !!
      !!              --> If ocean stress was really received:
      !!
      !!                  - transform the received ocean stress vector from the received
      !!                 referential and grid into an atmosphere-ocean stress in
      !!                 the (i,j) ocean referencial and at the ocean velocity point.
      !!                    The received stress are :
      !!                     - defined by 3 components (if cartesian coordinate)
      !!                            or by 2 components (if spherical)
      !!                     - oriented along geographical   coordinate (if eastward-northward)
      !!                            or  along the local grid coordinate (if local grid)
      !!                     - given at U- and V-point, resp.   if received on 2 grids
      !!                            or at T-point               if received on 1 grid
      !!                    Therefore and if necessary, they are successively
      !!                  processed in order to obtain them
      !!                     first  as  2 components on the sphere
      !!                     second as  2 components oriented along the local grid
      !!                     third  as  2 components on the U,V grid
      !!
      !!              -->
      !!
      !!              - In 'ocean only' case, non solar and solar ocean heat fluxes
      !!             and total ocean freshwater fluxes
      !!
      !! ** Method  :   receive all fields from the atmosphere and transform
      !!              them into ocean surface boundary condition fields
      !!
      !! ** Action  :   update  utau, vtau   ocean stress at T-point
      !!                        taum         wind stress module at T-point
      !!                        wndm         wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!                        qns          non solar heat fluxes including emp heat content    (ocean only case)
      !!                                     and the latent heat flux of solid precip. melting
      !!                        qsr          solar ocean heat fluxes   (ocean only case)
      !!                        emp          upward mass flux [evap. - precip. (- runoffs) (- calving)] (ocean only case)
      !!----------------------------------------------------------------------
      USE zdf_oce,  ONLY :   ln_zdfswm
      !
      INTEGER, INTENT(in) ::   kt          ! ocean model time step index
      INTEGER, INTENT(in) ::   k_fsbc      ! frequency of sbc (-> ice model) computation
      INTEGER, INTENT(in) ::   k_ice       ! ice management in the sbc (=0/1/2/3)
      INTEGER, INTENT(in) ::   Kbb, Kmm    ! ocean model time level indices
      !!
      LOGICAL  ::   llnewtx, llnewtau      ! update wind stress components and module??
      INTEGER  ::   ji, jj, jn             ! dummy loop indices
      INTEGER  ::   isec                   ! number of seconds since nit000 (assuming rdt did not change since nit000)
      REAL(wp) ::   zcumulneg, zcumulpos   ! temporary scalars
      REAL(wp) ::   zcoef                  ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3        ! drag coefficient
      REAL(wp) ::   zzx, zzy               ! temporary variables
      REAL(wp) ::   r1_grau                ! = 1.e0 / (grav * rho0)
      REAL(wp), DIMENSION(A2D(0)) ::   ztx, zty, zmsk, zemp
      REAL(wp), DIMENSION(A2D(0)) ::   zqns, zqsr, zcloud_fra
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
      !   cannot be done in the init phase when we use agrif as cpl_freq requires that oasis_enddef is done
         ncpl_qsr_freq = cpl_freq( 'O_QsrOce', nmodsbc ) + cpl_freq( 'O_QsrMix', nmodsbc ) + &
                     &   cpl_freq( 'I_QsrOce', nmodsbc ) + cpl_freq( 'I_QsrMix', nmodsbc )
         IF( ln_dm2dc .AND. ncpl_qsr_freq /= 86400 )   &
            &   CALL ctl_stop( 'sbc_cpl_rcv: diurnal cycle reconstruction (ln_dm2dc) needs daily couping for solar radiation' )

         IF ( ln_wave .AND. nn_components == 0 ) THEN
            ncpl_qsr_freq = 1;
            WRITE(numout,*) 'ncpl_qsr_freq is set to 1 when coupling NEMO with wave (without SAS) '
         ENDIF
      ENDIF
      !
      IF( ln_mixcpl )   zmsk (A2D(0)) = 1. - xcplmask (A2D(0),1,0)
      !
      !                                                      ! ======================================================= !
      !                                                      ! Receive all the atmos. fields (including ice information)
      !                                                      ! ======================================================= !
      isec = ( kt - nit000 ) * NINT( rn_Dt )                      ! date of exchanges
      DO jn = 1, jprcv                                          ! received fields sent by the atmosphere
         IF( srcv(nmodsbc)%fld(jn)%laction )   CALL cpl_rcv( nmodsbc, jn, isec, frcv(jn)%z3, nrcvinfo(jn), xcplmask(A2D(0),1:1,1:nn_cplmodel) )
      END DO

      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_otx1)%laction ) THEN                      !  ocean stress components  !
         !                                                   ! ========================= !
         ! define frcv(jpr_otx1)%z3(A2D(0),1) and frcv(jpr_oty1)%z3(A2D(0),1): stress at U/V point along model grid
         ! => need to be done only when we receive the field
         IF(  nrcvinfo(jpr_otx1) == OASIS_Rcv ) THEN
            !
            IF( TRIM( sn_rcv_tau%clvref ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               !
               CALL geo2oce( frcv(jpr_otx1)%z3(A2D(0),1), frcv(jpr_oty1)%z3(A2D(0),1), frcv(jpr_otz1)%z3(A2D(0),1), 'T', ztx, zty )
               frcv(jpr_otx1)%z3(A2D(0),1) = ztx(A2D(0))   ! overwrite 1st comp. on the 1st grid
               frcv(jpr_oty1)%z3(A2D(0),1) = zty(A2D(0))   ! overwrite 2nd comp. on the 1st grid
               !
            ENDIF
            !
            IF( TRIM( sn_rcv_tau%clvor ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(jpr_otx1)%z3(A2D(0),1), frcv(jpr_oty1)%z3(A2D(0),1), 'T', 'en->i', ztx )
               CALL rot_rep( frcv(jpr_otx1)%z3(A2D(0),1), frcv(jpr_oty1)%z3(A2D(0),1), 'T', 'en->j', zty )
               frcv(jpr_otx1)%z3(A2D(0),1) = ztx(A2D(0))      ! overwrite 1st component on the 1st grid
               frcv(jpr_oty1)%z3(A2D(0),1) = zty(A2D(0))      ! overwrite 2nd component on the 2nd grid
            ENDIF
            !
            llnewtx = .TRUE.
         ELSE
            llnewtx = .FALSE.
         ENDIF
         !                                                   ! ========================= !
      ELSE                                                   !   No dynamical coupling   !
         !                                                   ! ========================= !
         frcv(jpr_otx1)%z3(A2D(0),1) = 0.e0                               ! here simply set to zero
         frcv(jpr_oty1)%z3(A2D(0),1) = 0.e0                               ! an external read in a file can be added instead
         llnewtx = .TRUE.
         !
      ENDIF
      !                                                      ! ========================= !
      !                                                      !    wind stress module     !   (taum)
      !                                                      ! ========================= !
      IF( .NOT. srcv(nmodsbc)%fld(jpr_taum)%laction ) THEN                    ! compute wind stress module from its components if not received
         ! => need to be done only when otx1 was changed
         IF( llnewtx ) THEN
            DO_2D( 0, 0, 0, 0 )
               zzx = frcv(jpr_otx1)%z3(ji,jj,1)
               zzy = frcv(jpr_oty1)%z3(ji,jj,1) 
               frcv(jpr_taum)%z3(ji,jj,1) = 0.5 * SQRT( zzx * zzx + zzy * zzy )
            END_2D
            llnewtau = .TRUE.
         ELSE
            llnewtau = .FALSE.
         ENDIF
      ELSE
         llnewtau = nrcvinfo(jpr_taum) == OASIS_Rcv
         ! Stress module can be negative when received (interpolation problem)
         IF( llnewtau ) THEN
            frcv(jpr_taum)%z3(A2D(0),1) = MAX( 0._wp, frcv(jpr_taum)%z3(A2D(0),1) )
         ENDIF
      ENDIF
      !
      !                                                      ! ========================= !
      !                                                      !      10 m wind speed      !   (wndm)
      !                                                      ! ========================= !
      IF( .NOT. srcv(nmodsbc)%fld(jpr_w10m)%laction ) THEN                    ! compute wind spreed from wind stress module if not received
         ! => need to be done only when taumod was changed
         IF( llnewtau ) THEN
            zcoef = 1. / ( zrhoa * zcdrag )
            DO_2D( 0, 0, 0, 0 )
               frcv(jpr_w10m)%z3(ji,jj,1) = SQRT( frcv(jpr_taum)%z3(ji,jj,1) * zcoef )
            END_2D
         ENDIF
      ENDIF
!!$      !                                                      ! ========================= !
!!$      SELECT CASE( TRIM( sn_rcv_clouds%cldes ) )             !       cloud fraction      !
!!$      !                                                      ! ========================= !
!!$      cloud_fra(A2D(0)) = frcv(jpr_clfra)*z3(A2D(0),1)
!!$      END SELECT
!!$
      zcloud_fra(A2D(0)) = pp_cldf   ! should be real cloud fraction instead (as in the bulk) but needs to be read from atm.
      IF( ln_mixcpl ) THEN
         cloud_fra(A2D(0)) = cloud_fra(A2D(0)) * xcplmask(A2D(0),1,0) + zcloud_fra(A2D(0))* zmsk(A2D(0))
      ELSE
         cloud_fra(A2D(0)) = zcloud_fra(A2D(0))
      ENDIF
      !                                                      ! ========================= !
      ! u(v)tau and taum will be modified by ice model
      ! -> need to be reset before each call of the ice/fsbc
      IF( MOD( kt-1, k_fsbc ) == 0 ) THEN
         !
         IF( ln_mixcpl ) THEN
            utau(A2D(0)) = utau(A2D(0)) * xcplmask(A2D(0),1,0) + frcv(jpr_otx1)%z3(A2D(0),1) * zmsk(A2D(0))
            vtau(A2D(0)) = vtau(A2D(0)) * xcplmask(A2D(0),1,0) + frcv(jpr_oty1)%z3(A2D(0),1) * zmsk(A2D(0))
            taum(A2D(0)) = taum(A2D(0)) * xcplmask(A2D(0),1,0) + frcv(jpr_taum)%z3(A2D(0),1) * zmsk(A2D(0))
            wndm(A2D(0)) = wndm(A2D(0)) * xcplmask(A2D(0),1,0) + frcv(jpr_w10m)%z3(A2D(0),1) * zmsk(A2D(0))
         ELSE
            utau(A2D(0)) = frcv(jpr_otx1)%z3(A2D(0),1)
            vtau(A2D(0)) = frcv(jpr_oty1)%z3(A2D(0),1)
            taum(A2D(0)) = frcv(jpr_taum)%z3(A2D(0),1)
            wndm(A2D(0)) = frcv(jpr_w10m)%z3(A2D(0),1)
         ENDIF
         !
         CALL lbc_lnk( 'sbccpl', utau, 'T', -1.0_wp, vtau, 'T', -1.0_wp, ldfull = .TRUE. )
         !
         CALL iom_put( "taum_oce", taum )   ! output wind stress module
      ENDIF

      !                                                      ! ================== !
      !                                                      ! atmosph. CO2 (ppm) !
      !                                                      ! ================== !
      IF( srcv(nmodsbc)%fld(jpr_co2)%laction )   atm_co2(A2D(0)) = frcv(jpr_co2)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      ! Mean Sea Level Pressure   !   (taum)
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_mslp)%laction ) THEN                    ! UKMO SHELF effect of atmospheric pressure on SSH
          IF( kt /= nit000 )   ssh_ibb(A2D(0)) = ssh_ib(A2D(0))    !* Swap of ssh_ib fields

          r1_grau = 1.e0 / (grav * rho0)               !* constant for optimization
          ssh_ib(A2D(0)) = - ( frcv(jpr_mslp)%z3(A2D(0),1) - rpref ) * r1_grau    ! equivalent ssh (inverse barometer)
          apr   (A2D(0)) =     frcv(jpr_mslp)%z3(A2D(0),1)                        !atmospheric pressure
          CALL lbc_lnk( 'sbccpl', ssh_ib, 'T', 1.0_wp, apr, 'T', 1.0_wp, ldfull = .TRUE. )
          
          IF( kt == nit000 ) ssh_ibb(1:jpi,1:jpj) = ssh_ib(1:jpi,1:jpj)  ! correct this later (read from restart if possible)
      ENDIF
      !
      IF( ln_sdw .OR. ln_ice_wav ) THEN   ! Either Stokes drift correction or wave-ice interaction activated
      !                                                      ! ========================= !
      !                                                      !  Significant wave height  !
      !                                                      ! ========================= !
         IF( srcv(nmodsbc)%fld(jpr_hsig)%laction ) THEN
            hsw(A2D(0)) = frcv(jpr_hsig)%z3(A2D(0),1)
            CALL lbc_lnk( 'sbccpl', hsw, 'T', 1.0_wp, ldfull = .TRUE. )
         ENDIF
      ENDIF
      !
      IF( ln_ice_wav ) THEN   ! Wave-ice interaction activated
      !                                                      ! ========================= !
      !                                                      !    Wave peak frequency    !
      !                                                      ! ========================= !
         IF( srcv(nmodsbc)%fld(jpr_wpf)%laction ) THEN
            wpf(A2D(0)) = frcv(jpr_wpf)%z3(A2D(0),1)
            CALL lbc_lnk( 'sbccpl', wpf, 'T', 1.0_wp, ldfull = .TRUE. )
         ENDIF
      ENDIF
      !
      IF( ln_sdw ) THEN  ! Stokes Drift correction activated
      !                                                      ! ========================= !
      !                                                      !       Stokes drift        !
      !                                                      ! ========================= !
         IF( srcv(nmodsbc)%fld(jpr_sdrftx)%laction ) THEN
            ut0sd(A2D(0)) = frcv(jpr_sdrftx)%z3(A2D(0),1)
            vt0sd(A2D(0)) = frcv(jpr_sdrfty)%z3(A2D(0),1)
            CALL lbc_lnk( 'sbccpl', ut0sd, 'T', -1.0_wp, vt0sd, 'T', -1.0_wp, ldfull = .TRUE. )
         ENDIF
      !
      !                                                      ! ========================= !
      !                                                      !      Wave mean period     !
      !                                                      ! ========================= !
         IF( srcv(nmodsbc)%fld(jpr_wper)%laction ) THEN
            wmp(A2D(0)) = frcv(jpr_wper)%z3(A2D(0),1)
            CALL lbc_lnk( 'sbccpl', wmp, 'T', 1.0_wp, ldfull = .TRUE. )
         ENDIF
      !                                                      ! ========================= !
      !                                                      !    Vertical mixing Qiao   !
      !                                                      ! ========================= !
         IF( srcv(nmodsbc)%fld(jpr_wnum)%laction .AND. ln_zdfswm ) wnum(A2D(0)) = frcv(jpr_wnum)%z3(A2D(0),1)

         ! Calculate the 3D Stokes drift both in coupled and not fully uncoupled mode
         IF( srcv(nmodsbc)%fld(jpr_sdrftx)%laction .OR. srcv(nmodsbc)%fld(jpr_sdrfty)%laction .OR. &
             srcv(nmodsbc)%fld(jpr_wper)%laction .OR. srcv(nmodsbc)%fld(jpr_hsig)%laction )   THEN
            CALL sbc_stokes( Kmm )
         ENDIF
      ENDIF
      !                                                      ! ========================= !
      !                                                      ! Stress adsorbed by waves  !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_wstrf)%laction .AND. ln_tauoc )  tauoc_wave(A2D(0)) = frcv(jpr_wstrf)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      !   Wave drag coefficient   !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_wdrag)%laction .AND. ln_cdgw )   cdn_wave(A2D(0)) = frcv(jpr_wdrag)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      !   Chranock coefficient    !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_charn)%laction .AND. ln_charn )  charn(A2D(0)) = frcv(jpr_charn)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      ! net wave-supported stress !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_tawx)%laction .AND. ln_taw )     tawx(A2D(0)) = frcv(jpr_tawx)%z3(A2D(0),1)
      IF( srcv(nmodsbc)%fld(jpr_tawy)%laction .AND. ln_taw )     tawy(A2D(0)) = frcv(jpr_tawy)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      !wave to ocean momentum flux!
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_twox)%laction .AND. ln_taw )     twox(A2D(0)) = frcv(jpr_twox)%z3(A2D(0),1)
      IF( srcv(nmodsbc)%fld(jpr_twoy)%laction .AND. ln_taw )     twoy(A2D(0)) = frcv(jpr_twoy)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      !    wave TKE flux at sfc   !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_phioc)%laction .AND. ln_phioc )     phioc(A2D(0)) = frcv(jpr_phioc)%z3(A2D(0),1)
      !
      !                                                      ! ========================= !
      !                                                      !      Bernoulli head       !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_bhd)%laction .AND. ln_bern_srfc ) THEN
         bhd_wave(A2D(0)) = frcv(jpr_bhd)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', bhd_wave, 'T', 1.0_wp, ldfull = .TRUE. )
      ENDIF
      !                                                      ! ========================= !
      !                                                      !      Stokes transport     !
      !                                                      ! ========================= !
      IF( srcv(nmodsbc)%fld(jpr_tusd)%laction .AND. ln_breivikFV_2016 ) THEN
         tusd(A2D(0)) = frcv(jpr_tusd)%z3(A2D(0),1)
         tvsd(A2D(0)) = frcv(jpr_tvsd)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', tusd, 'T', -1.0_wp, tvsd, 'T', -1.0_wp, ldfull = .TRUE. )
      ENDIF
      !
      !  Fields received by SAS when OASIS coupling
      !  (arrays no more filled at sbcssm stage)
      !                                                      ! ================== !
      !                                                      !        SSS         !
      !                                                      ! ================== !
      IF( srcv(nmodsbc)%fld(jpr_soce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         sss_m(A2D(0)) = frcv(jpr_soce)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', sss_m, 'T', 1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'sss_m', sss_m )
      ENDIF
      !
      !                                                      ! ================== !
      !                                                      !        SST         !
      !                                                      ! ================== !
      IF( srcv(nmodsbc)%fld(jpr_toce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         IF( srcv(nmodsbc)%fld(jpr_soce)%laction .AND. l_useCT ) THEN     ! make sure that sst_m is the potential temperature
            CALL eos_pt_from_ct( frcv(jpr_toce)%z3(A2D(0),1), sss_m(A2D(0)), sst_m(A2D(0)), kbnd=0 )
         ELSE
            sst_m(A2D(0)) = frcv(jpr_toce)%z3(A2D(0),1)
         ENDIF
         CALL lbc_lnk( 'sbccpl', sst_m, 'T', 1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'sst_m', sst_m )
      ENDIF
      !                                                      ! ================== !
      !                                                      !        SSH         !
      !                                                      ! ================== !
      IF( srcv(nmodsbc)%fld(jpr_ssh )%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssh_m(A2D(0)) = frcv(jpr_ssh )%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', ssh_m, 'T', 1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'ssh_m', ssh_m )
      ENDIF
      !                                                      ! ================== !
      !                                                      !  surface currents  !
      !                                                      ! ================== !
      IF( srcv(nmodsbc)%fld(jpr_ocx1)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssu_m(A2D(0)) = frcv(jpr_ocx1)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', ssu_m, 'U', -1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'ssu_m', ssu_m )
         uu(1:jpi,1:jpj,1,Kbb) = ssu_m(1:jpi,1:jpj)          ! will be used in icestp in the call of ice_update_tau
         uu(1:jpi,1:jpj,1,Kmm) = ssu_m(1:jpi,1:jpj)          ! will be used in sbc_cpl_snd if atmosphere coupling
      ENDIF
      IF( srcv(nmodsbc)%fld(jpr_ocy1)%laction ) THEN
         ssv_m(A2D(0)) = frcv(jpr_ocy1)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', ssv_m, 'V', -1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'ssv_m', ssv_m )
         vv(1:jpi,1:jpj,1,Kbb) = ssv_m(1:jpi,1:jpj)          ! will be used in icestp in the call of ice_update_tau
         vv(1:jpi,1:jpj,1,Kmm) = ssv_m(1:jpi,1:jpj)          ! will be used in sbc_cpl_snd if atmosphere coupling
      ENDIF
      !                                                      ! ======================== !
      !                                                      !  first T level thickness !
      !                                                      ! ======================== !
      IF( srcv(nmodsbc)%fld(jpr_e3t1st )%laction ) THEN                   ! received by sas in case of opa <-> sas coupling
         e3t_m(A2D(0)) = frcv(jpr_e3t1st )%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', e3t_m, 'T', 1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'e3t_m', e3t_m )
      ENDIF
      !                                                      ! ================================ !
      !                                                      !  fraction of solar net radiation !
      !                                                      ! ================================ !
      IF( srcv(nmodsbc)%fld(jpr_fraqsr)%laction ) THEN                    ! received by sas in case of opa <-> sas coupling
         frq_m(A2D(0)) = frcv(jpr_fraqsr)%z3(A2D(0),1)
         CALL lbc_lnk( 'sbccpl', frq_m, 'T', 1.0_wp, ldfull = .TRUE. )
         CALL iom_put(  'frq_m', frq_m )
      ENDIF

      !                                                      ! ========================= !
      IF( k_ice <= 1 .AND. MOD( kt-1, k_fsbc ) == 0 ) THEN   !  heat & freshwater fluxes ! (Ocean only case)
         !                                                   ! ========================= !
         !
         !                                                       ! total freshwater fluxes over the ocean (emp)
         IF( srcv(nmodsbc)%fld(jpr_oemp)%laction .OR. srcv(nmodsbc)%fld(jpr_rain)%laction ) THEN
            SELECT CASE( TRIM( sn_rcv_emp%cldes ) )                                    ! evaporation - precipitation
            CASE( 'conservative' )
               zemp(A2D(0)) = frcv(jpr_tevp)%z3(A2D(0),1) - ( frcv(jpr_rain)%z3(A2D(0),1) + frcv(jpr_snow)%z3(A2D(0),1) )
            CASE( 'oce only', 'oce and ice' )
               zemp(A2D(0)) = frcv(jpr_oemp)%z3(A2D(0),1)
            CASE default
               CALL ctl_stop( 'sbc_cpl_rcv: wrong definition of sn_rcv_emp%cldes' )
            END SELECT
         ELSE
            zemp(A2D(0)) = 0._wp
         ENDIF
         !
         !                                                        ! runoffs and calving (added in emp)
         IF( srcv(nmodsbc)%fld(jpr_rnf)%laction )     rnf(A2D(0)) = frcv(jpr_rnf)%z3(A2D(0),1)
         IF( srcv(nmodsbc)%fld(jpr_cal)%laction )     zemp(A2D(0)) = zemp(A2D(0)) - frcv(jpr_cal)%z3(A2D(0),1)

         IF( srcv(nmodsbc)%fld(jpr_icb)%laction )  THEN
             fwficb(A2D(0)) = frcv(jpr_icb)%z3(A2D(0),1)
             rnf(A2D(0)) = rnf(A2D(0)) + fwficb(A2D(0))   ! iceberg added to runfofs
         ENDIF
         !
         ! ice shelf fwf
         IF( srcv(nmodsbc)%fld(jpr_isf)%laction )  THEN
            fwfisf_oasis(A2D(0)) = frcv(jpr_isf)%z3(A2D(0),1)  ! fresh water flux from the isf to the ocean ( > 0 = melting )
         END IF

         IF( ln_mixcpl ) THEN   ;   emp(A2D(0)) = emp(A2D(0)) * xcplmask(A2D(0),1,0) + zemp(A2D(0)) * zmsk(A2D(0))
         ELSE                   ;   emp(A2D(0)) =                                    zemp(A2D(0))
         ENDIF
         !
         IF(ln_rnf) THEN
            CALL lbc_lnk( 'sbccpl', emp, 'T', 1.0_wp, rnf, 'T', 1.0_wp, ldfull = .TRUE. )  
         ELSE
            CALL lbc_lnk( 'sbccpl', emp, 'T', 1.0_wp, ldfull = .TRUE. )  
         ENDIF
         !
         !                                                       ! non solar heat flux over the ocean (qns)
         IF(      srcv(nmodsbc)%fld(jpr_qnsoce)%laction ) THEN   ;   zqns(A2D(0)) = frcv(jpr_qnsoce)%z3(A2D(0),1)
         ELSE IF( srcv(nmodsbc)%fld(jpr_qnsmix)%laction ) THEN   ;   zqns(A2D(0)) = frcv(jpr_qnsmix)%z3(A2D(0),1)
         ELSE                                       ;   zqns(A2D(0)) = 0._wp
         ENDIF
         ! update qns over the free ocean with:
         IF( nn_components /= jp_iam_oce ) THEN
            zqns(A2D(0)) =  zqns(A2D(0)) - zemp(A2D(0)) * sst_m(A2D(0)) * rcp      ! remove heat content due to mass flux (assumed to be at SST)
            IF( srcv(nmodsbc)%fld(jpr_snow  )%laction ) THEN
               zqns(A2D(0)) = zqns(A2D(0)) - frcv(jpr_snow)%z3(A2D(0),1) * rLfus   ! energy for melting solid precipitation over the free ocean
            ENDIF
         ENDIF
         !
         IF( srcv(nmodsbc)%fld(jpr_icb)%laction )  zqns(A2D(0)) = zqns(A2D(0)) - frcv(jpr_icb)%z3(A2D(0),1) * rLfus ! remove heat content associated to iceberg melting
         !
         IF( ln_mixcpl ) THEN   ;   qns(A2D(0)) = qns(A2D(0)) * xcplmask(A2D(0),1,0) + zqns(A2D(0)) * zmsk(A2D(0))
         ELSE                   ;   qns(A2D(0)) =                                    zqns(A2D(0))
         ENDIF

         !                                                       ! solar flux over the ocean          (qsr)
         IF     ( srcv(nmodsbc)%fld(jpr_qsroce)%laction ) THEN   ;   zqsr(A2D(0)) = frcv(jpr_qsroce)%z3(A2D(0),1)
         ELSE IF( srcv(nmodsbc)%fld(jpr_qsrmix)%laction ) then   ;   zqsr(A2D(0)) = frcv(jpr_qsrmix)%z3(A2D(0),1)
         ELSE                                       ;   zqsr(A2D(0)) = 0._wp
         ENDIF
         IF( ln_dm2dc .AND. ln_cpl )   zqsr(A2D(0)) = sbc_dcy( zqsr )   ! modify qsr to include the diurnal cycle
         IF( ln_mixcpl ) THEN   ;   qsr(A2D(0)) = qsr(A2D(0)) * xcplmask(A2D(0),1,0) + zqsr(A2D(0)) * zmsk(A2D(0))
         ELSE                   ;   qsr(A2D(0)) =                                    zqsr(A2D(0))
         ENDIF
         !
         ! salt flux over the ocean (received by opa in case of opa <-> sas coupling)
         IF( srcv(nmodsbc)%fld(jpr_sflx )%laction )   sfx(A2D(0)) = frcv(jpr_sflx  )%z3(A2D(0),1)
         ! Ice cover  (received by opa in case of opa <-> sas coupling)
         IF( srcv(nmodsbc)%fld(jpr_fice )%laction ) THEN
            fr_i(A2D(0)) = frcv(jpr_fice )%z3(A2D(0),1)
            CALL lbc_lnk( 'sbccpl', fr_i, 'T', 1.0_wp, ldfull = .TRUE. )
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_cpl_rcv


   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_tau  ***
      !!
      !! ** Purpose :   provide the stress over sea-ice in coupled mode
      !!
      !! ** Method  :   transform the received stress from the atmosphere into
      !!             an atmosphere-ice stress in the (i,j) ocean referencial
      !!             and at the velocity point of the sea-ice model:
      !!                'C'-grid : i- (j-) components given at U- (V-) point
      !!
      !!                The received stress are :
      !!                 - defined by 3 components (if cartesian coordinate)
      !!                        or by 2 components (if spherical)
      !!                 - oriented along geographical   coordinate (if eastward-northward)
      !!                        or  along the local grid coordinate (if local grid)
      !!                 - given at U- and V-point, resp.   if received on 2 grids
      !!                        or at a same point (T or I) if received on 1 grid
      !!                Therefore and if necessary, they are successively
      !!             processed in order to obtain them
      !!                 first  as  2 components on the sphere
      !!                 second as  2 components oriented along the local grid
      !!                 third  as  2 components on the ice grid point
      !!
      !!                Except in 'oce and ice' case, only one vector stress field
      !!             is received. It has already been processed in sbc_cpl_rcv
      !!             so that it is now defined as (i,j) components given at U-
      !!             and V-points, respectively.
      !!
      !! ** Action  :   return ptau_i, ptau_j, the stress over the ice
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout), DIMENSION(A2D(nn_hls)) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(inout), DIMENSION(A2D(nn_hls)) ::   p_tauj   !                   at T-point
      !!
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER ::   itx      ! index of taux over ice
      REAL(wp)                    ::   zztmp1, zztmp2
      REAL(wp), DIMENSION(A2D(0)) ::   ztx, zty
      !!----------------------------------------------------------------------
      !
#if defined key_si3 || defined key_cice
      !
      IF( srcv(nmodsbc)%fld(jpr_itx1)%laction ) THEN   ;   itx =  jpr_itx1
      ELSE                                ;   itx =  jpr_otx1
      ENDIF

      ! do something only if we just received the stress from atmosphere
      IF(  nrcvinfo(itx) == OASIS_Rcv ) THEN
         !                                                      ! ======================= !
         IF( srcv(nmodsbc)%fld(jpr_itx1)%laction ) THEN                      !   ice stress received   !
            !                                                   ! ======================= !
            !
            IF( TRIM( sn_rcv_tau%clvref ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               CALL geo2oce(  frcv(jpr_itx1)%z3(A2D(0),1), frcv(jpr_ity1)%z3(A2D(0),1), frcv(jpr_itz1)%z3(A2D(0),1), 'T', ztx, zty )
               frcv(jpr_itx1)%z3(A2D(0),1) = ztx(A2D(0))   ! overwrite 1st comp. on the 1st grid
               frcv(jpr_ity1)%z3(A2D(0),1) = zty(A2D(0))   ! overwrite 2nd comp. on the 1st grid
               !
            ENDIF
            !
            IF( TRIM( sn_rcv_tau%clvor ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(jpr_itx1)%z3(A2D(0),1), frcv(jpr_ity1)%z3(A2D(0),1), 'T', 'en->i', ztx )
               CALL rot_rep( frcv(jpr_itx1)%z3(A2D(0),1), frcv(jpr_ity1)%z3(A2D(0),1), 'T', 'en->j', zty )
               frcv(jpr_itx1)%z3(A2D(0),1) = ztx(A2D(0))      ! overwrite 1st component on the 1st grid
               frcv(jpr_ity1)%z3(A2D(0),1) = zty(A2D(0))      ! overwrite 2nd component on the 1st grid
            ENDIF
            !                                                   ! ======================= !
         ELSE                                                   !     use ocean stress    !
            !                                                   ! ======================= !
            frcv(jpr_itx1)%z3(A2D(0),1) = frcv(jpr_otx1)%z3(A2D(0),1)
            frcv(jpr_ity1)%z3(A2D(0),1) = frcv(jpr_oty1)%z3(A2D(0),1)
            !
         ENDIF
         !                                                      ! ======================= !
         !                                                      !     put on ice grid     !
         !                                                      ! ======================= !
         p_taui(A2D(0)) = frcv(jpr_itx1)%z3(A2D(0),1)
         p_tauj(A2D(0)) = frcv(jpr_ity1)%z3(A2D(0),1)

      ENDIF
      !
#endif
      !
   END SUBROUTINE sbc_cpl_ice_tau


   SUBROUTINE sbc_cpl_ice_flx( kt, picefr, palbi, psst, pist, phs, phi )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_flx  ***
      !!
      !! ** Purpose :   provide the heat and freshwater fluxes of the ocean-ice system
      !!
      !! ** Method  :   transform the fields received from the atmosphere into
      !!             surface heat and fresh water boundary condition for the
      !!             ice-ocean system. The following fields are provided:
      !!               * total non solar, solar and freshwater fluxes (qns_tot,
      !!             qsr_tot and emp_tot) (total means weighted ice-ocean flux)
      !!             NB: emp_tot include runoffs and calving.
      !!               * fluxes over ice (qns_ice, qsr_ice, emp_ice) where
      !!             emp_ice = sublimation - solid precipitation as liquid
      !!             precipitation are re-routed directly to the ocean and
      !!             calving directly enter the ocean (runoffs are read but included in trasbc.F90)
      !!               * solid precipitation (sprecip), used to add to qns_tot
      !!             the heat lost associated to melting solid precipitation
      !!             over the ocean fraction.
      !!               * heat content of rain, snow and evap can also be provided,
      !!             otherwise heat flux associated with these mass flux are
      !!             guessed (qemp_oce, qemp_ice)
      !!
      !!             - the fluxes have been separated from the stress as
      !!               (a) they are updated at each ice time step compare to
      !!               an update at each coupled time step for the stress, and
      !!               (b) the conservative computation of the fluxes over the
      !!               sea-ice area requires the knowledge of the ice fraction
      !!               after the ice advection and before the ice thermodynamics,
      !!               so that the stress is updated before the ice dynamics
      !!               while the fluxes are updated after it.
      !!
      !! ** Details
      !!             qns_tot = (1-a) * qns_oce + a * qns_ice               => provided
      !!                     + qemp_oce + qemp_ice                         => recalculated and added up to qns
      !!
      !!             qsr_tot = (1-a) * qsr_oce + a * qsr_ice               => provided
      !!
      !!             emp_tot = emp_oce + emp_ice                           => calving is provided and added to emp_tot (and emp_oce).
      !!                                                                      runoff (which includes rivers+icebergs) and iceshelf
      !!                                                                      are provided but not included in emp here. Only runoff will
      !!                                                                      be included in emp in other parts of NEMO code
      !!
      !! ** Note : In case of the ice-atm coupling with conduction fluxes (such as Jules interface for the Met-Office),
      !!              qsr_ice and qns_ice are not provided and they are not supposed to be used in the ice code.
      !!              However, by precaution we also "fake" qns_ice and qsr_ice this way:
      !!              qns_ice = qml_ice + qcn_ice ??
      !!              qsr_ice = qtr_ice_top ??
      !!
      !! ** Action  :   update at each nf_ice time step:
      !!                   qns_tot, qsr_tot  non-solar and solar total heat fluxes
      !!                   qns_ice, qsr_ice  non-solar and solar heat fluxes over the ice
      !!                   emp_tot           total evaporation - precipitation(liquid and solid) (-calving)
      !!                   emp_ice           ice sublimation - solid precipitation over the ice
      !!                   dqns_ice          d(non-solar heat flux)/d(Temperature) over the ice
      !!                   sprecip           solid precipitation over the ocean
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)                                          ::   kt         ! ocean model time step index (only for a_i_last_couple)
      REAL(wp), INTENT(in)   , DIMENSION(A2D(0     )    )           ::   picefr     ! ice fraction                [0 to 1]
      !                                                        !!   optional arguments, used only in 'mixed oce-ice' case or for Met-Office coupling
      REAL(wp), INTENT(in)   , DIMENSION(A2D(0     ),jpl), OPTIONAL ::   palbi      ! all skies ice albedo
      REAL(wp), INTENT(in)   , DIMENSION(A2D(nn_hls)    ), OPTIONAL ::   psst       ! sea surface temperature     [Celsius]
      REAL(wp), INTENT(inout), DIMENSION(A2D(nn_hls),jpl), OPTIONAL ::   pist       ! ice surface temperature     [Kelvin] => inout for Met-Office
      REAL(wp), INTENT(in)   , DIMENSION(A2D(nn_hls),jpl), OPTIONAL ::   phs        ! snow depth                  [m]
      REAL(wp), INTENT(in)   , DIMENSION(A2D(nn_hls),jpl), OPTIONAL ::   phi        ! ice thickness               [m]
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop index
      REAL(wp), DIMENSION(A2D(0))     ::   zcptn, zcptrain, zcptsnw, ziceld, zmsk, zsnw
      REAL(wp), DIMENSION(A2D(0))     ::   zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip  , zevap_oce, zdevap_ice
      REAL(wp), DIMENSION(A2D(0))     ::   zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice
      REAL(wp), DIMENSION(A2D(0))     ::   zevap_ice_total
      REAL(wp), DIMENSION(A2D(0))     ::   ztri
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice, zevap_ice, zqtr_ice_top
      REAL(wp), DIMENSION(A2D(0),jpl) ::   ztsu
      REAL(wp) :: ztmp
      !!----------------------------------------------------------------------
      !
#if defined key_si3 || defined key_cice
      !
      IF( kt == nit000 ) THEN
         ! allocate ice fractions from last coupling time here and not in sbc_cpl_init because of jpl
         IF( .NOT.ALLOCATED(a_i_last_couple) )   ALLOCATE( a_i_last_couple(jpi,jpj,jpl) )
         ! initialize to a_i for the 1st time step
         a_i_last_couple(:,:,:) = a_i(:,:,:)
      ENDIF
      !
      IF( ln_mixcpl )   zmsk(A2D(0)) = 1. - xcplmask(A2D(0),1,0)
      ziceld(A2D(0)) = 1._wp - picefr(A2D(0))
      zcptn (A2D(0)) = rcp * sst_m(A2D(0))
      !
      !                                                      ! ========================= !
      !                                                      !    freshwater budget      !   (emp_tot)
      !                                                      ! ========================= !
      !
      !                                                           ! solid Precipitation                                (sprecip)
      !                                                           ! liquid + solid Precipitation                       (tprecip)
      !                                                           ! total Evaporation - total Precipitation            (emp_tot)
      !                                                           ! sublimation - solid precipitation (cell average)   (emp_ice)
      SELECT CASE( TRIM( sn_rcv_emp%cldes ) )
      CASE( 'conservative' )   ! received fields: jpr_rain, jpr_snow, jpr_ievp, jpr_tevp
         zsprecip(A2D(0)) = frcv(jpr_snow)%z3(A2D(0),1)                  ! May need to ensure positive here
         ztprecip(A2D(0)) = frcv(jpr_rain)%z3(A2D(0),1) + zsprecip(A2D(0))  ! May need to ensure positive here
         zemp_tot(A2D(0)) = frcv(jpr_tevp)%z3(A2D(0),1) - ztprecip(A2D(0))
      CASE( 'oce and ice'   )   ! received fields: jpr_sbpr, jpr_semp, jpr_oemp, jpr_ievp
         zemp_tot(A2D(0)) = ziceld(A2D(0)) * frcv(jpr_oemp)%z3(A2D(0),1) + picefr(A2D(0)) * frcv(jpr_sbpr)%z3(A2D(0),1)
         zemp_ice(A2D(0)) = frcv(jpr_semp)%z3(A2D(0),1) * picefr(A2D(0))
         zsprecip(A2D(0)) = frcv(jpr_ievp)%z3(A2D(0),1) - frcv(jpr_semp)%z3(A2D(0),1)
         ztprecip(A2D(0)) = frcv(jpr_semp)%z3(A2D(0),1) - frcv(jpr_sbpr)%z3(A2D(0),1) + zsprecip(A2D(0))
      CASE( 'none'      )       ! Not available as for now: needs additional coding below when computing zevap_oce
      !                         ! since fields received are not defined with none option
         CALL ctl_stop('STOP', 'sbccpl/sbc_cpl_ice_flx: some fields are not defined. Change sn_rcv_emp value in namelist namsbc_cpl')
      END SELECT

      ! --- evaporation over ice (kg/m2/s) --- !
      IF( ln_scale_ice_flux ) THEN ! typically met-office requirements
         IF( sn_rcv_emp%clcat == 'yes' ) THEN
            WHERE( a_i(A2D(0),:) > 1.e-10 ) ; zevap_ice(A2D(0),:) = frcv(jpr_ievp)%z3(A2D(0),:) * &
               &                                                 a_i_last_couple(A2D(0),:) / a_i(A2D(0),:)
            ELSEWHERE                       ; zevap_ice(A2D(0),:) = 0._wp
            END WHERE
            WHERE( picefr(A2D(0)) > 1.e-10 )   ; zevap_ice_total(A2D(0)) = SUM( zevap_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 ) / picefr(A2D(0))
            ELSEWHERE                       ; zevap_ice_total(A2D(0)) = 0._wp
            END WHERE
         ELSE
            WHERE( picefr(A2D(0)) > 1.e-10 )   ; zevap_ice(A2D(0),1) = frcv(jpr_ievp)%z3(A2D(0),1) * &
               &                                                 SUM( a_i_last_couple(A2D(0),:), dim=3 ) / picefr(A2D(0))
            ELSEWHERE                       ; zevap_ice(A2D(0),1) = 0._wp
            END WHERE
            zevap_ice_total(A2D(0)) = zevap_ice(A2D(0),1)
            DO jl = 2, jpl
               zevap_ice(A2D(0),jl) = zevap_ice(A2D(0),1)
            ENDDO
         ENDIF
      ELSE
         IF( sn_rcv_emp%clcat == 'yes' ) THEN
            zevap_ice(A2D(0),1:jpl) = frcv(jpr_ievp)%z3(A2D(0),1:jpl)
            WHERE( picefr(A2D(0)) > 1.e-10 ) ; zevap_ice_total(A2D(0)) = SUM( zevap_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 ) / picefr(A2D(0))
            ELSEWHERE                     ; zevap_ice_total(A2D(0)) = 0._wp
            END WHERE
         ELSE
            zevap_ice(A2D(0),1) = frcv(jpr_ievp)%z3(A2D(0),1)
            zevap_ice_total(A2D(0)) = zevap_ice(A2D(0),1)
            DO jl = 2, jpl
               zevap_ice(A2D(0),jl) = zevap_ice(A2D(0),1)
            ENDDO
         ENDIF
      ENDIF

      IF( TRIM( sn_rcv_emp%cldes ) == 'conservative' ) THEN
         ! For conservative case zemp_ice has not been defined yet. Do it now.
         zemp_ice(A2D(0)) = zevap_ice_total(A2D(0)) * picefr(A2D(0)) - frcv(jpr_snow)%z3(A2D(0),1) * picefr(A2D(0))
      ENDIF

      ! zsnw = snow fraction over ice after wind blowing (=picefr if no blowing)
      zsnw(A2D(0)) = 0._wp   ;   CALL ice_var_snwblow( ziceld, zsnw )

      ! --- evaporation minus precipitation corrected (because of wind blowing on snow) --- !
      zemp_ice(A2D(0)) = zemp_ice(A2D(0)) + zsprecip(A2D(0)) * ( picefr(A2D(0)) - zsnw(A2D(0)) )  ! emp_ice = A * sublimation - zsnw * sprecip
      zemp_oce(A2D(0)) = zemp_tot(A2D(0)) - zemp_ice(A2D(0))                                ! emp_oce = emp_tot - emp_ice

      ! --- evaporation over ocean (used later for qemp) --- !
      zevap_oce(A2D(0)) = frcv(jpr_tevp)%z3(A2D(0),1) - zevap_ice_total(A2D(0)) * picefr(A2D(0))

      ! since the sensitivity of evap to temperature (devap/dT) is not prescribed by the atmosphere, we set it to 0
      ! therefore, sublimation is not redistributed over the ice categories when no subgrid scale fluxes are provided by atm.
      zdevap_ice(A2D(0)) = 0._wp

      ! --- Continental fluxes --- !
      IF( srcv(nmodsbc)%fld(jpr_rnf)%laction ) THEN   ! runoffs (included in emp later on)
         rnf(A2D(0)) = frcv(jpr_rnf)%z3(A2D(0),1)
      ENDIF
      IF( srcv(nmodsbc)%fld(jpr_cal)%laction ) THEN   ! calving (put in emp_tot and emp_oce)
         zemp_tot(A2D(0)) = zemp_tot(A2D(0)) - frcv(jpr_cal)%z3(A2D(0),1)
         zemp_oce(A2D(0)) = zemp_oce(A2D(0)) - frcv(jpr_cal)%z3(A2D(0),1)
      ENDIF
      IF( srcv(nmodsbc)%fld(jpr_icb)%laction ) THEN   ! iceberg added to runoffs
         fwficb(A2D(0)) = frcv(jpr_icb)%z3(A2D(0),1)
         rnf(A2D(0)) = rnf(A2D(0)) + fwficb(A2D(0))
      ENDIF
      IF( srcv(nmodsbc)%fld(jpr_isf)%laction ) THEN   ! iceshelf (fwfisf > 0 mean melting)
        fwfisf_oasis(A2D(0)) = frcv(jpr_isf)%z3(A2D(0),1)
      ENDIF

      IF( ln_mixcpl ) THEN
         emp_tot(A2D(0)) = emp_tot(A2D(0)) * xcplmask(A2D(0),1,0) + zemp_tot(A2D(0)) * zmsk(A2D(0))
         emp_ice(A2D(0)) = emp_ice(A2D(0)) * xcplmask(A2D(0),1,0) + zemp_ice(A2D(0)) * zmsk(A2D(0))
         emp_oce(A2D(0)) = emp_oce(A2D(0)) * xcplmask(A2D(0),1,0) + zemp_oce(A2D(0)) * zmsk(A2D(0))
         sprecip(A2D(0)) = sprecip(A2D(0)) * xcplmask(A2D(0),1,0) + zsprecip(A2D(0)) * zmsk(A2D(0))
         tprecip(A2D(0)) = tprecip(A2D(0)) * xcplmask(A2D(0),1,0) + ztprecip(A2D(0)) * zmsk(A2D(0))
         DO jl = 1, jpl
            evap_ice (A2D(0),jl) = evap_ice (A2D(0),jl) * xcplmask(A2D(0),1,0) + zevap_ice (A2D(0),jl) * zmsk(A2D(0))
            devap_ice(A2D(0),jl) = devap_ice(A2D(0),jl) * xcplmask(A2D(0),1,0) + zdevap_ice(A2D(0))    * zmsk(A2D(0))
         END DO
      ELSE
         emp_tot (A2D(0))   = zemp_tot (A2D(0))
         emp_ice (A2D(0))   = zemp_ice (A2D(0))
         emp_oce (A2D(0))   = zemp_oce (A2D(0))
         sprecip (A2D(0))   = zsprecip (A2D(0))
         tprecip (A2D(0))   = ztprecip (A2D(0))
         evap_ice(A2D(0),:) = zevap_ice(A2D(0),:)
         DO jl = 1, jpl
            devap_ice(A2D(0),jl) = zdevap_ice(A2D(0))
         END DO
      ENDIF

!! for CICE ??
!!$      zsnw(A2D(0)) = picefr(A2D(0))
!!$      ! --- Continental fluxes --- !
!!$      IF( srcv(nmodsbc)%fld(jpr_rnf)%laction ) THEN   ! runoffs (included in emp later on)
!!$         rnf(A2D(0)) = frcv(jpr_rnf)%z3(A2D(0),1)
!!$      ENDIF
!!$      IF( srcv(nmodsbc)%fld(jpr_cal)%laction ) THEN   ! calving (put in emp_tot)
!!$         zemp_tot(A2D(0)) = zemp_tot(A2D(0)) - frcv(jpr_cal)%z3(A2D(0),1)
!!$      ENDIF
!!$      IF( srcv(nmodsbc)%fld(jpr_icb)%laction ) THEN   ! iceberg added to runoffs
!!$         fwficb(A2D(0)) = frcv(jpr_icb)%z3(A2D(0),1)
!!$         rnf(A2D(0))    = rnf(A2D(0)) + fwficb(A2D(0))
!!$      ENDIF
!!$      IF( srcv(nmodsbc)%fld(jpr_isf)%laction ) THEN   ! iceshelf (fwfisf >0 mean melting)
!!$        fwfisf_oasis(A2D(0)) = frcv(jpr_isf)%z3(A2D(0),1)
!!$      ENDIF
!!$      !
!!$      IF( ln_mixcpl ) THEN
!!$         emp_tot(A2D(0)) = emp_tot(A2D(0)) * xcplmask(A2D(0),0) + zemp_tot(A2D(0)) * zmsk(A2D(0))
!!$         emp_ice(A2D(0)) = emp_ice(A2D(0)) * xcplmask(A2D(0),0) + zemp_ice(A2D(0)) * zmsk(A2D(0))
!!$         sprecip(A2D(0)) = sprecip(A2D(0)) * xcplmask(A2D(0),0) + zsprecip(A2D(0)) * zmsk(A2D(0))
!!$         tprecip(A2D(0)) = tprecip(A2D(0)) * xcplmask(A2D(0),0) + ztprecip(A2D(0)) * zmsk(A2D(0))
!!$      ELSE
!!$         emp_tot(A2D(0)) =                                  zemp_tot(A2D(0))
!!$         emp_ice(A2D(0)) =                                  zemp_ice(A2D(0))
!!$         sprecip(A2D(0)) =                                  zsprecip(A2D(0))
!!$         tprecip(A2D(0)) =                                  ztprecip(A2D(0))
!!$      ENDIF
      !
      ! outputs
      IF( srcv(nmodsbc)%fld(jpr_cal)%laction )    CALL iom_put( 'calving_cea' , frcv(jpr_cal)%z3(A2D(0),1) * tmask(A2D(0),1)             )  ! calving
      IF( srcv(nmodsbc)%fld(jpr_icb)%laction )    CALL iom_put( 'iceberg_cea' , frcv(jpr_icb)%z3(A2D(0),1) * tmask(A2D(0),1)             )  ! icebergs
      IF( iom_use('snowpre') )       CALL iom_put( 'snowpre'     , sprecip(A2D(0))                                          )  ! Snow
      IF( iom_use('precip') )        CALL iom_put( 'precip'      , tprecip(A2D(0))                                          )  ! total  precipitation
      IF( iom_use('rain') )          CALL iom_put( 'rain'        , tprecip(A2D(0)) - sprecip(A2D(0))                           )  ! liquid precipitation
      IF( iom_use('snow_ao_cea') )   CALL iom_put( 'snow_ao_cea' , sprecip(A2D(0)) * ( 1._wp - zsnw(A2D(0)) )                  )  ! Snow over ice-free ocean  (cell average)
      IF( iom_use('snow_ai_cea') )   CALL iom_put( 'snow_ai_cea' , sprecip(A2D(0)) *           zsnw(A2D(0))                    )  ! Snow over sea-ice         (cell average)
      IF( iom_use('rain_ao_cea') )   CALL iom_put( 'rain_ao_cea' , ( tprecip(A2D(0)) - sprecip(A2D(0)) ) * ziceld(A2D(0))         )  ! liquid precipitation over ocean (cell average)
      IF( iom_use('subl_ai_cea') )   CALL iom_put( 'subl_ai_cea' , zevap_ice_total(A2D(0)) * picefr(A2D(0)) * smask0(A2D(0))      )  ! Sublimation over sea-ice (cell average)
      IF( iom_use('evap_ao_cea') )   CALL iom_put( 'evap_ao_cea' , ( frcv(jpr_tevp)%z3(A2D(0),1)  &
         &                                                         - zevap_ice_total(A2D(0)) * picefr(A2D(0)) ) * smask0(A2D(0))  )  ! ice-free oce evap (cell average)
      ! note: runoff output is done in sbcrnf (which includes icebergs too) and iceshelf output is done in sbcisf
!!      IF( srcv(nmodsbc)%fld(jpr_rnf)%laction )   CALL iom_put( 'runoffs' , rnf(A2D(0)) * tmask(A2D(0),1)                                 )  ! runoff
      IF( srcv(nmodsbc)%fld(jpr_isf)%laction )    CALL iom_put( 'iceshelf_cea', frcv(jpr_isf)%z3(A2D(0),1) * smask0(A2D(0))                 )  ! iceshelf
      !
      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_iceflx%cldes ) )             !  ice topmelt and botmelt  !
      !                                                      ! ========================= !
      CASE( 'coupled' )
         IF( ln_scale_ice_flux ) THEN
            WHERE( a_i(A2D(0),:) > 1.e-10_wp )
               qml_ice(A2D(0),:) = frcv(jpr_topm)%z3(A2D(0),:) * a_i_last_couple(A2D(0),:) / a_i(A2D(0),:)
               qcn_ice(A2D(0),:) = frcv(jpr_botm)%z3(A2D(0),:) * a_i_last_couple(A2D(0),:) / a_i(A2D(0),:)
            ELSEWHERE
               qml_ice(A2D(0),:) = 0.0_wp
               qcn_ice(A2D(0),:) = 0.0_wp
            END WHERE
         ELSE
            qml_ice(A2D(0),:) = frcv(jpr_topm)%z3(A2D(0),:)
            qcn_ice(A2D(0),:) = frcv(jpr_botm)%z3(A2D(0),:)
         ENDIF
      END SELECT
      !
      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_qns%cldes ) )                !   non solar heat fluxes   !   (qns)
      !                                                      ! ========================= !
      CASE( 'oce only' )         ! the required field is directly provided
         ! Get the sea ice non solar heat flux from conductive, melting and sublimation fluxes
         IF( TRIM(sn_rcv_iceflx%cldes) == 'coupled' ) THEN
            zqns_ice(A2D(0),:) = qml_ice(A2D(0),:) + qcn_ice(A2D(0),:)
         ELSE
            zqns_ice(A2D(0),:) = 0._wp
         ENDIF
         ! Calculate the total non solar heat flux. The ocean only non solar heat flux (zqns_oce) will be recalculated after this CASE
         ! statement to be consistent with other coupling methods even though .zqns_oce = frcv(jpr_qnsoce)%z3(A2D(0),1)
         zqns_tot(A2D(0)) = frcv(jpr_qnsoce)%z3(A2D(0),1) + SUM( zqns_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 )
      CASE( 'conservative' )     ! the required fields are directly provided
         zqns_tot(A2D(0)) = frcv(jpr_qnsmix)%z3(A2D(0),1)
         IF( TRIM(sn_rcv_qns%clcat) == 'yes' ) THEN
            zqns_ice(A2D(0),1:jpl) = frcv(jpr_qnsice)%z3(A2D(0),1:jpl)
         ELSE
            DO jl = 1, jpl
               zqns_ice(A2D(0),jl) = frcv(jpr_qnsice)%z3(A2D(0),1) ! Set all category values equal
            END DO
         ENDIF
      CASE( 'oce and ice' )      ! the total flux is computed from ocean and ice fluxes
         zqns_tot(A2D(0)) =  ziceld(A2D(0)) * frcv(jpr_qnsoce)%z3(A2D(0),1)
         IF( TRIM(sn_rcv_qns%clcat) == 'yes' ) THEN
            DO jl=1,jpl
               zqns_tot(A2D(0)   ) = zqns_tot(A2D(0)) + a_i(A2D(0),jl) * frcv(jpr_qnsice)%z3(A2D(0),jl)
               zqns_ice(A2D(0),jl) = frcv(jpr_qnsice)%z3(A2D(0),jl)
            ENDDO
         ELSE
            zqns_tot(A2D(0)) = zqns_tot(A2D(0)) + picefr(A2D(0)) * frcv(jpr_qnsice)%z3(A2D(0),1)
            DO jl = 1, jpl
               zqns_ice(A2D(0),jl) = frcv(jpr_qnsice)%z3(A2D(0),1)
            END DO
         ENDIF
      CASE( 'mixed oce-ice' )    ! the ice flux is cumputed from the total flux, the SST and ice informations
! ** NEED TO SORT OUT HOW THIS SHOULD WORK IN THE MULTI-CATEGORY CASE - CURRENTLY NOT ALLOWED WHEN INTERFACE INITIALISED **
         zqns_tot(A2D(0)  ) = frcv(jpr_qnsmix)%z3(A2D(0),1)
         IF ( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            DO jl = 1, jpl
               zqns_ice(A2D(0),jl) = frcv(jpr_qnsmix)%z3(A2D(0),jl)    &
                  &             + frcv(jpr_dqnsdt)%z3(A2D(0),jl) * ( pist(A2D(0),jl) - ( ( rt0 + psst(A2D(0)) ) * ziceld(A2D(0))   &
                  &                                             + pist(A2D(0),jl) * picefr(A2D(0)) ) )
            END DO
         ELSE
            DO jl = 1, jpl
               zqns_ice(A2D(0),jl) = frcv(jpr_qnsmix)%z3(A2D(0), 1)    &
                  &             + frcv(jpr_dqnsdt)%z3(A2D(0), 1) * ( pist(A2D(0),jl) - ( ( rt0 + psst(A2D(0)) ) * ziceld(A2D(0))   &
                  &                                             + pist(A2D(0),jl) * picefr(A2D(0)) ) )
            END DO
         ENDIF
      END SELECT
      !
      ! --- calving (removed from qns_tot) --- !
      IF( srcv(nmodsbc)%fld(jpr_cal)%laction )   zqns_tot(A2D(0)) = zqns_tot(A2D(0)) - frcv(jpr_cal)%z3(A2D(0),1) * rLfus  ! remove latent heat of calving
                                                                                                        ! we suppose it melts at 0deg, though it should be temp. of surrounding ocean
      ! --- iceberg (removed from qns_tot) --- !
      IF( srcv(nmodsbc)%fld(jpr_icb)%laction )   zqns_tot(A2D(0)) = zqns_tot(A2D(0)) - frcv(jpr_icb)%z3(A2D(0),1) * rLfus  ! remove latent heat of iceberg melting

      ! --- non solar flux over ocean --- !
      !         note: ziceld cannot be = 0 since we limit the ice concentration to amax
      zqns_oce = 0._wp
      WHERE( ziceld /= 0._wp )   zqns_oce(A2D(0)) = ( zqns_tot(A2D(0)) - SUM( a_i(A2D(0),:) * zqns_ice(A2D(0),:), dim=3 ) ) / ziceld(A2D(0))

      ! Heat content per unit mass of snow (J/kg)
      WHERE( SUM( a_i(A2D(0),:), dim=3 ) > 1.e-10 ) ; zcptsnw(A2D(0)) = rcpi * SUM( (tn_ice(A2D(0),:) - rt0) * a_i(A2D(0),:), dim=3 ) &
         &                                                                / SUM( a_i(A2D(0),:), dim=3 )
      ELSEWHERE                                     ; zcptsnw(A2D(0)) = zcptn(A2D(0))
      ENDWHERE
      ! Heat content per unit mass of rain (J/kg)
      zcptrain(A2D(0)) = rcp * ( SUM( (tn_ice(A2D(0),:) - rt0) * a_i(A2D(0),:), dim=3 ) + sst_m(A2D(0)) * ziceld(A2D(0)) )

      ! --- enthalpy of snow precip over ice in J/m3 (to be used in 1D-thermo) --- !
      zqprec_ice(A2D(0)) = rhos * ( zcptsnw(A2D(0)) - rLfus )

      ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) --- !
      DO jl = 1, jpl
         zqevap_ice(A2D(0),jl) = 0._wp ! should be -evap * ( ( Tice - rt0 ) * rcpi ) but atm. does not take it into account
      END DO

      ! --- heat flux associated with emp (W/m2) --- !
      zqemp_oce(A2D(0)) = -  zevap_oce(A2D(0))                                      *   zcptn   (A2D(0))   &        ! evap
         &             + ( ztprecip(A2D(0)) - zsprecip(A2D(0)) )                    *   zcptrain(A2D(0))   &        ! liquid precip
         &             +   zsprecip(A2D(0))                   * ( 1._wp - zsnw ) * ( zcptsnw (A2D(0)) - rLfus )  ! solid precip over ocean + snow melting
      zqemp_ice(A2D(0)) =     zsprecip(A2D(0))                   * zsnw             * ( zcptsnw (A2D(0)) - rLfus )  ! solid precip over ice (qevap_ice=0 since atm. does not take it into account)
!!    zqemp_ice(A2D(0)) = -   frcv(jpr_ievp)%z3(A2D(0),1)        * picefr(A2D(0))      *   zcptsnw (A2D(0))   &        ! ice evap
!!       &             +   zsprecip(A2D(0))                   * zsnw             * zqprec_ice(A2D(0)) * r1_rhos  ! solid precip over ice

      ! --- total non solar flux (including evap/precip) --- !
      zqns_tot(A2D(0)) = zqns_tot(A2D(0)) + zqemp_ice(A2D(0)) + zqemp_oce(A2D(0))

      ! --- in case both coupled/forced are active, we must mix values --- !
      IF( ln_mixcpl ) THEN
         qns_tot(A2D(0)) = qns_tot(A2D(0)) * xcplmask(A2D(0),1,0) + zqns_tot(A2D(0))* zmsk(A2D(0))
         qns_oce(A2D(0)) = qns_oce(A2D(0)) * xcplmask(A2D(0),1,0) + zqns_oce(A2D(0))* zmsk(A2D(0))
         DO jl=1,jpl
            qns_ice  (A2D(0),jl) = qns_ice  (A2D(0),jl) * xcplmask(A2D(0),1,0) +  zqns_ice  (A2D(0),jl)* zmsk(A2D(0))
            qevap_ice(A2D(0),jl) = qevap_ice(A2D(0),jl) * xcplmask(A2D(0),1,0) +  zqevap_ice(A2D(0),jl)* zmsk(A2D(0))
         ENDDO
         qprec_ice(A2D(0)) = qprec_ice(A2D(0)) * xcplmask(A2D(0),1,0) + zqprec_ice(A2D(0))* zmsk(A2D(0))
         qemp_oce (A2D(0)) =  qemp_oce(A2D(0)) * xcplmask(A2D(0),1,0) +  zqemp_oce(A2D(0))* zmsk(A2D(0))
         qemp_ice (A2D(0)) =  qemp_ice(A2D(0)) * xcplmask(A2D(0),1,0) +  zqemp_ice(A2D(0))* zmsk(A2D(0))
      ELSE
         qns_tot  (A2D(0)  ) = zqns_tot  (A2D(0)  )
         qns_oce  (A2D(0)  ) = zqns_oce  (A2D(0)  )
         qns_ice  (A2D(0),:) = zqns_ice  (A2D(0),:)
         qevap_ice(A2D(0),:) = zqevap_ice(A2D(0),:)
         qprec_ice(A2D(0)  ) = zqprec_ice(A2D(0)  )
         qemp_oce (A2D(0)  ) = zqemp_oce (A2D(0)  )
         qemp_ice (A2D(0)  ) = zqemp_ice (A2D(0)  )
      ENDIF

!! for CICE ??
!!$      ! --- non solar flux over ocean --- !
!!$      zcptsnw (A2D(0)) = zcptn(A2D(0))
!!$      zcptrain(A2D(0)) = zcptn(A2D(0))
!!$
!!$      ! clem: this formulation is certainly wrong... but better than it was...
!!$      zqns_tot(A2D(0)) = zqns_tot(A2D(0))                             &          ! zqns_tot update over free ocean with:
!!$         &          - (  ziceld(A2D(0)) * zsprecip(A2D(0)) * rLfus )  &          ! remove the latent heat flux of solid precip. melting
!!$         &          - (  zemp_tot(A2D(0))                          &          ! remove the heat content of mass flux (assumed to be at SST)
!!$         &             - zemp_ice(A2D(0)) ) * zcptn(A2D(0))
!!$
!!$     IF( ln_mixcpl ) THEN
!!$         qns_tot(A2D(0)) = qns(A2D(0)) * ziceld(A2D(0)) + SUM( qns_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 )   ! total flux from blk
!!$         qns_tot(A2D(0)) = qns_tot(A2D(0)) * xcplmask(A2D(0),0) +  zqns_tot(A2D(0))* zmsk(A2D(0))
!!$         DO jl=1,jpl
!!$            qns_ice(A2D(0),jl) = qns_ice(A2D(0),jl) * xcplmask(A2D(0),0) +  zqns_ice(A2D(0),jl)* zmsk(A2D(0))
!!$         ENDDO
!!$      ELSE
!!$         qns_tot(A2D(0)  ) = zqns_tot(A2D(0)  )
!!$         qns_ice(A2D(0),:) = zqns_ice(A2D(0),:)
!!$      ENDIF

      ! outputs
      IF ( srcv(nmodsbc)%fld(jpr_cal)%laction ) CALL iom_put('hflx_cal_cea' , - frcv(jpr_cal)%z3(A2D(0),1) * rLfus ) ! latent heat from calving
      IF ( srcv(nmodsbc)%fld(jpr_icb)%laction ) CALL iom_put('hflx_icb_cea' , - frcv(jpr_icb)%z3(A2D(0),1) * rLfus ) ! latent heat from icebergs melting
      IF (        iom_use('hflx_rain_cea') )    &                                                    ! heat flux from rain (cell average)
         &   CALL iom_put('hflx_rain_cea' , ( tprecip(A2D(0)) - sprecip(A2D(0)) ) * zcptrain(A2D(0)) )
      IF (        iom_use('hflx_evap_cea') )    &                                                    ! heat flux from evap (cell average)
         &   CALL iom_put('hflx_evap_cea' , ( frcv(jpr_tevp)%z3(A2D(0),1) - zevap_ice_total(A2D(0)) * picefr(A2D(0)) )  &
         &                                  * zcptn(A2D(0)) * smask0(A2D(0)) )
      IF (        iom_use('hflx_prec_cea') )    &                                                    ! heat flux from all precip (cell avg)
         &   CALL iom_put('hflx_prec_cea' ,    sprecip(A2D(0)) * ( zcptsnw(A2D(0)) - rLfus )  &
         &                                 + ( tprecip(A2D(0)) - sprecip(A2D(0)) ) * zcptrain(A2D(0)) )
      IF (        iom_use('hflx_snow_cea') )    &                                                    ! heat flux from snow (cell average)
         &   CALL iom_put('hflx_snow_cea'   , sprecip(A2D(0)) * ( zcptsnw(A2D(0)) - rLfus )  )
      IF (        iom_use('hflx_snow_ao_cea') ) &                                                    ! heat flux from snow (over ocean)
         &   CALL iom_put('hflx_snow_ao_cea', sprecip(A2D(0)) * ( zcptsnw(A2D(0)) - rLfus ) * ( 1._wp - zsnw(A2D(0)) ) )
      IF (        iom_use('hflx_snow_ai_cea') ) &                                                    ! heat flux from snow (over ice)
         &   CALL iom_put('hflx_snow_ai_cea', sprecip(A2D(0)) * ( zcptsnw(A2D(0)) - rLfus ) *  zsnw(A2D(0)) )
      IF(         iom_use('hflx_subl_cea') )    &                                                    ! heat flux from sublimation
         &   CALL iom_put('hflx_subl_cea' ,   SUM( qevap_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 ) * smask0(A2D(0)) )
      ! note: hflx for runoff and iceshelf are done in sbcrnf and sbcisf resp.
      !
      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_dqnsdt%cldes ) )             !          d(qns)/dt        !
      !                                                      ! ========================= !
      CASE ('coupled')
         IF( TRIM(sn_rcv_dqnsdt%clcat) == 'yes' ) THEN
            zdqns_ice(A2D(0),1:jpl) = frcv(jpr_dqnsdt)%z3(A2D(0),1:jpl)
         ELSE
            ! Set all category values equal for the moment
            DO jl=1,jpl
               zdqns_ice(A2D(0),jl) = frcv(jpr_dqnsdt)%z3(A2D(0),1)
            ENDDO
         ENDIF
      CASE( 'none' )
         zdqns_ice(A2D(0),:) = 0._wp
      END SELECT

      IF( ln_mixcpl ) THEN
         DO jl=1,jpl
            dqns_ice(A2D(0),jl) = dqns_ice(A2D(0),jl) * xcplmask(A2D(0),1,0) + zdqns_ice(A2D(0),jl) * zmsk(A2D(0))
         ENDDO
      ELSE
         dqns_ice(A2D(0),:) = zdqns_ice(A2D(0),:)
      ENDIF
      !                                                      ! ========================= !
      SELECT CASE( TRIM( sn_rcv_qsr%cldes ) )                !      solar heat fluxes    !   (qsr)
      !                                                      ! ========================= !
      CASE( 'oce only' )
         zqsr_tot(A2D(0)  ) = MAX( 0._wp , frcv(jpr_qsroce)%z3(A2D(0),1) )
         ! For the Met Office the only sea ice solar flux is the transmitted qsr which is added onto zqsr_ice
         ! further down. Therefore start zqsr_ice off at zero.
         zqsr_ice(A2D(0),:) = 0._wp
      CASE( 'conservative' )
         zqsr_tot(A2D(0)  ) = frcv(jpr_qsrmix)%z3(A2D(0),1)
         IF( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            zqsr_ice(A2D(0),1:jpl) = frcv(jpr_qsrice)%z3(A2D(0),1:jpl)
         ELSE
            ! Set all category values equal for the moment
            DO jl = 1, jpl
               zqsr_ice(A2D(0),jl) = frcv(jpr_qsrice)%z3(A2D(0),1)
            END DO
         ENDIF
      CASE( 'oce and ice' )
         zqsr_tot(A2D(0)  ) =  ziceld(A2D(0)) * frcv(jpr_qsroce)%z3(A2D(0),1)
         IF( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            DO jl = 1, jpl
               zqsr_tot(A2D(0)   ) = zqsr_tot(A2D(0)) + a_i(A2D(0),jl) * frcv(jpr_qsrice)%z3(A2D(0),jl)
               zqsr_ice(A2D(0),jl) = frcv(jpr_qsrice)%z3(A2D(0),jl)
            END DO
         ELSE
            zqsr_tot(A2D(0)) = zqsr_tot(A2D(0)) + picefr(A2D(0)) * frcv(jpr_qsrice)%z3(A2D(0),1)
            DO jl = 1, jpl
               zqsr_ice(A2D(0),jl) = frcv(jpr_qsrice)%z3(A2D(0),1)
            END DO
         ENDIF
      CASE( 'mixed oce-ice' )
         zqsr_tot(A2D(0)  ) = frcv(jpr_qsrmix)%z3(A2D(0),1)
! ** NEED TO SORT OUT HOW THIS SHOULD WORK IN THE MULTI-CATEGORY CASE - CURRENTLY NOT ALLOWED WHEN INTERFACE INITIALISED **
!       Create solar heat flux over ice using incoming solar heat flux and albedos
!       ( see OASIS3 user guide, 5th edition, p39 )
         IF ( TRIM(sn_rcv_qsr%clcat) == 'yes' ) THEN
            DO jl = 1, jpl
               zqsr_ice(A2D(0),jl) = frcv(jpr_qsrmix)%z3(A2D(0),jl) * ( 1.- palbi(A2D(0),jl) )   &
                  &               / (  1.- ( alb_oce_mix(A2D(0)   ) * ziceld(A2D(0))       &
                  &                        + palbi      (A2D(0),jl) * picefr(A2D(0)) ) )
            END DO
         ELSE
            DO jl = 1, jpl
               zqsr_ice(A2D(0),jl) = frcv(jpr_qsrmix)%z3(A2D(0), 1) * ( 1.- palbi(A2D(0),jl) )   &
                  &               / (  1.- ( alb_oce_mix(A2D(0)   ) * ziceld(A2D(0))       &
                  &                        + palbi      (A2D(0),jl) * picefr(A2D(0)) ) )
            END DO
         ENDIF
      CASE( 'none'      )       ! Not available as for now: needs additional coding
      !                         ! since fields received, here zqsr_tot,  are not defined with none option
         CALL ctl_stop('STOP', 'sbccpl/sbc_cpl_ice_flx: some fields are not defined. Change sn_rcv_qsr value in namelist namsbc_cpl')
      END SELECT
      IF( ln_dm2dc .AND. ln_cpl ) THEN   ! modify qsr to include the diurnal cycle
         zqsr_tot(A2D(0)) = sbc_dcy( zqsr_tot(A2D(0)) )
         DO jl = 1, jpl
            zqsr_ice(A2D(0),jl) = sbc_dcy( zqsr_ice(A2D(0),jl) )
         END DO
      ENDIF
      !                                                      ! ========================= !
      !                                                      !      Transmitted Qsr      !   [W/m2]
      !                                                      ! ========================= !
      IF( .NOT.ln_cndflx ) THEN                              !==  No conduction flux as surface forcing  ==!
         !
         IF( nn_qtrice == 0 ) THEN
            ! formulation derived from Grenfell and Maykut (1977), where transmission rate
            !    1) depends on cloudiness
            !       ! ===> used prescribed cloud fraction representative for polar oceans in summer (0.81)
            !       !      should be real cloud fraction instead (as in the bulk) but needs to be read from atm.
            !    2) is 0 when there is any snow
            !    3) tends to 1 for thin ice
            ztri(A2D(0)) = 0.18 * ( 1.0 - cloud_fra(A2D(0)) ) + 0.35 * cloud_fra(A2D(0))  ! surface transmission when hi>10cm
            DO jl = 1, jpl
               DO_2D( 0, 0, 0, 0 )
                  IF( phs(ji,jj,jl) <= 0._wp ) THEN
                     ztmp = ztri(ji,jj)                                       ! constant (ztri) when hi>10cm
                     IF( phi(ji,jj,jl) <  0.1_wp )   &                        ! linear decrease from hi=0 to 10cm
                        & ztmp = ztmp + ( 1._wp - ztmp ) * ( 1._wp - phi(ji,jj,jl) * 10._wp )
                     zqtr_ice_top(ji,jj,jl) = zqsr_ice(ji,jj,jl) * ztmp
                  ELSE                                                        ! zero when hs>0
                     zqtr_ice_top(ji,jj,jl) = 0._wp
                  ENDIF
               END_2D
            ENDDO
         ELSEIF( nn_qtrice == 1 ) THEN
            ! formulation is derived from the thesis of M. Lebrun (2019).
            !    It represents the best fit using several sets of observations
            !    It comes with snow conductivities adapted to freezing/melting conditions (see icethd_zdf_bl99.F90)
            zqtr_ice_top(A2D(0),:) = 0.3_wp * zqsr_ice(A2D(0),:)
         ENDIF
         !
      ELSEIF( ln_cndflx .AND. .NOT.ln_cndemulate ) THEN      !==  conduction flux as surface forcing  ==!
         !
         !
         SELECT CASE( TRIM( sn_rcv_qtrice%cldes ) )
            !
            !      ! ===> here we receive the qtr_ice_top array from the coupler
         CASE ('coupled')
            IF (ln_scale_ice_flux) THEN
               WHERE( a_i(A2D(0),:) > 1.e-10_wp )
                  zqtr_ice_top(A2D(0),:) = frcv(jpr_qtrice)%z3(A2D(0),:) * a_i_last_couple(A2D(0),:) / a_i(A2D(0),:)
               ELSEWHERE
                  zqtr_ice_top(A2D(0),:) = 0.0_wp
               ENDWHERE
            ELSE
               zqtr_ice_top(A2D(0),:) = frcv(jpr_qtrice)%z3(A2D(0),:)
            ENDIF
           
            ! Add retrieved transmitted solar radiation onto the ice and total solar radiation
            zqsr_ice(A2D(0),:) = zqsr_ice(A2D(0),:) + zqtr_ice_top(A2D(0),:)
            zqsr_tot(A2D(0))   = zqsr_tot(A2D(0)) + SUM( zqtr_ice_top(A2D(0),:) * a_i(A2D(0),:), dim=3 )
            
            !      if we are not getting this data from the coupler then assume zero (fully opaque ice)
         CASE ('none')
            zqtr_ice_top(A2D(0),:) = 0._wp
         END SELECT
            !
      ENDIF

      IF( ln_mixcpl ) THEN
         qsr_tot(A2D(0)) = qsr(A2D(0)) * ziceld(A2D(0)) + SUM( qsr_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 )   ! total flux from blk
         qsr_tot(A2D(0)) = qsr_tot(A2D(0)) * xcplmask(A2D(0),1,0) + zqsr_tot(A2D(0)) * zmsk(A2D(0))
         DO jl = 1, jpl
            qsr_ice    (A2D(0),jl) = qsr_ice    (A2D(0),jl) * xcplmask(A2D(0),1,0) + zqsr_ice    (A2D(0),jl) * zmsk(A2D(0))
            qtr_ice_top(A2D(0),jl) = qtr_ice_top(A2D(0),jl) * xcplmask(A2D(0),1,0) + zqtr_ice_top(A2D(0),jl) * zmsk(A2D(0))
         END DO
      ELSE
         qsr_tot    (A2D(0)  ) = zqsr_tot    (A2D(0)  )
         qsr_ice    (A2D(0),:) = zqsr_ice    (A2D(0),:)
         qtr_ice_top(A2D(0),:) = zqtr_ice_top(A2D(0),:)
      ENDIF
      
      ! --- solar flux over ocean --- !
      !         note: ziceld cannot be = 0 since we limit the ice concentration to amax
      zqsr_oce = 0._wp
      WHERE( ziceld /= 0._wp )  zqsr_oce(A2D(0)) = ( zqsr_tot(A2D(0)) - SUM( a_i(A2D(0),:) * zqsr_ice(A2D(0),:), dim=3 ) ) / ziceld(A2D(0))

      IF( ln_mixcpl ) THEN   ;   qsr_oce(A2D(0)) = qsr_oce(A2D(0)) * xcplmask(A2D(0),1,0) +  zqsr_oce(A2D(0))* zmsk(A2D(0))
      ELSE                   ;   qsr_oce(A2D(0)) = zqsr_oce(A2D(0))   ;   ENDIF

      !                                                      ! ================== !
      !                                                      !   ice skin temp.   !
      !                                                      ! ================== !
      ! needed by Met Office
      IF( srcv(nmodsbc)%fld(jpr_ts_ice)%laction ) THEN
         WHERE    ( frcv(jpr_ts_ice)%z3(A2D(0),:) > 0.0  )   ;   ztsu(A2D(0),:) =   0. + rt0
         ELSEWHERE( frcv(jpr_ts_ice)%z3(A2D(0),:) < -60. )   ;   ztsu(A2D(0),:) = -60. + rt0
         ELSEWHERE                                        ;   ztsu(A2D(0),:) = frcv(jpr_ts_ice)%z3(A2D(0),:) + rt0
         END WHERE
         !
         IF( ln_mixcpl ) THEN
            DO jl=1,jpl
               pist(A2D(0),jl) = pist(A2D(0),jl) * xcplmask(A2D(0),1,0) + ztsu(A2D(0),jl) * zmsk(A2D(0))
            ENDDO
         ELSE
            pist(A2D(0),:) = ztsu(A2D(0),:)
         ENDIF
         !
      ENDIF
      !
#endif
      !
   END SUBROUTINE sbc_cpl_ice_flx


   SUBROUTINE sbc_cpl_snd( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_snd  ***
      !!
      !! ** Purpose :   provide the ocean-ice informations to the atmosphere
      !!
      !! ** Method  :   send to the atmosphere through a call to cpl_snd
      !!              all the needed fields (as defined in sbc_cpl_init)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER, INTENT(in) ::   Kbb, Kmm    ! ocean model time level index
      !
      INTEGER ::   ji, jj, jl   ! dummy loop indices
      INTEGER ::   isec, info   ! local integer
      REAL(wp) ::   zumax, zvmax
      REAL(wp), DIMENSION(A2D(0))     ::   zat_i, zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1
      REAL(wp), DIMENSION(A2D(0),1)   ::   ze3t_i
      REAL(wp), DIMENSION(A2D(0),jpl) ::   ztmp3, ztmp4
      !!----------------------------------------------------------------------
      !
      isec = ( kt - nit000 ) * NINT( rn_Dt )        ! date of exchanges
      info = OASIS_idle

      zfr_l(A2D(0)) = 1.- fr_i(A2D(0))
      zat_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )
      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !   in Kelvin
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_toce)%laction .OR. ssnd(nmodsbc)%fld(jps_tice)%laction .OR. ssnd(nmodsbc)%fld(jps_tmix)%laction ) THEN

         IF( nn_components == jp_iam_oce ) THEN
            ztmp1(A2D(0)) = ts(A2D(0),1,jp_tem,Kmm)   ! send temperature as it is (potential or conservative) -> use of l_useCT on the received part
         ELSE
            ! we must send the surface potential temperature
            IF( l_useCT )  THEN ; CALL eos_pt_from_ct( ts(A2D(0),1,jp_tem,Kmm), ts(A2D(0),1,jp_sal,Kmm), ztmp1(A2D(0)), kbnd=0 )
            ELSE                ; ztmp1(A2D(0)) = ts(A2D(0),1,jp_tem,Kmm)
            ENDIF
            !
            SELECT CASE( sn_snd_temp%cldes)
            CASE( 'oce only'             )   ;   ztmp1(A2D(0)) =   ztmp1(A2D(0)) + rt0
            CASE( 'oce and ice'          )   ;   ztmp1(A2D(0)) =   ztmp1(A2D(0)) + rt0
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(A2D(0),1:jpl) = tn_ice(A2D(0),1:jpl)
               CASE( 'no' )
                  WHERE( zat_i(A2D(0)) /= 0. )
                     ztmp3(A2D(0),1) = SUM( tn_ice(A2D(0),:) * a_i(A2D(0),:), dim=3 ) / zat_i(A2D(0))
                  ELSEWHERE
                     ztmp3(A2D(0),1) = rt0
                  END WHERE
               CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'weighted oce and ice' )   ;   ztmp1(A2D(0)) = ( ztmp1(A2D(0)) + rt0 ) * zfr_l(A2D(0))
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(A2D(0),1:jpl) = tn_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
               CASE( 'no' )
                  ztmp3(A2D(0),:) = 0.0
                  DO jl=1,jpl
                     ztmp3(A2D(0),1) = ztmp3(A2D(0),1) + tn_ice(A2D(0),jl) * a_i(A2D(0),jl)
                  ENDDO
               CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'oce and weighted ice')    ;   ztmp1(A2D(0)) =   ts(A2D(0),1,jp_tem,Kmm) + rt0
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(A2D(0),1:jpl) = tn_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
               CASE( 'no' )
                  ztmp3(A2D(0),:) = 0.0
                  DO jl=1,jpl
                     ztmp3(A2D(0),1) = ztmp3(A2D(0),1) + tn_ice(A2D(0),jl) * a_i(A2D(0),jl)
                  ENDDO
               CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'mixed oce-ice'        )
               ztmp1(A2D(0)) = ( ztmp1(A2D(0)) + rt0 ) * zfr_l(A2D(0))
               DO jl=1,jpl
                  ztmp1(A2D(0)) = ztmp1(A2D(0)) + tn_ice(A2D(0),jl) * a_i(A2D(0),jl)
               ENDDO
            CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%cldes' )
            END SELECT
         ENDIF
         IF( ssnd(nmodsbc)%fld(jps_toce)%laction )   CALL cpl_snd( nmodsbc ,jps_toce, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info )
         IF( ssnd(nmodsbc)%fld(jps_tice)%laction )   CALL cpl_snd( nmodsbc ,jps_tice, isec, ztmp3, info )
         IF( ssnd(nmodsbc)%fld(jps_tmix)%laction )   CALL cpl_snd( nmodsbc ,jps_tmix, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      ! 1st layer ice/snow temp.  !
      !                                                      ! ------------------------- !
#if defined key_si3
      ! needed by  Met Office
      IF( ssnd(nmodsbc)%fld(jps_ttilyr)%laction) THEN
         SELECT CASE( sn_snd_ttilyr%cldes)
         CASE ('weighted ice')
            ztmp3(A2D(0),1:jpl) = t1_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
         CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_ttilyr%cldes' )
         END SELECT
         IF( ssnd(nmodsbc)%fld(jps_ttilyr)%laction )   CALL cpl_snd( nmodsbc ,jps_ttilyr, isec, ztmp3, info )
      ENDIF
#endif
      !                                                      ! ------------------------- !
      !                                                      !           Albedo          !
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_albice)%laction ) THEN                         ! ice
          SELECT CASE( sn_snd_alb%cldes )
          CASE( 'ice' )
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )
                ztmp3(A2D(0),1:jpl) = alb_ice(A2D(0),1:jpl)
             CASE( 'no' )
                WHERE( zat_i(A2D(0)) /= 0. )
                   ztmp1(A2D(0)) = SUM( alb_ice (A2D(0),1:jpl) * a_i(A2D(0),1:jpl), dim=3 ) / zat_i(A2D(0))
                ELSEWHERE
                   ztmp1(A2D(0)) = alb_oce_mix(A2D(0))
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%clcat' )
             END SELECT
          CASE( 'weighted ice' )   ;
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )
                ztmp3(A2D(0),1:jpl) = alb_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
             CASE( 'no' )
                WHERE( fr_i (A2D(0)) > 0. )
                   ztmp1(A2D(0)) = SUM ( alb_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl), dim=3 )
                ELSEWHERE
                   ztmp1(A2D(0)) = 0.
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_ice%clcat' )
             END SELECT
          CASE default      ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%cldes' )
         END SELECT

         SELECT CASE( sn_snd_alb%clcat )
            CASE( 'yes' )
               CALL cpl_snd( nmodsbc ,jps_albice, isec, ztmp3, info )      !-> MV this has never been checked in coupled mode
            CASE( 'no'  )
               CALL cpl_snd( nmodsbc ,jps_albice, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info )
         END SELECT
      ENDIF

      IF( ssnd(nmodsbc)%fld(jps_albmix)%laction ) THEN                         ! mixed ice-ocean
         ztmp1(A2D(0)) = alb_oce_mix(A2D(0)) * zfr_l(A2D(0))
         DO jl = 1, jpl
            ztmp1(A2D(0)) = ztmp1(A2D(0)) + alb_ice(A2D(0),jl) * a_i(A2D(0),jl)
         END DO
         CALL cpl_snd( nmodsbc ,jps_albmix, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness !
      !                                                      ! ------------------------- !
      ! Send ice fraction field to atmosphere
      IF( ssnd(nmodsbc)%fld(jps_fice)%laction ) THEN
         SELECT CASE( sn_snd_thick%clcat )
         CASE( 'yes' )   ;   ztmp3(A2D(0),1:jpl) =  a_i(A2D(0),1:jpl)
         CASE( 'no'  )   ;   ztmp3(A2D(0),1    ) = fr_i(A2D(0)      )
         CASE default    ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
         END SELECT
         CALL cpl_snd( nmodsbc ,jps_fice, isec, ztmp3, info )
      ENDIF

#if defined key_si3 || defined key_cice
      ! If this coupling was successful then save ice fraction for use between coupling points.
      ! This is needed for some calculations where the ice fraction at the last coupling point
      ! is needed.
      IF(  info == OASIS_Sent    .OR. info == OASIS_ToRest .OR. &
         & info == OASIS_SentOut .OR. info == OASIS_ToRestOut ) THEN
         IF ( sn_snd_thick%clcat == 'yes' ) THEN
           a_i_last_couple(A2D(0),1:jpl) = a_i(A2D(0),1:jpl)
         ENDIF
      ENDIF
#endif

      IF( ssnd(nmodsbc)%fld(jps_fice1)%laction ) THEN
         SELECT CASE( sn_snd_thick1%clcat )
         CASE( 'yes' )   ;   ztmp3(A2D(0),1:jpl) =  a_i(A2D(0),1:jpl)
         CASE( 'no'  )   ;   ztmp3(A2D(0),1    ) = fr_i(A2D(0)      )
         CASE default    ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick1%clcat' )
         END SELECT
         CALL cpl_snd( nmodsbc ,jps_fice1, isec, ztmp3, info )
      ENDIF

      ! Send ice fraction field to OCE (sent by SAS in SAS-OCE coupling)
      IF( ssnd(nmodsbc)%fld(jps_fice2)%laction ) THEN
         ztmp3(A2D(0),1) = fr_i(A2D(0))
         IF( ssnd(nmodsbc)%fld(jps_fice2)%laction )   CALL cpl_snd( nmodsbc, jps_fice2, isec, ztmp3, info )
      ENDIF

      ! Send ice and snow thickness field
      IF( ssnd(nmodsbc)%fld(jps_hice)%laction .OR. ssnd(nmodsbc)%fld(jps_hsnw)%laction ) THEN
         SELECT CASE( sn_snd_thick%cldes)
         CASE( 'none'                  )       ! nothing to do
         CASE( 'weighted ice and snow' )
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )
               ztmp3(A2D(0),1:jpl) =  h_i(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
               ztmp4(A2D(0),1:jpl) =  h_s(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
            CASE( 'no' )
               ztmp3(A2D(0),:) = 0.0   ;  ztmp4(A2D(0),:) = 0.0
               DO jl=1,jpl
                  ztmp3(A2D(0),1) = ztmp3(A2D(0),1) + h_i(A2D(0),jl) * a_i(A2D(0),jl)
                  ztmp4(A2D(0),1) = ztmp4(A2D(0),1) + h_s(A2D(0),jl) * a_i(A2D(0),jl)
               ENDDO
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE( 'ice and snow'         )
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )
               ztmp3(A2D(0),1:jpl) = h_i(A2D(0),1:jpl)
               ztmp4(A2D(0),1:jpl) = h_s(A2D(0),1:jpl)
            CASE( 'no' )
               WHERE( zat_i(A2D(0)) /= 0. )
                  ztmp3(A2D(0),1) = SUM( h_i(A2D(0),:) * a_i(A2D(0),:), dim=3 ) / zat_i(A2D(0))
                  ztmp4(A2D(0),1) = SUM( h_s(A2D(0),:) * a_i(A2D(0),:), dim=3 ) / zat_i(A2D(0))
               ELSEWHERE
                 ztmp3(A2D(0),1) = 0.
                 ztmp4(A2D(0),1) = 0.
               END WHERE
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%cldes' )
         END SELECT
         IF( ssnd(nmodsbc)%fld(jps_hice)%laction )   CALL cpl_snd( nmodsbc ,jps_hice, isec, ztmp3, info )
         IF( ssnd(nmodsbc)%fld(jps_hsnw)%laction )   CALL cpl_snd( nmodsbc ,jps_hsnw, isec, ztmp4, info )
      ENDIF

#if defined key_si3
      !                                                      ! ------------------------- !
      !                                                      !      Ice melt ponds       !
      !                                                      ! ------------------------- !
      ! needed by Met Office: 1) fraction of ponded ice 2) local/actual pond depth
      IF( ssnd(nmodsbc)%fld(jps_a_p)%laction .OR. ssnd(nmodsbc)%fld(jps_ht_p)%laction ) THEN
         SELECT CASE( sn_snd_mpnd%cldes)
         CASE( 'ice only' )
            SELECT CASE( sn_snd_mpnd%clcat )
            CASE( 'yes' )
               ztmp3(A2D(0),1:jpl) =  a_ip_eff(A2D(0),1:jpl)
               ztmp4(A2D(0),1:jpl) =  h_ip(A2D(0),1:jpl)
            CASE( 'no' )
               WHERE( zat_i(A2D(0)) /= 0. )
                  ztmp3(A2D(0),1) = SUM( a_i(A2D(0),:) * a_ip_eff(A2D(0),:), dim=3 ) / zat_i(A2D(0))
                  ztmp4(A2D(0),1) = SUM( a_i(A2D(0),:) * h_ip(A2D(0),:), dim=3 ) / zat_i(A2D(0))
               ELSEWHERE
                 ztmp3(A2D(0),1) = 0.
                 ztmp4(A2D(0),1) = 0.
               END WHERE
            CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_mpnd%clcat' )
            END SELECT
         CASE( 'weighted ice' )
            SELECT CASE( sn_snd_mpnd%clcat )
            CASE( 'yes' )
               ztmp3(A2D(0),1:jpl) =  a_ip_eff(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
               ztmp4(A2D(0),1:jpl) =  h_ip(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
            CASE( 'no' )
               ztmp3(A2D(0),:) = 0.
               ztmp4(A2D(0),:) = 0.
               DO jl=1,jpl
                 ztmp3(A2D(0),1) = ztmp3(A2D(0),1) + a_ip_eff(A2D(0),jl) * a_i(A2D(0),jl)
                 ztmp4(A2D(0),1) = ztmp4(A2D(0),1) + h_ip(A2D(0),jl) * a_i(A2D(0),jl)
               ENDDO
            CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_mpnd%clcat' )
            END SELECT
         CASE default      ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_mpnd%cldes' )
         END SELECT
         IF( ssnd(nmodsbc)%fld(jps_a_p)%laction  )   CALL cpl_snd( nmodsbc ,jps_a_p , isec, ztmp3, info )
         IF( ssnd(nmodsbc)%fld(jps_ht_p)%laction )   CALL cpl_snd( nmodsbc ,jps_ht_p, isec, ztmp4, info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      !     Ice conductivity      !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      IF( ssnd(nmodsbc)%fld(jps_kice)%laction ) THEN
         SELECT CASE( sn_snd_cond%cldes)
         CASE( 'weighted ice' )
            SELECT CASE( sn_snd_cond%clcat )
            CASE( 'yes' )
	       ztmp3(A2D(0),1:jpl) =  cnd_ice(A2D(0),1:jpl) * a_i(A2D(0),1:jpl)
            CASE( 'no' )
               ztmp3(A2D(0),:) = 0.0
               DO jl=1,jpl
                 ztmp3(A2D(0),1) = ztmp3(A2D(0),1) + cnd_ice(A2D(0),jl) * a_i(A2D(0),jl)
               ENDDO
            CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_cond%clcat' )
            END SELECT
         CASE( 'ice only' )
           ztmp3(A2D(0),1:jpl) = cnd_ice(A2D(0),1:jpl)
         CASE default      ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_cond%cldes' )
         END SELECT
         IF( ssnd(nmodsbc)%fld(jps_kice)%laction )   CALL cpl_snd( nmodsbc ,jps_kice, isec, ztmp3, info )
      ENDIF
#endif

      !                                                      ! ------------------------- !
      !                                                      !  CO2 flux from PISCES     !
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_co2)%laction .AND. l_co2cpl )   THEN
         ztmp1(A2D(0)) = oce_co2(A2D(0)) * 1000.  ! conversion in molC/m2/s
         CALL cpl_snd( nmodsbc, jps_co2, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ) , info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_ocx1)%laction ) THEN                      !      Surface current      !
         !                                                   ! ------------------------- !
         !
         !                                                  j     -----V---F
         ! surface velocity always sent from T point               !       |
         !                                                  j      |   T   U
         !                                                         |       |
         !                                                  j-1   -I-------|
         !                                                         |       |
         !                                                        i-1  i   i
         !!clem: make a new variable at T-point to replace uu and vv => uuT and vvT for instance
         IF( nn_components == jp_iam_oce ) THEN
            zotx1(A2D(0)) = uu(A2D(0),1,Kmm)
            zoty1(A2D(0)) = vv(A2D(0),1,Kmm)
            !!clem : should be demi sum, no? Or uuT and vvT
         ELSE
            SELECT CASE( TRIM( sn_snd_crt%cldes ) )
            CASE( 'oce only'             )      ! C-grid ==> T
               DO_2D( 0, 0, 0, 0 )
                  zotx1(ji,jj) = 0.5 * ( uu(ji,jj,1,Kmm) + uu(ji-1,jj  ,1,Kmm) )
                  zoty1(ji,jj) = 0.5 * ( vv(ji,jj,1,Kmm) + vv(ji  ,jj-1,1,Kmm) )
               END_2D
            CASE( 'oce and ice'          )      ! Ocean and Ice on C-grid ==> T
               DO_2D( 0, 0, 0, 0 )
                  zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) )
                  zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) )
                  zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  )     + u_ice(ji-1,jj    )     )
                  zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  )     + v_ice(ji  ,jj-1  )     )
               END_2D
            CASE( 'weighted oce and ice' )      ! Ocean and Ice on C-grid ==> T
               DO_2D( 0, 0, 0, 0 )
                  zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)
                  zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)
                  zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  )     + u_ice(ji-1,jj    )     ) *  fr_i(ji,jj)
                  zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  )     + v_ice(ji  ,jj-1  )     ) *  fr_i(ji,jj)
               END_2D
            CASE( 'mixed oce-ice'        )      ! Ocean and Ice on C-grid ==> T
               DO_2D( 0, 0, 0, 0 )
                  zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)   &
                     &         + 0.5 * ( u_ice(ji,jj  )     + u_ice(ji-1,jj    )     ) *  fr_i(ji,jj)
                  zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)   &
                     &         + 0.5 * ( v_ice(ji,jj  )     + v_ice(ji  ,jj-1  )     ) *  fr_i(ji,jj)
               END_2D
            END SELECT
            !
         ENDIF
         !
         !
         IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) THEN             ! Rotation of the components
            !                                                                     ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(nmodsbc)%fld(jps_ocx1)%clgrid, 'ij->e', ztmp1 )       ! 1st component
            CALL rot_rep( zotx1, zoty1, ssnd(nmodsbc)%fld(jps_ocx1)%clgrid, 'ij->n', ztmp2 )       ! 2nd component
            zotx1(A2D(0)) = ztmp1(A2D(0))                                                   ! overwrite the components
            zoty1(A2D(0)) = ztmp2(A2D(0))
            IF( ssnd(nmodsbc)%fld(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(nmodsbc)%fld(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component
               CALL rot_rep( zitx1, zity1, ssnd(nmodsbc)%fld(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component
               zitx1(A2D(0)) = ztmp1(A2D(0))                                                ! overwrite the components
               zity1(A2D(0)) = ztmp2(A2D(0))
            ENDIF
         ENDIF
         !
         ! spherical coordinates to cartesian -> 2 components to 3 components
         IF( TRIM( sn_snd_crt%clvref ) == 'cartesian' ) THEN
            ztmp1(A2D(0)) = zotx1(A2D(0))                     ! ocean currents
            ztmp2(A2D(0)) = zoty1(A2D(0))
            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
            !
            IF( ssnd(nmodsbc)%fld(jps_ivx1)%laction ) THEN           ! ice velocities
               ztmp1(A2D(0)) = zitx1(A2D(0))
               ztmp1(A2D(0)) = zity1(A2D(0))
               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
            ENDIF
         ENDIF
         !
         IF( ssnd(nmodsbc)%fld(jps_ocx1)%laction )   CALL cpl_snd( nmodsbc ,jps_ocx1, isec, RESHAPE ( zotx1, (/Ni_0,Nj_0,1/) ), info )   ! ocean x current 1st grid
         IF( ssnd(nmodsbc)%fld(jps_ocy1)%laction )   CALL cpl_snd( nmodsbc ,jps_ocy1, isec, RESHAPE ( zoty1, (/Ni_0,Nj_0,1/) ), info )   ! ocean y current 1st grid
         IF( ssnd(nmodsbc)%fld(jps_ocz1)%laction )   CALL cpl_snd( nmodsbc ,jps_ocz1, isec, RESHAPE ( zotz1, (/Ni_0,Nj_0,1/) ), info )   ! ocean z current 1st grid
         !
         IF( ssnd(nmodsbc)%fld(jps_ivx1)%laction )   CALL cpl_snd( nmodsbc ,jps_ivx1, isec, RESHAPE ( zitx1, (/Ni_0,Nj_0,1/) ), info )   ! ice   x current 1st grid
         IF( ssnd(nmodsbc)%fld(jps_ivy1)%laction )   CALL cpl_snd( nmodsbc ,jps_ivy1, isec, RESHAPE ( zity1, (/Ni_0,Nj_0,1/) ), info )   ! ice   y current 1st grid
         IF( ssnd(nmodsbc)%fld(jps_ivz1)%laction )   CALL cpl_snd( nmodsbc ,jps_ivz1, isec, RESHAPE ( zitz1, (/Ni_0,Nj_0,1/) ), info )   ! ice   z current 1st grid
         !
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      !  Surface current to waves !
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_ocxw)%laction .OR. ssnd(nmodsbc)%fld(jps_ocyw)%laction ) THEN
         !
         !                                                  j     -----V---F
         ! surface velocity always sent from T point               !       |
         !                                                  j      |   T   U
         !                                                         |       |
         !                                                  j-1   -I-------|
         !                                                         |       |
         !                                                        i-1  i   i
         !!clem: make a new variable at T-point to replace uu and vv => uuT and vvT for instance
         SELECT CASE( TRIM( sn_snd_crtw%cldes ) )
         CASE( 'oce only'             )      ! C-grid ==> T
            DO_2D( 0, 0, 0, 0 )
               zotx1(ji,jj) = 0.5 * ( uu(ji,jj,1,Kmm) + uu(ji-1,jj  ,1,Kmm) )
               zoty1(ji,jj) = 0.5 * ( vv(ji,jj,1,Kmm) + vv(ji , jj-1,1,Kmm) )
            END_2D
         CASE( 'weighted oce and ice' )      ! Ocean and Ice on C-grid ==> T
            DO_2D( 0, 0, 0, 0 )
               zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)
               zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)
               zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
               zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
            END_2D
         CASE( 'mixed oce-ice'        )      ! Ocean and Ice on C-grid ==> T
            DO_2D( 0, 0, 0, 0 )
               zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)   &
                  &         + 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
               zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)   &
                  &         + 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
            END_2D
         END SELECT
         !
         IF( TRIM( sn_snd_crtw%clvor ) == 'eastward-northward' ) THEN             ! Rotation of the components
            !                                                                        ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(nmodsbc)%fld(jps_ocxw)%clgrid, 'ij->e', ztmp1 )       ! 1st component
            CALL rot_rep( zotx1, zoty1, ssnd(nmodsbc)%fld(jps_ocxw)%clgrid, 'ij->n', ztmp2 )       ! 2nd component
            zotx1(A2D(0)) = ztmp1(A2D(0))                                                   ! overwrite the components
            zoty1(A2D(0)) = ztmp2(A2D(0))
            IF( ssnd(nmodsbc)%fld(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(nmodsbc)%fld(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component
               CALL rot_rep( zitx1, zity1, ssnd(nmodsbc)%fld(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component
               zitx1(A2D(0)) = ztmp1(A2D(0))                                                ! overwrite the components
               zity1(A2D(0)) = ztmp2(A2D(0))
            ENDIF
         ENDIF
         !
         !         ! spherical coordinates to cartesian -> 2 components to 3 components
         !         IF( TRIM( sn_snd_crtw%clvref ) == 'cartesian' ) THEN
         !            ztmp1(A2D(0)) = zotx1(A2D(0))                     ! ocean currents
         !            ztmp2(A2D(0)) = zoty1(A2D(0))
         !            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
         !            !
         !            IF( ssnd(nmodsbc)%fld(jps_ivx1)%laction ) THEN           ! ice velocities
         !               ztmp1(A2D(0)) = zitx1(A2D(0))
         !               ztmp1(A2D(0)) = zity1(A2D(0))
         !               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
         !            ENDIF
         !         ENDIF
         !
         IF( ssnd(nmodsbc)%fld(jps_ocxw)%laction )   CALL cpl_snd( nmodsbc ,jps_ocxw, isec, RESHAPE ( zotx1, (/Ni_0,Nj_0,1/) ), info )   ! ocean x current 1st grid
         IF( ssnd(nmodsbc)%fld(jps_ocyw)%laction )   CALL cpl_snd( nmodsbc ,jps_ocyw, isec, RESHAPE ( zoty1, (/Ni_0,Nj_0,1/) ), info )   ! ocean y current 1st grid
         !
      ENDIF
      !
      IF( ssnd(nmodsbc)%fld(jps_ficet)%laction ) THEN
         CALL cpl_snd( nmodsbc ,jps_ficet, isec, RESHAPE ( fr_i(A2D(0)), (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !   Water levels to waves   !
      !                                                      ! ------------------------- !
      IF( ssnd(nmodsbc)%fld(jps_wlev)%laction ) THEN
         IF( ln_apr_dyn ) THEN
            IF( kt /= nit000 ) THEN
               ztmp1(A2D(0)) = ssh(A2D(0),Kbb) - 0.5 * ( ssh_ib(A2D(0)) + ssh_ibb(A2D(0)) )
            ELSE
               ztmp1(A2D(0)) = ssh(A2D(0),Kbb)
            ENDIF
         ELSE
            ztmp1(A2D(0)) = ssh(A2D(0),Kmm)
         ENDIF
         CALL cpl_snd( nmodsbc ,jps_wlev  , isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !
      !  Fields sent by OCE to SAS when doing OCE<->SAS coupling
      !                                                        ! SSH
      IF( ssnd(nmodsbc)%fld(jps_ssh )%laction )  THEN
         !                          ! removed inverse barometer ssh when Patm
         !                          forcing is used (for sea-ice dynamics)
         IF( ln_apr_dyn ) THEN   ;   ztmp1(A2D(0)) = ssh(A2D(0),Kbb) - 0.5 * ( ssh_ib(A2D(0)) + ssh_ibb(A2D(0)) )
         ELSE                    ;   ztmp1(A2D(0)) = ssh(A2D(0),Kmm)
         ENDIF
         CALL cpl_snd( nmodsbc ,jps_ssh   , isec, RESHAPE ( ztmp1 , (/Ni_0,Nj_0,1/) ), info )

      ENDIF
      !                                                        ! SSS
      IF( ssnd(nmodsbc)%fld(jps_soce  )%laction )  THEN
         CALL cpl_snd( nmodsbc ,jps_soce  , isec, RESHAPE ( ts(A2D(0),1,jp_sal,Kmm), (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !                                                        ! first T level thickness
      ze3t_i(A2D(0),1) = e3t(Nis0:Nie0,Njs0:Nje0,1,Kmm)
      IF( ssnd(nmodsbc)%fld(jps_e3t1st )%laction )  THEN
         CALL cpl_snd( nmodsbc ,jps_e3t1st, isec, RESHAPE ( ze3t_i , (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !                                                        ! Qsr fraction
      IF( ssnd(nmodsbc)%fld(jps_fraqsr)%laction )  THEN
         CALL cpl_snd( nmodsbc ,jps_fraqsr, isec, RESHAPE ( fraqsr_1lev(A2D(0)) , (/Ni_0,Nj_0,1/) ), info )
      ENDIF
      !
      !  Fields sent by SAS to OCE when OASIS coupling
      !                                                        ! Solar heat flux
      IF( ssnd(nmodsbc)%fld(jps_qsroce)%laction )  CALL cpl_snd( nmodsbc, jps_qsroce, isec, RESHAPE ( qsr (A2D(0))   , (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_qnsoce)%laction )  CALL cpl_snd( nmodsbc, jps_qnsoce, isec, RESHAPE ( qns (A2D(0))   , (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_oemp  )%laction )  CALL cpl_snd( nmodsbc, jps_oemp  , isec, RESHAPE ( emp (A2D(0)), (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_sflx  )%laction )  CALL cpl_snd( nmodsbc, jps_sflx  , isec, RESHAPE ( sfx (A2D(0))   , (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_otx1  )%laction )  CALL cpl_snd( nmodsbc, jps_otx1  , isec, RESHAPE ( utau(A2D(0)), (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_oty1  )%laction )  CALL cpl_snd( nmodsbc, jps_oty1  , isec, RESHAPE ( vtau(A2D(0)), (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_rnf   )%laction )  CALL cpl_snd( nmodsbc, jps_rnf   , isec, RESHAPE ( rnf (A2D(0)), (/Ni_0,Nj_0,1/) ), info )
      IF( ssnd(nmodsbc)%fld(jps_taum  )%laction )  CALL cpl_snd( nmodsbc, jps_taum  , isec, RESHAPE ( taum(A2D(0))   , (/Ni_0,Nj_0,1/) ), info )

#if defined key_si3
      !                                                      ! ------------------------- !
      !                                                      ! Sea surface freezing temp !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      CALL eos_fzp( ts(A2D(0),1,jp_sal,Kmm), sstfrz(A2D(0)), kbnd=0 )
      ztmp1(A2D(0)) = sstfrz(A2D(0)) + rt0
      IF( ssnd(nmodsbc)%fld(jps_sstfrz)%laction )  CALL cpl_snd( nmodsbc, jps_sstfrz, isec, RESHAPE ( ztmp1, (/Ni_0,Nj_0,1/) ), info)
#endif
      !
   END SUBROUTINE sbc_cpl_snd

   !!======================================================================
END MODULE sbccpl
