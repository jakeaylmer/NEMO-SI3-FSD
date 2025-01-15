MODULE icethd_do
   !!======================================================================
   !!                       ***  MODULE icethd_do   ***
   !!   sea-ice: sea ice growth in the leads (open water)  
   !!======================================================================
   !! History :       !  2005-12  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd_do        : ice growth in open water (=lateral accretion of ice)
   !!   ice_thd_do_init   : initialization
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_oce
   USE dom_oce , ONLY : umask, vmask, smask0
   USE phycst         ! physical constants
   USE ice1D          ! sea-ice: thermodynamics variables
   USE ice            ! sea-ice: variables
   USE sbc_oce , ONLY : sss_m
   USE sbc_ice , ONLY : utau_ice, vtau_ice
   USE icetab         ! sea-ice: 2D <==> 1D
   USE icectl         ! sea-ice: conservation
   USE icevar  , ONLY : ice_var_vremap
   USE icethd_sal     ! sea-ice: salinity profiles

   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_do        ! called by ice_thd
   PUBLIC   ice_thd_frazil    ! called by ice_thd
   PUBLIC   ice_thd_do_init   ! called by ice_stp
   !
   !                             !!** namelist (namthd_do) **
   REAL(wp) ::   rn_hinew         !  thickness for new ice formation (m)
   LOGICAL  ::   ln_frazil        !  use of frazil ice collection as function of wind (T) or not (F)
   REAL(wp) ::   rn_maxfraz       !  maximum portion of frazil ice collecting at the ice bottom
   REAL(wp) ::   rn_vfraz         !  threshold drift speed for collection of bottom frazil ice
   REAL(wp) ::   rn_Cfraz         !  squeezing coefficient for collection of bottom frazil ice

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_do
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_thd_do  ***
      !!  
      !! ** Purpose : Computation of the evolution of the ice thickness and 
      !!              concentration as a function of the heat balance in the leads
      !!       
      !! ** Method  : Ice is formed in the open water when ocean looses heat
      !!              (heat budget of open water is negative) following
      !!
      !!       (dA/dt)acc = F[ (1-A)/(1-a) ] * [ Bl / (Li*h0) ]
      !!          where - h0 is the thickness of ice created in the lead
      !!                - a is a minimum fraction for leads
      !!                - F is a monotonic non-increasing function defined as:
      !!                  F(X)=( 1 - X**exld )**(1.0/exld)
      !!                - exld is the exponent closure rate (=2 default val.)
      !! 
      !! ** Action : - Adjustment of snow and ice thicknesses and heat
      !!                content in brine pockets
      !!             - Updating ice internal temperature
      !!             - Computation of variation of ice volume and mass
      !!             - Computation of a_i after lateral accretion and 
      !!               update h_s_1d, h_i_1d      
      !!------------------------------------------------------------------------
      INTEGER  ::   ii, ji, jj, jk, jl   ! dummy loop indices
      !
      REAL(wp) ::   ztmelts
      REAL(wp) ::   zdE
      REAL(wp) ::   zQm          ! enthalpy exchanged with the ocean (J/m2, >0 towards ocean)
      REAL(wp) ::   zEi          ! sea ice specific enthalpy (J/kg)
      REAL(wp) ::   zEw          ! seawater specific enthalpy (J/kg)
      REAL(wp) ::   zfmdt        ! mass flux x time step (kg/m2, >0 towards ocean)
      !
      INTEGER  ::   jcat        ! indexes of categories where new ice grows
      !
      REAL(wp) ::   zv_newfra
      REAL(wp) ::   zv_newice   ! volume of accreted ice
      REAL(wp) ::   za_newice   ! fractional area of accreted ice
      REAL(wp) ::   ze_newice   ! heat content of accreted ice
      REAL(wp) ::   zo_newice   ! age of accreted ice
      REAL(wp) ::   zdv_res     ! residual volume in case of excessive heat budget
      REAL(wp) ::   zda_res     ! residual area in case of excessive heat budget
      REAL(wp) ::   zv_frazb    ! accretion of frazil ice at the ice bottom
      !
      REAL(wp), DIMENSION(jpl) ::   zv_b    ! old volume of ice in category jl
      REAL(wp), DIMENSION(jpl) ::   za_b    ! old area of ice in category jl
      !
      REAL(wp), DIMENSION(jpij) ::   zs_newice     ! salinity of accreted ice
      REAL(wp), DIMENSION(jpij) ::   zh_newice     ! thickness of accreted ice
      REAL(wp), DIMENSION(jpij) ::   zfraz_frac_1d ! relative ice / frazil velocity (1D vector)
      !
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zh_i_old, ze_i_old, zs_i_old
      !!-----------------------------------------------------------------------!
      !
      IF( ln_timing    )   CALL timing_start('icethd_do')
      IF( ln_icediachk )   CALL ice_cons_hsm( 0, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft )
      IF( ln_icediachk )   CALL ice_cons2D  ( 0, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft )

      !------------------------------------------------------------------------------!
      ! 1) Compute thickness, salinity, enthalpy, age, area and volume of new ice
      !------------------------------------------------------------------------------!
      ! it occurs if cooling
      at_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )

      ! Identify grid points where new ice forms
#if defined key_si3_1D      
      DO_2D( 0, 0, 0, 0 )
         npti = 0   ;   nptidx(:) = 0
         IF ( qlead(ji,jj)  <  0._wp ) THEN
            npti = 1
            nptidx( npti ) = (jj - 1) * jpi + ji
         ENDIF
#else
      npti = 0   ;   nptidx(:) = 0
      DO_2D( 0, 0, 0, 0 )
         IF ( qlead(ji,jj)  <  0._wp ) THEN
            npti = npti + 1
            nptidx( npti ) = (jj - 1) * jpi + ji
         ENDIF
      END_2D
#endif
      ! Move from 2-D to 1-D vectors
      IF ( npti > 0 ) THEN

         CALL tab_2d_1d( npti, nptidx(1:npti), at_i_1d (1:npti)    , at_i        )
         CALL tab_3d_2d( npti, nptidx(1:npti), a_i_2d  (1:npti,:)  , a_i (:,:,:) )
         CALL tab_3d_2d( npti, nptidx(1:npti), v_i_2d  (1:npti,:)  , v_i (:,:,:) )
         CALL tab_3d_2d( npti, nptidx(1:npti), sv_i_2d (1:npti,:)  , sv_i(:,:,:) )
         CALL tab_4d_3d( npti, nptidx(1:npti), e_i_2d  (1:npti,:,:), e_i    )
         CALL tab_4d_3d( npti, nptidx(1:npti), szv_i_2d(1:npti,:,:), szv_i  )
         CALL tab_2d_1d( npti, nptidx(1:npti), qlead_1d     (1:npti), qlead      )
         CALL tab_2d_1d( npti, nptidx(1:npti), t_bo_1d      (1:npti), t_bo       )
         CALL tab_2d_1d( npti, nptidx(1:npti), sfx_opw_1d   (1:npti), sfx_opw    )
         CALL tab_2d_1d( npti, nptidx(1:npti), wfx_opw_1d   (1:npti), wfx_opw    )
         CALL tab_2d_1d( npti, nptidx(1:npti), zh_newice    (1:npti), ht_i_new   )
         CALL tab_2d_1d( npti, nptidx(1:npti), zfraz_frac_1d(1:npti), fraz_frac  )

         CALL tab_2d_1d( npti, nptidx(1:npti), hfx_thd_1d(1:npti), hfx_thd    )
         CALL tab_2d_1d( npti, nptidx(1:npti), hfx_opw_1d(1:npti), hfx_opw    )
         CALL tab_2d_1d( npti, nptidx(1:npti), rn_amax_1d(1:npti), rn_amax_2d )
         CALL tab_2d_1d( npti, nptidx(1:npti), sss_1d    (1:npti), sss_m      )

         ! Convert units for ice internal energy and salt content
         DO jl = 1, jpl
            DO jk = 1, nlay_i               
               WHERE( v_i_2d(1:npti,jl) > 0._wp )
                  e_i_2d  (1:npti,jk,jl) = e_i_2d  (1:npti,jk,jl) / v_i_2d(1:npti,jl) * REAL( nlay_i )
                  szv_i_2d(1:npti,jk,jl) = szv_i_2d(1:npti,jk,jl) / v_i_2d(1:npti,jl) * REAL( nlay_i )
               ELSEWHERE
                  e_i_2d  (1:npti,jk,jl) = 0._wp
                  szv_i_2d(1:npti,jk,jl) = 0._wp
               END WHERE
            END DO
         END DO

         ! --- Salinity of new ice --- ! 
         SELECT CASE ( nn_icesal )
         CASE ( 1 )                    ! Sice = constant 
            zs_newice(1:npti) = rn_icesal
         CASE ( 2 , 4 )                ! Sice = F(z,t) [Griewank and Notz 2013 ; Rees Jones and Worster 2014]
            zs_newice(1:npti) = rn_sinew * sss_1d(1:npti)
         CASE ( 3 )                    ! Sice = F(z) [multiyear ice]
            zs_newice(1:npti) =   2.3
         END SELECT
         !
         !                       ! ==================== !
         !                       ! Start main loop here !
         !                       ! ==================== !
         DO ii = 1, npti
            
            ! Keep old ice areas and volume in memory
            DO jl = 1, jpl
               zv_b(jl) = v_i_2d(ii,jl) 
               za_b(jl) = a_i_2d(ii,jl)
            ENDDO
            
            ! --- Heat content of new ice --- !
            ! We assume that new ice is formed at the seawater freezing point
            ztmelts   = - rTmlt * zs_newice(ii)                  ! Melting point (C)
            ze_newice =   rhoi * (  rcpi  * ( ztmelts - ( t_bo_1d(ii) - rt0 ) )                     &
               &                  + rLfus * ( 1.0 - ztmelts / MIN( t_bo_1d(ii) - rt0, -epsi10 ) )   &
               &                  - rcp   *         ztmelts )
            
            ! --- Age of new ice --- !
            zo_newice = 0._wp

            ! --- Volume of new ice --- !
            zEi           = - ze_newice * r1_rhoi                  ! specific enthalpy of forming ice [J/kg]

            zEw           = rcp * ( t_bo_1d(ii) - rt0 )            ! specific enthalpy of seawater at t_bo_1d [J/kg]
                                                                   ! clem: we suppose we are already at the freezing point (condition qlead<0 is satisfyied) 
                                                                   
            zdE           = zEi - zEw                              ! specific enthalpy difference [J/kg]
                                              
            zfmdt         = - qlead_1d(ii) / zdE                   ! Fm.dt [kg/m2] (<0) 
                                                                   ! clem: we use qlead instead of zqld (icethd) because we suppose we are at the freezing point   
            zv_newice     = - zfmdt * r1_rhoi

            zQm           = zfmdt * zEw                            ! heat to the ocean >0 associated with mass flux  

            ! Contribution to heat flux to the ocean [W.m-2], >0  
            hfx_thd_1d(ii) = hfx_thd_1d(ii) + zfmdt * zEw * r1_Dt_ice
            ! Total heat flux used in this process [W.m-2]  
            hfx_opw_1d(ii) = hfx_opw_1d(ii) - zfmdt * zdE * r1_Dt_ice
            ! mass flux
            wfx_opw_1d(ii) = wfx_opw_1d(ii) - zv_newice * rhoi * r1_Dt_ice
            ! salt flux
            sfx_opw_1d(ii) = sfx_opw_1d(ii) - zv_newice * rhoi * zs_newice(ii) * r1_Dt_ice
         
            ! A fraction fraz_frac of frazil ice is accreted at the ice bottom
            IF( at_i_1d(ii) > 0._wp ) THEN
               zv_frazb  =           zfraz_frac_1d(ii)   * zv_newice
               zv_newice = ( 1._wp - zfraz_frac_1d(ii) ) * zv_newice
            ELSE
               zv_frazb  = 0._wp
            ENDIF
            ! --- Area of new ice --- !
            za_newice = zv_newice / zh_newice(ii)


            ! --- Redistribute new ice area and volume into ice categories --- !

            ! --- lateral ice growth --- !
            ! If lateral ice growth gives an ice concentration > amax, then
            ! we keep the excessive volume in memory and attribute it later to bottom accretion
            IF ( za_newice >  MAX( 0._wp, rn_amax_1d(ii) - at_i_1d(ii) ) ) THEN ! max is for roundoff error
               zda_res   = za_newice - MAX( 0._wp, rn_amax_1d(ii) - at_i_1d(ii) )
               zdv_res   = zda_res * zh_newice(ii) 
               za_newice = MAX( 0._wp, za_newice - zda_res )
               zv_newice = MAX( 0._wp, zv_newice - zdv_res )
            ELSE
               zda_res = 0._wp
               zdv_res = 0._wp
            ENDIF

            ! find which category to fill
            at_i_1d(ii) = 0._wp
            DO jl = 1, jpl
               IF( zh_newice(ii) > hi_max(jl-1) .AND. zh_newice(ii) <= hi_max(jl) ) THEN
                  a_i_2d(ii,jl) = a_i_2d(ii,jl) + za_newice
                  v_i_2d(ii,jl) = v_i_2d(ii,jl) + zv_newice
                  jcat = jl
               ENDIF
               at_i_1d(ii) = at_i_1d(ii) + a_i_2d(ii,jl)
            END DO

            ! Heat content
            jl = jcat                                             ! categroy in which new ice is put
            IF( za_b(jl) > 0._wp ) THEN   
               e_i_2d  (ii,:,jl) = ( ze_newice     * zv_newice + e_i_2d  (ii,:,jl) * zv_b(jl) ) / MAX( v_i_2d(ii,jl), epsi20 )
               szv_i_2d(ii,:,jl) = ( zs_newice(ii) * zv_newice + szv_i_2d(ii,:,jl) * zv_b(jl) ) / MAX( v_i_2d(ii,jl), epsi20 )
            ELSE
               e_i_2d  (ii,:,jl) = ze_newice   
               szv_i_2d(ii,:,jl) = zs_newice(ii)   
            ENDIF
         
            ! --- bottom ice growth + ice enthalpy remapping --- !
            DO jl = 1, jpl
               
               ! for remapping
               zh_i_old(0:nlay_i+1) = 0._wp
               ze_i_old(0:nlay_i+1) = 0._wp
               zs_i_old(0:nlay_i+1) = 0._wp
               DO jk = 1, nlay_i
                  zh_i_old(jk) =                      v_i_2d(ii,jl) * r1_nlay_i
                  ze_i_old(jk) = e_i_2d  (ii,jk,jl) * v_i_2d(ii,jl) * r1_nlay_i
                  zs_i_old(jk) = szv_i_2d(ii,jk,jl) * v_i_2d(ii,jl) * r1_nlay_i
               END DO

               ! new volumes including lateral/bottom accretion + residual
               IF( at_i_1d(ii) >= epsi20 ) THEN
                  zv_newfra     = ( zdv_res + zv_frazb ) * a_i_2d(ii,jl) / MAX( at_i_1d(ii) , epsi20 )
               ELSE                  
                  zv_newfra     = 0._wp
                  a_i_2d(ii,jl) = 0._wp
               ENDIF
               v_i_2d(ii,jl) = v_i_2d(ii,jl) + zv_newfra
               ! for remapping
               zh_i_old(nlay_i+1) = zv_newfra
               ze_i_old(nlay_i+1) = ze_newice     * zv_newfra
               zs_i_old(nlay_i+1) = zs_newice(ii) * zv_newfra
           
               ! --- Update bulk salinity --- !
               sv_i_2d(ii,jl) = sv_i_2d(ii,jl) + zs_newice(ii) * ( v_i_2d(ii,jl) - zv_b(jl) )
              
               ! --- Ice enthalpy and salt remapping --- !
                                      CALL ice_var_vremap( zh_i_old, ze_i_old, e_i_2d  (ii,:,jl) ) 
               IF( nn_icesal == 4 )   CALL ice_var_vremap( zh_i_old, zs_i_old, szv_i_2d(ii,:,jl) ) 
               !
            END DO
            
         END DO ! npti
         !                       ! ================== !
         !                       ! End main loop here !
         !                       ! ================== !
         !
         ! Change units for e_i/szv_i
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               e_i_2d  (1:npti,jk,jl) = e_i_2d  (1:npti,jk,jl) * v_i_2d(1:npti,jl) * r1_nlay_i 
               szv_i_2d(1:npti,jk,jl) = szv_i_2d(1:npti,jk,jl) * v_i_2d(1:npti,jl) * r1_nlay_i 
            END DO
         END DO

         ! Move 2D vectors to 1D vectors 
         CALL tab_2d_3d( npti, nptidx(1:npti), a_i_2d  (1:npti,:)  , a_i (:,:,:) )
         CALL tab_2d_3d( npti, nptidx(1:npti), v_i_2d  (1:npti,:)  , v_i (:,:,:) )
         CALL tab_2d_3d( npti, nptidx(1:npti), sv_i_2d (1:npti,:)  , sv_i(:,:,:) )
         CALL tab_3d_4d( npti, nptidx(1:npti), e_i_2d  (1:npti,:,:), e_i   )
         CALL tab_3d_4d( npti, nptidx(1:npti), szv_i_2d(1:npti,:,:), szv_i )
         CALL tab_1d_2d( npti, nptidx(1:npti), sfx_opw_1d(1:npti), sfx_opw )
         CALL tab_1d_2d( npti, nptidx(1:npti), wfx_opw_1d(1:npti), wfx_opw )
         CALL tab_1d_2d( npti, nptidx(1:npti), hfx_thd_1d(1:npti), hfx_thd )
         CALL tab_1d_2d( npti, nptidx(1:npti), hfx_opw_1d(1:npti), hfx_opw )
         !
      ENDIF ! npti > 0
      !
#if defined key_si3_1D
      END_2D
#endif
      !
      ! the following fields need to be updated on the halos (done in icethd): a_i, v_i, sv_i, e_i 
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      IF( ln_timing    )   CALL timing_stop ('icethd_do')
      !
   END SUBROUTINE ice_thd_do


   SUBROUTINE ice_thd_frazil
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_frazil ***
      !!
      !! ** Purpose :   frazil ice collection thickness and fraction
      !!
      !! ** Inputs  :   u_ice, v_ice, utau_ice, vtau_ice
      !! ** Ouputs  :   ht_i_new, fraz_frac
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj             ! dummy loop indices
      INTEGER  ::   iter
      REAL(wp) ::   zvfrx, zvgx, ztaux, zf, ztenagm, zvfry, zvgy, ztauy, zvrel2, zfp, ztwogp
      REAL(wp), PARAMETER ::   zcai    = 1.4e-3_wp                       ! ice-air drag (clem: should be dependent on coupling/forcing used)
      REAL(wp), PARAMETER ::   zhicrit = 0.04_wp                         ! frazil ice thickness
      REAL(wp), PARAMETER ::   zsqcd   = 1.0_wp / SQRT( 1.3_wp * zcai )  ! 1/SQRT(airdensity*drag)
      REAL(wp), PARAMETER ::   zgamafr = 0.03_wp
      !!-----------------------------------------------------------------------
      !
      !---------------------------------------------------------!
      ! Collection thickness of ice formed in leads and polynyas
      !---------------------------------------------------------!    
      ! ht_i_new is the thickness of new ice formed in open water
      ! ht_i_new can be either prescribed (ln_frazil=F) or computed (ln_frazil=T)
      ! Frazil ice forms in open water, is transported by wind, accumulates at the edge of the consolidated ice edge
      ! where it forms aggregates of a specific thickness called collection thickness.
      !
      fraz_frac(:,:) = 0._wp
      !
      ! Default new ice thickness
      WHERE( qlead(:,:) < 0._wp ) ! cooling
         ht_i_new(:,:) = rn_hinew
      ELSEWHERE
         ht_i_new(:,:) = 0._wp
      END WHERE

      IF( ln_frazil ) THEN
         ztwogp  = 2._wp * rho0 / ( grav * 0.3_wp * ( rho0 - rhoi ) )  ! reduced grav
         !
         DO_2D( 0, 0, 0, 0 )
            IF ( qlead(ji,jj) < 0._wp ) THEN ! cooling
               ! -- Wind stress -- !
               ztaux = utau_ice(ji,jj) * smask0(ji,jj)
               ztauy = vtau_ice(ji,jj) * smask0(ji,jj)
               ! Square root of wind stress
               ztenagm = SQRT( SQRT( ztaux * ztaux + ztauy * ztauy ) )

               ! -- Frazil ice velocity -- !
               IF( ztenagm >= epsi10 ) THEN
                  zvfrx = zgamafr * zsqcd * ztaux / MAX( ztenagm, epsi10 )
                  zvfry = zgamafr * zsqcd * ztauy / MAX( ztenagm, epsi10 )
               ELSE
                  zvfrx = 0._wp
                  zvfry = 0._wp
               ENDIF
               ! -- Pack ice velocity -- !
               zvgx = ( u_ice(ji-1,jj  ) * umask(ji-1,jj  ,1)  + u_ice(ji,jj) * umask(ji,jj,1) ) * 0.5_wp
               zvgy = ( v_ice(ji  ,jj-1) * vmask(ji  ,jj-1,1)  + v_ice(ji,jj) * vmask(ji,jj,1) ) * 0.5_wp

               ! -- Relative frazil/pack ice velocity & fraction of frazil ice-- !
               IF( at_i(ji,jj) >= epsi10 ) THEN
                  zvrel2 = MAX( (zvfrx - zvgx)*(zvfrx - zvgx) + (zvfry - zvgy)*(zvfry - zvgy), 0.15_wp*0.15_wp )
                  fraz_frac(ji,jj) = ( TANH( rn_Cfraz * ( SQRT(zvrel2) - rn_vfraz ) ) + 1._wp ) * 0.5_wp * rn_maxfraz
               ELSE
                  zvrel2 = 0._wp
                  fraz_frac(ji,jj) = 0._wp
               ENDIF
               
               ! -- new ice thickness (iterative loop) -- !
               ht_i_new(ji,jj) = zhicrit +   ( zhicrit + 0.1_wp )    &
                  &                      / ( ( zhicrit + 0.1_wp ) * ( zhicrit + 0.1_wp ) -  zhicrit * zhicrit ) * ztwogp * zvrel2
               iter = 1
               DO WHILE ( iter < 20 ) 
                  zf  = ( ht_i_new(ji,jj) - zhicrit ) * ( ht_i_new(ji,jj) * ht_i_new(ji,jj) - zhicrit * zhicrit ) -   &
                     &    ht_i_new(ji,jj) * zhicrit * ztwogp * zvrel2
                  zfp = ( ht_i_new(ji,jj) - zhicrit ) * ( 3.0_wp * ht_i_new(ji,jj) + zhicrit ) - zhicrit * ztwogp * zvrel2

                  ht_i_new(ji,jj) = ht_i_new(ji,jj) - zf / MAX( zfp, epsi20 )
                  iter = iter + 1
               END DO
               !
               ! bound ht_i_new (though I don't see why it should be necessary)
               ht_i_new(ji,jj) = MAX( 0.01_wp, MIN( ht_i_new(ji,jj), hi_max(jpl) ) )
               !
            ELSE
               ht_i_new(ji,jj) = 0._wp
            ENDIF
            !
         END_2D
         ! 
      ENDIF
   END SUBROUTINE ice_thd_frazil

   SUBROUTINE ice_thd_do_init
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_do_init *** 
      !!                 
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice growth in the leads
      !!
      !! ** Method  :   Read the namthd_do namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_do
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer 
      !!
      NAMELIST/namthd_do/ rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_do)
      READ_NML_CFG(numnam_ice,namthd_do)
      IF(lwm) WRITE( numoni, namthd_do )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_do_init: Ice growth in open water'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_do:'
         WRITE(numout,*) '      ice thickness for lateral accretion                       rn_hinew   = ', rn_hinew
         WRITE(numout,*) '      Frazil ice thickness as a function of wind or not         ln_frazil  = ', ln_frazil
         WRITE(numout,*) '      Maximum proportion of frazil ice collecting at bottom     rn_maxfraz = ', rn_maxfraz
         WRITE(numout,*) '      Threshold relative drift speed for collection of frazil   rn_vfraz   = ', rn_vfraz
         WRITE(numout,*) '      Squeezing coefficient for collection of frazil            rn_Cfraz   = ', rn_Cfraz
      ENDIF
      !
      IF ( rn_hinew < rn_himin )   CALL ctl_stop( 'ice_thd_do_init : rn_hinew should be >= rn_himin' )
      !
   END SUBROUTINE ice_thd_do_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_do
