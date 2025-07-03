MODULE icewav
   !!======================================================================
   !!                       ***  MODULE icewav ***
   !!   sea-ice : ocean wave-ice interactions
   !!======================================================================
   !! History :  5.0  !  2025     (J.R. Aylmer)         Original code based
   !!                                                   on CPOM-CICE and
   !!                                                   CICE/Icepack
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_wav_newice  : calculate floe size of new ice from local wave properties
   !!   ice_wav_attn    : calculate attenuated wave spectrum under ice
   !!   ice_wav_calc    : calculate wave properties under ice
   !!   ice_wav_frac    : wave breakup of ice, main routine called by ice_stp
   !!   wav_frac_dist   : calculate distribution of ice fracture lengths due to waves
   !!   wav_spec_bret() : Bretschneider wave spectrum formula
   !!   ice_wav_init    : initialisation of wave-ice interaction module
   !!----------------------------------------------------------------------
   USE par_ice           ! SI3 parameters
   USE phycst , ONLY :   rpi, grav, rhoi
   USE sbc_oce, ONLY :   ln_wave, ln_wave_spec, nn_nwfreq   ! SBC: wave module
   USE sbcwave, ONLY :   hsw, wpf, wfreq, wdfreq, wspec     ! SBC: wave variables
   USE ice               ! sea-ice: variables
   USE icefsd , ONLY :   a_ifsd, nf_newice, floe_rl, floe_rc, floe_ru   ! floe size distribution parameters/variables
   USE icefsd , ONLY :   rDt_ice_fsd, fsd_cleanup                       ! floe size distribution functions/routines

   USE in_out_manager    ! I/O manager (needed for lwm and lwp logicals)
   USE iom               ! I/O manager library (needed for iom_put)
   USE lib_mpp           ! MPP library (needed for read_nml_substitute.h90)
   USE lbclnk            ! lateral boundary conditions (or mpp links)
   USE timing            ! Timing

   IMPLICIT NONE

   PUBLIC ::   ice_wav_newice   ! routine called by ice_thd_do
   PUBLIC ::   ice_wav_attn     ! routine called by sbc_wave
   PUBLIC ::   ice_wav_frac     ! routine called by ice_stp
   PUBLIC ::   ice_wav_init     ! routine called by ice_init

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   x1d     ! 1D subdomain for wave fracture
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   wknum   ! angular wave numbers corresponding to wfreq array (m-1)

   LOGICAL ::   l_frac_calc_spec   ! whether spectrum needs to be calculated in subroutine ice_wav_frac

   !                                     !!* namelist (namwav) *
   LOGICAL         ::   ln_ice_wav_rand   !: Use random phases for sea surface height in wave fracture calculation
   INTEGER         ::   nn_ice_wav_nx1d   !: Size of 1D subdomain for wave fracture calculation
   REAL(wp)        ::   rn_ice_wav_dx1d   !: Increment of 1D subdomain for wave fracture calculation (meters)
   INTEGER         ::   nn_ice_wav_rmin   !: Radius of smallest floes affected by wave fracture in units of rn_ice_wav_dx1d
   REAL(wp)        ::   rn_ice_wav_ecri   !: Critical strain at which ice fractures due to waves (dimensionless)
   !
   ! Note: other namwav parameters declared in par_ice as they are needed by module
   ! ----- sbcwave which cannot access this module (would create circular dependency)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"

CONTAINS


   SUBROUTINE ice_wav_newice( phsw, pwpf, kcat )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_newice ***
      !!
      !! ** Purpose :   Calculate floe size (category) of new ice grown in
      !!                open water dependent on local wave conditions
      !!
      !! ** Method  :   Maximum floe diameter is determined from:
      !!
      !!                   2*r_max = sqrt[ (2*C2*Lp^2) / (Wa*g*rhoi*pi^3) ]
      !!
      !!                where:
      !!
      !!                   C2   : tensile stress mode parameter (kg.m-1.s-2)
      !!                   Lp   : wavelength of peak wave energy density (m)
      !!                   Wa   : wave amplitude (m)
      !!                   g    : acceleration due to gravity (m.s-2)
      !!                   rhoi : density of sea ice (kg.m-3)
      !!
      !!                Following Roach et al. (2019), Wa is taken to be half
      !!                of the local significant wave height Hs, Lp is computed
      !!                from the local peak frequency fp using the dispersion
      !!                relation for deep water surface gravity waves:
      !!
      !!                   Lp = g / (2*pi*fp^2),
      !!
      !!                and C2 is a hard-coded constant, hence:
      !!
      !!                   r_max = [1/(2*pi^2*fp^2)] * sqrt[ (C2*g) / (pi*rhoi*Hs) ]
      !!
      !! ** Inputs  :   phsw : local significant wave height (m)
      !!                pwpf : local wave frequency of peak energy (Hz)
      !!
      !! ** Outputs :   kcat : integer index of floe size distribution category
      !!                       to which new ice formation is added
      !!
      !! ** Callers :   ice_thd_do -> [ice_wav_newice]
      !!
      !! ** References
      !!    ----------
      !!    Roach, L. A., Bitz, C. M., Horvat, C., & Dean, S. M. (2019).
      !!              Advances in modeling interactions between sea ice and ocean surface waves.
      !!              Journal of Advances in Modeling Earth Systems, 11, 4167-4181.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), INTENT(in)    ::   phsw   ! local significant wave height (m)
      REAL(wp), INTENT(in)    ::   pwpf   ! local wave peak frequency (Hz)
      INTEGER , INTENT(inout) ::   kcat   ! FSD category index for new ice growth
      !
      REAL(wp), PARAMETER ::   zC2 = .167_wp   ! tensile mode stress parameter (kg.m-1.s-2)
      REAL(wp)            ::   zrmax           ! maximum new ice floe size (m)
      INTEGER             ::   jf              ! dummy loop index
      !
      !!-------------------------------------------------------------------

      IF( (phsw < epsi06) .OR. (pwpf < epsi06) ) THEN
         !
         kcat = nf_newice   ! no waves present => set to default new floe size category
         !
      ELSE
         !
         ! --- Calculate maximum floe size from wave conditions, r_max:
         zrmax = SQRT( zC2 * grav / (rpi * rhoi * phsw) ) / (2._wp * rpi**2 * pwpf**2)
         !
         ! --- Find FSD category that r_max belongs to:
         !
         IF( zrmax < floe_rl(1) ) THEN
            ! Smaller than smallest 'resolved' floe size => put in smallest cat. anyway:
            kcat = 1
            !
         ELSE
            ! Find FSD category that r_max belongs to. Note that if r_max exceeds upper
            ! limit of largest floe size category, it goes into that category anyway:
            !
            DO jf = nn_nfsd, 1, -1
               IF( zrmax > floe_rl(jf) ) THEN
                  kcat = jf
                  EXIT   ! found kcat => stop iterating
               ENDIF
            ENDDO
            !
         ENDIF
         !
      ENDIF

   END SUBROUTINE ice_wav_newice


   SUBROUTINE ice_wav_attn
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_attn ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate attenuated wave spectrum and derived wave
      !!                properties in sea ice grid cells
      !!
      !! ** Method  :   For each ice-covered grid cell, locate the nearest equatorward
      !!                open-ocean grid cell along meridians. Assume waves propogate
      !!                from such open-ocean grid cells to the sea ice cell along such
      !!                meridians, with the wave energy being attenuated according to
      !!                the number of ice floes and mean ice thickness along the way.
      !!
      !!                This waves-in-ice attenuation scheme follows the method of
      !!                Roach et al. (2018) with related theory and expression for the
      !!                wave attenuation coefficients from Horvat and Tziperman (2015).
      !!
      !! ** Callers :   sbc_wave --> [ice_wav_attn]
      !! ** Calls   :                [ice_wav_attn] --> ice_wav_calc
      !! ** Invokes :                [ice_wav_attn] --> wav_spec_bret()
      !!
      !! ** Notes   :   This routine is only called when ln_ice_wav_attn=T. It either
      !!                uses the wave spectrum in the nearest open ocean if available
      !!                (ln_wave_spec=T, ln_ice_wav_spec=T), otherwise it estimates
      !!                the spectrum using the Bretschneider formula with the
      !!                significant wave height and peak frequency wave field inputs.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model.
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      !!-------------------------------------------------------------------
   END SUBROUTINE ice_wav_attn


   SUBROUTINE ice_wav_calc
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_calc ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate wave properties under sea ice from wave spectrum
      !!
      !! ** Method  :   Significant wave height, Hs, is given by (e.g., Holthuijsen 2007):
      !!
      !!                   Hs = 4 * sqrt{ int[ E(f)df ] }
      !!
      !!                where E(f) is the wave energy spectrum as a function of
      !!                frequency, f. Peak frequency, fp, satisfies:
      !!
      !!                   E(fp) = max[ E(f) ]
      !!
      !! ** Callers :   ice_wav_attn --> [ice_wav_calc]
      !!
      !! ** Note    :   This routine is only called when ln_ice_wav_attn=T.
      !!                It updates the global hsw and wpf arrays but only for
      !!                grid cells containing sea ice.
      !!
      !! ** References
      !!    ----------
      !!    Holthuijsen, L. H. (2007)
      !!              Waves in Oceanic and Coastal Waters.
      !!              Cambridge University Press, p70, ISBN 978-0-521-86028-4.
      !!-------------------------------------------------------------------
      !
      !!-------------------------------------------------------------------
   END SUBROUTINE ice_wav_calc


   SUBROUTINE ice_wav_frac( kt )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wav_frac  ***
      !!
      !! ** Purpose :   Evolve the floe size distribution (FSD) subject to fracture by ocean waves
      !!
      !! ** Method  :   1.   Calculate wave energy spectrum, E(f), if needed. This is done
      !!                     when ln_wave_spec=F OR ln_ice_wav_spec=F, in which case E(f) is
      !!                     estimated as a Bretschneider spectrum from local (under ice)
      !!                     significant wave height and peak frequency (from file or wave model).
      !!
      !!                     If in addition ln_ice_wav_attn=T, then this step is skipped as E(f)
      !!                     under ice will already have been calculated in subroutine ice_wav_attn.
      !!
      !!                     If ln_wave_spec=T AND ln_ice_wav_spec=T, this step is also skipped
      !!                     as these options mean the wave spectrum is being read in from a file
      !!                     or (more likely) from a wave model.
      !!
      !!                2.   Calculate sea surface height n(x) along 1D sub-grid domain x
      !!                     (conceptually aligned along direction of wave propagation):
      !!
      !!                        n(x) = SUM [ C(f) * COS( 2*pi*x/L(f) + phi(f) ) ]
      !!
      !!                     where C(f) = sqrt[2*E(f)df], L(f) is the wavelength of frequency
      !!                     component f, and phi(f) is its phase. The dispersion relation for
      !!                     surface deep water gravity waves, L(f) = g / (2*pi*f^2), where
      !!                     g is the acceleration due to gravity. The phase phi(f) is either
      !!                     fixed at pi for all f or chosen to be a random value between 0 and
      !!                     2*pi depending on ln_ice_wav_rand.
      !!
      !!                3.   Calculate distribution of ice fracture lengths in each grid
      !!                     cell (subroutine wav_frac_dist) from n(x).
      !!
      !!                4.   Calculate the tendency of the FSD:
      !!
      !!                        df(r,h)/dt = -Om(r,h) + iint [ Om(r',h')*Zt(r,h,r',h') dr'dh' ]
      !!
      !!                     where f(r,h) is the floe size-thickness distribution, Om(r,h) is the
      !!                     distribution of floes of size r and thickness h that are fractured by waves
      !!                     per unit ocean area per unit time, and Zt(r,h,r',h') is that for floes
      !!                     formed by fractures of floes of size r', thickness h' (Roach et al. 2018).
      !!                     Om and Zt are calculated from the fracture distribution (step 3):
      !!
      !!                        Om(r,h) = f(r,h) * int[ r'W(r') dr' ] / (D/2)
      !!
      !!                     where the integral is over all floe sizes up to (but excluding) r. Here
      !!                     there is an implicit dependence on wave group velocity and domain D size
      !!                     used to compute W(r), which introduces the rate (d/dt) and quantifies the
      !!                     fraction of ice reached by waves. This term comes from Horvat and Tziperman (2015)
      !!                     where the model is applied to a single grid box; in this context, the wave
      !!                     field is a local quantity and assumed to affect all ice in the grid cell,
      !!                     so this factor is set to 1 (per second). The other term:
      !!
      !!                        Zt(r,h,r',h') = rW(r) / int[ r'W(r') dr' ]   (h' == h) AND (r' > r)
      !!                                      = 0                            otherwise
      !!
      !!                     (see Horvat and Tziperman, 2015 and Roach et al. 2018). Note that the quantity
      !!                     returned by subroutine wav_frac_dist is rW(r) and normalised, so the factor
      !!                     of D/2 is not needed here (and cancels anyway in the equation for Zt).
      !!
      !!                5.   Evolve the FSD with tendency calculated in 4 using adaptive time stepping
      !!                     (Horvat and Tziperman, 2017).
      !!
      !! ** Callers :   ice_stp -> [ice_wav_frac]
      !! ** Calls   :              [ice_wav_frac] -> wav_frac_dist
      !! ** Invokes :              [ice_wav_frac] -> wav_spec_bret()   (ln_ice_wav_spec = F AND ln_ice_wav_attn = F)
      !!                           [ice_wav_frac] -> rDt_ice_fsd()
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Horvat, C., & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER , INTENT(in) ::   kt   ! ocean time step
      !
      REAL(wp), DIMENSION(nn_nwfreq)       ::   zCf            ! spectral coefficient (m)
      REAL(wp), DIMENSION(nn_nwfreq)       ::   zphi           ! spectral phase (rad)
      REAL(wp), DIMENSION(nn_ice_wav_nx1d) ::   zssh           ! sea surface height profile (m)
      REAL(wp), DIMENSION(nn_nfsd)         ::   zWfrac         ! wave fracture distribution
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd) ::   zZeta          ! term in FSD tendency equation
      REAL(wp), DIMENSION(nn_nfsd)         ::   zOmega         ! term in FSD tendency equation
      REAL(wp), DIMENSION(nn_nfsd)         ::   za_ifsd_tend   ! tendency of FSD due to wave fracture
      REAL(wp)                             ::   zfsd_res       ! correction term for area conservation
      !
      REAL(wp)  ::   zdt_sub     ! adaptive time step (s)
      REAL(wp)  ::   ztelapsed   ! to track time elapsed during adaptive time stepping (s)
      INTEGER   ::   isubt       ! to track number of adaptive time steps
      !
      INTEGER                             ::   ji, jj, jl, jf   ! dummy loop indices
      !
      REAL(wp), PARAMETER :: zat_i_min = .01_wp   ! minimum concentration for fracture to occur
      INTEGER , PARAMETER :: isubt_max = 100      ! maximum number of adaptive time steps before warning
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('ice_wav_frac')

      ! Calculate angular wave number (2 * pi / wavelength) assuming dispersion relation of
      ! surface deep water gravity waves (this array is a constant, so calculate once.
      ! Cannot do this in ice_wav_init as that is called before sbc_wave_init..)
      !
      IF( kt == nit000 ) THEN   ! at first time-step
         ALLOCATE( wknum(nn_nwfreq) )
         wknum(:) = 4._wp * rpi**2 * wfreq(:)**2 / grav
      ENDIF

      !-----------------!
      ! Begin main loop !
      !-----------------!
      DO_2D( 0, 0, 0, 0 )

         !==================================================!
         ! (1) Calculate wave spectrum, if needed           !
         !==================================================!
         !
         ! Condition depends on combination of various namelist flags; the net condition
         ! is saved in module variable l_frac_calc_spec calculated in ice_wav_init
         !
         IF( l_frac_calc_spec ) THEN
            wspec(ji,jj,:) = wav_spec_bret( hsw(ji,jj), wpf(ji,jj) )
         ENDIF

         ! Do not calculate fracture for: -total ice concentration below threshold
         !                                -negligible wave spectrum:
         IF( (at_i(ji,jj) > zat_i_min) .AND. (MAXVAL(wspec(ji,jj,:)) > epsi06) ) THEN
            !==================================================!
            ! (2) Calculate sea surface height in 1D subdomain !
            !==================================================!
            !
            ! Spectral coefficients in expansion of local sea surface height:
            zCf(:) = SQRT( 2._wp * wspec(ji,jj,:) * wdfreq(:) )

            ! Spectral phases [constant for now; possibly add (optional!) random phase later]:
            zphi(:) = rpi

            ! Calculate sea surface height along 1D subdomain (x1d):
            zssh(:) = 0._wp
            DO jf = 1, nn_nwfreq
               zssh(:) = zssh(:) + zCf(jf) * COS( zphi(jf) + wknum(jf) * x1d(:) )
            ENDDO

            !==================================================!
            ! (3) Calculate 'fracture histogram'               !
            !==================================================!
            !
            CALL wav_frac_dist( zssh(:), vt_i(ji,jj) / at_i(ji,jj), zWfrac(:) )

            ! Proceed only if some fractures can occur, implied by non-zero zWfrac.
            ! Fracturing quantified by zWfrac is applied to each ice thickness category
            ! if possible (enough ice to begin with), in proportion to its concentration
            !
            IF( MAXVAL(zWfrac(:)) > epsi10 ) THEN
               !--------------------------------------!
               ! Begin sub-loop: thickness categories !
               !--------------------------------------!
               DO jl = 1, jpl
                  !
                  CALL fsd_cleanup( a_ifsd(ji,jj,:,jl) )   ! necessary?
                  !
                  ! Conditions for wave fracture to occur in this ITD category:
                  !  - ice concentration (in this category) cannot be too small
                  !  - cannot have all ice in smallest floe size category
                  !  - total FSD cannot be zero (needed?)
                  !
                  IF(          (a_i(ji,jj,jl) > epsi06)                     &
                     &   .AND. (a_ifsd(ji,jj,1,jl) < 1._wp)                 &
                     &   .AND. (SUM(a_ifsd(ji,jj,:,jl)) >  epsi10) ) THEN
                     !
                     !==================================================!
                     ! (4) Calculate FSD tendency terms                 !
                     !==================================================!
                     !
                     ! zZeta(:,:) is independent of FSD, so calculate before adaptive time stepping
                     ! zOmega(:) depends on FSD and so must be (re-)calculated during adaptive time stepping
                     !
                     zZeta(:,:) = 0._wp
                     !
                     DO jf = 2, nn_nfsd
                        zZeta(jf,1:jf-1) = zWfrac(1:jf-1)
                     ENDDO
                     !
                     DO jf = 1, nn_nfsd
                        IF( SUM(zZeta(jf,:)) > 0._wp ) zZeta(jf,:) = zZeta(jf,:) / SUM(zZeta(jf,:))
                     ENDDO
                     !
                     !==================================================!
                     ! (5) Evolve the FSD with adaptive time stepping   !
                     !==================================================!
                     !
                     ! Initialise:
                     ztelapsed = 0._wp
                     isubt     = 0
                     !
                     DO WHILE (ztelapsed < rDt_ice)
                        !
                        ! Exit loop if all ice already in smallest floe size category:
                        IF( a_ifsd(ji,jj,1,jl) >= 1._wp - epsi10 ) EXIT
                        !
                        ! Calculate Omega term in wave fracture equation:
                        DO jf = 1, nn_nfsd
                           zOmega(jf) = a_ifsd(ji,jj,jf,jl) * SUM(zWfrac(1:jf-1))
                        ENDDO
                        !
                        ! Calculate FSD tendency due to wave fracture
                        ! (cannot combine with loop above as we need SUM over Omega):
                        DO jf = 1, nn_nfsd
                           za_ifsd_tend(jf) = SUM( zOmega(:) * zZeta(:,jf) ) - zOmega(jf)
                        ENDDO
                        !
                        WHERE( ABS(za_ifsd_tend) < epsi10 ) za_ifsd_tend = 0._wp
                        !
                        ! Compute adaptive timestep to increment FSD in each floe size
                        ! category, and make sure we do not overshoot actual time step:
                        zdt_sub = rDt_ice_fsd( a_ifsd(ji,jj,:,jl), za_ifsd_tend(:) )
                        zdt_sub = MIN(zdt_sub, rDt_ice - ztelapsed)
                        !
                        ! Update FSD and time elapsed:
                        a_ifsd(ji,jj,:,jl) = a_ifsd(ji,jj,:,jl) + zdt_sub * za_ifsd_tend(:)
                        ztelapsed          = ztelapsed + zdt_sub
                        isubt              = isubt + 1
                        !
                        IF( isubt == isubt_max ) THEN
                           CALL ctl_warn('ice_wav_frac not converging: ',               &
                              &          'reached maximum number of adaptive time steps')
                        ENDIF
                        !
                     ENDDO ! while loop
                     !
                     ! === Corrections and normalisation ===
                     !
                     ! The implementation of the wave fracture equation may lead to an FSD
                     ! in this category that no longer sums to exactly 1; it must, by
                     ! definition (a_ifsd is the 'per ITD category' FSD), and if it does
                     ! not, that represents loss or gain of sea ice area.
                     !
                     ! Since wave fracture physically cannot/should not lead to loss of sea
                     ! ice area (at least, not directly), correct for any residual FSD here
                     ! by adding it back to the smallest floe size category (if some area
                     ! is lost) or by taking it away from the largest floe size category
                     ! that has at least that residual available, if some area has been gained.
                     !
                     zfsd_res = SUM(a_ifsd(ji,jj,:,jl)) - 1._wp
                     !
                     IF( zfsd_res <= 0._wp ) THEN
                        a_ifsd(ji,jj,1,jl) = a_ifsd(ji,jj,1,jl) + ABS(zfsd_res)
                     ELSE
                        DO jf = nn_nfsd, 1, -1
                           IF( a_ifsd(ji,jj,jf,jl) > zfsd_res) THEN
                              a_ifsd(ji,jj,jf,jl) = a_ifsd(ji,jj,jf,jl) - ABS(zfsd_res)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                     !
                     CALL fsd_cleanup( a_ifsd(ji,jj,:,jl) )
                     !
                  ENDIF ! category jl can fracture
               ENDDO ! -- sub-loop (ice thickness categories)
            ENDIF ! ----- MAXVAL(zWfrac) > 0
         ENDIF ! -------- at_i(ji,jj) > zat_i_min and spectrum > 0
      END_2D ! ---------- main loop

      CALL lbc_lnk( 'icewav', a_ifsd, 'T', 1._wp )

      ! Control:
      IF( ln_timing )   CALL timing_stop('ice_wav_frac')

   END SUBROUTINE ice_wav_frac


   SUBROUTINE wav_frac_dist( pssh, ph_i, pWfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_dist ***
      !!
      !! ** Purpose :   Calculate distribution of fractured ice lengths from a given local
      !!                sea surface height (SSH) field, n(x), and mean ice thickness, h_i
      !!
      !! ** Method  :   Sea ice is subject to strain due to flexure by the varying SSH
      !!                associated with the local wave field. In this subroutine, ice of
      !!                thickness h_i is assumed/conceptualised to cover the whole 1D sub-
      !!                domain for which the input sea surface height is defined. Externally,
      !!                in subroutine ice_wav_frac, the result is applied in proportion to
      !!                the actual sea ice concentration in each thickness category.
      !!
      !!                Ice fractures at locations where the strain, e, exceeds a critical
      !!                threshold, e_crit:
      !!
      !!                    e = 0.5 * h_i * |d^2 n/dx^2| >= e_crit
      !!
      !!                For each extremum in SSH, the nearest neighbouring extrema either side
      !!                are located. The strain is then calculated across such triplets of
      !!                extrema in SSH that are either {min., max., min.} or {max., min., max.}
      !!                using a finite differencing approximation across the triplet. Ice breaks
      !!                at the central extremum if e >= e_crit. The distances between all such
      !!                breaking points along the 1D domain (x) determines the lengths of
      !!                fractured ice. The number distribution of these lengths is then binned
      !!                into the FSD category bins, and the result is called the 'fracture
      !!                distribution', W(r).
      !!
      !!                W(r) satisfies int[ rW(r) dr ] = D/2, where the integral is over all
      !!                floe sizes, since the sum of all fracture lengths (twice their radii)
      !!                must equal the domain size. The quantity returned by this subroutine is
      !!                rW(r)/(D/2), i.e., W(r) normalised and scaled by the floe size of each
      !!                category, since W(r) always appears multiplied by r in the terms of the
      !!                FSD tendency equation implemented in subroutine ice_wav_frac.
      !!
      !!                Theory and method summarised above is based on
      !!                Horvat and Tziperman (2015) and Roach et al. (2018).
      !!
      !! ** Inputs  :   pssh(nn_ice_wav_nx1d) :   sub-gridscale SSH associated with local wave field (m)
      !!                ph_i                  :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pWfrac(nn_nfsd)       :   fracture distribution (m-1; normalised; binned into
      !!                                          FSD categories; scaled by floe size in each category)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_dist]
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_ice_wav_nx1d), INTENT(in)    ::   pssh     ! sea surface height (m)
      REAL(wp)                            , INTENT(in)    ::   ph_i     ! grid cell mean ice thickness (m)
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pWfrac   ! wave fracture distribution (m-1)
      !
      INTEGER                              ::   jx, jy, jf              ! dummy loop indices
      INTEGER                              ::   ixlo, ixhi              ! indices of x1d to locate extrema
      INTEGER                              ::   ixfrac                  ! number of fracture points along x1d
      LOGICAL , DIMENSION(nn_ice_wav_nx1d) ::   llmin, llmax, llext     ! sea surface height is a min / is a max / is an extrema
      REAL(wp)                             ::   zdx, zdxlo, zdxhi       ! distances between x1d points in finite difference computation (m)
      REAL(wp)                             ::   zstrain                 ! strain experienced by sea ice due to wave field
      REAL(wp), DIMENSION(nn_ice_wav_nx1d) ::   zxfrac                  ! distances to points along x1d at which ice fractures
      REAL(wp), DIMENSION(nn_ice_wav_nx1d) ::   zdxfrac                 ! distances between fracture points along x1d (m)
      REAL(wp)                             ::   zfrac_rad               ! floe radius of a piece of fractured ice (m)
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('wav_frac_dist')

      ! Initialisation:
      llmin(:)   = .FALSE.
      llmax(:)   = .FALSE.
      ixfrac     = 1
      zxfrac (:) = 0._wp
      zdxfrac(:) = 0._wp
      pWfrac (:) = 0._wp

      ! Find local extrema in sea surface height, defined to be minima or maxima over
      ! a 'moving window' of (2*nn_ice_wav_rmin + 1) points in the 1D subdomain x1d:
      !
      DO jx = 1 + nn_ice_wav_rmin, nn_ice_wav_nx1d - nn_ice_wav_rmin
         llmax(jx) = ( MAXLOC( pssh(jx-nn_ice_wav_rmin:jx+nn_ice_wav_rmin), DIM=1 ) == nn_ice_wav_rmin )
         llmin(jx) = ( MINLOC( pssh(jx-nn_ice_wav_rmin:jx+nn_ice_wav_rmin), DIM=1 ) == nn_ice_wav_rmin )
         llext(jx) = (llmin(jx) .OR. llmax(jx))
      ENDDO

      ! Loop over all points again to identify series of three consecutive, alternating
      ! extrema {min., max., min.} or {max., min., max.}, from which calculate strain
      ! and hence determine whether ice fractures there or not.
      !
      ! Loop start/end indices correspond to first/last possible index that could
      ! possibly be at the centre of a triplet.
      !
      DO jx = 2 + nn_ice_wav_rmin, nn_ice_wav_nx1d - nn_ice_wav_rmin - 1
         !
         ! Reset values for next loop iteration. Note: re-using local integer variables
         ! ixlo and ixhi from above; now they are the indices of x1d corresponding to the
         ! nearest extrema on either side of the current extrema being considered):
         !
         ixlo = 0
         ixhi = 0
         !
         IF( llext(jx) ) THEN
            !
            ! Identify nearest extrema on the left [such that x1d(ixlo) < x1d(jx)]:
            !
            DO jy = jx-1, 1, -1
               IF( llext(jy) ) THEN
                  ixlo = jy
                  EXIT
               ENDIF
            ENDDO
            !
            ! Identify nearest extrema on the right [such that x1d(jx) < x1d(ixhi)]:
            !
            DO jy = jx+1, nn_ice_wav_nx1d
               IF( llext(jy) ) THEN
                  ixhi = jy
                  EXIT
               ENDIF
            ENDDO
            !
            ! If we have a series of three extrema, with the central one being current jx, then both
            ! ixlo and and ixhi will have changed from 0. If they are alternating {max., min., max.}
            ! or {min., max., min.}, then calculate strain at x1d(jx) and determine if ice fractures
            ! there. If it does, append x1d(jx) to zxfrac array and increment ixfrac.
            !
            IF( (ixlo > 0) .AND. (ixhi > 0) ) THEN
               !
               IF(       ( llmax(ixlo) .AND. llmin(jx) .AND. llmax(ixhi) )   &
                  & .OR. ( llmin(ixlo) .AND. llmax(jx) .AND. llmin(ixhi) )    ) THEN
                  !
                  ! Calculate second derivative of SSH w.r.t. x at index jx
                  !
                  ! Centred finite difference for second derivative, forward and backward differences
                  ! on the first/inner derivatives, using the points x1d(ixlo) < x1d(jx) < x1d(ixhi).
                  ! Some simplifying algebra results in the calculation below, where also multiplying
                  ! by half of the mean ice thickness gives the strain at x1d(jx).
                  !
                  zdxlo = x1d(jx  ) - x1d(ixlo)
                  zdx   = x1d(ixhi) - x1d(ixlo)
                  zdxhi = x1d(ixhi) - x1d(jx  )
                  !
                  ! Note: zdx* are all strictly > 0, since ixlo <= jx - 1 and ixhi >= jx + 1
                  !
                  zstrain = ABS( .5_wp * ph_i * ( pssh(ixhi) * zdxlo - pssh(jx) * zdx + pssh(ixlo) * zdxhi)   &
                     &                          / ( zdxlo * zdx * zdxhi ) )
                  !
                  ! Only need to know whether this strain exceeds the critical strain
                  ! If it does, save it as a fracture point in array zxfrac:
                  IF( zstrain >= rn_ice_wav_ecri ) THEN
                     zxfrac(ixfrac) = x1d(jx)
                     ixfrac = ixfrac + 1
                  ENDIF
                  !
               ENDIF ! pssh(jx)  is at the centre of a triplet of alternating min/max
            ENDIF ! -- pssh(jx)  is at the centre of a triplet of extrema
         ENDIF ! ----- llext(jx) [pssh(jx) is an extrema]
      ENDDO ! -------- jx        [loop of x1d points]

      ! Now have locations of strain points, zxfrac(1:ixfrac). The distances between such points, when
      ! converted to radii and binned into floe size categories, gives the fracture histogram.
      !
      ! In loop above, index ixfrac is used to populate zxfrac(:), so now, ixfrac - 1 = number of
      ! fracture points. So, if:
      !
      !    ixfrac == 1    ==> no fractures at all
      !    ixfrac == 2    ==> 1 fracture point
      !    ixfrac >= 3    ==> at least 2 fracture points
      !
      ! We only compute fracture lengths between fracture points; the end points, x = x1d(1) = 0 and
      ! x = x1d(nn_ice_wav_nx1d), do not count as fracture points (because we can never calculate the
      ! strain there). Therefore, only proceed if ixfrac is at least 3.
      !
      IF( ixfrac >= 3 ) THEN
         !
         DO jx = 2, ixfrac
            !
            zfrac_rad = .5_wp * (zxfrac(jx) - zxfrac(jx-1))   ! factor of 0.5 ==> radius of fractured ice
            !
            ! Populate appropriate floe size category in fracture histogram, pWfrac(:)
            ! Just add 1 for now to get relative proportions in each category; scale whole thing afterwards
            !
            DO jf = 1, nn_nfsd - 1
               IF( zfrac_rad < floe_ru(jf) ) THEN
                  pWfrac(jf) = pWfrac(jf) + 1._wp
                  EXIT
               ENDIF
            ENDDO
            !
            ! Separate check for largest fractures (even if it exceeds upper bound of largest
            ! floe size category, it goes into that category anyway; note similar for very small
            ! fractures accounted for in above loop anyway):
            IF( zfrac_rad >= floe_rl(nn_nfsd) ) pWfrac(nn_nfsd) = pWfrac(nn_nfsd) + 1._wp
            !
         ENDDO

         ! Scale fracture histogram with floe size of each category and normalise:
         DO jf = 1, nn_nfsd
            pWfrac(jf) = floe_rc(jf) * pWfrac(jf)
         ENDDO
         !
         IF( SUM(pWfrac(:)) > 0._wp ) pWfrac(:) = pWfrac(:) / SUM(pWfrac(:))
         !
      ENDIF

      ! Control:
      IF( ln_timing )   CALL timing_stop('wav_frac_dist')

   END SUBROUTINE wav_frac_dist


   FUNCTION wav_spec_bret( phsw_l, pwpf_l )
      !!-------------------------------------------------------------------
      !!                *** FUNCTION wav_spec_bret ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Estimate local wave energy spectrum from local wave properties
      !!                calculated as Bretschneider spectrum
      !!
      !! ** Method  :   SB(f) = (5/16) * Hs^2 * (fp^4 / f^5) * EXP[ -(5/4)*(fp/f)^4 ]
      !!
      !!                where Hs is significant wave height (m), fp is the frequency
      !!                of peak wave energy (Hz) and f is frequency (Hz). SB(f) is the
      !!                Bretschneider wave energy spectrum (a.k.a., power spectral density)
      !!                and has units of m2.Hz-1 (i.e., m2.s).
      !!
      !! ** Inputs  :   phsw_l :   local significant wave height (m)
      !!                pwpf_l :   local peak frequency (Hz)
      !!
      !! ** Outputs :   local wave energy spectrum (Bretschneider; m2.Hz-1)
      !!
      !! ** Notes   :   The approach of using the Bretschneider formula to estimate the wave spectrum
      !!                for wave-ice interactions follows Horvat and Tziperman (2015) and Roach et al. (2018)
      !!                although the functional form here differs. The above matches Williams et al. (2013, Eq. 21)
      !!                and Bateson et al. (2020, Eq. 5), just converting from their expressions for SB(w) in terms
      !!                of angular frequency, w = 2*pi*f, using SB(w)dw = SB(w)(dw/df)df = 2*pi*SB(w)df,
      !!                hence SB(f) = 2*pi*SB(w). This form can be traced back to Bretschneider (1959).
      !!
      !! ** Invokers:   ice_wav_frac  --> [wav_spec_bret()]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=F)
      !!                wav_attn_spec --> [wav_spec_bret()]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=T)
      !!
      !! ** References
      !!    ----------
      !!    Bateson, A. W., Feltham, D. L., Schroeder, D., Hosekova, L., Ridley, J. K., & Aksenov, Y. (2020).
      !!              Impact of sea ice floe size distribution on seasonal fragmentation and melt of Arctic sea ice.
      !!              The Cryosphere, 14, 403-428.
      !!    Bretschneider, C. L. (1959).
      !!              Wave variability and wave spectra for wind-generated gravity waves.
      !!              Technical Memorandum No. 118, Beach Erosion Board, U.S. Army Corps of Engineers, Washington, DC, USA.
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model.
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!    Williams, T. D., Bennetts, L. G., Squire, V. A., Dumon, D. & Bertino, L. (2013).
      !!              Wave-ice interactions in the marginal ice zone. Part 1: Theoretical foundations.
      !!              Ocean Modelling, 71, 81-91.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), INTENT(in)           ::   phsw_l          ! local significant wave height (m)
      REAL(wp), INTENT(in)           ::   pwpf_l          ! local peak frequency (Hz)
      !
      REAL(wp), DIMENSION(nn_nwfreq) ::   wav_spec_bret   ! local wave energy spectrum (m2.Hz-1)
      !
      !!-------------------------------------------------------------------

      wav_spec_bret(:) = .3125_wp * phsw_l**2 * pwpf_l**4 * EXP( -1.25_wp * ( pwpf_l / wfreq(:) )**4 ) / wfreq(:)**5

   END FUNCTION wav_spec_bret


   SUBROUTINE ice_wav_init
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wav_init  ***
      !!
      !! ** Purpose :   Initialise ice wave impacts module.
      !!
      !! ** Method  :   Namelist read.
      !!                Check flags suitably set and stop if not.
      !!                Calculate some module constants.
      !!
      !! ** Callers :   ice_init --> [ice_wav_init]
      !!
      !! ** Note    :   Must be called after ice_fsd_init so that FSD-related
      !!                namelist flags are read and set. Note that some flag
      !!                checking cannot be done until after namelist group
      !!                namsbc_wave is read, which occurs later. Such checks
      !!                are handled in subroutine sbc_wave_init directly.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER ::   ji            ! Dummy loop index
      INTEGER ::   ierr          ! Local integer output status for allocate
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !
      !!
      NAMELIST/namwav/ ln_ice_wav     , ln_ice_wav_attn, ln_ice_wav_spec, ln_ice_wav_rand,   &
         &             nn_ice_wav_nx1d, rn_ice_wav_dx1d, nn_ice_wav_rmin, rn_ice_wav_ecri
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice, namwav)
      READ_NML_CFG(numnam_ice, namwav)
      IF(lwm) WRITE(numoni, namwav)
      !
      IF(lwp) THEN   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_wav_init: parameters for wave-ice interactions'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namwav:'
         WRITE(numout,*) '      Wave-ice interactions active or not                        ln_ice_wav = ', ln_ice_wav
         WRITE(numout,*) '         Activate wave-in-ice attenuation scheme or not     ln_ice_wav_attn = ', ln_ice_wav_attn
         WRITE(numout,*) '         Read full wave energy spectrum or not              ln_ice_wav_spec = ', ln_ice_wav_spec
         WRITE(numout,*) '         Use random phases for wave breakup or not          ln_ice_wav_rand = ', ln_ice_wav_rand
         WRITE(numout,*) '         Size of 1D subdomain for wave breakup              nn_ice_wav_nx1d = ', nn_ice_wav_nx1d
         WRITE(numout,*) '         Increment of 1D subdomain for wave breakup (m)     rn_ice_wav_dx1d = ', rn_ice_wav_dx1d
         WRITE(numout,*) '         Radius of smallest floes affected by waves (dx1d)  nn_ice_wav_rmin = ', nn_ice_wav_rmin
         WRITE(numout,*) '         Critical strain at which ice breaks due to waves   rn_ice_wav_ecri = ', rn_ice_wav_ecri
      ENDIF

      IF( ln_ice_wav ) THEN

         ! Checks on flags that do not require knowning SBC wave module flags
         ! (those are handled in subroutine sbc_wave_init).
         !
         ! Wave-ice interactions module requires both wave inputs and FSD:
         IF( .NOT. ln_wave ) CALL ctl_stop('ice_wav_init: ln_ice_wav=T but SBC wave module inactive (ln_wave=F)')
         IF( .NOT. ln_fsd  ) CALL ctl_stop('ice_wav_init: ln_ice_wav=T but FSD inactive (ln_fsd=F)')

         ! Warn if both attenuation scheme and reading of full wave spectrum selected
         ! (possible to do so, but usually spectrum comes from a coupled wave model so
         ! should not need to attenuate waves under ice as that is done in the wave model)
         !
         ! (check that spectrum is actually read in is done in sbc_wave_init)
         !
         IF( ln_ice_wav_attn .AND. ln_ice_wav_spec )   &
            &   CALL ctl_warn('ice_wav_init: ln_ice_wav_attn=T but also using spectrum (ln_ice_wav_spec=T): intentional?')

         ! Allocate and define module constants
         !
         ! Sub-gridscale domain (1D axis in direction of wave propagation) for
         ! computation of wave fracture distribution in subroutine wav_frac_dist:
         ALLOCATE( x1d(nn_ice_wav_nx1d), STAT=ierr )
         !
         IF (ierr /= 0) CALL ctl_stop('ice_wav_init: could not allocate array: x1d')
         !
         x1d(1) = 0._wp   ! x1d = 0., dx, 2*dx, ...
         !
         DO ji = 2, nn_ice_wav_nx1d
            x1d(ji) = x1d(ji-1) + rn_ice_wav_dx1d
         ENDDO

         ! Constant flag for subroutine ice_wav_frac: conditions under which it needs to
         ! compute local wave spectrum (T), otherwise it is already available in wspec (F):
         !
         ! Key point is that attenuation scheme (ln_ice_wav_attn) calculates attenuated wave spectrum
         ! So if NOT reading spectrum from file/model, need to calculate spectrum only if that scheme
         ! is also NOT activated (first line of condition below).
         !
         ! If ARE reading spectrum, spectrum is thus available whether attenuation scheme is active
         ! or not (second line of condition below).
         !
         l_frac_calc_spec =      ( (.NOT. ln_ice_wav_spec) .AND. (.NOT. ln_ice_wav_attn) )   &
            &               .OR. ( ln_ice_wav_spec )

      ENDIF

   END SUBROUTINE ice_wav_init


#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewav
