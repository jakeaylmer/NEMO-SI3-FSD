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
   USE sbc_oce, ONLY :   ln_wave, ln_wave_spec   ! SBC: wave module
   USE sbcwave, ONLY :   hsw, wpf                ! SBC: wave variables
   USE ice               ! sea-ice: variables
   USE icefsd , ONLY :   a_ifsd, nf_newice, floe_rl   ! floe size distribution parameters/variables
   USE icefsd , ONLY :   rDt_ice_fsd                  ! floe size distribution functions/routines

   USE in_out_manager    ! I/O manager (needed for lwm and lwp logicals)
   USE iom               ! I/O manager library (needed for iom_put)
   USE lib_mpp           ! MPP library (needed for read_nml_substitute.h90)

   IMPLICIT NONE
   PRIVATE

   PUBLIC ::   ice_wav_newice   ! routine called by ice_thd_do
   PUBLIC ::   ice_wav_attn     ! routine called by sbc_wave
   PUBLIC ::   ice_wav_frac     ! routine called by ice_stp
   PUBLIC ::   ice_wav_init     ! routine called by ice_init

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   x1d   ! 1D subdomain for wave fracture

   !                                     !!* namelist (namwav) *
   LOGICAL         ::   ln_ice_wav_rand   !: Use random phases for sea surface height in wave fracture calculation
   INTEGER         ::   nn_ice_wav_nx1d   !: Size of 1D subdomain for wave fracture calculation
   REAL(wp)        ::   rn_ice_wav_dx1d   !: Increment of 1D subdomain for wave fracture calculation (meters)
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


   SUBROUTINE ice_wav_frac
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
      !!                3.   Calculate distribution W(r) of ice fracture lengths in each grid
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
      !!                     Om and Zt are calculated from the fracture distribution (step 3).
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
      !!    Horvat, C., & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      !!-------------------------------------------------------------------
   END SUBROUTINE ice_wav_frac


   SUBROUTINE wav_frac_dist( pssh, ph_i, pWfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_dist ***
      !!
      !! ** Purpose :   Calculate distribution of fractured ice lengths from a given local
      !!                sea surface height (SSH) field, n(x), and mean ice thickness, h_i
      !!
      !! ** Method  :   Sea ice is subject to strain due to flexure by the varying SSH
      !!                associated with the local wave field. Ice fractures at locations
      !!                where the strain exceeds a critical threshold, e_crit:
      !!
      !!                    e = 0.5 * h_i * |d^2 n/dx^2| > e_crit
      !!
      !!                For each extremum in SSH, the nearest neighbouring extrema either side
      !!                are located. The strain e is calculated from such triplets of extrema
      !!                in SSH that are either {min, max, min} or {max, min, max}, which
      !!                correspond to local maxima or minima in SSH. Ice breaks at locations for
      !!                which e > e_crit. The distances between all such breaking points along the
      !!                1D domain (x) determines the lengths of fractured ice. These lengths are
      !!                then binned into the FSD bins, called the 'fracture distribution'.
      !!                Theory and method follows Roach et al. (2018).
      !!
      !! ** Inputs  :   pssh(nn_ice_wav_nx1d) :   sub-gridscale SSH associated with local wave field (m)
      !!                ph_i                  :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pWfrac(nn_nfsd)       :   fracture distribution (m-1; normalised; binned into FSD categories)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_dist]
      !!
      !! ** References
      !!    ----------
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_ice_wav_nx1d), INTENT(in)  ::   pssh     ! sea surface height (m)
      REAL(wp)                            , INTENT(in)  ::   ph_i     ! grid cell mean ice thickness (m)
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(out) ::   pWfrac   ! wave fracture distribution (m-1)
      !
      !!-------------------------------------------------------------------
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
      REAL(wp),      INTENT(in) ::   phsw_l          ! local significant wave height (m)
      REAL(wp),      INTENT(in) ::   pwpf_l          ! local peak frequency (Hz)
      !
      REAL(wp), DIMENSION(25)   ::   wav_spec_bret   ! local wave energy spectrum (m2.Hz-1)
      !
      ! NOTE: dimension corresponds to number of frequency classes resolved by wave spectrum.
      ! ----- This will be a variable declared in sbcwave but is not set up yet, so the 25
      !       here is a placeholder, as functions cannot have deferred dimensions.
      !
      !!-------------------------------------------------------------------
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
      INTEGER ::   ii            ! Dummy loop index
      INTEGER ::   ierr          ! Local integer output status for allocate
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !
      !!
      NAMELIST/namwav/ ln_ice_wav     , ln_ice_wav_attn, ln_ice_wav_spec, ln_ice_wav_rand,   &
         &             nn_ice_wav_nx1d, rn_ice_wav_dx1d, rn_ice_wav_ecri
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
         DO ii = 2, nn_ice_wav_nx1d
            x1d(ii) = x1d(ii-1) + rn_ice_wav_dx1d
         ENDDO

      ENDIF

   END SUBROUTINE ice_wav_init


#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewav
