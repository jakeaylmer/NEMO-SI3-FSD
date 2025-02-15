MODULE icefsd
   !!======================================================================
   !!                       ***  MODULE icefsd ***
   !!   sea-ice : floe size distribution
   !!======================================================================
   !! History :  5.0  !  2024     (J.R. Aylmer)         Original code based
   !!                                                   on CPOM-CICE and
   !!                                                   CICE/Icepack
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_fsd_init : namelist read
   !!----------------------------------------------------------------------
   USE par_ice          ! SI3 parameters
   USE ice1D            ! sea-ice: thermodynamics variables
   USE ice              ! sea-ice: variables

   USE in_out_manager   ! I/O manager (needed for lwm and lwp logicals)
   USE iom              ! I/O manager library (needed for iom_put)
   USE lib_mpp          ! MPP library (needed for read_nml_substitute.h90)

   IMPLICIT NONE
   PRIVATE

   PUBLIC ::   ice_fsd_init               ! routine called by ice_stp
   PUBLIC ::   ice_fsd_wri                ! routine called by ice_stp
   PUBLIC ::   ice_fsd_cleanup            ! routine called by ice_dyn_adv_pra
   PUBLIC ::   ice_fsd_partition_newice   ! routine called by ice_thd_do
   PUBLIC ::   ice_fsd_add_newice         ! routine called by ice_thd_do
   PUBLIC ::   ice_fsd_thd_evolve         ! routine called by ice_thd_d{a,o}
   PUBLIC ::   fsd_peri_dens              ! function called by ice_thd_da

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   floe_rl   !: FSD floe radii, lower bounds of categories (m)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   floe_rc   !: FSD floe radii, centre       of categories (m)
   REAL(wp),         ALLOCATABLE, DIMENSION(:) ::   floe_ru   !: FSD floe radii, upper bounds of categories (m)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   floe_dr   !: FSD category widths (m)
   REAL(wp),         ALLOCATABLE, DIMENSION(:) ::   floe_ac   !: FSD floe areas, floes of radii floe_rc (m2)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   a_ifsd      !: FSD per ice thickness category
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_ifsd_2d   !: Reduced-dimension version of a_ifsd for thermodynamic routines
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   a_ifsd_1d   !: Reduced-dimension version of a_ifsd for thermodynamic routines

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"

CONTAINS


   SUBROUTINE ice_fsd_partition_newice( ki, pv_newice, pv_latgro, pda_latgro )
      !!-------------------------------------------------------------------
      !!            ***  ROUTINE ice_fsd_partition_newice  ***
      !!
      !! ** Purpose :   Partition total new ice volume into new ice formation
      !!                in open water and lateral growth of existing ice
      !!
      !! ** Method  :   Calculate lead area and lateral surface area of floes
      !!                following Horvat and Tziperman (2015):
      !!
      !!                A_lead = int [ f(r,h) * (2*r_lw/r + (r_lw/r)**2) ] dr dh
      !!                A_lat  = int [ f(r,h) * (2*h/r) ] dr dh
      !!
      !!                where A_lead = lead area (per unit ocean area)
      !!                      A_lat  = total area of vertical edges of floes
      !!                               (per unit ocean area)
      !!                      int    = integral over floe size r and thickness h
      !!                      f(r,h) = g(h)L(r,h) floe size-thickness distribution
      !!
      !!                The lead width r_lw is the annulus surrounding floes for
      !!                which freezing of existing floes occurs, distinguished
      !!                from 'open water area' in which new ice forms away from
      !!                existing ice. Hence, 'lead area' plus 'open water area'
      !!                equals one minus sea ice concentration. See Horvat and
      !!                Tziperman (2015), sect. 2.1 and Fig. 1 for details.
      !!                Following Roach et al. (2018), the smallest floe size
      !!                is used for r_lw.
      !!
      !!                The lateral growth volume (in the 'lead area') is then:
      !!
      !!                v_latgro = v_newice * A_lead / [ 1 + (at_i / A_lat) ]
      !!
      !!                where v_newice = total new ice growth
      !!                      at_i     = sea ice concentration
      !!
      !!                v_newice, initially calculated in ice_thd_do from which
      !!                this routine is called, is then updated by subtracting
      !!                v_latgro. In ice_thd_do, it is then used to add new ice
      !!                in open water (as usual since there is no lateral growth
      !!                by default, i.e., without FSD) and the total new ice
      !!                volume added, now (v_newice + v_latgro), is unchanged.
      !!
      !! ** Input   :   ki              : 1D thermodynamic array index
      !!                pv_newice       : volume of new ice to grow in total as
      !!                                  calculated in ice_thd_do, before
      !!                                  accounting for FSD (ln_fsd) and/or
      !!                                  frazil ice collection (ln_frazil)
      !!
      !! ** Output  :   pv_newice       : input updated by subtracting pv_latgro
      !!                                  (so it now represents new ice volume
      !!                                  grown in open water).
      !!                pv_latgro       : volume of ice to grow laterally on
      !!                                  existing floes in the lead area
      !!                pda_latgro(jpl) : a_i change due to lateral growth of
      !!                                  existing floes in each thickness cat.
      !!
      !! ** Note    :   no updates to ice concentration, volume, or FSD
      !!                prognostic variables are made in this routine. That is
      !!                done in the routines: ice_thd_do, ice_fsd_thd_evolve,
      !!                and ice_fsd_add_newice.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size distribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER                 , INTENT(in)    ::   ki           ! 1-D thermodynamic array index
      REAL(wp)                , INTENT(inout) ::   pv_newice    ! total new ice volume (from ice_thd_do)
      REAL(wp)                , INTENT(out)   ::   pv_latgro    ! lateral growth volume
      REAL(wp), DIMENSION(jpl), INTENT(out)   ::   pda_latgro   ! a_i change due to lateral growth
      !
      INTEGER  ::   jl, jf        ! dummy loop indices
      REAL(wp) ::   za_lead       ! lead area for open water growth (per unit ocean area)
      REAL(wp) ::   za_lat_surf   ! lateral surface area of floes (per unit ocean area)
      REAL(wp) ::   zh_i          ! ice thickness (m)
      REAL(wp) ::   zr_lw         ! width of lead region (m)
      !
      !!-------------------------------------------------------------------

      pv_latgro     = 0._wp   ! initialise
      pda_latgro(:) = 0._wp
      za_lead       = 0._wp
      za_lat_surf   = 0._wp
      zr_lw         = floe_rc(1)   ! smallest floe size for width of lead region

      ! --- Calculate za_lead and za_lat_surf (integrate/sum over thickness
      !     and floe size cats.):
      DO jl = 1, jpl

         ! need ice thickness (m) for za_lat_surf:
         IF ( a_i_2d(ki,jl) > 0._wp ) THEN
            zh_i = v_i_2d(ki,jl) / a_i_2d(ki,jl)
         ELSE
            zh_i = 0._wp
         ENDIF

         DO jf = 1, nn_nfsd
            !
            za_lead = za_lead + a_i_2d(ki,jl) * a_ifsd_2d(ki,jf,jl)   &
               &                * (2._wp * zr_lw / floe_rc(jf)        &
               &                   + zr_lw**2 / floe_rc(jf)**2)
            !
            za_lat_surf = za_lat_surf + a_ifsd_2d(ki,jf,jl) * a_i_2d(ki,jl)   &
               &                        * 2._wp * zh_i / floe_rc(jf)
            !
         ENDDO
      ENDDO

      ! --- Lead area cannot exceed open water fraction and must be > 0:
      za_lead = MAX( 0._wp, MIN( za_lead, 1._wp - at_i_1d(ki) ) )

      ! --- Calculate lateral growth volume, zv_latgro, and the change in a_i
      !     due to lateral growth, or leave both as 0 if no lateral growth:
      IF (za_lat_surf > epsi10) THEN

         pv_latgro = pv_newice * za_lead / (1._wp + at_i_1d(ki) / za_lat_surf)

         DO jl = 1, jpl
            DO jf = 1, nn_nfsd

               ! note lateral growth rate = zv_latgro / rDt_ice, but here we
               ! calculate the growth over time step which is just zv_latgro:
               pda_latgro(jl) = pda_latgro(jl)                                  &
                  &             + 2._wp * a_i_2d(ki,jl) * a_ifsd_2d(ki,jf,jl)   &
                  &                     * pv_latgro / floe_rc(jf)

            ENDDO
         ENDDO

         IF ( SUM(pda_latgro) >= za_lead ) THEN
            ! --- Cannot expand ice laterally beyond the lead region
            !     so normalise net lateral area growth to equal lead area:
            pda_latgro(:) = pda_latgro(:) / SUM(pda_latgro)
            pda_latgro(:) = pda_latgro(:) * za_lead
         ENDIF

      ENDIF

      ! --- Update volume of new ice to grow in open water:
      pv_newice = pv_newice - pv_latgro

   END SUBROUTINE ice_fsd_partition_newice


   SUBROUTINE ice_fsd_add_newice( ki, kl, pa_newice, pa_i_before )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_fsd_add_newice  ***
      !!
      !! ** Purpose :   Add new ice growth in open water (not lateral growth
      !!                of existing ice) to the floe size distribution in the
      !!                appropriate floe size and thickness categories
      !!
      !! ** Method  :   New ice is added to the smallest floe size category.
      !!                The floe size distribution for the specified thickness
      !!                category is updated accordingly.
      !!
      !! ** Input   :   ki          : 1D thermodynamic array index
      !!                kl          : thickness category new ice is added to
      !!                pa_newice   : area fraction of new ice formation
      !!                pa_i_before : ice concentration *after* lateral growth
      !!                              but *before* new ice growth at 1-D array
      !!                              index ki and in thickness category kl
      !!
      !! ** Note    :   This routine only updates the floe size distribution,
      !!                not ice concentration a_i, which is done in ice_thd_do
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER , INTENT(in) ::   ki            ! 1-D thermodynamic array index
      INTEGER , INTENT(in) ::   kl            ! thickness cat. new ice added to
      REAL(wp), INTENT(in) ::   pa_newice     ! area fraction of new ice
      REAL(wp), INTENT(in) ::   pa_i_before   ! a_i at 1-D index ki, in thickness
      !                                       ! cat. = kl, after lat. growth of
      !                                       ! existing ice but before addition
      !                                       ! of pa_newice
      !
      INTEGER ::   jf   ! dummy loop index
      !
      !!-------------------------------------------------------------------

      IF( pa_newice > 0._wp ) THEN
         IF( SUM(a_ifsd_2d(ki,:,kl)) > epsi10 ) THEN
            !
            ! --- Add new ice to smallest floe size category
            !
            ! The area fraction of ice in the smallest floe size category, r0,
            ! and thickness category to which new ice is added, h, is:
            !
            !    [ L(r0,h)g(h)drdh ]_before = a_ifsd_2d(ki,1,kl) * pa_i_before
            !
            ! before addition of pa_newice. Then, after addition of new ice:
            !
            !    [ L(r0,h)g(h)drdh ]_after = [ L(r0,h)g(h)drdh ]_before + pa_newice
            !
            ! g(h) is already updated in ice_thd_do, but L(r0,h) needs updating
            ! too, achieved by rearranging the above. This is why it is necessary
            ! to pass pa_i_before to this routine rather than just using a_i_2d.
            !
            a_ifsd_2d(ki,1,kl) = (a_ifsd_2d(ki,1,kl)*pa_i_before + pa_newice)   &
               &                 / (pa_i_before + pa_newice)

            ! --- Adjust other floe size categories
            !
            ! New ice area is only added to one floe size category, r0.
            ! So for the remaining floe size categories, r:
            !
            !    [ L(r,h)g(h)drdh ]_after = [ L(r,h)g(h)drdh ]_before
            !
            ! Since g(h)_before /= g(h)_after, L(r,h)_before /= L(r,h)_after.
            ! Rearranging gives L(r,h)_after and is thus updated:
            !
            DO jf = 2, nn_nfsd
               a_ifsd_2d(ki,jf,kl) = a_ifsd_2d(ki,jf,kl)*pa_i_before   &
                  &                  / (pa_i_before + pa_newice)
            ENDDO

         ELSE
            !
            ! --- Entirely new ice: put in smallest floe size category and
            !     specified thickness category:
            a_ifsd_2d(ki, 1        , kl) = 1._wp
            a_ifsd_2d(ki, 2:nn_nfsd, kl) = 0._wp

         ENDIF
      ENDIF

      CALL fsd_cleanup( a_ifsd_2d(ki,:,kl) )

   END SUBROUTINE ice_fsd_add_newice


   SUBROUTINE ice_fsd_thd_evolve( pa_ifsd, pG_r )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_fsd_thd_evolve  ***
      !!
      !! ** Purpose :   Evolve the floe size distribution subject to lateral
      !!                growth/melt
      !!
      !! ** Method  :   dL(r,h)/dt = -G_r * div_r(L) + (2/r) * G_r * L(r,h)
      !!
      !!                where L(r,h) is the floe size (r) distribution at
      !!                             thickness h
      !!                      div_r  is divergence in r-space
      !!                      G_r    is the lateral growth/melt rate, assumed
      !!                             to be independent of r and h, and G_r > 0
      !!                             implies growth
      !!
      !!                This equation is derived by Horvat and Tziperman (2015)
      !!                and adapted to the modified-areal floe size distribution
      !!                L(r,h) by Roach et al. (2018). This routine integrates
      !!                it forwards (for one thickness category) by one model
      !!                time step using adaptive time stepping (Horvat and
      !!                Tziperman, 2017). The adaptive time step is calculated
      !!                by function rDt_ice_fsd().
      !!
      !! ** Input   :   pa_ifsd(nn_nfsd) : floe size distribution at one grid
      !!                                   point and for one thickness category
      !!                pG_r             : lateral growth/melt rate in m/s
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
      !!              An emergent sea ice floe size distribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd   ! FSD at one location, one thickness cat.
      REAL(wp),                     INTENT(in)    ::   pG_r      ! lateral growth/melt rate (m/s)
      !
      REAL(wp), DIMENSION(nn_nfsd) ::   za_ifsd_tend      ! FSD tendency (left side of eq. above)
      REAL(wp), DIMENSION(nn_nfsd) ::   zdiv_fsd          ! divergence term in equation
      REAL(wp)                     ::   zfsd_cor          ! correction factor to ensure area conservation
      REAL(wp)                     ::   zdt_sub           ! adaptive time step (s)
      REAL(wp)                     ::   ztelapsed         ! time elapsed during adaptive time stepping
      INTEGER                      ::   isubt             ! to track number of adaptive time steps used
      INTEGER                      ::   jf                ! dummy loop index
      !
      INTEGER , PARAMETER          ::   isubt_max = 100   ! max. adaptive time steps before warning
      !!-------------------------------------------------------------------

      ! --- Start adaptive time stepping
      ztelapsed = 0._wp
      isubt     = 0

      DO WHILE (ztelapsed < rDt_ice)

         za_ifsd_tend(:) = 0._wp   ! initialise (or reset with loop iteration)
         zdiv_fsd    (:) = 0._wp

         ! --- Calculate the divergence term [div(L), without the -G_r factor]
         !     of the FSD thermodynamic tendency equation using the divergence
         !     theorem.
         !
         ! The divergence in floe category jf equals the net 'flux' of floes
         ! out of that category. Since only growth or melt occurs at once, the
         ! array indices are different depending on the sign of pG_r. In
         ! growth, floes move from smaller to larger floe size categories only,
         ! so the 'flux' of floes from category jf is directed into category
         ! jf+1, while for melt it is into category jf-1.
         !
         IF( pG_r > 0._wp ) THEN   ! lateral growth

            DO jf = 2, nn_nfsd-1
               zdiv_fsd(jf) = (   (pa_ifsd(jf  ) / floe_dr(jf  ) )     &
                  &             - (pa_ifsd(jf-1) / floe_dr(jf-1) ) )
            ENDDO

            ! Smallest category: no 'floe flux' from smaller category:
            zdiv_fsd(1) = pa_ifsd(1) / floe_dr(1)

            ! Largest category: no 'floe flux' leaving this category:
            zdiv_fsd(nn_nfsd) = -pa_ifsd(nn_nfsd-1) / floe_dr(nn_nfsd-1)

         ELSE   ! pG_r < 0; lateral melt
            !
            ! Note 'flux' of floes from category jf+1 goes into category jf
            ! which is a convergence in category jf, so need a minus sign
            ! for that term. But that minus sign is already provided by
            ! pG_r < 0, so need to add an extra minus sign to this and
            ! other 'flux' terms throughout so it cancels out later.
            !
            ! ToDo: may be clearer to write these fluxes with zG_r
            !       included in both growth and melt cases?
            !
            DO jf = 2, nn_nfsd-1
               zdiv_fsd(jf) = (   (pa_ifsd(jf+1) / floe_dr(jf+1) )     &
                  &             - (pa_ifsd(jf  ) / floe_dr(jf  ) ) )
            ENDDO

            ! Smallest category: no 'floe flux' leaving this category:
            zdiv_fsd(1) = pa_ifsd(2) / floe_dr(2)

            ! Largest category: no 'floe flux' from larger category:
            zdiv_fsd(nn_nfsd) = -pa_ifsd(nn_nfsd) / floe_dr(nn_nfsd)

         ENDIF

         ! --- Correction term
         !
         ! Sum over all floe size categories of the tendency equation must (in
         ! theory) be zero, because int(L dr) = 1 by definition, and so
         ! d/dt( int(L dr)) = 0. The divergence term also integrates to zero:
         ! indeed all elements of zdiv_fsd computed above cancel out when summed.
         ! Therefore, second term on RHS should sum to zero. In case of noise,
         ! which would manifest as spurious ice area, compute its integral,
         ! zfsd_cor, and subtract it from the actual tendency in each category
         ! weighted by that category's area fraction.
         !
         zfsd_cor = 2._wp * pG_r * SUM( pa_ifsd(:) / floe_rc(:) )

         ! --- Compute rate of change of FSD in each floe size category:
         DO jf = 1, nn_nfsd
            za_ifsd_tend(jf) = -pG_r * zdiv_fsd(jf)                                   &
               &               + 2._wp * pG_r * pa_ifsd(jf) * (1._wp / floe_rc(jf))   &
               &               - pa_ifsd(jf) * zfsd_cor
         ENDDO

         ! --- Compute adaptive timestep to increment FSD at this rate
         !     and make sure we do not overshoot actual time step:
         zdt_sub = rDt_ice_fsd( pa_ifsd(:), za_ifsd_tend(:) )
         zdt_sub = MIN(zdt_sub, rDt_ice - ztelapsed)

         ! --- Update FSD and elapsed time:
         pa_ifsd(:) = pa_ifsd(:) + zdt_sub * za_ifsd_tend(:)
         ztelapsed  = ztelapsed + zdt_sub
         isubt      = isubt + 1

         IF( isubt > isubt_max ) THEN
            CALL ctl_warn('ice_fsd_thd_evolve not converging: ',              &
               &          ' reached maximum number of adaptive time steps')
         ENDIF

      ENDDO

      CALL fsd_cleanup( pa_ifsd(:) )

   END SUBROUTINE ice_fsd_thd_evolve


   FUNCTION rDt_ice_fsd( pa_ifsd_init, pa_ifsd_tend )
      !!-------------------------------------------------------------------
      !!                   *** FUNCTION rDt_ice_fsd ***
      !!
      !! ** Purpose :   Calculate adaptive time step for evolving the floe
      !!                size distribution subject to lateral growth/melt
      !!
      !! ** Method  :   Calculate time step restrictions for incrementing the
      !!                current FSD at a specified rate, in each floe size
      !!                category. See Horvat and Tziperman (2017), Appendix A.
      !!
      !! ** Input   :   pa_ifsd_init(nn_nfsd) : current value of FSD
      !!                pa_ifsd_tend(nn_nfsd) : required tendency of FSD
      !!
      !! ** Output  :   rDt_ice_fsd         : maximum time step satisfying all
      !!                                      restrictions in each floe size
      !!                                      category and
      !!                                      0 < rDt_ice_fsd <= rDt_ice
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C. & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in) ::   pa_ifsd_init   ! current FSD
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in) ::   pa_ifsd_tend   ! required FSD tendency
      !
      REAL(wp), DIMENSION(nn_nfsd)             ::   zdt_restr    ! time step restrictions
      INTEGER                                  ::   jf           ! dummy loop index
      !
      REAL(wp) ::   rDt_ice_fsd   ! adaptive time step for FSD calculations
      !
      !!-------------------------------------------------------------------

      ! --- Calculate maximum possible time step in each floe category
      !     and save to zdt_restr
      !
      ! Afterwards we select the maximum possible time step, but it cannot be
      ! larger than the model time step so can safely use that as the initial/
      ! default value (in case of no tendency) of zdt_restr:
      zdt_restr(:) = rDt_ice

      DO jf = 1, nn_nfsd
         IF( pa_ifsd_tend(jf) > epsi10 ) THEN
            zdt_restr(jf) = (1._wp - pa_ifsd_init(jf)) / pa_ifsd_tend(jf)
         ENDIF
         IF( pa_ifsd_tend(jf) < -epsi10 ) THEN
            zdt_restr(jf) = pa_ifsd_init(jf) / ABS(pa_ifsd_tend(jf))
         ENDIF
      ENDDO

      rDt_ice_fsd = MIN(rDt_ice, MINVAL(zdt_restr))

   END FUNCTION rDt_ice_fsd


   FUNCTION fsd_peri_dens( pa_ifsd )
      !!-------------------------------------------------------------------
      !!                   *** ROUTINE ice_fsd_peri ***
      !!
      !! ** Purpose :   Calculate floe perimeter density from floe size
      !!                distribution
      !!
      !! ** Method  :   P = int[ (2/r) * F(r) dr ]
      !!
      !!                where F(r) = floe size distribution
      !!                      int  = integral over all floe sizes, r
      !!
      !! ** Note    :  Perimeter density is the total perimeter of an
      !!               ensemble of floes divided by the total sea ice area
      !!               (Bateson et al. 2022). Multiply result by sea ice
      !!               concentation to get floe perimeter per unit ocean
      !!               area.
      !!
      !! ** Input   :  pa_ifsd(nn_nfsd) : floe size distribution. Can be for one
      !!               ice thickness category [e.g., a_ifsd(ji,jj,:,jl) ] or for
      !!               all ice [e.g., sum over jl of a_i(ji,jj,jl) * a_ifsd(ji,jj,:,jl)].
      !!
      !! ** Output  :  Perimeter density [m.m-2]
      !!
      !! ** References
      !!    ----------
      !!    Bateson, A. W., Feltham, D. L., Schroeder, D. S., Wang, Y., Hwang, B., Ridley, J. K. & Aksenov, Y. (2022).
      !!              Sea ice floe size: its impact on pan-Arctic and local ice mass and required model complexity.
      !!              The Cryosphere, 16, 2565-2593.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in)  ::   pa_ifsd   ! floe size distribution
      REAL(wp)  :: fsd_peri_dens
      !
      INTEGER ::   jf   ! dummy loop index
      !
      !!-------------------------------------------------------------------

      fsd_peri_dens = 0._wp   ! initialise

      DO jf = 1, nn_nfsd
         fsd_peri_dens = fsd_peri_dens + 2._wp * pa_ifsd(jf) / floe_rc(jf)
      ENDDO

   END FUNCTION fsd_peri_dens


   FUNCTION fsd_leff_cat()
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_leff_cat  ***
      !!
      !! ** Purpose :   Calculates the effective floe size (diameter) per ice
      !!                thickness category.
      !!
      !! ** Method  :   Effective floe size (diameter) per ITD category is:
      !! 
      !!                   2 / integral[ (1/r) * L(r,h) * dr ]
      !!
      !!                where r is floe size (radius), L(r,h) is the modified-
      !!                areal FSTD, i.e., L(r,h)*dr = a_ifsd, and integral is
      !!                over all floe sizes (Bateson et al., 2022, Cryosphere,
      !!                doi:10.5194/tc-16-2565-2022)
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(A2D(0),jpl) ::   fsd_leff_cat   ! effective floe size in each thickness cat.
      !
      INTEGER ::   ji, jj, jl, jf   ! dummy variables for loop indices
      !
      !!-------------------------------------------------------------------

      fsd_leff_cat(:,:,:) = 0._wp   ! initial value

      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )

            ! Integral over floe categories [use radius at category centres
            ! and note that a_ifsd corresponds to L(r,h)*dr]:
            DO jf = 1, nn_nfsd
               fsd_leff_cat(ji,jj,jl) = fsd_leff_cat(ji,jj,jl)               &
                  &                     + a_ifsd(ji,jj,jf,jl) / floe_rc(jf)
            ENDDO

            ! 2.0 divided by above integral, except where integral is zero
            ! (or effectively zero):
            IF (fsd_leff_cat(ji,jj,jl) >= epsi06) THEN
               fsd_leff_cat(ji,jj,jl) = 2._wp / fsd_leff_cat(ji,jj,jl)
            ELSE
               fsd_leff_cat(ji,jj,jl) = 0._wp
            ENDIF

         END_2D
      ENDDO

   END FUNCTION fsd_leff_cat


   FUNCTION fsd_leff(p_leff_per_cat)
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_leff  ***
      !!
      !! ** Purpose :   Calculates the effective floe size (diameter) per
      !!                grid cell
      !!
      !! ** Method  :   Area-weighted average of effective floe size per
      !!                ice thickness category (calculated by the function
      !!                fsd_leff_cat(), the result of which should be input
      !!                to the present function). See Bateson et al. (2022,
      !!                Cryosphere, doi:10.5194/tc-16-2565-2022).
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(A2D(0),jpl), INTENT(in) ::   p_leff_per_cat   ! effective floe size per thickness cat.
      !
      REAL(wp), DIMENSION(A2D(0)) ::   fsd_leff   ! effective floe size (diameter, m)
      !
      INTEGER ::   ji, jj, jl     ! dummy variables for loop indices
      !
      !!-------------------------------------------------------------------

      fsd_leff(:,:) = 0._wp   ! initial value

      DO_2D( 0, 0, 0, 0 )
         ! Only calculate if ice is present (at_i > 0), otherwise leave as 0:
         IF (at_i(ji,jj) >= epsi06) THEN
            DO jl = 1, jpl
               fsd_leff(ji,jj) = fsd_leff(ji,jj)   &
                  &              + p_leff_per_cat(ji,jj,jl) * a_i(ji,jj,jl) / at_i(ji,jj)
            ENDDO
         ENDIF
      END_2D

   END FUNCTION fsd_leff


   SUBROUTINE ice_fsd_wri( kt )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_fsd_wri  ***
      !!
      !! ** Purpose :   Writes output fields related to the FSD.
      !!
      !! ** Method  :   Calculates metrics if requested for output and
      !!                writes using iom routines.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt        ! time step (not used?)
      !
      REAL(wp), DIMENSION(A2D(0))     ::   zleff_t   ! effective floe size, grid cell
      REAL(wp), DIMENSION(A2D(0))     ::   zmsk00    ! 0% conc. mask, grid cell
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zleff     ! effective floe size, each ITD category
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zmsk00c   ! 0% conc. mask, each ITD category
      !
      !!-------------------------------------------------------------------

      ! Copied from ice_wri generic subroutine; thresholds for outputs:
      zmsk00 (:,:)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= epsi06  )
      zmsk00c(:,:,:) = MERGE( 1._wp, 0._wp, a_i(A2D(0),:) >= epsi06  )

      ! Effective floe size (per ITD category and/or for grid cell)
      ! 
      ! (even if only for grid cell is required, it is still necessary to
      ! first calculate it per ice thickness category)
      ! 
      IF (iom_use( 'icefsdleff' ) .or. iom_use( 'icefsdleff_cat' )) THEN

         zleff = fsd_leff_cat()

         IF (iom_use( 'icefsdleff_cat' )) THEN
             CALL iom_put( 'icefsdleff_cat' , zleff(A2D(0),:) * zmsk00c)
         ENDIF

         IF (iom_use( 'icefsdleff' )) THEN
            zleff_t = fsd_leff(zleff)
            CALL iom_put( 'icefsdleff' , zleff_t(A2D(0)) * zmsk00 )
         ENDIF

      ENDIF

   END SUBROUTINE ice_fsd_wri


   SUBROUTINE fsd_cleanup(pa_ifsd_jl)
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE fsd_cleanup  ***
      !!
      !! ** Purpose :   Remove small/negative values and re-normalise floe
      !!                size distribution
      !! ** Input   :   FSD for one grid cell and one thickness category
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd_jl
      !
      REAL(wp) ::   ztotfrac   ! for normalisation
      INTEGER  ::   jf         ! dummy loop index
      !
      !!-------------------------------------------------------------------

      ! Remove negative and/or very small values in each FSD category.
      ! (note: icevar.F90 subroutine ice_var_zapsmall uses epsi10 = 1e-10)
      DO jf = 1, nn_nfsd
         IF (pa_ifsd_jl(jf) <= epsi10) THEN
            pa_ifsd_jl(jf) = 0._wp
         ENDIF
      ENDDO

      ! Compute total ice fraction for this grid cell and thickness category:
      ztotfrac = SUM(pa_ifsd_jl(:))

      IF (ztotfrac >= epsi10) THEN
         ! Ensure normalisation:
         DO jf = 1, nn_nfsd
            pa_ifsd_jl(jf) = pa_ifsd_jl(jf) / ztotfrac
         ENDDO
      ELSE
         ! Assume an ice-free grid cell and set to exactly zero:
         pa_ifsd_jl(:) = 0._wp
      ENDIF

   END SUBROUTINE fsd_cleanup


   SUBROUTINE ice_fsd_cleanup(pa_ifsd)
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_fsd_cleanup  ***
      !!
      !! ** Purpose :   Remove small/negative values and re-normalise floe
      !!                size distribution (wrapper of fsd_cleanup intended
      !!                for main FSD variable a_ifsd or intermediate variables
      !!                used in other calculations such as advection).
      !! 
      !! ** Input   :   4-D array representing FSD for all grid cells (2-D),
      !!                all thickness categories.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj,nn_nfsd,jpl), INTENT(inout) ::   pa_ifsd
      !
      INTEGER ::   ji, jj, jl   ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            CALL fsd_cleanup( pa_ifsd(ji,jj,:,jl) )
         END_2D
      ENDDO

   END SUBROUTINE ice_fsd_cleanup


   SUBROUTINE fsd_initbounds
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_init_bounds  ***
      !!
      !! ** Purpose :   Creates the FSD category boundaries and related arrays
      !!
      !! ** Method  :   Determines FSD category boundaries from namelist
      !!                parameter nn_nfsd. The actual boundaries are hard-
      !!                coded for specific values of nn_nfsd = 24, 16, 12, or
      !!                1, since this is how it is done in CICE/Icepack. But
      !!                this could be changed in the future.
      !!-------------------------------------------------------------------
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zlims   ! floe size category limits
      !
      INTEGER ::   ierr   ! allocate status return value
      !
      !!-------------------------------------------------------------------

      ! Floe size category boundaries are hard-coded based on number of
      ! categories specified (here nn_nfsd) in CICE/Icepack. Below are the
      ! same limits. Note that, except when nn_nfsd = 1, increasing nn_nfsd
      ! only adds larger floe size categories (i.e., smallest floe categories
      ! are the same):
      IF (nn_nfsd == 24) THEN

         ALLOCATE(zlims(25))

         zlims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp,   &
            &       2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp,   &
            &       1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp,   &
            &       4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp,   &
            &       9.45812834e+02_wp, 1.34354446e+03_wp, 1.82265364e+03_wp,   &
            &       2.47261361e+03_wp, 3.35434988e+03_wp, 4.55051413e+03_wp,   &
            &       6.17323164e+03_wp, 8.37461170e+03_wp, 1.13610059e+04_wp,   &
            &       1.54123510e+04_wp, 2.09084095e+04_wp, 2.83643675e+04_wp,   &
            &       3.84791270e+04_wp /)

      ELSEIF (nn_nfsd == 16) THEN

         ALLOCATE(zlims(17))

         zlims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp,   &
            &       2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp,   &
            &       1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp,   &
            &       4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp,   &
            &       9.45812834e+02_wp, 1.34354446e+03_wp, 1.82265364e+03_wp,   &
            &       2.47261361e+03_wp, 3.35434988e+03_wp /)

      ELSEIF (nn_nfsd == 12) THEN

         ALLOCATE(zlims(13))

         zlims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp,   &
            &       2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp,   &
            &       1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp,   &
            &       4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp,   &
            &       9.45812834e+02_wp /)

      ELSEIF (nn_nfsd == 1) THEN

         ALLOCATE(zlims(2))

         zlims = (/ 6.65000000e-02_wp, 3.0e+02_wp /)

      ELSE
         CALL ctl_stop('fsd_init_bounds: floe size categories not defined ',   &
            &          'for specified value of nn_nfsd')
      ENDIF

      ALLOCATE(floe_rl(nn_nfsd),   &
         &     floe_ru(nn_nfsd),   &
         &     floe_rc(nn_nfsd),   &
         &     floe_dr(nn_nfsd),   &
         &     floe_ac(nn_nfsd),   &
         &     STAT=ierr)

      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_init_bounds: could not allocate arrays')
      ENDIF

      floe_rl = zlims(1:nn_nfsd)
      floe_ru = zlims(2:nn_nfsd+1)
      floe_rc = 0.5_wp * (floe_ru + floe_rl)

      floe_dr = floe_ru - floe_rl

      floe_ac = 4._wp * rn_floeshape * floe_rc ** 2

      IF (ALLOCATED(zlims)) DEALLOCATE(zlims)

   END SUBROUTINE fsd_initbounds


   SUBROUTINE fsd_alloc
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE fsd_alloc ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Allocate floe size distribution variables
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER ::   ierr   ! ALLOCATE status return value
      !
      !!-------------------------------------------------------------------

      ALLOCATE(a_ifsd(jpi, jpj, nn_nfsd, jpl), STAT=ierr)

      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_alloc: could not allocate FSD array (a_ifsd)')
      ENDIF

      ! Allocate reduced-dimensions versions for thermodynamics.
      !
      ! Note: for other SI3 variables these are allocated by ice1D_alloc() (in
      ! ice1d.F90) which is called by ice_init (in icestp.F90) only. Seems no
      ! harm to have this here (this subroutine is only called if
      ! ln_fsd = .true.) but in the future it may be better to move this
      ! allocation there (and the above allocation of a_ifsd into ice.F90).
      ALLOCATE(a_ifsd_2d(jpij, nn_nfsd, jpl), STAT=ierr)

      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_alloc: could not allocate FSD array (a_ifsd_2d)')
      ENDIF

      ALLOCATE(a_ifsd_1d(jpij, nn_nfsd), STAT=ierr)

      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_alloc: could not allocate FSD array (a_ifsd_1d)')
      ENDIF

   END SUBROUTINE fsd_alloc


   SUBROUTINE fsd_init
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_init  ***
      !!
      !! ** Purpose :   Set initial values of floe size distribution variable
      !!                a_ifsd.
      !!
      !! ** Method  :   Either initialise to zero or to a power-law
      !!                distribution depending on ln_iceini.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), PARAMETER ::   zalpha = 2.1_wp   ! parameter from Perovich and Jones (2014)
      !
      REAL(wp)            ::   ztotfrac          ! for normalising
      INTEGER             ::   jf, jl            ! dummy variables for loop indices
      !
      !!-------------------------------------------------------------------

      ! This logical probably needs to have an "and no FSD info saved"
      ! (e.g., if restarting?)...
      IF (ln_iceini) THEN

         ! Following CICE/Icepack, initialise with a power law distribution
         ! using parameters given by Perovich and Jones (2014, JGR: Oceans,
         ! 119(12), 8767-8777, doi:10.1002/2014JC010136)
         ! 
         ! Fraction of sea ice in each floe size and thickness category is
         ! the same for all grid cells (even where there is no sea ice)
         ! initially

         ztotfrac = 0._wp

         ! Initial FSD is the same for each ice thickness category; calculate
         ! for first category:        
         DO jf = 1, nn_nfsd
            ! 
            !   ! Calculate power law FSD number distribution based on Perovich
            !   ! and Jones (2014) and convert to area fraction distribution:
            ! 
            a_ifsd(:,:,jf,1) = (2._wp * floe_rc(jf)) ** (-zalpha - 1._wp)   &
               &               * floe_ac(jf) * floe_dr(jf)

            ztotfrac = ztotfrac + a_ifsd(1,1,jf,1)
         ENDDO

         a_ifsd(:,:,:,1) = a_ifsd(:,:,:,1) / ztotfrac   ! normalise

         ! Assign same initial FSD to remaining thickness categories:
         DO jl = 2, jpl
            a_ifsd(:,:,:,jl) = a_ifsd(:,:,:,1)
         ENDDO

         IF(lwp) WRITE(numout,*) 'Initialised a_ifsd (to power law distribution)'

      ELSE

         !  Initialise FSD to zero for all categories, which allows the FSD to
         !  emerge from physical processes

         a_ifsd(:,:,:,:) = 0._wp

         IF(lwp) WRITE(numout,*) 'Initialised a_ifsd (0 everywhere and for all categories)'

      ENDIF

   END SUBROUTINE fsd_init


   SUBROUTINE ice_fsd_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_init   ***
      !! 
      !! ** Purpose :   Check whether FSD is to be activated, and if so carry
      !!                out initialisation of FSD, printing parameter values
      !!                to STDOUT, and call other FSD initialisation routines
      !! 
      !! ** Method  :   Read the namfsd namelist, call other initialisation
      !!                subroutines in the module if FSD is activated.
      !! 
      !! ** input   :   Namelist namfsd
      !!-------------------------------------------------------------------
      INTEGER ::   jf            ! Local loop index for FSD categories
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namfsd/ ln_fsd, nn_nfsd, rn_floeshape
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice, namfsd)
      READ_NML_CFG(numnam_ice, namfsd)
      IF(lwm) WRITE(numoni, namfsd)
      !
      IF(lwp) THEN   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_fsd_init: ice parameters for floe size distribution'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namfsd:'
         WRITE(numout,*) '      Floe size distribution activated or not    ln_fsd = ', ln_fsd
         WRITE(numout,*) '         Number of floe size categories         nn_nfsd = ', nn_nfsd
         WRITE(numout,*) '         Floe shape parameter              rn_floeshape = ', rn_floeshape
      ENDIF

      IF(ln_fsd) THEN
         CALL fsd_initbounds

         ! Writing the FSD bounds in CICE/Icepack is done within its analogue
         ! of the fsd_init_bounds subroutine. I think it makes more sense here,
         ! continuing from the above printing.

         IF(lwp) THEN   ! continue control print
            DO jf = 1, nn_nfsd
               WRITE(numout,*) floe_rl(jf), ' < fsd Cat ', jf, ' < ', floe_ru(jf)
            ENDDO
         ENDIF

         CALL fsd_alloc  ! could come before fsd_initbounds
         CALL fsd_init   ! must come after fsd_alloc

      ENDIF

   END SUBROUTINE ice_fsd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icefsd
