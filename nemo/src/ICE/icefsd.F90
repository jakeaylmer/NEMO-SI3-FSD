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
   USE par_ice         ! SI3 parameters
   
   USE in_out_manager  ! I/O manager (needed for lwm and lwp logicals)
   USE lib_mpp         ! MPP library (needed for read_nml_substitute.h90)
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC ::   ice_fsd_init   ! routine called by icestp.F90
   
   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   &
      floe_rad_u,    &  ! FSD categories upper bounds (floe radii in m)
      floe_area_c       ! FSD category floe area at centre of bounds (m2)
   
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   &
      ! Icepack declares these in a high-level driver module, not FSD module,
      ! and already sets the dimension to the number of FSD categories:
      ! 
      floe_rad_l,    &  ! FSD categories lower bounds (floe radii in m)
      floe_rad_c,    &  ! FSD size of floes in centre of categories (radii in m)
      floe_binwidth     ! FSD category bin widths (floe radii in m)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   &
      a_ifsd            ! FSD per ice thickness category (called modified-areal
      !                 ! FSD in Roach et al., 2018, JGR: Oceans)

   !! * Substitutions
#  include "read_nml_substitute.h90"

CONTAINS
   
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
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   &
         lims   ! floe size category limits: smallest category lower limit to
         !      ! largest category upper limit (i.e., will be size nn_nfsd+1)
      
      INTEGER ::   ierr   ! allocate status return value
      !
      !!-------------------------------------------------------------------
      
      ! Floe size category boundaries are hard-coded based on number of
      ! categories specified (here nn_nfsd) in CICE/Icepack. Below are the
      ! same limits. Note that, except when nn_nfsd = 1, increasing nn_nfsd
      ! only adds larger floe size categories (i.e., smallest floe categories
      ! are the same):
      IF (nn_nfsd == 24) THEN
         
         ALLOCATE(lims(25))
         
         lims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp, &
                   2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp, &
                   1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp, &
                   4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp, &
                   9.45812834e+02_wp, 1.34354446e+03_wp, 1.82265364e+03_wp, &
                   2.47261361e+03_wp, 3.35434988e+03_wp, 4.55051413e+03_wp, &
                   6.17323164e+03_wp, 8.37461170e+03_wp, 1.13610059e+04_wp, &
                   1.54123510e+04_wp, 2.09084095e+04_wp, 2.83643675e+04_wp, &
                   3.84791270e+04_wp /)
      
      ELSEIF (nn_nfsd == 16) THEN
         
         ALLOCATE(lims(17))
         
         lims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp, &
                   2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp, &
                   1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp, &
                   4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp, &
                   9.45812834e+02_wp, 1.34354446e+03_wp, 1.82265364e+03_wp, &
                   2.47261361e+03_wp, 3.35434988e+03_wp /)
      
      ELSEIF (nn_nfsd == 12) THEN
         
         ALLOCATE(lims(13))
         
         lims = (/ 6.65000000e-02_wp, 5.31030847e+00_wp, 1.42865861e+01_wp, &
                   2.90576686e+01_wp, 5.24122136e+01_wp, 8.78691405e+01_wp, &
                   1.39518470e+02_wp, 2.11635752e+02_wp, 3.08037274e+02_wp, &
                   4.31203059e+02_wp, 5.81277225e+02_wp, 7.55141047e+02_wp, &
                   9.45812834e+02_wp /)
      
      ELSEIF (nn_nfsd == 1) THEN
         
         ALLOCATE(lims(2))
         
         lims = (/ 6.65000000e-02_wp, 3.0e+02_wp /)
         
      ELSE
         CALL ctl_stop('fsd_init_bounds: floe size categories not defined ', &
                       'for specified value of nn_nfsd')
      ENDIF
      
      ALLOCATE(floe_rad_l(nn_nfsd),    &
         &     floe_rad_u(nn_nfsd),    &
         &     floe_rad_c(nn_nfsd),    &
         &     floe_binwidth(nn_nfsd), &
         &     floe_area_c(nn_nfsd),   &
         &     STAT=ierr)
      
      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_init_bounds: could not allocate arrays')
      ENDIF
      
      floe_rad_l = lims(1:nn_nfsd)
      floe_rad_u = lims(2:nn_nfsd+1)
      floe_rad_c = 0.5_wp * (floe_rad_u + floe_rad_l)
      
      floe_binwidth = floe_rad_u - floe_rad_l
      
      floe_area_c = 4.0_wp * rn_floeshape * floe_rad_c ** 2
      
      IF (ALLOCATED(lims)) DEALLOCATE(lims)   ! we do not need lims any more
      
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
      
   END SUBROUTINE fsd_alloc
   
   
   SUBROUTINE fsd_init
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_init  ***
      !!
      !! ** Purpose :   Allocate and initialise floe size distribution
      !!                variables.
      !!
      !! ** Method  :   Allocate arrays based in FSD category bounds.
      !! 
      !! ** Input   :   ??
      !!-------------------------------------------------------------------
      !
      REAL(wp), PARAMETER ::   &
         alpha = 2.1_wp   ! parameter from Perovich and Jones (2014) 
      
      REAL(wp) ::   &
         totfrac          ! for normalising

      INTEGER ::   &
         j_itd_cat,   &   ! dummy variable for loop over ITD categories
         j_fsd_cat,   &   ! dummy variable for loop over FSD categories
         ierr             ! ALLOCATE status return value
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
         
         totfrac = 0.0_wp
         
         ! Initial FSD is the same for each ice thickness category; calculate
         ! for first category:        
         DO j_fsd_cat = 1, nn_nfsd
            ! 
            !   ! Calculate power law FSD number distribution based on Perovich
            !   ! and Jones (2014) and convert to area fraction distribution:
            ! 
            a_ifsd(:,:,j_fsd_cat,1) = (2.0_wp * floe_rad_c(j_fsd_cat))   &
               &                         ** (-alpha - 1.0_wp)            &
               &                      * floe_area_c(j_fsd_cat)           &
               &                      * floe_binwidth(j_fsd_cat)
            
            totfrac = totfrac + a_ifsd(1,1,j_fsd_cat,1)
         ENDDO
         
         a_ifsd(:,:,:,1) = a_ifsd(:,:,:,1) / totfrac   ! normalise
         
         ! Assign same initial FSD to remaining thickness categories:
         DO j_itd_cat = 2, jpl
            a_ifsd(:,:,:,j_itd_cat) = a_ifsd(:,:,:,1)
         ENDDO
         
         IF(lwp) WRITE(numout,*) 'Initialised a_ifsd (to power law distribution)'
      
      ELSE
         
         !  Initialise FSD to zero for all categories, which allows the FSD to
         !  emerge from physical processes
            
         a_ifsd(:,:,:,:) = 0.0_wp
         
         IF(lwp) WRITE(numout,*) 'Initialised a_ifsd (0 everywhere and for all categories)'
         
      ENDIF
      
   END SUBROUTINE fsd_init
   
   
   SUBROUTINE ice_fsd_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_init   ***
      !! 
      !! ** Purpose :   Check whether FSD is to be activated, and if so carry
      !!                out initialisation of FSD, printing parameter values
      !!                to stdout.
      !! 
      !! ** Method  :   Read the namfsd namelist, call other initialisation
      !!                subroutines in the module if FSD is activated.
      !! 
      !! ** input   :   Namelist namfsd
      !!-------------------------------------------------------------------
      INTEGER ::   j_fsd_cat     ! Local loop index
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
      !
      IF(ln_fsd) THEN
         CALL fsd_initbounds
      !  
      !  ! Writing the FSD bounds in CICE/Icepack is done within its analogue
      !  ! of the fsd_init_bounds subroutine. I think it makes more sense here,
      !  ! continuing from the above printing.
      !  
         IF(lwp) THEN   ! continue control print
            DO j_fsd_cat = 1, nn_nfsd
               WRITE(numout,*) floe_rad_l(j_fsd_cat), ' < fsd Cat ', &
                  &            j_fsd_cat, ' < ', floe_rad_u(j_fsd_cat)
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
