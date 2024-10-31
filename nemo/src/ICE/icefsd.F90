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
   
   REAL(wp), DIMENSION(:), ALLOCATABLE ::   &
      floe_rad_u        ! FSD categories upper bounds (floe radii in m)
   
   REAL(wp), DIMENSION(:), ALLOCATABLE, PUBLIC ::   &
      ! Icepack declares these in a high-level driver module, not FSD module,
      ! and already sets the dimension to the number of FSD categories:
      ! 
      floe_rad_l,    &  ! FSD categories lower bounds (floe radii in m)
      floe_rad_c,    &  ! FSD size of floes in centre of categories (radii in m)
      floe_binwidth     ! FSD category bin widths (floe radii in m)

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
      REAL(wp), DIMENSION(:), ALLOCATABLE ::   &
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
         &     STAT=ierr)
      
      IF (ierr /= 0) THEN
         CALL ctl_stop('fsd_init_bounds: could not allocate arrays')
      ENDIF
      
      floe_rad_l = lims(1:nn_nfsd)
      floe_rad_u = lims(2:nn_nfsd+1)
      floe_rad_c = 0.5_wp * (floe_rad_u + floe_rad_l)
      
      floe_binwidth = floe_rad_u - floe_rad_l
      
      IF (ALLOCATED(lims)) DEALLOCATE(lims)   ! we do not need lims any more
      
   END SUBROUTINE fsd_initbounds
   
   
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
      NAMELIST/namfsd/ ln_fsd, nn_nfsd
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
      
      ENDIF
      
   END SUBROUTINE ice_fsd_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================
END MODULE icefsd
