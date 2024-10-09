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
   
   PUBLIC   ice_fsd_init    ! routine called by icestp.F90
   
   !! * Substitutions
#  include "read_nml_substitute.h90"

CONTAINS
   
   SUBROUTINE ice_fsd_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_init   ***
      !! 
      !! ** Purpose :   Write whether FSD is activated to stdout
      !! 
      !! ** Method  :   Read the namfsd namelist
      !! 
      !! ** input   :   Namelist namfsd
      !!-------------------------------------------------------------------
      INTEGER  ::  ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namfsd/ ln_fsd
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namfsd)
      READ_NML_CFG(numnam_ice,namfsd)
      IF(lwm) WRITE ( numoni, namfsd )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*) 
         WRITE(numout,*) 'ice_fsd_init: ice parameters for floe size distribution'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namfsd:'
         WRITE(numout,*) '      Floe size distribution activated or not    ln_fsd = ', ln_fsd
      ENDIF
      !
   END SUBROUTINE ice_fsd_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================
END MODULE icefsd
