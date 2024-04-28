MODULE ldftra_split
   !!======================================================================
   !!                       ***  MODULE  ldftra_split  ***
   !! Ocean physics:  splitting out a large and small-scale field in the tracer
   !!=====================================================================
   !! History :  4.2  !  2024-04  (J. Mak) adapt code for NEMO 4.2
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_split      : filtering of fields based on diffusion
   !!   ldf_tra_split_init : initialization, namelist read, and parameters 
   !!                        control
   !!----------------------------------------------------------------------
   USE oce            ! ocean: dynamics and active tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! domain: ocean
   USE ldfslp         ! lateral physics: slope of iso-neutral surfaces
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_split        ! routine called in step module
   PUBLIC   ldf_tra_split_init   ! routine called in opa module
   
   !                                 !!** Namelist  namldf_tra_split  **
   REAL(wp) ::   rn_L                   !  filtering length-scale                               [m]
   INTEGER  ::   nn_period              !  number of time-steps between field splitting / large-scale field update
   REAL(wp) ::   rn_gamma               !  (dx dependent) preconditioning param for convergence [-]
   REAL(wp) ::   rn_it_tol              !  iteration tolerance with global supremum norm        [tracer units]
   INTEGER  ::   nn_it_max              !  max number of iteration                              [-]
   LOGICAL  ::   ln_diag                !  switch on ldftra_split diagnostics
   !
   
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: ldfslp.F90 15062 2021-06-28 11:19:48Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE ldf_tra_split( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ldf_tra_split  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Notes   :   
      !!
      !! ** Method  :   
      !!
      !! ** Action  : 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jn, jm, jt   ! dummy loop arguments
      !
      REAL(wp)                         ::   zlap, zres
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztu, ztv, zaheeu, zaheev
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z_orig, z_prev, z_now
      !
      IF( ln_timing )  CALL timing_start('ldf_tra_split')
      !
      IF( kt == nit000 .AND. lwp) THEN       !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_split : Tracer field splitting procedure'
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      
      ! only update every so often
      IF( MOD( kt-1, nn_period ) == 0 ) THEN
         !                                !==  Initialization of metric arrays used for all tracers  ==!
         DO_3D( 1, 0, 1, 0, 1, jpkm1 )
            zaheeu(ji,jj,jk) = (rn_L * rn_L) * e2_e1u(ji,jj) * umask(ji,jj,jk)
            zaheev(ji,jj,jk) = (rn_L * rn_L) * e1_e2v(ji,jj) * vmask(ji,jj,jk)
         END_3D
         !
         ! 1)                             !==  wslpi ==!
         ! jm: probably should make this a subroutine, being lazy here
         ! make a copy of the field to be filtered: 
         z_orig(:,:,:) = wslpi(:,:,:)
         !
         DO jm = 1, 2
            !
            z_prev(:,:,:) = z_orig(:,:,:)
            z_now (:,:,:) = 0._wp
            !
            ! initialise stop check params for iteration
            jt = 0
            zres = 1.e3
            !
            DO WHILE ( ( zres > rn_it_tol ) .AND. ( jt < nn_it_max ) )
               !               
               DO_3D( 1, 0, 1, 0, 1, jpkm1 )  !== First derivative (gradient)    ==!
                  ztu(ji,jj,jk) = zaheeu(ji,jj,jk) * ( z_prev(ji+1,jj  ,jk) - z_prev(ji,jj,jk) )
                  ztv(ji,jj,jk) = zaheev(ji,jj,jk) * ( z_prev(ji  ,jj+1,jk) - z_prev(ji,jj,jk) )
               END_3D
               !
               DO_3D( 0, 0, 0, 0, 1, jpkm1 )  !== Second derivative (divergence) ==! TODO: check against traldf_lap_blp
                  zlap = ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk)     &
                     &   + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)     &
                     &   ) / e1e2t(ji,jj)
                  z_now(ji,jj,jk) = (1.0 / rn_gamma) * z_orig(ji,jj,jk)                                                    &
                     &            + (1.0 / rn_gamma) * ( (rn_gamma - 1.0)  * z_prev(ji,jj,jk) + zlap ) * tmask(ji,jj,jk)
               END_3D
               !
               ! calculate global sup residual, stop if need be
               zres = MAXVAL(  ABS( z_now(:,:,:) - z_prev(:,:,:) )  )
               CALL mpp_max( 'ldftra_split', zres )
               IF ( zres > 1.e6 ) THEN 
                  CALL ctl_stop( 'ldftra_split: solver blowing up, turn up rn_gamma or turn down rn_L' )
                  IF(lwp)WRITE(numout,*)'    power = ', jm,' wslpi global sup residual = ', zres, 'it = ', jt
                  EXIT
               END IF
               !
               ! updates (halo, field, iteration)
               !
               CALL lbc_lnk( 'ldftra_split', z_now, 'T', 1.)
               z_prev(:,:,:) = z_now(:,:,:)
               jt = jt + 1
               !
            END DO
            !
            ! end of iteration, overwrite the orig field
            z_orig(:,:,:) = z_now(:,:,:)
            !
            ! extra diagnostics
            IF(lwp .AND. ln_diag)WRITE(numout,*)'ldftra_split iteration: power = ', jm,' wslpi global sup residual = ', zres, 'it = ', jt
            !
            IF ( zres > rn_it_tol ) THEN
               IF(lwp)WRITE(numout,*)'ldftra_split iteration: max_it reached at wslpi but glob sup not below tol, might want to be careful...'
            END IF
            !
         END DO
         ! 
         ! output the filtered field
         wslpi_l(:,:,:) = z_now(:,:,:)
         !
         ! 2)                             !==  wslpj ==!
         ! jm: probably should make this a subroutine, being lazy here
         ! make a copy of the field to be filtered: 
         z_orig(:,:,:) = wslpj(:,:,:)
         !
         DO jm = 1, 2
            !
            z_prev(:,:,:) = z_orig(:,:,:)
            z_now (:,:,:) = 0._wp
            !
            ! initialise stop check params for iteration
            jt = 0
            zres = 1.e3
            !
            DO WHILE ( ( zres > rn_it_tol ) .AND. ( jt < nn_it_max ) )
               !               
               DO_3D( 1, 0, 1, 0, 1, jpkm1 )  !== First derivative (gradient)    ==!
                  ztu(ji,jj,jk) = zaheeu(ji,jj,jk) * ( z_prev(ji+1,jj  ,jk) - z_prev(ji,jj,jk) )
                  ztv(ji,jj,jk) = zaheev(ji,jj,jk) * ( z_prev(ji  ,jj+1,jk) - z_prev(ji,jj,jk) )
               END_3D
               !
               DO_3D( 0, 0, 0, 0, 1, jpkm1 )  !== Second derivative (divergence) ==! TODO: check against traldf_lap_blp
                  zlap = ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk)     &
                     &   + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)     &
                     &   ) / e1e2t(ji,jj)
                  z_now(ji,jj,jk) = (1.0 / rn_gamma) * z_orig(ji,jj,jk)                                                    &
                     &            + (1.0 / rn_gamma) * ( (rn_gamma - 1.0)  * z_prev(ji,jj,jk) + zlap ) * tmask(ji,jj,jk)
               END_3D
               !
               ! calculate global sup residual, stop if need be
               zres = MAXVAL(  ABS( z_now(:,:,:) - z_prev(:,:,:) )  )
               CALL mpp_max( 'ldftra_split', zres )
               IF ( zres > 1.e6 ) THEN 
                  CALL ctl_stop( 'ldftra_split: solver blowing up, turn up rn_gamma or turn down rn_L' )
                  IF(lwp)WRITE(numout,*)'    power = ', jm,' wslpj global sup residual = ', zres, 'it = ', jt
                  EXIT
               END IF
               !
               ! updates (halo, field, iteration)
               !
               CALL lbc_lnk( 'ldftra_split', z_now, 'T', 1.)
               z_prev(:,:,:) = z_now(:,:,:)
               jt = jt + 1
               !
            END DO
            !
            ! end of iteration, overwrite the orig field
            z_orig(:,:,:) = z_now(:,:,:)
            !
            ! extra diagnostics
            IF(lwp .AND. ln_diag)WRITE(numout,*)'ldftra_split iteration: power = ', jm,' wslpj global sup residual = ', zres, 'it = ', jt
            !
            IF ( zres > rn_it_tol ) THEN
               IF(lwp)WRITE(numout,*)'ldftra_split iteration: max_it reached at wslpj but glob sup not below tol, might want to be careful...'
            END IF
            !
         END DO
         ! 
         ! output the filtered field
         wslpj_l(:,:,:) = z_now(:,:,:)
         !
      END IF
      !
      IF( ln_timing )  CALL timing_stop('ldf_tra_split')
      !
   END SUBROUTINE ldf_tra_split


   SUBROUTINE ldf_tra_split_init
      
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_split_init  ***
      !!                     
      !! ** Purpose :   Initialization of the ldf_tra_split
      !!
      !! ** Method  :   Read the namldf_tra_split namelist and check the parameters
      !!              called at the first timestep (nit000)
      !!
      !! ** input   :   Namlist namldf_tra_split
      !!----------------------------------------------------------------------
      INTEGER ::   ios, inum, ierr
      INTEGER ::   ji, jj
      !!
      NAMELIST/namldf_tra_split/  rn_L, nn_period, rn_gamma, rn_it_tol, nn_it_max,   &
         &                        ln_diag
      !!----------------------------------------------------------------------
      !
      !
      READ  ( numnam_cfg, namldf_tra_split, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namldf_tra_split in configuration namelist' )
      IF(lwm) WRITE ( numond, namldf_tra_split )
      !
      IF(lwp) THEN                    !* Control print
         WRITE(numout,*)
         WRITE(numout,*)    'ldf_tra_split_init : Tracer field splitting procedure for lateral diffusion'
         WRITE(numout,*)    '~~~~~~~~~~~~'
         WRITE(numout,*)    '   Namelist namldf_tra_split : set the field splitting parameters'
         WRITE(numout,*)    '      filtering length scale          rn_L      = ', rn_L
         WRITE(numout,*)    '      steps between updating ts_large nn_period = ', nn_period
         WRITE(numout,*)    '      preconditioning factor for stabilisation  = ', rn_gamma
         WRITE(numout,*)    '      iteration tolerance (supremum norm)       = ', rn_it_tol
         WRITE(numout,*)    '      max number of iterations                  = ', nn_it_max
         WRITE(numout,*)    '      solver diagnostic flag                    = ', ln_diag
      ENDIF
      !
      !
   END SUBROUTINE ldf_tra_split_init

   !!======================================================================
END MODULE ldftra_split
