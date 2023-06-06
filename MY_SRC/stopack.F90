#define NEMO_V4
#define NEMO_V42

! Back-Portability to NEMO 3.4 is achieved through (aggressive) precompiler directives:
! Uncomment the first line for NEMO v3.4

#if defined NEMO_V34

! Change name of trend indexes
#define  jptra_xad jptra_trd_xad
#define  jptra_yad jptra_trd_yad
#define  jptra_zad jptra_trd_zad
#define  jptra_ldf jptra_trd_ldf
#define  jptra_zdf jptra_trd_zdf
#define  jptra_bbc jptra_trd_bbc
#define  jptra_bbl jptra_trd_bbl
#define  jptra_npc jptra_trd_npc
#define  jptra_dmp jptra_trd_dmp
#define  jptra_qsr jptra_trd_qsr
#define  jptra_nsr jptra_trd_nsr
#define  jptra_atf jptra_trd_atf

#define  jptra_sad  -999
#define  jptra_zdfp -999
#define  jptra_evd  -999

#define  jpdyn_hpg jpdyn_trd_hpg
#define  jpdyn_keg jpdyn_trd_keg
#define  jpdyn_rvo jpdyn_trd_rvo
#define  jpdyn_pvo jpdyn_trd_pvo
#define  jpdyn_ldf jpdyn_trd_ldf
#define  jpdyn_had jpdyn_trd_had
#define  jpdyn_zad jpdyn_trd_zad
#define  jpdyn_zdf jpdyn_trd_zdf
#define  jpdyn_spg jpdyn_trd_spg
#define  jpdyn_dat jpdyn_trd_dat
#define  jpdyn_swf jpdyn_trd_swf
#define  jpdyn_bfr jpdyn_trd_bfr

#define jpdyn_spgflt  -999
#define jpdyn_spgexp  -999
#define jpdyn_atf     -999
#define jpdyn_tau     -999
#define jpdyn_bfri    -999

! Change name of namelist units
#define numnam_ref numnam
#define numnam_cfg numnam
#define lwm        lwp   
#define numond     numout

#define wmask      tmask  

#endif

#if defined NEMO_V4
#define _LBCNAME_ 'stopack',
#define fse3w e3w_n
#define fse3t e3t_n
#define fse3u e3u_n
#define fse3v e3v_n
#define mbathy mbkt
#define jpdyn_spgexp  jpdyn_spg
#define jpdyn_spgflt  -999
#define umask_i tmask_i
#define vmask_i tmask_i
#define nstock nn_stock
#define lk_vvl .NOT.ln_linssh
#define key_dynldf_c3d
#define key_traldf_c3d
#else
#define _LBCNAME_ 
#endif

MODULE stopack
   !!======================================================================
   !!                       ***  MODULE  stopack ***
   !!        Calculate and Apply sotchastic physics perturbations
   !!======================================================================
   !! History :  1.0  !  2018-02  (A. Storto), Original SPPT code
   !!                                          for NEMO 3.6
   !!            2.0  !  2019-05  (A. Storto), upgrades and updates:
   !!                                          (SPP, SKEB and sea-ice)
   !!----------------------------------------------------------------------
   !!   
   !!   stopack       : Generate stochastic physics perturbations
   !!   
   !!                   Method
   !!                   ======
   !!                   The module allows users to activate:
   !!   		- SPPT (Stochastically perturbed parameterization
   !!			  tendencies )scheme for user-selected trends for
   !!  			  tracers, momentum and sea-ice 
   !!   		- SPP (Schastically perturbed parameters) scheme
   !!			  for some (namelist) parameters
   !!   		- SKEB (Stochastic Energy backscatter)
   !!			  backscatter energy dissipated numerically or
   !! 			  through deep convection.
   !!   
   !!   
   !!                   Acknowledgements: C3S funded ERGO project
   !!   
   !!----------------------------------------------------------------------
   USE par_kind
   USE timing          ! Timing
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
#ifdef NEMO_V34
   USE trdmod_oce
   USE restart
#else
   USE trd_oce
#endif
   USE iom             ! I/O manager library
   USE lib_mpp         ! MPP library
   USE ldfdyn
#if defined key_lim3
   USE ice             ! LIM-3 variables
#endif
#if defined key_lim2
   USE ice_2           ! LIM-2 variables
#endif
#if defined NEMO_V4
#if defined key_si3
   USE ice            ! sea-ice: variables
   USE ice1D
#endif
#endif
#if ! defined NEMO_V4
   USE wrk_nemo
#endif
   USE diaptr
   USE zdf_oce         
   USE phycst
   USE ioipsl
   USE fldread        ! read input fields

   IMPLICIT NONE
   PRIVATE

   ! Accessible routines
   PUBLIC tra_sppt_collect
   PUBLIC dyn_sppt_collect
   PUBLIC tra_sppt_apply 
   PUBLIC dyn_sppt_apply 
   PUBLIC stopack_rst
   PUBLIC stopack_init
   PUBLIC stopack_pert
#if defined key_lim2 || defined key_lim3 || defined key_si3
   PUBLIC sppt_ice
#endif
   PUBLIC spp_aht
   PUBLIC spp_ahm
   PUBLIC spp_gen
   PUBLIC skeb_comp
   PUBLIC skeb_apply
   PUBLIC dyn_tausppt

   ! Main logical switch
   LOGICAL, PUBLIC :: ln_stopack

   ! Debug option, "hidden" namelist parameter
   LOGICAL, SAVE :: ln_stopack_debug = .FALSE.

   ! Internal switches for SPPT
   LOGICAL, PUBLIC, SAVE :: ln_sppt_tra = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_dyn = .FALSE.
   LOGICAL, PUBLIC, SAVE :: ln_sppt_ice = .FALSE.

   ! SPPT Options
   INTEGER :: sppt_filter_pass, sppt_step, nn_vwei, &
   & nn_sppt_step_bound, nn_rndm_freq, nn_deftau
   REAL(wp), SAVE :: rn_sppt_tau, rn_sppt_bound, rn_distcoast, rn_sppt_stdev
   REAL(wp), SAVE :: rn_skeb_tau, rn_skeb_stdev
   REAL(wp), SAVE :: rn_spp_tau, rn_spp_stdev
   INTEGER :: skeb_filter_pass, spp_filter_pass
 
   ! SPPT Logical switches for individual tendencies
   LOGICAL :: ln_sppt_taumap, ln_stopack_restart, ln_distcoast, &
   ln_sppt_traxad, ln_sppt_trayad, ln_sppt_trazad, ln_sppt_trasad, ln_sppt_traldf, &
   ln_sppt_trazdf, ln_sppt_trazdfp,ln_sppt_traevd, ln_sppt_trabbc, ln_sppt_trabbl, &
   ln_sppt_tranpc, ln_sppt_tradmp, ln_sppt_traqsr, ln_sppt_transr, ln_sppt_traatf    
   LOGICAL :: &
   ln_sppt_dynhpg, ln_sppt_dynspg, ln_sppt_dynkeg, ln_sppt_dynrvo, ln_sppt_dynpvo, ln_sppt_dynzad,&
   ln_sppt_dynldf, ln_sppt_dynzdf, ln_sppt_dynbfr, ln_sppt_dynatf, ln_sppt_dyntau, ln_sppt_glocon
   LOGICAL, PUBLIC :: ln_sppt_tau
   LOGICAL, PUBLIC :: &
   ln_sppt_icehdf, ln_sppt_icelat, ln_sppt_icezdf

   ! Arrays
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   gauss_n_2d
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   gauss_n_2d_p
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   gauss_n_2d_k
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   coeff_bfr

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gauss_n, zdc
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gauss_nu, gauss_nv
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   spptt, sppts, spptu, spptv
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   gauss_b, sppt_tau, sppt_a, sppt_b
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   skeb_a, skeb_b
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   spp_a, spp_b
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   gauss_bp, gauss_bk
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   g2d_save
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:    ) ::   gauss_w
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zicewrk
   REAL(wp),              SAVE                   ::   flt_fac,flt_fac_k,flt_fac_p
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   spp_tau
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   skeb_tau

   INTEGER, PARAMETER, PUBLIC :: jk_spp_alb    = 1
   INTEGER, PARAMETER, PUBLIC :: jk_spp_rhg    = 2
   INTEGER, PARAMETER, PUBLIC :: jk_spp_relw   = 3
   INTEGER, PARAMETER, PUBLIC :: jk_spp_dqdt   = 4
   INTEGER, PARAMETER, PUBLIC :: jk_spp_deds   = 5
   INTEGER, PARAMETER, PUBLIC :: jk_spp_arnf   = 6
   INTEGER, PARAMETER, PUBLIC :: jk_spp_geot   = 7
   INTEGER, PARAMETER, PUBLIC :: jk_spp_qsi0   = 8
   INTEGER, PARAMETER, PUBLIC :: jk_spp_bfr    = 9
   INTEGER, PARAMETER, PUBLIC :: jk_spp_aevd   = 10 
   INTEGER, PARAMETER, PUBLIC :: jk_spp_avt    = 11
   INTEGER, PARAMETER, PUBLIC :: jk_spp_avm    = 12
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tkelc  = 13
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tkedf  = 14
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tkeds  = 15
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tkebb  = 16
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tkefr  = 17

   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahtu   = 18
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahtv   = 19
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahtw   = 20
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahtt   = 21

   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahubbl = 22
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahvbbl = 23

   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahm1   = 24
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahm2   = 25
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahm3   = 26
   INTEGER, PARAMETER, PUBLIC :: jk_spp_ahm4   = 27

   INTEGER, PARAMETER, PUBLIC :: jk_spp_blkd   = 28
   INTEGER, PARAMETER, PUBLIC :: jk_spp_blkh   = 29
   INTEGER, PARAMETER, PUBLIC :: jk_spp_blke   = 30
   INTEGER, PARAMETER, PUBLIC :: jk_spp_tdmp   = 31

   INTEGER, PARAMETER, PUBLIC :: jk_skeb_dnum  = 41
   INTEGER, PARAMETER, PUBLIC :: jk_skeb_dcon  = 42
   INTEGER, PARAMETER, PUBLIC :: jk_skeb_deke  = 43
   INTEGER, PARAMETER, PUBLIC :: jk_skeb_tot   = 44

   INTEGER, PARAMETER, PUBLIC :: jk_sppt_tem   = 51
   INTEGER, PARAMETER, PUBLIC :: jk_sppt_sal   = 52
   INTEGER, PARAMETER, PUBLIC :: jk_sppt_uvl   = 53
   INTEGER, PARAMETER, PUBLIC :: jk_sppt_vvl   = 54

   INTEGER, PARAMETER, PUBLIC :: jk_spp        = 31
   INTEGER, PARAMETER, PUBLIC :: jk_stpk_tot   = 54
   REAL(wp), SAVE :: rn_mmax( jk_stpk_tot ) = -9.E+12_wp
   REAL(wp), SAVE :: rn_mmin( jk_stpk_tot ) =  9.E+12_wp

   INTEGER, SAVE :: numdiag      = 601
   INTEGER, SAVE :: numrep        = 602
   INTEGER, SAVE :: lkt
   
   ! Randome generator seed
   INTEGER, SAVE   :: nn_stopack_seed(4)
   LOGICAL,   SAVE, PUBLIC :: ln_stopack_repr = .TRUE.
   LOGICAL,   SAVE :: ln_spp

   REAL(wp), SAVE :: rn_uv_infl = 1._wp
   REAL(wp), SAVE :: rn_ice_infl= 1._wp
   ! SPP switches
   INTEGER, PUBLIC, SAVE :: nn_spp_bfr,nn_spp_dqdt,nn_spp_dedt,nn_spp_avt,nn_spp_avm,nn_spp_qsi0,nn_spp_relw, nn_spp_arnf,nn_spp_geot,nn_spp_aevd,nn_spp_ahubbl,nn_spp_ahvbbl
   INTEGER, PUBLIC, SAVE :: nn_spp_tkelc,nn_spp_tkedf,nn_spp_tkeds,nn_spp_tkebb,nn_spp_tkefr
   INTEGER, PUBLIC, SAVE :: nn_spp_ahtu, nn_spp_ahtv, nn_spp_ahtw, nn_spp_ahtt
   INTEGER, PUBLIC, SAVE :: nn_spp_ahm1, nn_spp_ahm2, nn_spp_tdmp
   ! SPP Sea-ice switches
   INTEGER, PUBLIC, SAVE :: nn_spp_icealb,nn_spp_icestr
   INTEGER, PUBLIC, SAVE :: nn_spp_blkd,nn_spp_blkh,nn_spp_blke

   ! SPP parameters
   REAL(wp), PUBLIC, SAVE :: rn_bfr_sd, rn_dqdt_sd, rn_dedt_sd, rn_avt_sd, rn_avm_sd, rn_qsi0_sd, rn_relw_sd, rn_arnf_sd, rn_geot_sd, rn_aevd_sd, rn_ahubbl_sd, rn_ahvbbl_sd
   REAL(wp), PUBLIC, SAVE :: rn_tkelc_sd,rn_tkedf_sd,rn_tkeds_sd,rn_tkebb_sd,rn_tkefr_sd
   REAL(wp), PUBLIC, SAVE :: rn_ahtu_sd, rn_ahtv_sd, rn_ahtw_sd, rn_ahtt_sd
   REAL(wp), PUBLIC, SAVE :: rn_ahm1_sd, rn_ahm2_sd, rn_tdmp_sd
   REAL(wp), PUBLIC, SAVE :: rn_icestr_sd, rn_icealb_sd
   REAL(wp), PUBLIC, SAVE :: rn_blkd_sd, rn_blkh_sd, rn_blke_sd

   ! Internal switches for SPP
   LOGICAL,   SAVE :: ln_spp_own_gauss
   INTEGER         :: nn_spp

   LOGICAL,  SAVE :: ln_stopack_diags

   LOGICAL,  SAVE :: ln_skeb_own_gauss
   LOGICAL,  SAVE, PUBLIC :: ln_skeb, ln_skeb_tune
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: dpsiv, dpsiu
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dnum , dcon, dekeh, dekev, deke, ut, vt
   INTEGER,  SAVE :: jpri,jprj
   INTEGER,  SAVE :: nn_dcom_freq = 1, nn_uvdta
   REAL(wp), SAVE :: rn_skeb, rn_kh, rn_kc
   REAL(wp), ALLOCATABLE, DIMENSION(:,:), SAVE :: rn_kh2, rn_kc2
   REAL(wp), ALLOCATABLE, SAVE :: stun(:,:,:)
   LOGICAL,  SAVE, PUBLIC :: ln_dpsiv = .FALSE.
   LOGICAL,  SAVE, PUBLIC :: ln_dpsiu = .FALSE.
   LOGICAL,  PUBLIC :: ln_skeb_apply = .TRUE.
   REAL(wp), SAVE :: rn_beta_num, rn_beta_con, rn_beta_eke, rn_heke, rn_veke
   INTEGER,  SAVE :: nn_dconv = 1, ktun
   INTEGER,  SAVE, PUBLIC :: nn_skst  = 2
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_uc 
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_vc 

   ! Default/initial seeds
   INTEGER(KIND=i8) :: x=1234567890987654321_8
   INTEGER(KIND=i8) :: y=362436362436362436_8
   INTEGER(KIND=i8) :: z=1066149217761810_8
   INTEGER(KIND=i8) :: w=123456123456123456_8
   ! Parameters to generate real random variates
   REAL(KIND=wp), PARAMETER :: huge64=9223372036854775808.0  ! +1
   REAL(KIND=wp), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0
   ! Variables to store 2 Gaussian random numbers with current index (ig)
   INTEGER(KIND=i8), SAVE :: ig=1
   REAL(KIND=wp), SAVE :: gran1, gran2

   !! * Substitutions
#if ! defined NEMO_V4
#  include "domzgr_substitute.h90"
#endif
#if defined NEMO_V42
#include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
#define ln_ctl .FALSE.
#define nmember 1
#define rau0 rho0
#define rn_rdt rn_Dt
#define e3w_n(ji,jj,jk) e3w(ji,jj,jk,Kmm)
#else
#include "vectopt_loop_substitute.h90"
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdytra.F90 4292 2013-11-20 16:28:04Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

#if defined NEMO_V34
   SUBROUTINE tra_sppt_collect( ptrdx, ptrdy, ktrd )
#else
   SUBROUTINE tra_sppt_collect( ptrdx, ptrdy, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sppt_collect ***
      !!
      !! ** Purpose :   Collect tracer tendencies (additive)
      !!                This function is called by the tendency diagnostics 
      !!                module
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
#endif
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptrdx   ! Temperature 
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptrdy   ! Salinity
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*) 
#if defined NEMO_V34
         WRITE(numout,*) ' CALL to tra_sppt_collect : ', ktrd
#else
         WRITE(numout,*) ' CALL to tra_sppt_collect : ', ktrd, kt
#endif
         WRITE(numout,*) ' spptt MAXVAL before      : ', MAXVAL( ABS( spptt ) )
      ENDIF

      IF( ktrd .eq. jptra_xad .AND. ln_sppt_traxad ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_yad .AND. ln_sppt_trayad ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_zad .AND. ln_sppt_trazad ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_sad .AND. ln_sppt_trasad ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_ldf .AND. ln_sppt_traldf ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_zdf .AND. ln_sppt_trazdf ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_zdfp.AND. ln_sppt_trazdfp) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_evd .AND. ln_sppt_traevd ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_bbc .AND. ln_sppt_trabbc ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_bbl .AND. ln_sppt_trabbl ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_npc .AND. ln_sppt_tranpc ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_dmp .AND. ln_sppt_tradmp ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_qsr .AND. ln_sppt_traqsr ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_nsr .AND. ln_sppt_transr ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF
      IF( ktrd .eq. jptra_atf .AND. ln_sppt_traatf ) THEN
           spptt = spptt + ptrdx ; sppts = sppts + ptrdy
      ENDIF

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*) ' spptt MAXVAL after       : ', MAXVAL( ABS( spptt ) )
         WRITE(numout,*)
      ENDIF

   END SUBROUTINE tra_sppt_collect

#if defined NEMO_V34
   SUBROUTINE dyn_sppt_collect( ptrdx, ptrdy, ktrd )
#else
   SUBROUTINE dyn_sppt_collect( ptrdx, ptrdy, ktrd, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_sppt_collect ***
      !!
      !! ** Purpose :   Collect momentum tendencies (additive)
      !!                This function is called by the tendency diagnostics 
      !!                module
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
#endif
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptrdx   ! Temperature 
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptrdy   ! Salinity
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*)
#if defined NEMO_V34
         WRITE(numout,*) ' CALL to dyn_sppt_collect : ', ktrd
#else
         WRITE(numout,*) ' CALL to dyn_sppt_collect : ', ktrd, kt
#endif
         WRITE(numout,*) ' spptu MAXVAL before      : ', MAXVAL( ABS( spptu ) )
      ENDIF

      IF( ktrd .eq. jpdyn_hpg .AND. ln_sppt_dynhpg ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_spg .AND. ln_sppt_dynspg ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_spgflt .AND. ln_sppt_dynspg ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_spgexp .AND. ln_sppt_dynspg ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_keg .AND. ln_sppt_dynkeg ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_rvo .AND. ln_sppt_dynrvo ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_pvo .AND. ln_sppt_dynpvo ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_zad .AND. ln_sppt_dynzad ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_ldf .AND. ln_sppt_dynldf ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_zdf .AND. ln_sppt_dynzdf ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_bfr .AND. ln_sppt_dynbfr ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_atf .AND. ln_sppt_dynatf ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_tau .AND. ln_sppt_dyntau ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF
      IF( ktrd .eq. jpdyn_bfri.AND. ln_sppt_dynbfr ) THEN
           spptu = spptu + ptrdx ; spptv = spptv + ptrdy
      ENDIF

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*) ' spptu MAXVAL after       : ', MAXVAL( ABS( spptu ) )
         WRITE(numout,*)
      ENDIF

   END SUBROUTINE dyn_sppt_collect

   SUBROUTINE stopack_pert( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE stopack_pert ***
      !!
      !! ** Purpose :   Update perturbation every nn_rndm_freq timestep
      !!                Calculate for SPPT, eventually for SPP and SKEB
      !!                If they are activeated and with different error scales
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt
      INTEGER :: ji,jj

      IF( ln_stopack_diags ) lkt = kt

      IF( MOD( kt - 1, nn_rndm_freq ) == 0 .OR. kt == nit000 ) THEN

       IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Calculating perturbation at timestep ', kt
         WRITE(numout,*)
       ENDIF

       ! Generate red noise
       CALL gaussian_ar1_field ( gauss_n , sppt_filter_pass, gauss_w , gauss_b,&
        sppt_a, sppt_b, gauss_n_2d,flt_fac,1)

       ! Generate red noise for SPP, SKEB, if required
       IF( ln_spp_own_gauss ) CALL gaussian_ar1_field ( gauss_n , spp_filter_pass, gauss_w , gauss_bp,&
        spp_a, spp_b, gauss_n_2d_p,flt_fac_p,2)
       IF( ln_skeb_own_gauss ) CALL gaussian_ar1_field ( gauss_n ,skeb_filter_pass, gauss_w , gauss_bk,&
        skeb_a, skeb_b, gauss_n_2d_k,flt_fac_k,3)

      ! Staggering on U-,V- grids
      ! (later should account also for possibly different bottom levels)
       DO jj=1,jpj-1
         DO ji=1,jpi-1
          gauss_nu(ji,jj,:) = 0.5_wp * (gauss_n(ji,jj,:)+gauss_n(ji+1,jj,:))
          gauss_nv(ji,jj,:) = 0.5_wp * (gauss_n(ji,jj,:)+gauss_n(ji,jj+1,:))
         ENDDO
       ENDDO
       CALL lbc_lnk(_LBCNAME_ gauss_nu, 'U', -1._wp)
       CALL lbc_lnk(_LBCNAME_ gauss_nv, 'V', -1._wp)

       IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Random noise MAXVAL SPPT : ', MAXVAL( ABS( gauss_n_2d) )
         IF( ln_spp_own_gauss ) &
         & WRITE(numout,*) ' Random noise MAXVAL SPP  : ', MAXVAL( ABS( gauss_n_2d_p) )
         IF( ln_skeb_own_gauss ) &
         & WRITE(numout,*) ' Random noise MAXVAL SKEB : ', MAXVAL( ABS( gauss_n_2d_k) )
         WRITE(numout,*)
       ENDIF

      ENDIF

#ifdef key_iomput
      ! Output the perturbation field
      CALL iom_put( "sppt_ar1" , gauss_n )
      IF(ln_spp_own_gauss)  CALL iom_put( "spp_ar1"  , gauss_n_2d_p)
      IF(ln_skeb_own_gauss) CALL iom_put( "skeb_ar1" , gauss_n_2d_k)
#endif

   END SUBROUTINE stopack_pert

   SUBROUTINE dyn_tausppt( kt, tau, tau_b, tav, tav_b )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_tausppt ***
      !!
      !! ** Purpose :   Apply collinear perturbation to wind stress
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: kt   ! Start (1) or End (2) of ice routine
      REAL(wp), INTENT(INOUT), DIMENSION(jpi,jpj) :: tau, tau_b, tav, tav_b

      tau(:,:) = tau(:,:) + (tau(:,:)-tau_b(:,:))*gauss_n_2d(:,:)*umask(:,:,1)
      tav(:,:) = tav(:,:) + (tav(:,:)-tav_b(:,:))*gauss_n_2d(:,:)*vmask(:,:,1)

      IF( kt .le. 100 ) THEN
              IF(lwp) WRITE(numout,*)
              IF(lwp) WRITE(numout,*) ' dyn_tausppt at kt =',kt
              IF(lwp) WRITE(numout,*) ' MIN / MAX utau Perturbation : ',&
              & MINVAL ( (tau(:,:)-tau_b(:,:))*gauss_n_2d(:,:), MASK = ( umask(:,:,1) .GT. 0.5_wp) ), &
              & MAXVAL ( (tau(:,:)-tau_b(:,:))*gauss_n_2d(:,:), MASK = ( umask(:,:,1) .GT. 0.5_wp) )
              IF(lwp) WRITE(numout,*) ' MIN / MAX vtau Perturbation : ',&
              & MINVAL ( (tav(:,:)-tav_b(:,:))*gauss_n_2d(:,:), MASK = ( umask(:,:,1) .GT. 0.5_wp) ), &
              & MAXVAL ( (tav(:,:)-tav_b(:,:))*gauss_n_2d(:,:), MASK = ( umask(:,:,1) .GT. 0.5_wp) )
              IF(lwp) WRITE(numout,*)
      ENDIF

    END SUBROUTINE dyn_tausppt

#ifdef NEMO_V4
#ifdef key_si3
   SUBROUTINE sppt_ice( kstep, kl  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_ice ***
      !!
      !! ** Purpose :   Apply collinear perturbation to ice fields
      !!                For specific processes coded in IC3
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: kstep         ! Start (1) or End (2) of ice routine
      INTEGER, INTENT(IN) :: kl            ! Category
      INTEGER :: jmt                       ! Number of sea-ice variables (depends on LIM version, process)
      jmt=6
      IF( kstep == 1 ) THEN        ! Store the before values
        IF ( .NOT. ALLOCATED ( zicewrk) ) ALLOCATE ( zicewrk(jpi,jpj,jmt) )
        zicewrk(:,:,1) = at_i(:,:)
        zicewrk(:,:,2) = a_i (:,:,kl)
        zicewrk(:,:,3) = h_i (:,:,kl)
        zicewrk(:,:,4) = h_s (:,:,kl)
        zicewrk(:,:,5) = t_su(:,:,kl)
        zicewrk(:,:,6) = s_i (:,:,kl)
      ELSEIF ( kstep == 2 ) THEN
        at_i(:,:)     = at_i(:,:)   + rn_ice_infl * ( at_i(:,:)   - zicewrk(:,:,1) ) * gauss_n_2d(:,:)
        a_i(:,:,kl)   = a_i(:,:,kl) + rn_ice_infl * ( a_i(:,:,kl) - zicewrk(:,:,2) ) * gauss_n_2d(:,:)
        h_i(:,:,kl)   = h_i(:,:,kl) + rn_ice_infl * ( h_i(:,:,kl) - zicewrk(:,:,3) ) * gauss_n_2d(:,:)
        h_s(:,:,kl)   = h_s(:,:,kl) + rn_ice_infl * ( h_s(:,:,kl) - zicewrk(:,:,4) ) * gauss_n_2d(:,:)
        t_su(:,:,kl)  = t_su(:,:,kl)+ rn_ice_infl * ( t_su(:,:,kl)- zicewrk(:,:,5) ) * gauss_n_2d(:,:)
        s_i(:,:,kl)   = s_i(:,:,kl) + rn_ice_infl * ( s_i(:,:,kl) - zicewrk(:,:,6) ) * gauss_n_2d(:,:)
      ELSE
        CALL ctl_stop ( ' Step in sppt_ice is not valid')
      ENDIF
   END SUBROUTINE sppt_ice
#endif
#else
#if defined key_lim2
   SUBROUTINE sppt_ice( kstep, kopt, z1,z2,z3,z4,z5,z6,z7 )
#else
   SUBROUTINE sppt_ice( kstep, kopt)
#endif
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_ice ***
      !!
      !! ** Purpose :   Apply collinear perturbation to ice fields
      !!                For specific processes coded in LIM2/LIM3 
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: kstep         ! Start (1) or End (2) of ice routine
      INTEGER, INTENT(IN) :: kopt          ! Option for sea-ice routine
#if defined key_lim2
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(INOUT) :: z1,z2,z3,z4,z5,z6,z7
#endif
#if defined key_lim3 || defined key_lim2
      INTEGER :: jmt                       ! Number of sea-ice variables (depends on LIM version, process)
      INTEGER :: jm, jl, jk

#if defined key_lim3
      jmt = 3*jpl+jpl*nlay_i
      IF( kopt == 2 ) jmt=jmt+3*jpl+1
#else
      IF( kopt == 1 ) jmt=6
      IF( kopt == 2 ) jmt=7
      IF( kopt == 3 ) jmt=8
#endif
      IF( kopt == 4 ) jmt=2

      IF( kstep == 1 ) THEN        ! Store the before values
       IF ( ALLOCATED ( zicewrk) ) DEALLOCATE ( zicewrk )
       ALLOCATE ( zicewrk(jpi,jpj,jmt) )
       jm=1
#if defined key_lim3
       DO jl = 1, jpl
          zicewrk(:,:,jm) = a_i(:,:,jl) ; jm=jm+1
          zicewrk(:,:,jm) = v_i(:,:,jl) ; jm=jm+1
          zicewrk(:,:,jm) =smv_i(:,:,jl); jm=jm+1
          DO jk = 1, nlay_i
             zicewrk(:,:,jm) = e_i(:,:,jk,jl) ; jm=jm+1
          END DO
       END DO
       IF( kopt .EQ. 2 ) THEN
         DO jl = 1, jpl
           zicewrk(:,:,jm) =oa_i(:,:,jl) ; jm=jm+1
           zicewrk(:,:,jm) = e_s(:,:,1,jl) ; jm=jm+1
           zicewrk(:,:,jm) = v_s(:,:,jl) ; jm=jm+1
         ENDDO
         zicewrk(:,:,jm) =ato_i(:,:); jm=jm+1
       ENDIF
#else
       IF( kopt .EQ. 1 ) THEN
         zicewrk(:,:,jm) = frld   ; jm=jm+1
         zicewrk(:,:,jm) = hsnif  ; jm=jm+1
         zicewrk(:,:,jm) = hicif  ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,1)   ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,2)   ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,3)   ; jm=jm+1
       ENDIF
       IF( kopt .EQ. 2 ) THEN
         IF(.NOT. PRESENT(z1) ) CALL ctl_stop( 'sppt icehdf problem, step 1')
         zicewrk(:,:,jm) = z1     ; jm=jm+1
         zicewrk(:,:,jm) = z2     ; jm=jm+1
         zicewrk(:,:,jm) = z3     ; jm=jm+1
         zicewrk(:,:,jm) = z4     ; jm=jm+1
         zicewrk(:,:,jm) = z5     ; jm=jm+1
         zicewrk(:,:,jm) = z6     ; jm=jm+1
         zicewrk(:,:,jm) = z7   
       ENDIF
       IF( kopt .EQ. 3 ) THEN
         zicewrk(:,:,jm) = frld   ; jm=jm+1
         zicewrk(:,:,jm) = hsnif  ; jm=jm+1
         zicewrk(:,:,jm) = hicif  ; jm=jm+1
         zicewrk(:,:,jm) = sist   ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,1)   ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,2)   ; jm=jm+1
         zicewrk(:,:,jm) = tbif(:,:,3)   ; jm=jm+1
       ENDIF
#endif
       IF ( kopt .EQ. 4 ) THEN
         zicewrk(:,:,jm) = u_ice  ; jm=jm+1
         zicewrk(:,:,jm) = v_ice  ; jm=jm+1
       ENDIF
      ELSEIF ( kstep == 2 ) THEN   ! Add collinear perturbation
          jm=1
#if defined key_lim3
          DO jl = 1, jpl
            a_i(:,:,jl) = a_i(:,:,jl) + (a_i(:,:,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_a_i(:,:,jl), 'T', 1. )
            v_i(:,:,jl) = v_i(:,:,jl) + (v_i(:,:,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_v_i(:,:,jl), 'T', 1. )
           smv_i(:,:,jl)=smv_i(:,:,jl)+ (smv_i(:,:,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_smv_i(:,:,jl), 'T', 1. )
            DO jk = 1, nlay_i
               e_i(:,:,jk,jl) = e_i(:,:,jk,jl)+(e_i(:,:,jk,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_e_i(:,:,jk,jl), 'T', 1. )
            END DO
          END DO
          IF( kopt .EQ. 2 ) THEN
            DO jl = 1, jpl
              oa_i(:,:,jl)=oa_i(:,:,jl)+(oa_i(:,:,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_oa_i(:,:,jl), 'T', 1. )
               e_s(:,:,1,jl)= e_s(:,:,1,jl)+( e_s(:,:,1,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_e_s(:,:,1,jl), 'T', 1. )
               v_s(:,:,jl)= v_s(:,:,jl)+( v_s(:,:,jl)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_v_s(:,:,jl), 'T', 1. )
            ENDDO
            ato_i(:,:)=ato_i(:,:)+(ato_i(:,:)-zicewrk(:,:,jm))*gauss_n_2d ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ato_i(:,:), 'T', 1. )
          ENDIF
          DEALLOCATE ( zicewrk )
#else
          IF( kopt .EQ. 1 ) THEN
            frld   = frld   + gauss_n_2d * ( frld   - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_frld, 'T', 1. )
            hsnif  = hsnif  + gauss_n_2d * ( hsnif  - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_hsnif, 'T', 1. )
            hicif  = hicif  + gauss_n_2d * ( hicif  - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_hicif, 'T', 1. )
        tbif(:,:,1)=tbif(:,:,1) + gauss_n_2d * ( tbif(:,:,1) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,1), 'T', 1. )
        tbif(:,:,2)=tbif(:,:,2) + gauss_n_2d * ( tbif(:,:,2) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,2), 'T', 1. )
        tbif(:,:,3)=tbif(:,:,3) + gauss_n_2d * ( tbif(:,:,3) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,3), 'T', 1. )
          ENDIF
          IF( kopt .EQ. 2 ) THEN
            IF(.NOT. PRESENT(z1) ) CALL ctl_stop( 'sppt icehdf problem, step 2')
            z1 = z1 + gauss_n_2d * ( z1 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z1, 'T', 1. )
            z2 = z2 + gauss_n_2d * ( z2 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z2, 'T', 1. )
            z3 = z3 + gauss_n_2d * ( z3 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z3, 'T', 1. )
            z4 = z4 + gauss_n_2d * ( z4 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z4, 'T', 1. )
            z5 = z5 + gauss_n_2d * ( z5 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z5, 'T', 1. )
            z6 = z6 + gauss_n_2d * ( z6 - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_ z6, 'T', 1. )
            z7 = z7 + gauss_n_2d * ( z7 - zicewrk(:,:,jm) )           ; CALL lbc_lnk(_LBCNAME_ z7, 'T', 1. )
          ENDIF
          IF( kopt .EQ. 3 ) THEN
            frld   = frld   + gauss_n_2d * ( frld   - zicewrk(:,:,jm) ) ; jm=jm+1  ; CALL lbc_lnk(_LBCNAME_frld, 'T', 1. )
            hsnif  = hsnif  + gauss_n_2d * ( hsnif  - zicewrk(:,:,jm) ) ; jm=jm+1  ; CALL lbc_lnk(_LBCNAME_hsnif, 'T', 1. )
            hicif  = hicif  + gauss_n_2d * ( hicif  - zicewrk(:,:,jm) ) ; jm=jm+1  ; CALL lbc_lnk(_LBCNAME_hicif, 'T', 1. )
            sist   = sist   + gauss_n_2d * ( sist   - zicewrk(:,:,jm) ) ; jm=jm+1  ; CALL lbc_lnk(_LBCNAME_sist, 'T', 1. )
        tbif(:,:,1)=tbif(:,:,1) + gauss_n_2d * ( tbif(:,:,1) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,1), 'T', 1. )
        tbif(:,:,2)=tbif(:,:,2) + gauss_n_2d * ( tbif(:,:,2) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,2), 'T', 1. )
        tbif(:,:,3)=tbif(:,:,3) + gauss_n_2d * ( tbif(:,:,3) - zicewrk(:,:,jm) ) ; jm=jm+1 ; CALL lbc_lnk(_LBCNAME_tbif(:,:,3), 'T', 1. )
          ENDIF
#endif
   ! EVP or VP rheology common to LIM2-EVP and LIM3
          IF ( kopt .EQ. 4 ) THEN
            u_ice = u_ice + (u_ice - zicewrk(:,:,jm) ) * gauss_n_2d ; jm=jm+1
            v_ice = v_ice + (v_ice - zicewrk(:,:,jm) ) * gauss_n_2d ; jm=jm+1
#if defined key_lim3 || (  defined key_lim2 && ! defined key_lim2_vp )
            CALL lbc_lnk(_LBCNAME_ u_ice, 'U', -1. )
            CALL lbc_lnk(_LBCNAME_ v_ice, 'V', -1. )
#endif         
#if defined key_lim2   &&   defined key_lim2_vp
            CALL lbc_lnk(_LBCNAME_ u_ice(:,1:jpj), 'I', -1. )
            CALL lbc_lnk(_LBCNAME_ v_ice(:,1:jpj), 'I', -1. )
#endif         
          ENDIF
          DEALLOCATE ( zicewrk )
      ENDIF
#endif
   END SUBROUTINE sppt_ice
#endif

   SUBROUTINE spp_gen  ( kt, coeff, nn_type, rn_sd, kspp, klev )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_gen ***
      !!
      !! ** Purpose :   Perturbing parameters (generic function)  
      !!                Given a value of standard deviation, the 2D parameter
      !!                coeff is perturbed following an additive Normal,
      !!                a lognormal with mean the original coeff,
      !!                a lognormal with median the original coeff,
      !!		or the previous functions but with value bounded [0.1]
      !!----------------------------------------------------------------------
   INTEGER, INTENT( in ) :: kt
   REAL(wp), INTENT( inout ), DIMENSION(jpi,jpj) :: coeff
   INTEGER, INTENT( in ) ::  nn_type
   REAL(wp), INTENT( in ) :: rn_sd
   INTEGER, INTENT( in ) ::  kspp
   INTEGER, INTENT( in ), OPTIONAL ::  klev
   REAL(wp), POINTER, DIMENSION(:,:) ::   gauss
   REAL(wp) :: zsd,xme,mm
   CHARACTER (LEN=99) :: cstrng
   INTEGER :: jklev

   IF ( nn_spp .eq. 0 ) RETURN

   IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Applying SPP perturbation on at timestep: ', kt
         WRITE(numout,*) ' Parameter : ',kspp,' nn_type/rn_sdr:',nn_type, rn_sd
   ENDIF

#if defined NEMO_V4
   ALLOCATE( gauss(jpi,jpj) )
#else
   CALL wrk_alloc(jpi,jpj,gauss)
#endif

   IF( ln_spp_own_gauss ) THEN
      gauss = gauss_n_2d_p
   ELSE
      gauss = rn_spp_stdev * gauss_n_2d / rn_sppt_stdev
   ENDIF

   IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*) ' MAXVAL of perturbations : ',MAXVAL( ABS(gauss) )
         WRITE(numout,*) ' MAXVAL of param before  : ',MAXVAL( ABS(coeff) )
   ENDIF

   IF( nn_type == 1 ) THEN
       gauss = gauss * rn_sd
       coeff = coeff * ( 1._wp + gauss )
   ELSEIF ( nn_type  == 2 ) THEN
       zsd = rn_sd
       xme = -0.5_wp*zsd*zsd
       gauss = gauss * zsd + xme
       coeff = exp(gauss) * coeff
   ELSEIF ( nn_type == 3 ) THEN
       zsd = rn_sd
       xme = 0._wp
       gauss = gauss * zsd + xme
       coeff = exp(gauss) * coeff
   ELSEIF( nn_type == 4 ) THEN
       gauss = gauss * rn_sd
       coeff = coeff * ( 1._wp + gauss )
       WHERE( coeff > 1._wp ) coeff = 1._wp
       WHERE( coeff < 0._wp ) coeff = 0._wp
   ELSEIF ( nn_type  == 5 ) THEN
       zsd = rn_sd
       xme = -0.5_wp*zsd*zsd
       gauss = gauss * zsd + xme
       coeff = exp(gauss) * coeff
       WHERE( coeff > 1._wp ) coeff = 1._wp
       WHERE( coeff < 0._wp ) coeff = 0._wp
   ELSEIF ( nn_type == 6 ) THEN
       zsd = rn_sd
       xme = 0._wp
       gauss = gauss * zsd + xme
       coeff = exp(gauss) * coeff
       WHERE( coeff > 1._wp ) coeff = 1._wp
       WHERE( coeff < 0._wp ) coeff = 0._wp
   ELSE
       CALL ctl_stop( 'spp wrong option for nn_type')
   ENDIF

   IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*) ' MAXVAL of param after   : ',MAXVAL( ABS(coeff) )
         WRITE(numout,*)
   ENDIF

#ifdef key_iomput
   write(cstrng,'(A,I2.2)') 'spp_par',kspp
   IF (iom_use(TRIM(cstrng)) ) CALL iom_put( TRIM(cstrng) , coeff )
#endif

   IF( ln_stopack_diags ) THEN
     IF(PRESENT(klev)) THEN
       jklev = klev
     ELSE
       jklev = 0 
     ENDIF
     CALL spp_stats(kt,kspp,jklev,coeff)
   ENDIF

#if defined NEMO_V4
   DEALLOCATE( gauss )
#else
   CALL wrk_dealloc(jpi,jpj,gauss)
#endif

   END SUBROUTINE

   SUBROUTINE spp_stats(mt,kp,kl,rcf)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_stats ***
      !!
      !! ** Purpose :   Compute and print basic SPP statistics
      !!----------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: mt,kp,kl
   REAL(wp), INTENT(IN) :: rcf(jpi,jpj)
   REAL(wp) :: mi,ma
   CHARACTER(LEN=16) :: cstr = '                '
   SELECT CASE ( kp ) 
     CASE( jk_spp_alb   ) 
       cstr = 'ALBEDO      '
     CASE( jk_spp_rhg   ) 
       cstr = 'ICE RHEOLOGY' 
     CASE( jk_spp_relw  ) 
       cstr = 'RELATIVE WND' 
     CASE( jk_spp_dqdt  ) 
       cstr = 'SST RELAXAT.' 
     CASE( jk_spp_deds  ) 
       cstr = 'SSS RELAXAT.' 
     CASE( jk_spp_arnf  ) 
       cstr = 'RIVER MIXING' 
     CASE( jk_spp_geot  ) 
       cstr = 'GEOTHERM.FLX' 
     CASE( jk_spp_qsi0  ) 
       cstr = 'SOLAR EXTIN.' 
     CASE( jk_spp_bfr   ) 
       cstr = 'BOTTOM FRICT' 
     CASE( jk_spp_aevd  ) 
       cstr = 'EDDY VISCDIF' 
     CASE( jk_spp_avt   ) 
       cstr = 'VERT. DIFFUS'
     CASE( jk_spp_avm   ) 
       cstr = 'VERT. VISCOS'
     CASE( jk_spp_tkelc ) 
       cstr = 'TKE LANGMUIR'
     CASE( jk_spp_tkedf ) 
       cstr = 'TKE RN_EDIFF' 
     CASE( jk_spp_tkeds ) 
       cstr = 'TKE RN_EDISS' 
     CASE( jk_spp_tkebb ) 
       cstr = 'TKE RN_EBB  '
     CASE( jk_spp_tkefr ) 
       cstr = 'TKE RN_EFR  '
     CASE( jk_spp_ahtu  ) 
       cstr = 'TRALDF AHTU '
     CASE( jk_spp_ahtv  ) 
       cstr = 'TRALDF AHTV '
     CASE( jk_spp_ahtw  ) 
       cstr = 'TRALDF AHTW '
     CASE( jk_spp_ahtt  ) 
       cstr = 'TRALDF AHTT '
     CASE( jk_spp_ahubbl )
       cstr = 'BBL DIFFUSU '
     CASE( jk_spp_ahvbbl )
       cstr = 'BBL DIFFUSV '
     CASE( jk_spp_ahm1 )
       cstr = 'DYNLDF AHM1 '
     CASE( jk_spp_ahm2 )
       cstr = 'DYNLDF AHM2 '
     CASE( jk_spp_ahm3 )
       cstr = 'DYNLDF AHM3 '
     CASE( jk_spp_ahm4 )
       cstr = 'DYNLDF AHM4 '
     CASE( jk_spp_blkd )
       cstr = 'SBCBLK CD   '
     CASE( jk_spp_blkh )
       cstr = 'SBCBLK CH   '
     CASE( jk_spp_blke )
       cstr = 'SBCBLK CE   '
     CASE DEFAULT
       CALL ctl_stop('Unrecognized SPP parameter: add it or turn off diagnostics')
   END SELECT
   mi = MINVAL(rcf, MASK = (tmask(:,:,1) > 0._wp ))
   ma = MAXVAL(rcf, MASK = (tmask(:,:,1) > 0._wp ))
   IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
   IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
   IF(lwp) THEN
      IF ( kl > 0 ) write(cstr(14:16),'(I3.3)') kl
      WRITE(numdiag,9300) lkt,cstr,mi,ma
   ENDIF
9300  FORMAT(' it :', i8, ' ', A16, ' Min: ',d10.3 , ' Max: ',d10.3)
   rn_mmin ( kp ) =  MIN ( mi, rn_mmin ( kp ) )
   rn_mmax ( kp ) =  MAX ( ma, rn_mmax ( kp ) )
   END SUBROUTINE spp_stats

   SUBROUTINE stopack_report
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_stats ***
      !!
      !! ** Purpose :  Report at the end of the run
      !!----------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: jp, numrep
   CHARACTER(LEN=16) :: cstr = '                '
   REAL(wp) :: zmul

   CALL ctl_opn(numrep, 'stopack.report', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )

   WRITE(numrep,'(A)') ' REPORT from STOPACK stochastic physics package'
   WRITE(numrep,'(A)')

   IF ( nn_spp > 0 ) THEN

   WRITE(numrep,'(A)') ' SPP  part'

   DO jp=1,jk_spp
     SELECT CASE ( jp ) 
       CASE( jk_spp_alb   ) 
         cstr = 'ALBEDO      '
       CASE( jk_spp_rhg   ) 
         cstr = 'ICE RHEOLOGY' 
       CASE( jk_spp_relw  ) 
         cstr = 'RELATIVE WND' 
       CASE( jk_spp_dqdt  ) 
         cstr = 'SST RELAXAT.' 
       CASE( jk_spp_deds  ) 
         cstr = 'SSS RELAXAT.' 
       CASE( jk_spp_arnf  ) 
         cstr = 'RIVER MIXING' 
       CASE( jk_spp_geot  ) 
         cstr = 'GEOTHERM.FLX' 
       CASE( jk_spp_qsi0  ) 
         cstr = 'SOLAR EXTIN.' 
       CASE( jk_spp_bfr   ) 
         cstr = 'BOTTOM FRICT' 
       CASE( jk_spp_aevd  ) 
         cstr = 'EDDY VISCDIF' 
       CASE( jk_spp_avt   ) 
         cstr = 'VERT. DIFFUS'
       CASE( jk_spp_avm   ) 
         cstr = 'VERT. VISCOS'
       CASE( jk_spp_tkelc ) 
         cstr = 'TKE LANGMUIR'
       CASE( jk_spp_tkedf ) 
         cstr = 'TKE RN_EDIFF' 
       CASE( jk_spp_tkeds ) 
         cstr = 'TKE RN_EDISS' 
       CASE( jk_spp_tkebb ) 
         cstr = 'TKE RN_EBB  '
       CASE( jk_spp_tkefr ) 
         cstr = 'TKE RN_EFR  '
       CASE( jk_spp_ahtu  ) 
         cstr = 'TRALDF AHTU '
       CASE( jk_spp_ahtv  ) 
         cstr = 'TRALDF AHTV '
       CASE( jk_spp_ahtw  ) 
         cstr = 'TRALDF AHTW '
       CASE( jk_spp_ahtt  ) 
         cstr = 'TRALDF AHTT '
       CASE( jk_spp_ahubbl )
         cstr = 'BBL DIFFUSU '
       CASE( jk_spp_ahvbbl )
         cstr = 'BBL DIFFUSV '
       CASE( jk_spp_ahm1 )
         cstr = 'DYNLDF AHM1 '
       CASE( jk_spp_ahm2 )
         cstr = 'DYNLDF AHM2 '
       CASE( jk_spp_ahm3 )
         cstr = 'DYNLDF AHM3 '
       CASE( jk_spp_ahm4 )
         cstr = 'DYNLDF AHM4 '
       CASE( jk_spp_blkd )
         cstr = 'SBCBLK CD   '
       CASE( jk_spp_blkh )
         cstr = 'SBCBLK CH   '
       CASE( jk_spp_blke )
         cstr = 'SBCBLK CE   '
       CASE DEFAULT
         CALL ctl_stop('Unrecognized SPP parameter: add it or turn off diagnostics')
     END SELECT
     IF ( rn_mmax(jp) .GT. 0._wp ) THEN
       WRITE(numrep,'(A,A,A,D10.3)') ' Minimum of values for parameter ', trim(cstr), ' = ', rn_mmin(jp)
       WRITE(numrep,'(A,A,A,D10.3)') ' Maximum of values for parameter ', trim(cstr), ' = ', rn_mmax(jp)
     ENDIF
   ENDDO

   ENDIF

   IF ( ln_sppt_dyn .OR. ln_sppt_tra ) THEN

   IF(sppt_step .eq. 2 ) THEN
        zmul = rdt
      ELSEIF (sppt_step .eq. 0 .or. sppt_step .eq. 1 ) THEN
        zmul = 1._wp
   ENDIF
   rn_mmax(jk_sppt_tem:jk_sppt_vvl) = rn_mmax(jk_sppt_tem:jk_sppt_vvl) * zmul
   rn_mmin(jk_sppt_tem:jk_sppt_vvl) = rn_mmin(jk_sppt_tem:jk_sppt_vvl) * zmul

   WRITE(numrep,'(A)')
   WRITE(numrep,'(A)') ' SPPT part'
   WRITE(numrep,'(A,D10.3,A)') ' ( zmul =',zmul,')'

   IF ( ln_sppt_tra ) THEN

     WRITE(numrep,'(A,D10.3)') ' Minimum of values for TEM  ', rn_mmin(jk_sppt_tem)
     WRITE(numrep,'(A,D10.3)') ' Maximum of values for TEM  ', rn_mmax(jk_sppt_tem)
     IF( rn_mmax(jk_sppt_tem) > 0.5 ) WRITE(numrep,'(A)' ) ' Larger than 0.5, might be too big  '
  
     WRITE(numrep,'(A,D10.3)') ' Minimum of values for SAL  ', rn_mmin(jk_sppt_sal)
     WRITE(numrep,'(A,D10.3)') ' Maximum of values for SAL  ', rn_mmax(jk_sppt_sal)
     IF( rn_mmax(jk_sppt_sal) > 0.2 ) WRITE(numrep,'(A)' ) ' Larger than 0.2, might be too big  '

   ENDIF

   IF ( ln_sppt_dyn ) THEN

     WRITE(numrep,'(A,D10.3)') ' Minimum of values for UVEL ', rn_mmin(jk_sppt_uvl)
     WRITE(numrep,'(A,D10.3)') ' Maximum of values for UVEL ', rn_mmax(jk_sppt_uvl)
     IF( rn_mmax(jk_sppt_uvl) > 0.1 ) WRITE(numrep,'(A)' ) ' Larger than 0.1, might be too big  '
  
     WRITE(numrep,'(A,D10.3)') ' Minimum of values for VVEL ', rn_mmin(jk_sppt_vvl)
     WRITE(numrep,'(A,D10.3)') ' Maximum of values for VVEL ', rn_mmax(jk_sppt_vvl)
     IF( rn_mmax(jk_sppt_vvl) > 0.1 ) WRITE(numrep,'(A)' ) ' Larger than 0.1, might be too big  '

   ENDIF

   ENDIF

   IF ( ln_skeb ) THEN

   WRITE(numrep,'(A)')
   WRITE(numrep,'(A)') ' SKEB part'
   WRITE(numrep,'(A,D10.3)') ' Maximum of absolute values for NUM  ', rn_mmax(jk_skeb_dnum)
   WRITE(numrep,'(A)'          ) ' (Perturbation from numerical dissipation)'
   IF( rn_mmax(jk_skeb_dnum) < 0.5e-4 ) WRITE(numrep,'(A)' ) ' Smaller than 0.5e-4, might be too small'
   IF( rn_mmax(jk_skeb_dnum) > 0.2e-2 ) WRITE(numrep,'(A)' ) ' Larger  than 0.2e-2, might be too big  '
   WRITE(numrep,'(A)')
   WRITE(numrep,'(A,D10.3)') ' Maximum of absolute values for CONV ', rn_mmax(jk_skeb_dcon)
   WRITE(numrep,'(A)'          ) ' (Perturbation from convection dissipation)'
   IF( rn_mmax(jk_skeb_dcon) < 0.5e-4 ) WRITE(numrep,'(A)' ) ' Smaller than 0.5e-4, might be too small'
   IF( rn_mmax(jk_skeb_dcon) > 0.2e-2 ) WRITE(numrep,'(A)' ) ' Larger  than 0.2e-2, might be too big  '
   WRITE(numrep,'(A)')
   WRITE(numrep,'(A,D10.3)') ' Maximum of absolute values for EKE  ', rn_mmax(jk_skeb_deke)
   WRITE(numrep,'(A)'          ) ' (Perturbation from mesoscale dissipation)'
   IF( rn_mmax(jk_skeb_deke) < 0.5e-4 ) WRITE(numrep,'(A)' ) ' Smaller than 0.5e-4, might be too small'
   IF( rn_mmax(jk_skeb_deke) > 0.2e-2 ) WRITE(numrep,'(A)' ) ' Larger  than 0.2e-2, might be too big  '
   WRITE(numrep,'(A)')
   WRITE(numrep,'(A,D10.3)') ' Maximum of absolute values for TOTAL', rn_mmax(jk_skeb_tot)
   WRITE(numrep,'(A)'          ) ' (Perturbation from total energy dissipation)'
   IF( rn_mmax(jk_skeb_tot) < 0.5e-4 ) WRITE(numrep,'(A)' ) ' Smaller than 0.5e-4, might be too small'
   IF( rn_mmax(jk_skeb_tot) > 0.2e-2 ) WRITE(numrep,'(A)' ) ' Larger  than 0.2e-2, might be too big  '

   ENDIF

   CLOSE(numrep)
   CLOSE(numdiag)

   END SUBROUTINE stopack_report

   SUBROUTINE spp_aht ( kt, coeff,nn_type, rn_sd, kspp  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_aht ***
      !!
      !! ** Purpose :   Perturbing diffusivity parameters
      !!                As spp_gen, but specifically designed for diffusivity
      !!                where coeff can be 2D or 3D
      !!----------------------------------------------------------------------
   INTEGER, INTENT( in ) :: kt
#if defined key_traldf_c3d
   REAL(wp), INTENT( inout ), DIMENSION(jpi,jpj,jpk) :: coeff
#elif defined key_traldf_c2d
   REAL(wp), INTENT( inout ), DIMENSION(jpi,jpj) :: coeff
#elif defined key_traldf_c1d
   REAL(wp), INTENT( inout ), DIMENSION(jpk) :: coeff
#else
   REAL(wp), INTENT( inout ) :: coeff
#endif
   INTEGER, INTENT( in ) ::  kspp
   INTEGER, INTENT( in ) ::  nn_type
   REAL(wp), INTENT( in ) :: rn_sd
   REAL(wp), POINTER, DIMENSION(:,:) ::   gauss
   REAL(wp) :: zsd,xme
   INTEGER :: jk

   IF ( nn_spp .eq. 0 ) RETURN

#if defined NEMO_V4
   ALLOCATE( gauss(jpi,jpj) )
#else
   CALL wrk_alloc(jpi,jpj,gauss)
#endif

   IF( ln_spp_own_gauss ) THEN
      gauss = gauss_n_2d_p
   ELSE
      gauss = rn_spp_stdev * gauss_n_2d / rn_sppt_stdev
   ENDIF

   IF( nn_type == 1 ) THEN
       gauss = gauss * rn_sd
#if defined key_traldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = coeff(:,:,jk) * ( 1._wp + gauss )
       ENDDO
#elif defined key_traldf_c2d
       coeff = coeff * ( 1._wp + gauss )
#endif
   ELSEIF ( nn_type == 2 ) THEN
       zsd = rn_sd
       xme = -0.5_wp*zsd*zsd
       gauss = gauss * zsd + xme
#if defined key_traldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = exp(gauss) * coeff(:,:,jk)
       ENDDO
#elif defined key_traldf_c2d
       coeff = exp(gauss) * coeff
#endif
   ELSEIF ( nn_type == 3 ) THEN
       zsd = rn_sd
       xme = 0._wp
       gauss = gauss * zsd + xme
#if defined key_traldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = exp(gauss) * coeff(:,:,jk)
       ENDDO
#elif defined key_traldf_c2d
       coeff = exp(gauss) * coeff
#endif
   ELSE
       CALL ctl_stop( 'spp aht wrong option')
   ENDIF

   IF( ln_stopack_diags ) THEN
     jk=0
     CALL spp_stats(kt,kspp,jk,coeff)
   ENDIF

#if defined NEMO_V4
   DEALLOCATE( gauss )
#else
   CALL wrk_dealloc(jpi,jpj,gauss)
#endif

   END SUBROUTINE

   SUBROUTINE spp_ahm ( kt, coeff,nn_type, rn_sd, kspp  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE spp_ahm ***
      !!
      !! ** Purpose :   Perturbing viscosity parameters
      !!                As spp_gen, but specifically designed for viscosity
      !!                where coeff can be 2D or 3D
      !!----------------------------------------------------------------------
   INTEGER, INTENT( in ) :: kt
#if defined key_dynldf_c3d
   REAL(wp), INTENT( inout ), DIMENSION(jpi,jpj,jpk) :: coeff
#elif defined key_dynldf_c2d
   REAL(wp), INTENT( inout ), DIMENSION(jpi,jpj) :: coeff
#elif defined key_dynldf_c1d
   REAL(wp), INTENT( inout ), DIMENSION(jpk) :: coeff
#else
   REAL(wp), INTENT( inout ) :: coeff
#endif
   INTEGER, INTENT( in ) ::  kspp
   INTEGER, INTENT( in ) ::  nn_type
   REAL(wp), INTENT( in ) :: rn_sd
   REAL(wp), POINTER, DIMENSION(:,:) ::   gauss
   REAL(wp) :: zsd,xme
   INTEGER :: jk

   IF ( nn_spp .eq. 0 ) RETURN

#if defined NEMO_V4
   ALLOCATE( gauss(jpi,jpj) )
#else
   CALL wrk_alloc(jpi,jpj,gauss)
#endif

   IF( ln_spp_own_gauss ) THEN
      gauss = gauss_n_2d_p
   ELSE
      gauss = rn_spp_stdev * gauss_n_2d / rn_sppt_stdev
   ENDIF

   IF( nn_type == 1 ) THEN
       gauss = gauss * rn_sd
#if defined key_dynldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = coeff(:,:,jk) * ( 1._wp + gauss )
       ENDDO
#elif defined key_dynldf_c2d
       coeff = coeff * ( 1._wp + gauss )
#endif
   ELSEIF ( nn_type == 2 ) THEN
       zsd = rn_sd
       xme = -0.5_wp*zsd*zsd
       gauss = gauss * zsd + xme
#if defined key_dynldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = exp(gauss) * coeff(:,:,jk)
       ENDDO
#endif
#if defined key_dynldf_c2d
       coeff = exp(gauss) * coeff
#endif
   ELSEIF ( nn_type == 3 ) THEN
       zsd = rn_sd
       xme = 0._wp
       gauss = gauss * zsd + xme
#if defined key_dynldf_c3d
       DO jk=1,jpk
         coeff(:,:,jk) = exp(gauss) * coeff(:,:,jk)
       ENDDO
#elif defined key_dynldf_c2d
       coeff = exp(gauss) * coeff
#endif
   ELSE
       CALL ctl_stop( 'spp ahm wrong option')
   ENDIF

   IF( ln_stopack_diags ) THEN
     jk=0
     CALL spp_stats(kt,kspp,jk,coeff)
   ENDIF

#if defined NEMO_V4
   DEALLOCATE( gauss )
#else
   CALL wrk_dealloc(jpi,jpj,gauss)
#endif

   END SUBROUTINE

#if defined NEMO_V42
   SUBROUTINE tra_sppt_apply ( kt , ks, Kbb, Kmm, pts, Krhs )
#else
   SUBROUTINE tra_sppt_apply ( kt , ks )
#endif
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sppt_apply ***
      !!
      !! ** Purpose :   Apply perturbation to tracers
      !! 		Step 0/1 for tendencies, step 2 for variables
      !! 		after timestepping
      !!----------------------------------------------------------------------

      INTEGER, INTENT( in ) :: kt, ks     ! Main time step counter
#if defined NEMO_V42
      INTEGER, INTENT( in ) :: Kbb, Kmm, Krhs  ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
#endif
      REAL(wp) :: mi,ma
#if defined NEMO_V34
      INTEGER  :: itrst
      itrst = mod(kt,nn_trd)
#else
      INTEGER, PARAMETER  :: itrst = 0
#endif
      CHARACTER (LEN=99) :: cstrng

      IF( ks .NE. sppt_step ) RETURN

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Applying SPPT perturbation on T/S at timestep: ', kt
         WRITE(numout,*) ' Step: ',ks,' rdt:',rdt
         WRITE(numout,*) ' MAXVAL of perturbations for T/S: ',MAXVAL( ABS(gauss_n) )
         WRITE(numout,*) ' MAXVAL of tendencies    for T  : ',MAXVAL( ABS(spptt) ), MAXVAL( ABS(spptt) )*rdt
         WRITE(numout,*) ' MAXVAL of tendencies    for S  : ',MAXVAL( ABS(sppts) ), MAXVAL( ABS(sppts) )*rdt
         WRITE(numout,*)
      ENDIF

      ! Modify the tendencies
      IF(sppt_step .eq. 2 ) THEN
        ! At the end of the timestepping, after array swapping
#ifdef NEMO_V42
        DO_3D( 0, 0, 0, 0, 1, jpk )
          pts(ji,jj,jk,jp_tem,Kmm) = pts(ji,jj,jk,jp_tem,Kmm) + gauss_n(ji,jj,jk) * spptt(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          pts(ji,jj,jk,jp_sal,Kmm) = pts(ji,jj,jk,jp_sal,Kmm) + gauss_n(ji,jj,jk) * sppts(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          pts(ji,jj,jk,jp_tem,Kbb) = pts(ji,jj,jk,jp_tem,Kbb) + gauss_n(ji,jj,jk) * spptt(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          pts(ji,jj,jk,jp_sal,Kbb) = pts(ji,jj,jk,jp_sal,Kbb) + gauss_n(ji,jj,jk) * sppts(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
        END_3D
#else
        tsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem) + gauss_n * spptt * rdt / REAL ( itrst + 1, wp )
        tsn(:,:,:,jp_sal) = tsn(:,:,:,jp_sal) + gauss_n * sppts * rdt / REAL ( itrst + 1, wp )
        tsb(:,:,:,jp_tem) = tsb(:,:,:,jp_tem) + gauss_n * spptt * rdt / REAL ( itrst + 1, wp )
        tsb(:,:,:,jp_sal) = tsb(:,:,:,jp_sal) + gauss_n * sppts * rdt / REAL ( itrst + 1, wp )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsn(:,:,:,jp_tem), 'T' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsn(:,:,:,jp_sal), 'T' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsb(:,:,:,jp_tem), 'T' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsb(:,:,:,jp_sal), 'T' )
        CALL lbc_lnk(_LBCNAME_ tsb(:,:,:,jp_tem) , 'T', 1._wp)
        CALL lbc_lnk(_LBCNAME_ tsb(:,:,:,jp_sal) , 'T', 1._wp)
        CALL lbc_lnk(_LBCNAME_ tsn(:,:,:,jp_tem) , 'T', 1._wp)
        CALL lbc_lnk(_LBCNAME_ tsn(:,:,:,jp_sal) , 'T', 1._wp)
#endif
      ELSEIF (sppt_step .eq. 0 .or. sppt_step .eq. 1 ) THEN
        ! At the beginning / before vertical diffusion
#ifdef NEMO_V42
        DO_3D( 0, 0, 0, 0, 1, jpk )
           pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs) + gauss_n(ji,jj,jk) * spptt(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
           pts(ji,jj,jk,jp_sal,Krhs) = pts(ji,jj,jk,jp_sal,Krhs) + gauss_n(ji,jj,jk) * sppts(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
        END_3D
#else
        tsa(:,:,:,jp_tem) = tsa(:,:,:,jp_tem) + gauss_n * spptt / REAL ( itrst + 1, wp )
        tsa(:,:,:,jp_sal) = tsa(:,:,:,jp_sal) + gauss_n * sppts / REAL ( itrst + 1, wp )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsa(:,:,:,jp_tem), 'T' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( tsa(:,:,:,jp_sal), 'T' )
        CALL lbc_lnk(_LBCNAME_ tsa(:,:,:,jp_tem) , 'T', 1._wp)
        CALL lbc_lnk(_LBCNAME_ tsa(:,:,:,jp_sal) , 'T', 1._wp)
#endif
      ENDIF

#ifdef key_iomput
   cstrng='sppt_tem'
   IF (iom_use(TRIM(cstrng)) ) CALL iom_put( TRIM(cstrng) , gauss_n * spptt )
   cstrng='sppt_sal'
   IF (iom_use(TRIM(cstrng)) ) CALL iom_put( TRIM(cstrng) , gauss_n * sppts )
#endif

      IF( ln_stopack_diags ) THEN
        mi = MINVAL(spptt / REAL ( itrst + 1, wp ), MASK = (tmask > 0._wp ) )
        ma = MAXVAL(spptt / REAL ( itrst + 1, wp ), MASK = (tmask > 0._wp ) )
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
        IF(lwp) WRITE(numdiag,9301) lkt,'SPPT TEMPERATURE' ,mi,ma
        rn_mmin ( jk_sppt_tem ) =  MIN ( rn_mmin ( jk_sppt_tem ), mi )
        rn_mmax ( jk_sppt_tem ) =  MAX ( rn_mmax ( jk_sppt_tem ), ma )

        mi = MINVAL(sppts / REAL ( itrst + 1, wp ), MASK = (tmask > 0._wp ) )
        ma = MAXVAL(sppts / REAL ( itrst + 1, wp ), MASK = (tmask > 0._wp ) )
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
        IF(lwp) WRITE(numdiag,9301) lkt,'SPPT SALINITY   ' ,mi,ma
9301  FORMAT(' it :', i8, ' ', A16, ' Min: ',d10.3 , ' Max: ',d10.3)
        rn_mmin ( jk_sppt_sal ) =  MIN ( rn_mmin ( jk_sppt_sal ), mi )
        rn_mmax ( jk_sppt_sal ) =  MAX ( rn_mmax ( jk_sppt_sal ), ma )
      ENDIF

      ! Reset the tendencies. For NEMO3.4 only when nn_trd using persistent tendencies
      ! For later versions, always do that
#if defined NEMO_V34
      IF( ( itrst == (nn_trd-1 ) .OR. kt == (nitend-1) ) )   THEN
         spptt = 0._wp  ; sppts = 0._wp
      ENDIF
#else
      spptt = 0._wp  ; sppts = 0._wp
#endif

   END SUBROUTINE tra_sppt_apply

#if defined NEMO_V42
   SUBROUTINE dyn_sppt_apply ( kt , ks, Kbb, Kmm, puu, pvv, Krhs )
#else
   SUBROUTINE dyn_sppt_apply ( kt , ks)
#endif
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_sppt_apply ***
      !!
      !! ** Purpose :   Apply perturbation to momentum
      !! 		Step 0/1 for tendencies, step 2 for variables
      !! 		after timestepping
      !!----------------------------------------------------------------------

      INTEGER, INTENT( in ) :: kt, ks     ! Main time step counter
      REAL(wp) :: mi,ma
#if defined NEMO_V42
      INTEGER, INTENT( in ) :: Kbb, Kmm, Krhs  ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) :: puu, pvv
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
#endif
#if defined NEMO_V34
      INTEGER  :: itrst
      itrst = mod(kt,nn_trd)
#else
      INTEGER, PARAMETER  :: itrst = 0
#endif
      CHARACTER (LEN=99) :: cstrng

      IF( ks .NE. sppt_step ) RETURN

      IF ( ln_stopack_debug .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Applying SPPT perturbation on U/V at timestep: ', kt
         WRITE(numout,*) ' Step: ',ks,' rdt:',rdt
         WRITE(numout,*) ' MAXVAL of perturbations for U: ',MAXVAL( ABS(spptu) ), MAXVAL( ABS(spptu) )*rdt
         WRITE(numout,*) ' MAXVAL of perturbations for V: ',MAXVAL( ABS(spptv) ), MAXVAL( ABS(spptv) )*rdt
         WRITE(numout,*)
      ENDIF

      ! Modify the tendencies
      IF(sppt_step .eq. 2 ) THEN
#ifdef NEMO_V42
        DO_3D( 0, 0, 0, 0, 1, jpk )
          puu(ji,jj,jk,Kmm) = puu(ji,jj,jk,Kmm) + rn_uv_infl * gauss_nu(ji,jj,jk)* spptu(ji,jj,jk) * umask(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          pvv(ji,jj,jk,Kmm) = pvv(ji,jj,jk,Kmm) + rn_uv_infl * gauss_nv(ji,jj,jk)* spptv(ji,jj,jk) * vmask(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          puu(ji,jj,jk,Kbb) = puu(ji,jj,jk,Kbb) + rn_uv_infl * gauss_nu(ji,jj,jk)* spptu(ji,jj,jk) * umask(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
          pvv(ji,jj,jk,Kbb) = pvv(ji,jj,jk,Kbb) + rn_uv_infl * gauss_nv(ji,jj,jk)* spptv(ji,jj,jk) * vmask(ji,jj,jk) * rdt / REAL ( itrst + 1, wp )
        END_3D
#else
        ! At the end of the timestepping, after array swapping
        un(:,:,:) = un(:,:,:) + rn_uv_infl * gauss_nu* umask * spptu * rdt / REAL ( itrst + 1, wp )
        vn(:,:,:) = vn(:,:,:) + rn_uv_infl * gauss_nv* vmask * spptv * rdt / REAL ( itrst + 1, wp )
        ub(:,:,:) = ub(:,:,:) + rn_uv_infl * gauss_nu* umask * spptu * rdt / REAL ( itrst + 1, wp )
        vb(:,:,:) = vb(:,:,:) + rn_uv_infl * gauss_nv* vmask * spptv * rdt / REAL ( itrst + 1, wp )
        IF( ln_sppt_glocon ) CALL sppt_glocon( ub, 'U' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( vb, 'V' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( un, 'U' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( vn, 'V' )
        CALL lbc_lnk(_LBCNAME_ un , 'U', -1._wp)
        CALL lbc_lnk(_LBCNAME_ vn , 'V', -1._wp)
        CALL lbc_lnk(_LBCNAME_ ub , 'U', -1._wp)
        CALL lbc_lnk(_LBCNAME_ vb , 'V', -1._wp)
#endif
      ELSEIF (sppt_step .eq. 0 .or. sppt_step .eq. 1 ) THEN
        ! At the beginning / before vertical diffusion
#ifdef NEMO_V42
        DO_3D( 0, 0, 0, 0, 1, jpk )
           puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + rn_uv_infl * gauss_nu(ji,jj,jk)* spptu(ji,jj,jk) / REAL ( itrst + 1, wp )
           pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + rn_uv_infl * gauss_nv(ji,jj,jk)* spptv(ji,jj,jk) / REAL ( itrst + 1, wp )
        END_3D
#else
        ua(:,:,:) = ua(:,:,:) + rn_uv_infl * gauss_nu* spptu / REAL ( itrst + 1, wp )
        va(:,:,:) = va(:,:,:) + rn_uv_infl * gauss_nv* spptv / REAL ( itrst + 1, wp )
        IF( ln_sppt_glocon ) CALL sppt_glocon( ua, 'U' )
        IF( ln_sppt_glocon ) CALL sppt_glocon( va, 'V' )
        CALL lbc_lnk(_LBCNAME_ ua , 'U', -1._wp)
        CALL lbc_lnk(_LBCNAME_ va , 'V', -1._wp)
#endif
      ENDIF

#ifdef key_iomput
   cstrng='sppt_uvl'
   IF (iom_use(TRIM(cstrng)) ) CALL iom_put( TRIM(cstrng) , gauss_nu * spptu )
   cstrng='sppt_vvl'
   IF (iom_use(TRIM(cstrng)) ) CALL iom_put( TRIM(cstrng) , gauss_nv * spptv )
#endif

      IF( ln_stopack_diags ) THEN
        mi = MINVAL(spptu, MASK = (umask > 0._wp ) )
        ma = MAXVAL(spptu, MASK = (umask > 0._wp ) )
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
        IF(lwp) WRITE(numdiag,9301) lkt,'SPPT ZONAL CURR.' ,mi,ma
        rn_mmin ( jk_sppt_uvl ) =  MIN ( rn_mmin ( jk_sppt_uvl ), mi )
        rn_mmax ( jk_sppt_uvl ) =  MAX ( rn_mmax ( jk_sppt_uvl ), ma )

        mi = MINVAL(spptv, MASK = (vmask > 0._wp ) )
        ma = MAXVAL(spptv, MASK = (vmask > 0._wp ) )
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
        IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
        IF(lwp) WRITE(numdiag,9301) lkt,'SPPT MERID CURR.' ,mi,ma
9301  FORMAT(' it :', i8, ' ', A16, ' Min: ',d10.3 , ' Max: ',d10.3)
        rn_mmin ( jk_sppt_vvl ) =  MIN ( rn_mmin ( jk_sppt_vvl ), mi )
        rn_mmax ( jk_sppt_vvl ) =  MAX ( rn_mmax ( jk_sppt_vvl ), ma )
      ENDIF

      ! Reset the tendencies. For NEMO3.4 only when nn_trd using persistent tendencies
      ! For later versions, always do that
#if defined NEMO_V34
      IF( ( itrst == (nn_trd-1 ) .OR. kt == (nitend-1) ) )   THEN
         spptu = 0._wp  ; spptv = 0._wp
      ENDIF
#else
      spptu = 0._wp  ; spptv = 0._wp
#endif

   END SUBROUTINE dyn_sppt_apply

   SUBROUTINE sppt_glocon ( zts, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_glocon ***
      !!
      !! ** Purpose :   Apply global conservation constraint
      !!                Note: not coded for vvl (code will stop at initializ.)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(INOUT) :: zts(jpi,jpj,jpk)
      CHARACTER(len=1), INTENT(IN) :: cd_type
      INTEGER :: ji,jj,jk
      REAL(wp) :: zv, vl
      ! Calculate volume tendencies and renormalize
      ! Note: sshn should be staggered before being used.
#ifndef NEMO_V42
      SELECT CASE ( cd_type )
               CASE ( 'T' )  
                jk=1
                zv = SUM( tmask_i(:,:)*tmask(:,:,jk)*e1t(:,:)*e2t(:,:)*sshn(:,:)*zts(:,:,1) )
                vl = SUM( tmask_i(:,:)*tmask(:,:,jk)*e1t(:,:)*e2t(:,:)*sshn(:,:) )
                DO jk = 1, jpkm1
                   DO jj=1,jpj
                      DO ji=1,jpi
                        zv = zv + SUM( tmask_i(:,:)*tmask(:,:,jk)*e1t(:,:)*e2t(:,:)*fse3t(:,:,jk)*zts(:,:,jk) )
                        vl = vl + SUM( tmask_i(:,:)*tmask(:,:,jk)*e1t(:,:)*e2t(:,:)*fse3t(:,:,jk) )
                      END DO
                   END DO
                END DO
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ zv)
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ vl)
                zv = zv / vl
                zts = zts - zv*tmask
               CASE ( 'U' )
                jk=1
#if defined NEMO_V34
                zv = SUM( umask(:,:,jk)*e1u(:,:)*e2u(:,:)*sshn(:,:)*zts(:,:,1) )
                vl = SUM( umask(:,:,jk)*e1u(:,:)*e2u(:,:)*sshn(:,:) )
#else
                zv = SUM( umask_i(:,:)*umask(:,:,jk)*e1u(:,:)*e2u(:,:)*sshn(:,:)*zts(:,:,1) )
                vl = SUM( umask_i(:,:)*umask(:,:,jk)*e1u(:,:)*e2u(:,:)*sshn(:,:) )
#endif
                DO jk = 1, jpkm1
                   DO jj=1,jpj
                      DO ji=1,jpi
#if defined NEMO_V34
                        zv = zv + SUM( umask(:,:,jk)*e1u(:,:)*e2u(:,:)*fse3u(:,:,jk)*zts(:,:,jk) )
                        vl = vl + SUM( umask(:,:,jk)*e1u(:,:)*e2u(:,:)*fse3u(:,:,jk) )
#else
                        zv = zv + SUM( umask_i(:,:)*umask(:,:,jk)*e1u(:,:)*e2u(:,:)*fse3u(:,:,jk)*zts(:,:,jk) )
                        vl = vl + SUM( umask_i(:,:)*umask(:,:,jk)*e1u(:,:)*e2u(:,:)*fse3u(:,:,jk) )
#endif
                      END DO
                   END DO
                END DO
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ zv)
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ vl)
                zv = zv / vl
                zts = zts - zv*umask
               CASE ( 'V' )
                jk=1
#if defined NEMO_V34
                zv = SUM( vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*sshn(:,:)*zts(:,:,1) )
                vl = SUM( vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*sshn(:,:) )
#else
                zv = SUM( vmask_i(:,:)*vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*sshn(:,:)*zts(:,:,1) )
                vl = SUM( vmask_i(:,:)*vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*sshn(:,:) )
#endif
                DO jk = 1, jpkm1
                   DO jj=1,jpj
                      DO ji=1,jpi
#if defined NEMO_V34
                        zv = zv + SUM( vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*fse3v(:,:,jk)*zts(:,:,jk) )
                        vl = vl + SUM( vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*fse3v(:,:,jk) )
#else
                        zv = zv + SUM( vmask_i(:,:)*vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*fse3v(:,:,jk)*zts(:,:,jk) )
                        vl = vl + SUM( vmask_i(:,:)*vmask(:,:,jk)*e1v(:,:)*e2v(:,:)*fse3v(:,:,jk) )
#endif
                      END DO
                   END DO
                END DO
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ zv)
                IF(lk_mpp) CALL mpp_sum(_LBCNAME_ vl)
                zv = zv / vl
                zts = zts - zv*vmask
              ENDSELECT
#endif

   END SUBROUTINE sppt_glocon

   SUBROUTINE gaussian_ar1_field ( gn, nk, wei, gb, a, b,gn0, fltf, istep)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gaussian_ar1_field ***
      !!
      !! ** Purpose :   Generate correlated perturbation field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: a, b
      REAL(wp), DIMENSION(jpk)    , INTENT(INOUT) :: wei
      REAL(wp), DIMENSION(jpi,jpj), INTENT(INOUT) :: gb
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(INOUT) :: gn
      REAL(wp), DIMENSION(jpi,jpj    ), INTENT(OUT) :: gn0
      INTEGER, INTENT(IN) :: nk,istep
      REAL(wp), INTENT(IN) :: fltf
      REAL(wp) :: g2d(jpi,jpj)
      INTEGER :: jf, jk, ji, jj, nbl

      ! Random noise on 2d field
      IF ( istep == 1 ) THEN
         CALL sppt_rand2d( g2d ) 
         CALL lbc_lnk(_LBCNAME_ g2d , 'T', 1._wp)
         g2d_save = g2d
         IF( nn_sppt_step_bound .EQ. 1 ) THEN
           WHERE ( g2d < -ABS(rn_sppt_bound) ) g2d = -ABS(rn_sppt_bound)
           WHERE ( g2d >  ABS(rn_sppt_bound) ) g2d =  ABS(rn_sppt_bound)
         ENDIF
      ELSEIF ( istep == 2 ) THEN
         g2d = rn_spp_stdev * g2d_save / rn_sppt_stdev
      ELSEIF ( istep == 3 ) THEN
         g2d = rn_skeb_stdev * g2d_save / rn_sppt_stdev
      ENDIF
   
      ! Laplacian filter and re-normalization
      DO jf = 1, nk
         CALL sppt_flt( g2d )
      END DO
      g2d = g2d * fltf

#ifdef key_iommput
      ! Output the random field
      IF(istep==1) THEN
        CALL iom_put( "sppt_ran" , g2d )
      ELSEIF (istep==2) THEN
        CALL iom_put( "spp_ran" , g2d )
      ELSEIF (istep==3) THEN
        CALL iom_put( "skeb_ran" , g2d )
      ENDIF
#endif
   
      ! AR(1) process and array swap
      g2d = a*gb + b*g2d
      gb = g2d
      gn0 = g2d * zdc(:,:,1)

      IF (istep==2 .or. istep==3 .or. COUNT( (/ln_sppt_tra,ln_sppt_dyn/) ) .EQ. 0 ) RETURN

      ! From 2- to 3-d and vertical weigth
      IF(nn_vwei .eq. 2) THEN
       DO jj=1,jpj
        DO ji=1,jpi
         nbl = mbathy(ji,jj)
         wei(:) = 0._wp
         IF(nbl>1) THEN
           wei(1:nbl) = 1._wp
           wei(1) = 0._wp
           wei(2) = 0.5_wp
           wei(nbl-1) = 0.5_wp
           wei(nbl) = 0._wp
           DO jk=1,jpk
             gn(ji,jj,jk) = g2d(ji,jj) * wei(jk) * zdc(ji,jj,jk)
           ENDDO
         ENDIF
        ENDDO
       ENDDO
      ELSEIF(nn_vwei .eq. 3) THEN
       DO jj=1,jpj
        DO ji=1,jpi
         nbl = mbathy(ji,jj)
         wei(:) = 0._wp
         IF(nbl>1) THEN
           wei(1:nbl) = 1._wp
           wei(nbl-1) = 0.5_wp
           wei(nbl) = 0._wp
           DO jk=1,jpk
             gn(ji,jj,jk) = g2d(ji,jj) * wei(jk) * zdc(ji,jj,jk)
           ENDDO
         ENDIF
        ENDDO
       ENDDO
       ELSE
        DO jk=1,jpk
          gn(:,:,jk) = g2d * wei(jk) * zdc(:,:,jk)
        ENDDO
      ENDIF
      
      ! Bound
      IF( nn_sppt_step_bound .EQ. 2 ) THEN
        WHERE ( gn < -ABS(rn_sppt_bound) ) gn = -ABS(rn_sppt_bound)
        WHERE ( gn >  ABS(rn_sppt_bound) ) gn =  ABS(rn_sppt_bound)
        WHERE ( gn0< -ABS(rn_sppt_bound) ) gn0= -ABS(rn_sppt_bound)
        WHERE ( gn0>  ABS(rn_sppt_bound) ) gn0=  ABS(rn_sppt_bound)
      ENDIF

   END SUBROUTINE gaussian_ar1_field
   !
   SUBROUTINE stopack_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE stopack_init ***
      !!
      !! ** Purpose :   Initialize stopack
      !!----------------------------------------------------------------------
      !!
      INTEGER  :: ios, inum, ierror
      INTEGER  :: nn_sppt_tra, nn_sppt_dyn, nn_spp_aht, nn_sppt_ice
      INTEGER  :: ivals(8)
      INTEGER(KIND=8)     ::   ziseed(4)           ! RNG seeds in integer type
      REAL(KIND=8)        ::   zrseed(4)           ! RNG seeds in real type
      REAL(KIND=8)        ::   rdt_sppt
      TYPE(FLD_N) ::   sn_uc , sn_vc         ! informations about the fields to be read
      CHARACTER(len=100) ::  cn_dir
      !!
      NAMELIST/namstopack/ ln_stopack, ln_sppt_taumap, rn_sppt_tau, &
      ln_stopack_restart, sppt_filter_pass, rn_sppt_bound, sppt_step, nn_deftau,rn_sppt_stdev,&
      ln_distcoast, rn_distcoast, nn_vwei, nn_stopack_seed, nn_sppt_step_bound, nn_rndm_freq, &
      ln_sppt_glocon, ln_stopack_repr, rn_uv_infl, rn_ice_infl, &
      rn_icestr_sd, rn_icealb_sd, &
      ln_sppt_traxad, ln_sppt_trayad, ln_sppt_trazad, ln_sppt_trasad, ln_sppt_traldf, &
      ln_sppt_trazdf, ln_sppt_trazdfp,ln_sppt_traevd, ln_sppt_trabbc, ln_sppt_trabbl, &
      ln_sppt_tranpc, ln_sppt_tradmp, ln_sppt_traqsr, ln_sppt_transr, ln_sppt_traatf ,&
      ln_sppt_dynhpg, ln_sppt_dynspg, ln_sppt_dynkeg, ln_sppt_dynrvo, ln_sppt_dynpvo, &
      ln_sppt_dynzad, ln_sppt_dynldf, ln_sppt_dynzdf, ln_sppt_dynbfr, ln_sppt_dynatf, ln_sppt_dyntau,&
      ln_sppt_icehdf, ln_sppt_icelat, ln_sppt_icezdf, ln_sppt_tau, & 
      nn_spp_icealb, nn_spp_icestr,&
      spp_filter_pass,rn_spp_stdev,rn_spp_tau,&
      nn_spp_bfr,nn_spp_dqdt,nn_spp_dedt,nn_spp_avt,nn_spp_avm,nn_spp_qsi0,&
      nn_spp_ahtu,nn_spp_ahtv,nn_spp_ahm1,nn_spp_ahm2,&
      nn_spp_blkd,nn_spp_blkh,nn_spp_blke,&
      nn_spp_ahtw,nn_spp_ahtt,nn_spp_relw,nn_spp_arnf,nn_spp_geot,nn_spp_aevd,nn_spp_ahubbl,nn_spp_ahvbbl,&
      nn_spp_tkelc,nn_spp_tkedf,nn_spp_tkeds,nn_spp_tkebb,nn_spp_tkefr,&
      rn_bfr_sd,rn_dqdt_sd,rn_dedt_sd,rn_avt_sd,rn_avm_sd,rn_qsi0_sd,rn_ahtu_sd,rn_ahtv_sd,&
      rn_ahm1_sd,rn_ahm2_sd,&
      rn_ahtw_sd,rn_ahtt_sd, rn_relw_sd, rn_arnf_sd,rn_geot_sd, rn_aevd_sd,rn_ahubbl_sd,rn_ahvbbl_sd,&
      rn_tkelc_sd,rn_tkedf_sd,rn_tkeds_sd,rn_tkebb_sd,rn_tkefr_sd,&
      rn_blkd_sd,rn_blkh_sd,rn_blke_sd,&
      nn_spp_tdmp,rn_tdmp_sd,&
      ln_skeb,rn_skeb,nn_dcom_freq,rn_kh,rn_kc,ln_stopack_diags,skeb_filter_pass,rn_skeb_stdev,rn_skeb_tau,&
      nn_dconv,rn_beta_num, rn_beta_con, ln_stopack_debug, ln_skeb_tune, ln_spp, rn_beta_eke, &
      nn_uvdta, sn_uc, sn_vc, cn_dir, rn_heke, rn_veke, ln_skeb_apply, nn_skst

      ! Default values
      rn_sppt_bound = 3._wp
      ln_stopack = .false.
      ln_sppt_taumap = .false.
      rn_sppt_tau =  86400._wp * 5._wp
      sppt_filter_pass = 30
      nn_vwei = 0
      ln_distcoast = .false.
      ln_sppt_glocon = .false.
      rn_distcoast = 10._wp
      ln_stopack_restart = .true.
      nn_stopack_seed = (/1,2,3,narea/)
      nn_rndm_freq = 1
      nn_sppt_step_bound = 2
      nn_deftau = 1
      nn_skst   = 2
      ln_stopack_debug = .FALSE.
      ln_sppt_traxad       = .false.  ! Switch for x advection
      ln_sppt_trayad       = .false.  ! Switch for y advection
      ln_sppt_trazad       = .false.  ! Switch for z advection
      ln_sppt_trasad       = .false.  ! Switch for z advection  (s- case)
      ln_sppt_traldf       = .false.  ! Switch for lateral diffusion
      ln_sppt_trazdf       = .false.  ! Switch for vertical diffusion
      ln_sppt_trazdfp      = .false.  ! Switch for pure vertical diffusion
      ln_sppt_traevd       = .false.  ! Switch for enhanced vertical diffusion
      ln_sppt_trabbc       = .false.  ! Switch for bottom boundary condition
      ln_sppt_trabbl       = .false.  ! Switch for bottom boundary layer
      ln_sppt_tranpc       = .false.  ! Switch for non-penetrative convection
      ln_sppt_tradmp       = .false.  ! Switch for tracer damping
      ln_sppt_traqsr       = .false.  ! Switch for solar radiation
      ln_sppt_transr       = .false.  ! Switch for non-solar radiation / freshwater flux
      ln_sppt_traatf       = .false.  ! Switch for Asselin time-filter

      rn_uv_infl           = 1._wp
      ln_sppt_dynhpg       = .false.  ! Switch for hydrost. press. grad.
      ln_sppt_dynspg       = .false.  ! Switch for surface  press. grad.
      ln_sppt_dynkeg       = .false.  ! Switch for horiz. advcetion
      ln_sppt_dynrvo       = .false.  ! Switch for Relative vorticity
      ln_sppt_dynpvo       = .false.  ! Switch for planetary vortic.
      ln_sppt_dynzad       = .false.  ! Switch for vertical advection
      ln_sppt_dynldf       = .false.  ! Switch for lateral viscosity
      ln_sppt_dynzdf       = .false.  ! Switch for vertical viscosity
      ln_sppt_dynbfr       = .false.  ! Switch for bottom friction
      ln_sppt_dynatf       = .false.  ! Switch for Asselin filter
      ln_sppt_dyntau       = .false.  ! Switch for wind stress

      ln_sppt_tau          = .false.  ! Switch for wind stress

      rn_ice_infl          = 1._wp
      ln_sppt_icehdf       = .false.
      ln_sppt_icelat       = .false.
      ln_sppt_icezdf       = .false.

      ln_spp = .true.

      nn_spp_icestr  = 0
      nn_spp_icealb  = 0
      rn_icestr_sd = 0.30_wp
      rn_icealb_sd = 0.30_wp

      nn_spp_bfr  =0
      nn_spp_dqdt =0
      nn_spp_arnf =0
      nn_spp_aevd =0
      nn_spp_geot =0
      nn_spp_dedt =0
      nn_spp_avt  =0
      nn_spp_avm  =0
      nn_spp_qsi0 =0
      nn_spp_relw =0
      nn_spp_tkelc =0
      nn_spp_tkedf =0
      nn_spp_tkeds =0
      nn_spp_tkebb =0
      nn_spp_tkefr =0
      nn_spp_blkd =0
      nn_spp_blkh =0
      nn_spp_blke =0
      nn_spp_tdmp =0

      ln_skeb = .false.
      ln_stopack_diags = .false.
      rn_skeb_stdev = 1.0_wp
      skeb_filter_pass = 50
      rn_skeb_tau      = 50
      ln_skeb_tune = .false.
      nn_uvdta   = 1
      rn_beta_num = 1._wp
      rn_beta_con = 1._wp
      rn_beta_eke = 1._wp
      rn_heke     = 1._wp
      rn_veke     = 1._wp

#ifdef NEMO_V34
      REWIND( numnam )            
      READ  ( numnam, namstopack )
#else
#ifndef NEMO_V42
      REWIND( numnam_ref ) 
#endif
      READ  ( numnam_ref, namstopack, IOSTAT = ios, ERR = 901)
#ifdef NEMO_V4
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namstopack in reference namelist')
#else
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namstopack in reference namelist', lwp )
#endif

#ifndef NEMO_V42
      REWIND( numnam_cfg )  
#endif
      READ  ( numnam_cfg, namstopack, IOSTAT = ios, ERR = 902 )
#ifdef NEMO_V4
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namstopack in configuration namelist')
#else
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namstopack in configuration namelist', lwp )
#endif
      IF(lwm) WRITE ( numond, namstopack )
#endif

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'init_stopack : Stochastic physics package'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namstopack: '
         WRITE(numout,*)
         WRITE(numout,*) '       Switch on STOPACK                                     ln_stopack : ',ln_stopack
         WRITE(numout,*) '       Verbose diagnostics for STOPACK                  ln_stopack_diags:', ln_stopack_diags
         WRITE(numout,*)
         WRITE(numout,*) '  ***  Random generation'
         WRITE(numout,*) '       Read from restart file previous perturbation  ln_stopack_restart :', ln_stopack_restart
         WRITE(numout,*) '       Frequency of calls of the SPPT perturbation field   nn_rndm_freq :', nn_rndm_freq
         WRITE(numout,*) '       Reproducibility of random generation             ln_stopack_repr :', ln_stopack_repr
         WRITE(numout,*) '       Seed for random number generator (no restart)    nn_stopack_seed :', nn_stopack_seed(:)
         WRITE(numout,*)
         WRITE(numout,*) '  ***  SPPT scheme '
         WRITE(numout,*) '       SPPT step (0: beginning; 1: before ZVD; 2: after ZVD)  sppt_step : ',sppt_step
         WRITE(numout,*) '       Use map of decorr. time scale                     ln_sppt_taumap :', ln_sppt_taumap
         WRITE(numout,*) '       If ln_sppt_taumap FALSE, use this constant (in days) rn_sppt_tau :', rn_sppt_tau
         WRITE(numout,*) '       Number of filter passes to correlate in space   sppt_filter_pass :', sppt_filter_pass
         WRITE(numout,*) '       Standard deviation of the white noise              rn_sppt_stdev :', rn_sppt_stdev
         WRITE(numout,*) '       Apply global conservation constraints             ln_sppt_glocon :', ln_sppt_glocon
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on tracers:'
         WRITE(numout,*) '       Switch for x advection                            ln_sppt_traxad :', ln_sppt_traxad
         WRITE(numout,*) '       Switch for y advection                            ln_sppt_trayad :', ln_sppt_trayad
         WRITE(numout,*) '       Switch for z advection                            ln_sppt_trazad :', ln_sppt_trazad
         WRITE(numout,*) '       Switch for z advection  (s- case)                 ln_sppt_trasad :', ln_sppt_trasad
         WRITE(numout,*) '       Switch for lateral diffusion                      ln_sppt_traldf :', ln_sppt_traldf
         WRITE(numout,*) '       Switch for vertical diffusion                     ln_sppt_trazdf :', ln_sppt_trazdf
         WRITE(numout,*) '       Switch for pure vertical diffusion               ln_sppt_trazdfp :', ln_sppt_trazdfp
         WRITE(numout,*) '       Switch for enhanced vertical diffusion            ln_sppt_traevd :', ln_sppt_traevd
         WRITE(numout,*) '       Switch for bottom boundary condition              ln_sppt_trabbc :', ln_sppt_trabbc
         WRITE(numout,*) '       Switch for bottom boundary layer                  ln_sppt_trabbl :', ln_sppt_trabbl
         WRITE(numout,*) '       Switch for non-penetrative convection             ln_sppt_tranpc :', ln_sppt_tranpc
         WRITE(numout,*) '       Switch for tracer damping                         ln_sppt_tradmp :', ln_sppt_tradmp
         WRITE(numout,*) '       Switch for solar radiation                        ln_sppt_traqsr :', ln_sppt_traqsr
         WRITE(numout,*) '       Switch for non-solar rad. / freshwater flx flux   ln_sppt_transr :', ln_sppt_transr
         WRITE(numout,*) '       Switch for Asselin time-filter                    ln_sppt_traatf :', ln_sppt_traatf
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on dynamics:'
         WRITE(numout,*) '       Inflation coefficient for SPPT on DYN             rn_uv_infl     :', rn_uv_infl
         WRITE(numout,*) '       Switch for horiz advection                        ln_sppt_dynkeg :', ln_sppt_dynkeg
         WRITE(numout,*) '       Switch for z advection                            ln_sppt_dynzad :', ln_sppt_dynzad
         WRITE(numout,*) '       Switch for wind stress                            ln_sppt_dyntau :', ln_sppt_dyntau
         WRITE(numout,*) '       Switch for lateral diffusion                      ln_sppt_dynldf :', ln_sppt_dynldf
         WRITE(numout,*) '       Switch for vertical diffusion                     ln_sppt_dynzdf :', ln_sppt_dynzdf
         WRITE(numout,*) '       Switch for planet. vorticity                      ln_sppt_dynpvo :', ln_sppt_dynpvo
         WRITE(numout,*) '       Switch for relative vorticity                     ln_sppt_dynrvo :', ln_sppt_dynrvo
         WRITE(numout,*) '       Switch for hydrost. press. grad.                  ln_sppt_dynhpg :', ln_sppt_dynhpg
         WRITE(numout,*) '       Switch for surface  press. grad.                  ln_sppt_dynspg :', ln_sppt_dynspg
         WRITE(numout,*) '       Switch for bottom friction                        ln_sppt_dynbfr :', ln_sppt_dynbfr
         WRITE(numout,*) '       Switch for Asselin time-filter                    ln_sppt_dynatf :', ln_sppt_dynatf
         WRITE(numout,*)
         WRITE(numout,*) '       Switch for tau perturbation                       ln_sppt_tau    :', ln_sppt_tau    
         WRITE(numout,*)
         WRITE(numout,*) '       Perturbation on sea-ice:'
         WRITE(numout,*) '       Inflation coefficient for SPPT on ICE             rn_ice_infl    :', rn_ice_infl
         WRITE(numout,*) '       Switch for sea-ice diffusivity                    ln_sppt_icehdf :', ln_sppt_icehdf
         WRITE(numout,*) '       Switch for sea-ice lateral accretion              ln_sppt_icelat :', ln_sppt_icelat
         WRITE(numout,*) '       Switch for sea-ice vertical thermodyn.            ln_sppt_icezdf :', ln_sppt_icezdf
         WRITE(numout,*)
         WRITE(numout,*) '       Bound for gaussian random numbers                  rn_sppt_bound :', rn_sppt_bound
         WRITE(numout,*) '       Bound Step (0: nobound; 1: Gaussian; 2: Pert) nn_sppt_step_bound :', nn_sppt_step_bound
         WRITE(numout,*) '       Definition of Tau (1: days; 2: timesteps)              nn_deftau :', nn_deftau
         WRITE(numout,*) '       Smoothing of perturbation close to coast            ln_distcoast :', ln_distcoast
         WRITE(numout,*) '       Spatial scale of the smoothing near coasts (m)      rn_distcoast :', rn_distcoast
         WRITE(numout,*) '       Type of vertical weight:                                 nn_vwei :', nn_vwei
         WRITE(numout,*) '                            0 : No weight '
         WRITE(numout,*) '                            1 : First/Last level smoothing '
         WRITE(numout,*) '                            2 : Top/Bottom smoothing'
         WRITE(numout,*) '                            3 : Bottom smoothing'
         WRITE(numout,*)
         WRITE(numout,*) '  ***  SPP schemes (0: off ; 1: normal; 2:lognormal, mean as unpert.; 3: lognormal, median as unpert.'
         WRITE(numout,*) '                   (the meaning of standard dev. depends on the distribution and parametr.)'
         WRITE(numout,*)
         WRITE(numout,*) '       General switch for the SPP scheme                   ln_spp       :', ln_spp
         WRITE(numout,*)
         WRITE(numout,*) '       Number of passes for spatial filter (AR1 field)     spp_filter_pass:', spp_filter_pass
         WRITE(numout,*) '       Standard deviation of random generator (AR1 field)  rn_spp_stdev :', rn_spp_stdev      
         WRITE(numout,*) '       Decorr. time scale                     (AR1 field)  rn_spp_tau   :', rn_spp_tau
         WRITE(numout,*)
         WRITE(numout,*) '       SPP for bottom friction coeff                       nn_spp_bfr   :', nn_spp_bfr  
         WRITE(numout,*) '                                            STDEV          rn_bfr_sd    :', rn_bfr_sd   
         WRITE(numout,*) '       SPP for SST relaxation  coeff                       nn_spp_dqdt  :', nn_spp_dqdt 
         WRITE(numout,*) '                                            STDEV          rn_dqdt_sd   :', rn_dqdt_sd   
         WRITE(numout,*) '       SPP for SSS relaxation  coeff                       nn_spp_dedt  :', nn_spp_dedt 
         WRITE(numout,*) '                                            STDEV          rn_dedt_sd   :', rn_dedt_sd   
         WRITE(numout,*) '       SPP for vertical tra mixing coeff (only TKE, GLS)   nn_spp_avt   :', nn_spp_avt  
         WRITE(numout,*) '                                            STDEV          rn_avt_sd    :', rn_avt_sd   
         WRITE(numout,*) '       SPP for vertical dyn mixing coeff (only TKE, GLS)   nn_spp_avm   :', nn_spp_avm  
         WRITE(numout,*) '                                            STDEV          rn_avm_sd    :', rn_avm_sd   
         WRITE(numout,*) '       SPP for solar penetration scheme  (only RGB)        nn_spp_qsi0  :', nn_spp_qsi0
         WRITE(numout,*) '                                            STDEV          rn_qsi0_sd   :', rn_qsi0_sd   
         WRITE(numout,*) '       SPP for horiz. diffusivity  U                       nn_spp_ahtu  :', nn_spp_ahtu
         WRITE(numout,*) '                                            STDEV          rn_ahtu_sd   :', rn_ahtu_sd   
         WRITE(numout,*) '       SPP for horiz. diffusivity  V                       nn_spp_ahtv  :', nn_spp_ahtv
         WRITE(numout,*) '                                            STDEV          rn_ahtv_sd   :', rn_ahtv_sd   
         WRITE(numout,*) '       SPP for horiz. diffusivity  W                       nn_spp_ahtw  :', nn_spp_ahtw
         WRITE(numout,*) '                                            STDEV          rn_ahtw_sd   :', rn_ahtw_sd   
         WRITE(numout,*) '       SPP for horiz. diffusivity  T                       nn_spp_ahtt  :', nn_spp_ahtt
         WRITE(numout,*) '                                            STDEV          rn_ahtt_sd   :', rn_ahtt_sd   
         WRITE(numout,*) '       SPP for horiz. viscosity (1/3)                      nn_spp_ahm1  :', nn_spp_ahm1
         WRITE(numout,*) '                                            STDEV          rn_ahm1_sd   :', rn_ahm1_sd   
         WRITE(numout,*) '       SPP for horiz. viscosity (2/4)                      nn_spp_ahm2  :', nn_spp_ahm2
         WRITE(numout,*) '                                            STDEV          rn_ahm2_sd   :', rn_ahm2_sd   
         WRITE(numout,*) '       SPP for relative wind factor                        nn_spp_relw  :', nn_spp_relw
         WRITE(numout,*) '       (use 4, 5, 6 for nn_spp_relw to have options 1, 2, 3 with limits bounded to [0,1]'
         WRITE(numout,*) '                                            STDEV          rn_relw_sd   :', rn_relw_sd   
         WRITE(numout,*) '       SPP for mixing close to river mouth                 nn_spp_arnf  :', nn_spp_arnf
         WRITE(numout,*) '                                            STDEV          rn_arnf_sd   :', rn_arnf_sd   
         WRITE(numout,*) '       SPP for geothermal heating                          nn_spp_geot  :', nn_spp_geot
         WRITE(numout,*) '                                            STDEV          rn_geot_sd   :', rn_geot_sd   
         WRITE(numout,*) '       SPP for enhanced vertical diffusion                 nn_spp_aevd  :', nn_spp_aevd
         WRITE(numout,*) '                                            STDEV          rn_aevd_sd   :', rn_aevd_sd   
         WRITE(numout,*) '       SPP for TKE rn_lc    Langmuir cell coefficient      nn_spp_tkelc :', nn_spp_tkelc
         WRITE(numout,*) '                                            STDEV          rn_tkelc_sd  :', rn_tkelc_sd   
         WRITE(numout,*) '       SPP for TKE rn_ediff Eddy diff. coefficient         nn_spp_tkedf :', nn_spp_tkedf
         WRITE(numout,*) '                                            STDEV          rn_tkedf_sd  :', rn_tkedf_sd   
         WRITE(numout,*) '       SPP for TKE rn_ediss Kolmogoroff dissipation coeff. nn_spp_tkeds :', nn_spp_tkeds
         WRITE(numout,*) '                                            STDEV          rn_tkeds_sd  :', rn_tkeds_sd   
         WRITE(numout,*) '       SPP for TKE rn_ebb   Surface input of tke           nn_spp_tkebb :', nn_spp_tkebb
         WRITE(numout,*) '                                            STDEV          rn_tkebb_sd  :', rn_tkebb_sd   
         WRITE(numout,*) '       SPP for TKE rn_efr   Fraction of srf TKE below ML   nn_spp_tkefr :', nn_spp_tkefr
         WRITE(numout,*) '                                            STDEV          rn_tkefr_sd  :', rn_tkefr_sd   
         WRITE(numout,*) '       SPP for BBL U  diffusivity                          nn_spp_ahubbl:', nn_spp_ahubbl
         WRITE(numout,*) '                                            STDEV          rn_ahubbl_sd :', rn_ahubbl_sd
         WRITE(numout,*) '       SPP for BBL V  diffusivity                          nn_spp_ahvbbl:', nn_spp_ahvbbl
         WRITE(numout,*) '                                            STDEV          rn_ahvbbl_sd :', rn_ahvbbl_sd
         WRITE(numout,*)
         WRITE(numout,*) '       SPP for Air-sea drag coefficient (CORE BULK)        nn_spp_blkd  :', nn_spp_blkd
         WRITE(numout,*) '                                            STDEV          rn_blkd_sd   :', rn_blkd_sd
         WRITE(numout,*) '       SPP for Air-sea sensible heat transfer (CORE BULK)  nn_spp_blkh  :', nn_spp_blkh
         WRITE(numout,*) '                                            STDEV          rn_blkh_sd   :', rn_blkh_sd
         WRITE(numout,*) '       SPP for Air-sea evaporation transfer (CORE BULK)    nn_spp_blke  :', nn_spp_blke
         WRITE(numout,*) '                                            STDEV          rn_blke_sd   :', rn_blke_sd
         WRITE(numout,*) '       SPP for tracer damping                              nn_spp_tdmp  :', nn_spp_tdmp
         WRITE(numout,*) '                                            STDEV          rn_tdmp_sd   :', rn_tdmp_sd
         WRITE(numout,*)

         WRITE(numout,*)
         WRITE(numout,*) '  ***  SPP schemes for sea-ice '
         WRITE(numout,*) '       Albedo                                              nn_spp_icealb:', nn_spp_icealb
         WRITE(numout,*) '            St. dev. for ice albedo                        rn_icealb_sd :', rn_icealb_sd
         WRITE(numout,*) '       Ice Strength                                        nn_spp_icestr:', nn_spp_icestr
         WRITE(numout,*) '            St. dev. for ice strength                      rn_icestr_sd :', rn_icestr_sd
         WRITE(numout,*)
         WRITE(numout,*) ' SKEB Perturbation scheme '
         WRITE(numout,*) '       SKEB switch                                         ln_skeb      :', ln_skeb    
         WRITE(numout,*) '       SKEB ratio of backscattered energy                  rn_skeb      :', rn_skeb    
         WRITE(numout,*) '       Frequency update for dissipation mask               nn_dcom_freq :', nn_dcom_freq
         WRITE(numout,*) '       Numerical dissipation factor (resolut. dependent)   rn_kh        :', rn_kh
         WRITE(numout,*) '       Number of passes for spatial filter (AR1 field)     skeb_filter_pass:', skeb_filter_pass
         WRITE(numout,*) '       Standard deviation of random generator (AR1 field)  rn_skeb_stdev:', rn_skeb_stdev      
         WRITE(numout,*) '       Decorr. time scale                     (AR1 field)  rn_skeb_tau  :', rn_skeb_tau
         WRITE(numout,*) '       Option of convection energy dissipation             nn_dconv     :', nn_dconv
         WRITE(numout,*) '       Convection dissipation factor (resolut. dependent)  rn_kc        :', rn_kc
         WRITE(numout,*) '       Multiplier for numerical dissipation                rn_beta_num  :', rn_beta_num
         WRITE(numout,*) '       Multiplier for convection dissipation               rn_beta_con  :', rn_beta_con
         WRITE(numout,*) '       Multiplier for EKE dissipation                      rn_beta_eke  :', rn_beta_eke
         WRITE(numout,*) '       Read external U/V climatology for EKE dissipation   nn_uv_dta    :', nn_uvdta
         WRITE(numout,*) '       Tuning configuration for SKEB                       ln_skeb_tune :', ln_skeb_tune
         WRITE(numout,*) '       Step of application of SKEB                         nn_skst      :', nn_skst     

         WRITE(numout,*)
      ENDIF

      IF( .NOT. ln_stopack ) THEN
         IF(lwp) WRITE(numout,*) '     STOPACK is switched off'
         nn_spp = 0
         ln_spp = .false.
         ln_sppt_tra = .false.
         ln_sppt_dyn = .false.
         ln_sppt_ice = .false.
         ln_skeb = .false.
         RETURN
      ENDIF

      IF( MOD( nitend - nit000 + 1, nn_rndm_freq) /= 0 .OR.   &
          ( MOD( nstock, nn_rndm_freq) /= 0 .AND. nstock .GT. 0 ) ) THEN
         WRITE(ctmp1,*) 'experiment length (', nitend - nit000 + 1, ') or nstock (', nstock,   &
            &           ' is NOT a multiple of nn_rndm_freq (', nn_rndm_freq, ')'
         CALL ctl_stop( ctmp1, 'Impossible to properly setup STOPACK restart' )
      ENDIF

#ifdef NEMO_V4
      IF( ln_sppt_icehdf ) THEN
         CALL ctl_stop( 'ln_sppt_icehdf not available in NEMO4')
      ENDIF
#endif

      nn_sppt_tra = COUNT( (/ln_sppt_traxad, ln_sppt_trayad, ln_sppt_trazad, &
                             ln_sppt_trasad, ln_sppt_traldf, ln_sppt_trazdf, &
                            ln_sppt_trazdfp, ln_sppt_traevd, ln_sppt_trabbc, &
                             ln_sppt_trabbl, ln_sppt_tranpc, ln_sppt_tradmp, &
                             ln_sppt_traqsr, ln_sppt_transr, ln_sppt_traatf/) )
      nn_sppt_dyn = COUNT( (/ln_sppt_dynhpg, ln_sppt_dynspg, ln_sppt_dynkeg, &
                             ln_sppt_dynrvo, ln_sppt_dynpvo, ln_sppt_dynzad, &
                             ln_sppt_dynldf, ln_sppt_dynzdf, ln_sppt_dynbfr, &
                             ln_sppt_dynatf, ln_sppt_dyntau, ln_sppt_tau /) )
      nn_sppt_ice = COUNT( (/ln_sppt_icehdf, ln_sppt_icelat, ln_sppt_icezdf/) )

      ln_sppt_tra = ( nn_sppt_tra > 0 )
      ln_sppt_dyn = ( nn_sppt_dyn > 0 )
      ln_sppt_ice = ( nn_sppt_ice > 0 )

      nn_spp = nn_spp_bfr+nn_spp_dqdt+nn_spp_dedt+nn_spp_avt+nn_spp_avm+nn_spp_qsi0+&
      & nn_spp_ahtu+nn_spp_ahtv+nn_spp_ahtw+nn_spp_ahtt+nn_spp_relw+nn_spp_arnf+nn_spp_geot+nn_spp_aevd+&
      & nn_spp_tkelc+nn_spp_tkedf+nn_spp_tkeds+nn_spp_tkebb+nn_spp_tkefr+nn_spp_ahubbl+nn_spp_ahvbbl+&
      & nn_spp_ahm1+nn_spp_ahm2+nn_spp_blkd+nn_spp_blkh+nn_spp_blke+nn_spp_icealb+nn_spp_icestr+nn_spp_tdmp

      IF(.not. ln_spp ) THEN
          ! Overwrite default values to force switching off SPP
          nn_spp = 0
          nn_spp_bfr = 0
          nn_spp_dqdt = 0
          nn_spp_dedt = 0
          nn_spp_avt = 0
          nn_spp_avm = 0
          nn_spp_qsi0 = 0
          nn_spp_ahtu = 0
          nn_spp_ahtv = 0
          nn_spp_ahtw = 0
          nn_spp_ahtt = 0
          nn_spp_relw = 0
          nn_spp_arnf = 0
          nn_spp_geot = 0
          nn_spp_aevd = 0
          nn_spp_tkelc = 0
          nn_spp_tkedf = 0
          nn_spp_tkeds = 0
          nn_spp_tkebb = 0
          nn_spp_tkefr = 0
          nn_spp_ahubbl = 0
          nn_spp_ahvbbl = 0
          nn_spp_ahm1 = 0
          nn_spp_ahm2 = 0
          nn_spp_blkd = 0
          nn_spp_blkh = 0
          nn_spp_blke = 0
          nn_spp_icealb = 0
          nn_spp_icestr = 0
          nn_spp_tdmp = 0
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*) '    SPPT on tracers     : ',ln_sppt_tra
         WRITE(numout,*) '    SPPT on dynamics    : ',ln_sppt_dyn
         WRITE(numout,*) '    SPPT on sea-ice     : ',ln_sppt_ice
         WRITE(numout,*) '    SPP perturbed params: ',nn_spp
         WRITE(numout,*) '    SKEB scheme         : ',ln_skeb
      ENDIF

#ifdef NEMO_V34

#ifndef key_trdtra
      IF( ln_sppt_tra ) &
         & CALL ctl_stop( ' SPPT on tracers on .AND. .NOT. key_trdtra ', &
         &                ' SPPT cannot work without tracer tendency computation')
#endif
#ifndef key_trddyn
      IF( ln_sppt_dyn ) &
         & CALL ctl_stop( ' SPPT on momentum on .AND. .NOT. key_trddyn ', &
         &                ' SPPT cannot work without dynamic tendency computation')
#endif

#else
      IF( ln_sppt_tra .AND. .NOT.  ln_tra_trd ) &
         & CALL ctl_stop( ' SPPT on tracers on .AND. .NOT. ln_tra_trd ', &
         &                ' SPPT cannot work without tracer tendency computation')

      IF( ln_sppt_dyn .AND. .NOT.  ln_dyn_trd ) &
         & CALL ctl_stop( ' SPPT on momentum on .AND. .NOT. ln_dyn_trd ', &
         &                ' SPPT cannot work without dynamic tendency computation')
#endif

      IF( ln_sppt_glocon .AND. lk_vvl ) &
         & CALL ctl_stop( ' ln_sppt_glocon .AND. lk_vvl ', &
         &                ' SPPT conservation not coded yet for VVL')

      IF( nn_deftau .NE. 1 .AND. nn_deftau .NE. 2 ) &
         & CALL ctl_stop( ' nn_deftau must be 1 or 2 ')

      nn_spp_aht = nn_spp_ahtu + nn_spp_ahtv + nn_spp_ahtw + nn_spp_ahtt
      IF(nn_spp_aht .GT. 0) THEN
#if defined key_traldf_c3d
        IF(lwp) WRITE(numout,*) 'SPP : diffusivity perturbation with 3D coefficients'
#elif defined key_traldf_c2d
        IF(lwp) WRITE(numout,*) 'SPP : diffusivity perturbation with 2D coefficients'
#else
        CALL ctl_stop( 'SPP : diffusivity perturbation requires key_traldf_c3d or key_traldf_c2d')
#endif
      ENDIF

      IF(nn_spp_ahm1+nn_spp_ahm2 .GT. 0) THEN
#if defined key_dynldf_c3d
        IF(lwp) WRITE(numout,*) 'SPP : viscosity perturbation with 3D coefficients'
#elif defined key_dynldf_c2d
        IF(lwp) WRITE(numout,*) 'SPP : viscosity perturbation with 2D coefficients'
#else
        CALL ctl_stop( 'SPP : viscosity perturbation requires key_dynldf_c3d or key_dynldf_c2d')
#endif
      ENDIF

#ifdef NEMO_V34
      ! In case of SPPT and NEMO3.4, make sure tendency frequency calculation is sub-daily
      IF( REAL( nn_trd ) * rdt / 3600._wp .GT. 6._wp .AND. (nn_sppt_tra+nn_sppt_dyn) .GT. 0) THEN
          IF( REAL( nn_trd ) * rdt / 3600._wp .GT. 12._wp ) THEN
             CALL ctl_stop( 'SPPT : tendency terms persist for more than 12 hours')
          ELSE
             CALL ctl_warn( 'SPPT : tendency terms persist for more than  6 hours')
          ENDIF
      ENDIF
#endif

      ln_skeb_own_gauss = .FALSE.
      IF( ln_skeb ) THEN
        IF( skeb_filter_pass > 0 .AND. skeb_filter_pass .NE. sppt_filter_pass ) ln_skeb_own_gauss = .TRUE.
        IF( rn_skeb_tau > 0._wp .AND. rn_skeb_tau .NE. rn_sppt_tau ) ln_skeb_own_gauss = .TRUE.
      ENDIF

      ln_spp_own_gauss = .FALSE.
      IF( nn_spp > 0 ) THEN
        IF( spp_filter_pass > 0 .AND. spp_filter_pass .NE. sppt_filter_pass ) ln_spp_own_gauss = .TRUE.
        IF( rn_spp_tau > 0._wp .AND. rn_spp_tau .NE. rn_sppt_tau ) ln_spp_own_gauss = .TRUE.
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*) '    SPP with own scales : ',ln_spp_own_gauss
         WRITE(numout,*) '    SKEB with own scales: ',ln_skeb_own_gauss
      ENDIF

      IF( sppt_tra_alloc() /= 0) CALL ctl_stop( 'STOP', 'sppt_alloc: unable to allocate arrays' )

      spptt = 0._wp  ; sppts = 0._wp
      spptu = 0._wp  ; spptv = 0._wp

      ! Find filter attenuation factor
   
      flt_fac = sppt_fltfac( sppt_filter_pass )
      rdt_sppt = nn_rndm_freq * rn_rdt
   
      IF( ln_sppt_taumap ) THEN
         CALL iom_open ( 'sppt_tau_map', inum )
         CALL iom_get  ( inum, jpdom_global, 'tau', sppt_tau )
         CALL iom_close( inum )
      ELSE
         sppt_tau(:,:) = rn_sppt_tau
      ENDIF

      IF( nn_deftau .EQ. 1 ) THEN
         ! Decorrelation time-scale defined in days
         sppt_tau(:,:) = exp( -rdt_sppt / (sppt_tau(:,:)*86400._wp) )
      ELSE
         ! Decorrelation time-scale defined in timesteps
         sppt_tau(:,:) = exp( -REAL( nn_rndm_freq, wp) / sppt_tau(:,:) )
      ENDIF

      IF( ln_distcoast ) THEN
        CALL iom_open('dist.coast', inum )
        CALL iom_get(inum,jpdom_global,'Tcoast',zdc)
        CALL iom_close( inum )
        zdc = 1._wp - exp(-zdc/rn_distcoast)
        CALL lbc_lnk(_LBCNAME_ zdc , 'T', 1._wp)
      ELSE
        zdc = 1._wp
      ENDIF

      IF( ln_skeb ) THEN
        ALLOCATE(rn_kh2 ( jpi,jpj ), rn_kc2 ( jpi,jpj ) )
        rn_kh2(:,:) = rn_kh
        rn_kc2(:,:) = rn_kc
        IF ( rn_kh .LT. 0._wp ) THEN
             IF(lwp) WRITE(numout,*) '    Using 2D varying rn_kh from file skeb.nc'
             CALL iom_open('skeb', inum )
             CALL iom_get(inum,jpdom_global,'rn_kh',rn_kh2)
             CALL iom_close( inum )
        ENDIF
        IF ( rn_kc .LT. 0._wp ) THEN
             IF(lwp) WRITE(numout,*) '    Using 2D varying rn_kc from file skeb.nc'
             CALL iom_open('skeb', inum )
             CALL iom_get(inum,jpdom_global,'rn_kc',rn_kc2)
             CALL iom_close( inum )
        ENDIF
        IF ( ln_skeb_tune ) THEN
             IF(lwp) WRITE(numout,*) '    SKEB: tuning diagnostics activated'
             rn_kh2 (:,:) = 1._wp
             rn_kc2 (:,:) = 1._wp
             ktun = 0
             ALLOCATE( stun ( jpi,jpj,3 ) )
             stun (:,:,:) = 0._wp
        ENDIF
      ENDIF

      ! Initialize
      sppt_a = sppt_tau
      sppt_b = sqrt ( 1._wp - sppt_tau*sppt_tau )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '    **** SPPT SCHEME'
         WRITE(numout,*) '    Definit. time-scale : ',nn_deftau
         WRITE(numout,*) '     Decorr. time-scale : ',rn_sppt_tau
         WRITE(numout,*) '        Mean AR1 a term : ',SUM(sppt_a)/REAL(jpi*jpj)
         WRITE(numout,*) '        Mean AR1 b term : ',SUM(sppt_b)/REAL(jpi*jpj)
         WRITE(numout,*)
      ENDIF

      gauss_b = 0._wp
      ! Weigths
      gauss_w(:)    = 1.0_wp 
      IF( nn_vwei .eq. 1 ) THEN
        gauss_w(1)    = 0.0_wp
        gauss_w(2)    = 0.5_wp
        gauss_w(jpk)  = 0.0_wp
        gauss_w(jpkm1)= 0.5_wp
      ENDIF

      IF( ln_spp_own_gauss ) THEN
        IF( spp_alloc() /= 0) CALL ctl_stop( 'STOP', 'spp_alloc: unable to allocate arrays' )
        flt_fac_p = sppt_fltfac( spp_filter_pass )
        spp_tau (:, :) = rn_spp_tau
        IF( nn_deftau .EQ. 1 ) THEN
          spp_tau(:,:) = spp_tau(:,:) * 86400._wp
          spp_tau(:,:) = exp( -rdt_sppt / (spp_tau) )
        ELSE
          spp_tau(:,:) = exp( -1._wp / spp_tau )
        ENDIF
        spp_a = spp_tau
        spp_b = sqrt ( 1._wp - spp_tau*spp_tau )
        gauss_bp = 0._wp
        IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '    **** SPP  SCHEME'
          WRITE(numout,*) '    Definit. time-scale : ',nn_deftau
          WRITE(numout,*) '     Decorr. time-scale : ',rn_spp_tau
          WRITE(numout,*) '        Mean AR1 a term : ',SUM(spp_a)/REAL(jpi*jpj)
          WRITE(numout,*) '        Mean AR1 b term : ',SUM(spp_b)/REAL(jpi*jpj)
          WRITE(numout,*)
        ENDIF
      ENDIF

      IF( ln_skeb_own_gauss ) THEN
        IF( skeb_alloc() /= 0) CALL ctl_stop( 'STOP', 'skeb_alloc: unable to allocate arrays' )
        flt_fac_k = sppt_fltfac( skeb_filter_pass )
        skeb_tau (:, :) = rn_skeb_tau
        IF( nn_deftau .EQ. 1 ) THEN
          skeb_tau(:,:) = skeb_tau(:,:) * 86400._wp
          skeb_tau(:,:) = exp( -rdt_sppt / (skeb_tau) )
        ELSE
          skeb_tau(:,:) = exp( -1._wp / skeb_tau )
        ENDIF
        skeb_a = skeb_tau
        skeb_b = sqrt ( 1._wp - skeb_tau*skeb_tau )
        gauss_bk = 0._wp
        IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '    **** SKEB  SCHEME'
          WRITE(numout,*) '    Definit. time-scale : ',nn_deftau
          WRITE(numout,*) '     Decorr. time-scale : ',rn_skeb_tau
          WRITE(numout,*) '        Mean AR1 a term : ',SUM(skeb_a)/REAL(jpi*jpj)
          WRITE(numout,*) '        Mean AR1 b term : ',SUM(skeb_b)/REAL(jpi*jpj)
          WRITE(numout,*)
        ENDIF
      ENDIF

      IF( ln_skeb_own_gauss .OR. ln_spp_own_gauss ) &
      & ALLOCATE ( g2d_save(jpi,jpj) )

      CALL stopack_rst( nit000, 'READ' )

      IF(lwp .and. ln_stopack_diags) &
      CALL ctl_opn(numdiag, 'stopack.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )

      IF( ln_skeb .AND. nn_uvdta ) THEN
              ALLOCATE( sf_uc(1), sf_vc(1), STAT=ierror )
              IF( ierror > 0 ) THEN
                CALL ctl_stop( 'stopack_init: unable to allocate sf_uvc structure' )   ;   RETURN
              ENDIF
              ALLOCATE( sf_uc(1)%fnow(jpi,jpj,jpk), sf_vc(1)%fnow(jpi,jpj,jpk)   )
              IF( sn_uc%ln_tint )   ALLOCATE( sf_uc(1)%fdta(jpi,jpj,jpk,2) )
              IF( sn_vc%ln_tint )   ALLOCATE( sf_vc(1)%fdta(jpi,jpj,jpk,2) )
              CALL fld_fill( sf_uc, (/ sn_uc /), cn_dir, 'stopack_init',   &
              &           'U-Velocity climatology', 'namstopack' , no_print )
              CALL fld_fill( sf_vc, (/ sn_vc /), cn_dir, 'stopack_init',   &
              &           'V-Velocity climatology', 'namstopack' , no_print )
      ENDIF

      IF(lwp) WRITE(numout,*)  '    STOPACK INIT - END'
   
   END SUBROUTINE stopack_init
   !
   SUBROUTINE stopack_rst( kt, cdrw )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE stopack_rst ***
      !!
      !! ** Purpose :   Read/write from/to restarsts seed and perturbation field
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      INTEGER :: id1, jseed
      CHARACTER(LEN=10)   ::   clseed='spsd0_0000'
      INTEGER(KIND=8)     ::   ziseed(4)           ! RNG seeds in integer type 
      INTEGER(KIND=8)     ::   ivals(8) 
      REAL(wp)            ::   zrseed4(4)           ! RNG seeds in integer type
      REAL(wp)            ::   zrseed2d(jpi,jpj)
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise
         !                                   ! ---------------
         IF(lwp) WRITE(numout,*) ' *** stopack-rst ***'
         IF( ln_rstart ) THEN !* Read the restart file
            id1 = iom_varid( numror, 'sppt_b', ldstop = .FALSE. )
            IF(.NOT.ln_stopack_restart) id1=0 ! Overwrite if not needed
            !
            IF( id1 > 0 ) THEN
                CALL iom_get( numror, jpdom_auto, 'sppt_b', gauss_b )
                DO jseed = 1 , 4
                   WRITE(clseed(5:5) ,'(i1.1)') jseed
                   CALL iom_get( numror, jpdom_auto, clseed(1:5) , zrseed2d )
                   zrseed4(jseed) = zrseed2d(2,2)
                   ziseed(jseed) = TRANSFER( zrseed4(jseed) , ziseed(jseed) )
               END DO
               IF (lwp) THEN
                  WRITE(numout,*) 'Read ziseed4 = ',zrseed4(:)
                  WRITE(numout,*) 'Read ziseed  = ',ziseed(:)
               ENDIF
               CALL kiss_seed( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
            ELSE
               IF( lwp ) WRITE(numout,*) ' ===>>>> : previous run without sppt_b'
               IF( ln_stopack_restart ) CALL ctl_stop( ' ln_stopack_restart TRUE :',&
               & ' variable not found in restart file ')
               IF(lwp) WRITE(numout,*) ' ===>>>> : Initialisation of sppt_b to 0.'
               gauss_b = 0._wp
               IF ( .NOT. ln_stopack_repr ) THEN
                  CALL date_and_time(VALUES=ivals)
                  DO jseed=1,4
                    nn_stopack_seed(jseed) = nn_stopack_seed(jseed) + INT(ivals(4+jseed))
                  ENDDO
               ENDIF
               DO jseed=1,4
                 ziseed(jseed) = TRANSFER(nn_stopack_seed(jseed)+narea, ziseed(jseed))
               ENDDO
               IF ( nmember .gt. 0 ) THEN
                   DO jseed=1,4
                       ziseed(jseed) = ziseed(jseed) + (nmember*(10.**jseed))
                   ENDDO
               ENDIF
               IF(lwp) WRITE(numout,*) ' ===>>>> Seed, nn_stopack_seed+narea',nn_stopack_seed+narea
               IF(lwp) WRITE(numout,*) ' ===>>>> Seed, ziseed',ziseed
               CALL kiss_seed( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
            ENDIF
            IF( ln_spp_own_gauss ) THEN
              id1 = iom_varid( numror, 'spp_b', ldstop = .FALSE. )
              IF( id1 > 0 ) THEN
                 CALL iom_get( numror, jpdom_auto, 'spp_b', gauss_bp )
              ELSE
                 IF(lwp) WRITE(numout,*) ' ===>>>> : Initialisation of spp_b to 0.'
                 gauss_bp = 0._wp
              ENDIF
            ENDIF
            IF( ln_skeb_own_gauss ) THEN
              id1 = iom_varid( numror, 'skeb_b', ldstop = .FALSE. )
              IF( id1 > 0 ) THEN
                 CALL iom_get( numror, jpdom_auto, 'skeb_b', gauss_bk )
              ELSE
                 IF(lwp) WRITE(numout,*) ' ===>>>> : Initialisation of skeb_b to 0.'
                 gauss_bk = 0._wp
              ENDIF
            ENDIF
         ELSE
            gauss_b = 0._wp
            gauss_bp = 0._wp
            gauss_bk = 0._wp
            IF(lwp) WRITE(numout,*) ' ===>>>> : Initialisation of STOPACK to 0.'
            IF ( .NOT. ln_stopack_repr ) THEN
                CALL date_and_time(VALUES=ivals)
                DO jseed=1,4
                   nn_stopack_seed(jseed) = nn_stopack_seed(jseed) + INT(ivals(4+jseed))
                ENDDO
            ENDIF
            DO jseed=1,4
               ziseed(jseed) = TRANSFER(nn_stopack_seed(jseed)+narea, ziseed(jseed))
            ENDDO
            IF(lwp) WRITE(numout,*) ' ===>>>> Seed, nn_stopack_seed+narea',nn_stopack_seed+narea
            IF(lwp) WRITE(numout,*) ' ===>>>> Seed, ziseed',ziseed
            CALL kiss_seed( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- stopack-rst ----'
         CALL iom_rstput( kt, nitrst, numrow, 'sppt_b', gauss_b )
         IF(ln_spp_own_gauss) CALL iom_rstput( kt, nitrst, numrow, 'spp_b', gauss_bp )
         IF(ln_skeb_own_gauss) CALL iom_rstput( kt, nitrst, numrow, 'skeb_b', gauss_bk )
         CALL kiss_state( ziseed(1) , ziseed(2) , ziseed(3) , ziseed(4) )
         DO jseed=1,4
            zrseed4(jseed) = TRANSFER(ziseed(jseed), zrseed4(jseed))
         ENDDO
         IF (lwp) THEN
            WRITE(numout,*) 'Write ziseed4 = ',zrseed4(:)
            WRITE(numout,*) 'Write ziseed  = ',ziseed(:)
         ENDIF
         DO jseed = 1 , 4
            WRITE(clseed(5:5) ,'(i1.1)') jseed
            zrseed2d(:,:) = zrseed4(jseed)
            CALL iom_rstput( kt, nitrst, numrow, clseed(1:5) , zrseed2d )
         END DO
         !
      ENDIF
      !
   END SUBROUTINE stopack_rst
   !
   INTEGER FUNCTION sppt_tra_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION sppt_tra_alloc  ***
      !!---------------------------------------------------------------------
      !
      ALLOCATE( spptt(jpi,jpj,jpk) , sppts(jpi,jpj,jpk) , gauss_n(jpi,jpj,jpk) ,& 
      gauss_nu(jpi,jpj,jpk) , gauss_nv(jpi,jpj,jpk) , & 
      spptu(jpi,jpj,jpk) , spptv(jpi,jpj,jpk) , gauss_n_2d(jpi,jpj) ,&
      gauss_b (jpi,jpj), sppt_tau(jpi,jpj), sppt_a(jpi,jpj), sppt_b(jpi,jpj), gauss_w(jpk),&
      zdc(jpi,jpj,jpk), STAT= sppt_tra_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum (_LBCNAME_ sppt_tra_alloc )
      IF( sppt_tra_alloc /= 0)   CALL ctl_warn('sppt_tra_alloc: failed to allocate arrays')

      IF(nn_spp_bfr>0) THEN
        ALLOCATE(coeff_bfr(jpi,jpj), STAT= sppt_tra_alloc )
        IF( lk_mpp             )   CALL mpp_sum (_LBCNAME_ sppt_tra_alloc )
        IF( sppt_tra_alloc /= 0)   CALL ctl_warn('sppt_tra_alloc: failed to allocate arrays')
        coeff_bfr = 0._wp
      ENDIF
      gauss_b = 0._wp
      gauss_n = 0._wp
      gauss_nu= 0._wp
      gauss_nv= 0._wp
      gauss_n_2d = 0._wp

   END FUNCTION sppt_tra_alloc

   INTEGER FUNCTION skeb_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION skeb_alloc  ***
      !!---------------------------------------------------------------------
      !
      ALLOCATE( gauss_n_2d_k(jpi,jpj) , gauss_bk (jpi,jpj),skeb_tau(jpi,jpj),&
      skeb_a(jpi,jpj), skeb_b(jpi,jpj), STAT= skeb_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum (_LBCNAME_ skeb_alloc )
      IF( skeb_alloc /= 0)   CALL ctl_warn('sppt_tra_alloc: failed to allocate arrays')
      gauss_n_2d_k = 0._wp
      gauss_bk = 0._wp

   END FUNCTION skeb_alloc

   INTEGER FUNCTION spp_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION spp_alloc  ***
      !!---------------------------------------------------------------------
      !
      ALLOCATE( gauss_n_2d_p(jpi,jpj), gauss_bp (jpi,jpj),spp_tau(jpi,jpj), &
      spp_a(jpi,jpj), spp_b(jpi,jpj), STAT= spp_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum (_LBCNAME_ spp_alloc )
      IF( spp_alloc /= 0)   CALL ctl_warn('sppt_tra_alloc: failed to allocate arrays')
      gauss_n_2d_p = 0._wp
      gauss_bp = 0._wp

   END FUNCTION spp_alloc

   SUBROUTINE sppt_flt( psto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_flt  ***
      !!
      !! ** Purpose :   apply horizontal Laplacian filter to input array
      !!                Adapted from STO package
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   psto
      REAL(wp), DIMENSION(jpi,jpj)                ::   pstm
      !!
      INTEGER  :: ji, jj

      pstm = psto
      DO jj = 2, jpj-1
         DO ji = 2, jpi-1
            psto(ji,jj) = 0.5_wp * pstm(ji,jj) + 0.125_wp * &
                              &  ( pstm(ji-1,jj) + pstm(ji+1,jj) +  &
                              &    pstm(ji,jj-1) + pstm(ji,jj+1) )
         END DO
      END DO
      CALL lbc_lnk(_LBCNAME_ psto , 'T', 1._wp )

   END SUBROUTINE sppt_flt

   SUBROUTINE sppt_rand2d( psto )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sppt_rand2d ***
      !!
      !! ** Purpose :   fill input array with white Gaussian noise
      !!                Adapted from STO package
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out)           ::   psto
      !!
      INTEGER  :: ji, jj
      REAL(KIND=wp) :: gran   ! Gaussian random number (forced KIND=8 as in kiss_gaussian)

      DO jj = 1, jpj
         DO ji = 1, jpi
            CALL kiss_gaussian( gran )
            psto(ji,jj) = gran*rn_sppt_stdev
         END DO
      END DO

   END SUBROUTINE sppt_rand2d

   FUNCTION sppt_fltfac( kpasses )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION sppt_fltfac  ***
      !!
      !! ** Purpose :   compute factor to restore standard deviation
      !!                as a function of the number of passes
      !!                of the Laplacian filter
      !!                From STO package
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kpasses
      REAL(wp) :: sppt_fltfac
      !!
      INTEGER :: jpasses, ji, jj, jflti, jfltj
      INTEGER, DIMENSION(-1:1,-1:1) :: pflt0
      ! 16-b reals to avoid overflows
      INTEGER,  PARAMETER ::   qp = SELECTED_REAL_KIND(33,4931)
      REAL(qp), DIMENSION(:,:), ALLOCATABLE :: pfltb
      REAL(qp), DIMENSION(:,:), ALLOCATABLE :: pflta
      REAL(qp) :: ratio
      REAL(qp) :: aux0, aux1
      pflt0(-1,-1) = 0 ; pflt0(-1,0) = 1 ; pflt0(-1,1) = 0
      pflt0( 0,-1) = 1 ; pflt0( 0,0) = 4 ; pflt0( 0,1) = 1
      pflt0( 1,-1) = 0 ; pflt0( 1,0) = 1 ; pflt0( 1,1) = 0
      ALLOCATE(pfltb(-kpasses-1:kpasses+1,-kpasses-1:kpasses+1))
      ALLOCATE(pflta(-kpasses-1:kpasses+1,-kpasses-1:kpasses+1))
      pfltb(:,:) = 0
      pfltb(0,0) = 1
      DO jpasses = 1, kpasses
        pflta(:,:) = 0._qp
        DO jflti= -1, 1
        DO jfltj= -1, 1
          DO ji= -kpasses, kpasses
          DO jj= -kpasses, kpasses
            pflta(ji,jj) = pflta(ji,jj) + pfltb(ji+jflti,jj+jfltj) * pflt0(jflti,jfltj)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        pfltb(:,:) = pflta(:,:)
      ENDDO
      ratio = SUM(pfltb(:,:))
      aux0  = SUM(pfltb(:,:)*pfltb(:,:))
      aux1  = ratio*ratio
      ratio = aux1 / aux0
      ratio = SQRT(ratio)
      DEALLOCATE(pfltb,pflta)
      sppt_fltfac = REAL(ratio, wp)
   END FUNCTION sppt_fltfac

   SUBROUTINE skeb_comp( lt,  Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_comp  ***
      !!
      !! ** Purpose :   Computation of energy dissipation terms
      !!                This is a wrapper to the enrgy terms computation
      !!                routines
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: lt
#ifdef NEMO_V42
   INTEGER, INTENT(IN)  :: Kbb, Kmm, Krhs
#endif

#ifdef NEMO_V42
      CALL skeb_dnum ( lt )
      CALL skeb_dcon ( lt )
      CALL skeb_deke ( lt, Kbb, Kmm, Krhs)
#else
      CALL skeb_dnum ( lt )
      CALL skeb_dcon ( lt )
      CALL skeb_deke ( lt )
#endif
   END SUBROUTINE skeb_comp

#ifndef NEMO_V4
#ifndef NEMO_V42
   SUBROUTINE bsfcomp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bsfcomp  ***
      !!
      !! ** Purpose :   Online computation of barotropic streamfunction
      !!                For diagnostics purposes
      !!                This routine is not explicitly called anywhere
      !!                and utilizes low-level MPI routines
      !!----------------------------------------------------------------------
      USE MPI
      !
      INTEGER, INTENT(IN)  :: kt
      !
      REAL(wp) :: dtrpv(jpi,jpj)
      INTEGER  :: jk,ji,jj
      !
#if defined key_mpp_mpi
      INTEGER, DIMENSION(1) ::   ish
      INTEGER, DIMENSION(2) ::   ish2
      INTEGER               ::   ijp,tag,ierr,jp,isend,irecv
      REAL(wp), DIMENSION(jpj) ::   zwrk    ! mask flux array at V-point
      REAL(wp), DIMENSION(jpi) ::   zwr2    ! mask flux array at V-point
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
      IF(kt .EQ. nit000 ) THEN
#ifdef NEMO_V34
        ln_dpsiu = .FALSE.
        ln_dpsiv = .FALSE.
#else
        IF (iom_use('bsft') ) ln_dpsiu = .TRUE.
        IF (iom_use('bsftv') ) ln_dpsiv = .TRUE.
#endif
        IF( ln_dpsiv ) ALLOCATE ( dpsiv(jpi,jpj) )
        IF( ln_dpsiu ) ALLOCATE ( dpsiu(jpi,jpj) )
        IF( ln_dpsiv .OR. ln_dpsiu ) THEN
          jpri = nint ( REAL(nimpp) / REAL(jpi) ) + 1
          jprj = nint ( REAL(njmpp) / REAL(jpj) ) + 1
        ENDIF
      ENDIF

      IF( ln_dpsiv ) THEN
       dtrpv = 0._wp
       DO jk = 1,jpkm1
           dtrpv = dtrpv + fse3v(:,:,jk) * e1v(:,:) * vn(:,:,jk)
       ENDDO
       dpsiv(1,:)=0._wp
       DO ji = 2,jpi
          dpsiv(ji,:) = dpsiv(ji-1,:) + dtrpv(ji,:)
       END DO
      ENDIF

      IF( ln_dpsiu ) THEN
       dtrpv = 0._wp
       DO jk = 1,jpkm1
           dtrpv = dtrpv + fse3u(:,:,jk) * e2u(:,:) * un(:,:,jk)
       ENDDO
       dpsiu(:,1)=0._wp
       DO jj = 2,jpj
          dpsiu(:,jj) = dpsiu(:,jj-1) + dtrpv(:,jj)
       END DO
      ENDIF

#if defined key_mpp_mpi
      IF ( ln_dpsiv ) THEN
       DO jp=1,jpni-1
         IF( jpri == jp ) THEN ! SEND TO EAST 
          zwrk(1:jpj) = dpsiv(jpi-1,:)
          tag=2000+narea
          CALL mpi_isend(zwrk, jpj, mpi_double_precision, noea, tag, mpi_comm_opa, isend,ierr)
         ELSEIF ( jpri == jp+1 ) THEN ! RECEIVE FROM WEST
          CALL mpi_irecv(zwrk, jpj, mpi_double_precision, nowe, mpi_any_tag, mpi_comm_opa, irecv,ierr)
         ENDIF
         IF(jpri == jp) CALL mpi_wait(isend, istatus, ierr)
         IF(jpri == jp+1 ) THEN
          CALL mpi_wait(irecv, istatus, ierr)
          DO ji=1,jpi
            dpsiv(ji,:) = dpsiv(ji,:) + zwrk(:)
          ENDDO
         ENDIF
       ENDDO
      ENDIF

      IF ( ln_dpsiv ) THEN
       DO jp=1,jpnj-1
         IF( jprj == jp ) THEN ! SEND TO NORTH
          zwr2(1:jpi) = dpsiu(:,jpj-1)
          tag=3000+narea
          CALL mpi_isend(zwr2, jpi, mpi_double_precision, nono, tag, mpi_comm_opa, isend,ierr)
         ELSEIF ( jprj == jp+1 ) THEN ! RECEIVE FROM SOUTH
          CALL mpi_irecv(zwr2, jpi, mpi_double_precision, noso, mpi_any_tag, mpi_comm_opa, irecv,ierr)
         ENDIF
         IF(jprj == jp) CALL mpi_wait(isend, istatus, ierr)
         IF(jprj == jp+1 ) THEN
          CALL mpi_wait(irecv, istatus, ierr)
          DO ji=1,jpj
            dpsiu(:,ji) = dpsiu(:,ji) + zwr2(:)
          ENDDO
         ENDIF
       ENDDO
      ENDIF
#endif
      IF (ln_dpsiu ) THEN
          CALL lbc_lnk(_LBCNAME_ dpsiu,'T',1._wp)
#ifdef key_iomput
          CALL iom_put( "bsft" , dpsiu )
#endif
      ENDIF
      IF (ln_dpsiv ) THEN
          CALL lbc_lnk(_LBCNAME_ dpsiv,'T',1._wp)
#ifdef key_iomput
          CALL iom_put( "bsftv" , dpsiv )
#endif
      ENDIF

      IF(kt .EQ. nitend ) THEN
        IF (ln_dpsiv ) DEALLOCATE ( dpsiv )
        IF (ln_dpsiu ) DEALLOCATE ( dpsiu )
      ENDIF

   END SUBROUTINE
#endif
#endif

   SUBROUTINE skeb_dnum ( mt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_dnum  ***
      !!
      !! ** Purpose :   Computation of numerical energy dissipation
      !!                For later use in the SKEB scheme
      !!----------------------------------------------------------------------
   INTEGER, INTENT(IN)  :: mt
   REAL(wp) :: ds,dt,dtot
   INTEGER :: ji,jj,jk
   
   IF ( mt .eq. nit000 ) THEN
        ALLOCATE ( dnum(jpi,jpj,jpk) )
        dnum (:,:,: ) = 0._wp
   ENDIF

#ifndef NEMO_V42
   IF( mt .eq. nit000 .OR. MOD( mt - 1, nn_dcom_freq ) == 0 ) THEN

     DO jk=1,jpkm1
      DO jj=2,jpj
       DO ji=2,jpi
          ! Shear
          ds = (vn(ji,jj,jk)-vn(ji-1,jj,jk))*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj) + &
               (un(ji,jj,jk)-un(ji,jj-1,jk))*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
          ! Tension
          dt = (vn(ji,jj,jk)-vn(ji-1,jj,jk))*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e2v(ji,jj) + &
               (un(ji,jj,jk)-un(ji,jj-1,jk))*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e1u(ji,jj)
  
          dtot = sqrt( ds*ds + dt*dt ) * tmask(ji,jj,jk)
          dnum(ji,jj,jk) = dtot*dtot*dtot*rn_kh2(ji,jj)*rn_kh2(ji,jj)*e1t(ji,jj)*e2t(ji,jj)
       ENDDO
      ENDDO
     ENDDO
   
     CALL lbc_lnk(_LBCNAME_ dnum,'T',1._wp)

   ENDIF
#endif

   END SUBROUTINE 

#ifdef NEMO_V42
   SUBROUTINE skeb_deke ( mt,   Kbb, Kmm, Krhs )
#else
   SUBROUTINE skeb_deke ( mt )
#endif
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_eked  ***
      !!
      !! ** Purpose :   Computation of eddy dissipation
      !!                For later use in the SKEB scheme
      !!----------------------------------------------------------------------
   INTEGER, INTENT(IN)  :: mt
#ifdef NEMO_V42
   INTEGER, INTENT(IN)  :: Kbb, Kmm, Krhs
#endif
   REAL(wp) :: uTp,uTm,dudz,vTp,vTm,dvdz,Av,dudx,dvdy,Ah,mi,ma
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: tmpe
   INTEGER :: ji,jj,jk

   IF ( mt .eq. nit000 ) THEN
        ALLOCATE ( deke(jpi,jpj,jpk), dekeh(jpi,jpj,jpk), dekev(jpi,jpj,jpk) )
        ALLOCATE ( ut(jpi,jpj,jpk), vt(jpi,jpj,jpk) )
        deke (:,:,: ) = 0._wp
        dekeh(:,:,: ) = 0._wp
        dekev(:,:,: ) = 0._wp
        ut(:,:,: ) = 0._wp
        vt(:,:,: ) = 0._wp
   ENDIF

   IF( mt .eq. nit000 .OR. MOD( mt - 1, nn_dcom_freq ) == 0 ) THEN

     IF( nn_uvdta ) THEN
             CALL fld_read( mt, nn_dcom_freq, sf_uc )
             CALL fld_read( mt, nn_dcom_freq, sf_vc )
#ifdef NEMO_V42
             ut = ( uu(:,:,:,Kmm) - sf_uc(1)%fnow(:,:,:) ) * umask
             vt = ( vv(:,:,:,Kmm) - sf_vc(1)%fnow(:,:,:) ) * vmask
#else
             ut = ( un - sf_uc(1)%fnow(:,:,:) ) * umask
             vt = ( vn - sf_vc(1)%fnow(:,:,:) ) * vmask
#endif
     ELSE
#ifdef NEMO_V42
             CALL ctl_stop ( ' nn_uvdta != 1 not supported in NEMO 4.2 ' )
#else
             ut = ( un - ub ) * umask
             vt = ( vn - vb ) * vmask
#endif
     ENDIF 

     ALLOCATE ( tmpe(jpi,jpj,jpk) )

#ifdef NEMO_V42
 dekev(:,:,:) = 0._wp
#else
     ! Vertical part, on grid w
     DO jk=2,jpk
      DO jj=2,jpj
       DO ji=2,jpi

        uTp = 0.5_wp*(ut(ji,jj,jk-1)+ut(ji-1,jj,jk-1))  * tmask(ji,jj,jk-1) * umask(ji-1,jj,jk-1)
        uTm = 0.5_wp*(ut(ji,jj,jk)  +ut(ji-1,jj,jk))    * tmask(ji,jj,jk)   * umask(ji-1,jj,jk)
        dudz = (uTp-uTm) / e3w_n(ji,jj,jk)              * wmask(ji,jj,jk)
        vTp = 0.5_wp*(vt(ji,jj,jk-1)+vt(ji,jj-1,jk-1))  * tmask(ji,jj,jk-1) * vmask(ji,jj-1,jk-1)
        vTm = 0.5_wp*(vt(ji,jj,jk)  +vt(ji,jj-1,jk))    * tmask(ji,jj,jk)   * vmask(ji,jj-1,jk)
        dvdz = (vTp-vTm) / e3w_n(ji,jj,jk)              * wmask(ji,jj,jk)
        Av = avm(ji,jj,jk) 
        tmpe(ji,jj,jk) = rau0*(-Av)*(dudz*dudz+dvdz*dvdz) 

       ENDDO
      ENDDO
     ENDDO
     
     ! on grid W
     dekev(:,:,1) = 0._wp
     DO jk=2,jpkm1
        dekev(:,:,jk) = 0.5_wp*( tmpe(:,:,jk) + tmpe(:,:,jk+1) ) * tmask(:,:,jk)
     ENDDO
     CALL lbc_lnk(_LBCNAME_ dekev,'T',1._wp)
     DEALLOCATE ( tmpe )
#endif

     DO jk=2,jpkm1
      DO jj=2,jpj
       DO ji=2,jpi

        dudx= (ut(ji,jj,jk)-ut(ji-1,jj,jk))/e1t(ji,jj) * tmask(ji,jj,jk) * umask(ji-1,jj,jk)
        dvdy= (vt(ji,jj,jk)-vt(ji,jj-1,jk))/e2t(ji,jj) * tmask(ji,jj,jk) * vmask(ji,jj-1,jk)
        Ah  = ahmt(ji,jj,jk)

        dekeh(ji,jj,jk) = rau0*(-Ah)*(dudx*dudx+dvdy*dvdy)

       ENDDO
      ENDDO
     ENDDO
     CALL lbc_lnk(_LBCNAME_ dekeh,'T',1._wp)

     deke(:,:,:) = ( rn_heke*dekeh(:,:,:) + rn_veke*dekev(:,:,:) ) * tmask(:,:,:)

     IF( ln_stopack_diags ) THEN
       mi = MAXVAL(ABS(rn_heke*dekeh))
       ma = MAXVAL(ABS(rn_veke*dekev))
       IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ mi)
       IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
       IF(lwp) WRITE(numdiag,9304) mt,'DEKE ABS ENERGY' ,mi,ma
9304 FORMAT(' it :', i8, ' ', A16, ' Ekh: ',d10.3 , ' Ekv: ',d10.3)
     ENDIF

   ENDIF

   END SUBROUTINE 
 
   SUBROUTINE skeb_dcon ( mt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_dcon  ***
      !!
      !! ** Purpose :   Computation of convective energy dissipation
      !!                For later use in the SKEB scheme
      !!                Two formulation are implemented, based on buoyancy or
      !!                convective mass flux formulations, respectively
      !!----------------------------------------------------------------------
   INTEGER, INTENT(IN)  :: mt
   REAL(wp) :: zz, mf1,mf2
   INTEGER  :: ji,jj,jk

   IF ( mt .eq. nit000 ) THEN
        ALLOCATE ( dcon(jpi,jpj,jpk) )
        dcon (:,:,: ) = 0._wp
   ENDIF

#ifndef NEMO_V42
   IF( mt .eq. nit000 .OR. MOD( mt - 1, nn_dcom_freq ) == 0 ) THEN

    IF(nn_dconv  .eq. 1 ) THEN

     DO jk=2,jpkm1
      DO jj=1,jpj
       DO ji=1,jpi

           zz = - grav*avt(ji,jj,jk) * ( rhd(ji,jj,jk)-rhd(ji,jj,jk-1) ) * wmask(ji,jj,jk) * tmask(ji,jj,jk) * tmask(ji,jj,jk-1) &
              & / ( rau0 * fse3w(ji,jj,jk) ) 

           dcon(ji,jj,jk) = rn_kc2(ji,jj)*rn_kc2(ji,jj)*zz*e1t(ji,jj)*e2t(ji,jj)*rau0 / fse3w(ji,jj,jk)

       ENDDO
      ENDDO
     ENDDO

    ELSEIF (nn_dconv .eq. 2 ) THEN

     DO jk=2,jpkm1
      DO jj=1,jpj
       DO ji=1,jpi

           mf1 = wn(ji,jj,jk+1)*e1t(ji,jj)*e2t(ji,jj)
           mf2 = wn(ji,jj,jk-1)*e1t(ji,jj)*e2t(ji,jj)
           dcon(ji,jj,jk) = rn_kc2(ji,jj) * (mf1-mf2)*(mf1-mf2) * tmask(ji,jj,jk+1) / (fse3w(ji,jj,jk)*rau0*rau0)

       ENDDO
      ENDDO
     ENDDO

    ENDIF

     CALL lbc_lnk(_LBCNAME_ dcon,'T',1._wp)

   ENDIF
#endif

   END SUBROUTINE

#ifdef NEMO_V42
   SUBROUTINE skeb_apply ( mt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE skeb_apply  ***
      !!
      !! ** Purpose :   Application of SKEB perturbation
      !!                Convective and Numerical energy dissipation are
      !!                multiplied by the beta terms
      !!----------------------------------------------------------------------
   INTEGER, INTENT(IN)  :: Kbb, Kmm, Krhs
#else
   SUBROUTINE skeb_apply ( mt )
#endif
   INTEGER, INTENT(IN)  :: mt
   INTEGER :: ji,jj,jk
   REAL(wp), DIMENSION(jpi,jpj,jpk) :: ustar,vstar,psi
   REAL(wp), DIMENSION(jpi,jpj    ) :: ustar2,vstar2,zinthu,zinthv
   REAL(wp) :: mi,ma,m2

   ! Do not apply perturbations at the first timestep
   IF( mt .eq. nit000 ) RETURN

   IF( ln_stopack_diags .OR. ln_skeb_tune) THEN

     psi = 0._wp
     IF(ln_skeb_own_gauss) THEN
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb * sqrt( rn_beta_num * dnum(:,:,jk) ) * gauss_n_2d_k(:,:) 
       ENDDO
     ELSE
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb_stdev * rn_skeb * sqrt( rn_beta_num * dnum(:,:,jk) ) &
         & * gauss_n_2d(:,:) / rn_sppt_stdev
       ENDDO
     ENDIF
     ustar = 0._wp
     vstar = 0._wp
     DO jk=1,jpkm1
      DO jj=2,jpj
       DO ji=2,jpi
          ustar(ji,jj,jk) =   ( psi(ji,jj,jk)-psi(ji,jj-1,jk) )*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
          vstar(ji,jj,jk) = - ( psi(ji,jj,jk)-psi(ji-1,jj,jk) )*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj)
       ENDDO
      ENDDO
     ENDDO
     mi = MAXVAL(ABS(ustar))
     ma = MAXVAL(ABS(vstar))
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ mi)
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
     IF(lwp) WRITE(numdiag,9304) lkt,'DNUM ABS CURRENT' ,mi,ma

     IF( ln_skeb_tune ) THEN
         ktun=ktun+1
         DO jj=2,jpj
           DO ji=2,jpi
              mi = MAXVAL(ABS(ustar(ji,jj,:)))
              ma = MAXVAL(ABS(vstar(ji,jj,:)))
              stun(ji,jj,1) = stun(ji,jj,1) + max(mi,ma)
           ENDDO
         ENDDO
     ENDIF

     rn_mmax ( jk_skeb_dnum ) =  MAX ( MAX( mi, ma), rn_mmax ( jk_skeb_dnum ) )

     psi = 0._wp
     IF(ln_skeb_own_gauss) THEN
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb * sqrt( rn_beta_con * dcon(:,:,jk) ) * gauss_n_2d_k(:,:) 
       ENDDO
     ELSE
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb_stdev * rn_skeb * sqrt( rn_beta_con * dcon(:,:,jk) ) &
         & * gauss_n_2d(:,:) / rn_sppt_stdev
       ENDDO
     ENDIF
     ustar = 0._wp
     vstar = 0._wp
     DO jk=1,jpkm1
      DO jj=2,jpj
       DO ji=2,jpi
          ustar(ji,jj,jk) =   ( psi(ji,jj,jk)-psi(ji,jj-1,jk) )*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
          vstar(ji,jj,jk) = - ( psi(ji,jj,jk)-psi(ji-1,jj,jk) )*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj)
       ENDDO
      ENDDO
     ENDDO
     mi = MAXVAL(ABS(ustar))
     ma = MAXVAL(ABS(vstar))
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ mi)
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
     IF(lwp) WRITE(numdiag,9304) lkt,'DCON ABS CURRENT' ,mi,ma

     IF( ln_skeb_tune ) THEN
         DO jj=2,jpj
           DO ji=2,jpi
              mi = MAXVAL(ABS(ustar(ji,jj,:)))
              ma = MAXVAL(ABS(vstar(ji,jj,:)))
              stun(ji,jj,2) = stun(ji,jj,2) + max(mi,ma)
           ENDDO
         ENDDO
     ENDIF

     rn_mmax ( jk_skeb_dcon ) =  MAX ( MAX( mi, ma), rn_mmax ( jk_skeb_dcon ) )

     psi = 0._wp
     IF(ln_skeb_own_gauss) THEN
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb * sqrt( rn_beta_eke * deke(:,:,jk) ) * gauss_n_2d_k(:,:) 
       ENDDO
     ELSE
       DO jk=1,jpkm1
         psi(:,:,jk) = rn_skeb_stdev * rn_skeb * sqrt( rn_beta_eke * deke(:,:,jk) ) &
         & * gauss_n_2d(:,:) / rn_sppt_stdev
       ENDDO
     ENDIF
     ustar = 0._wp
     vstar = 0._wp
     DO jk=1,jpkm1
      DO jj=2,jpj
       DO ji=2,jpi
          ustar(ji,jj,jk) =   ( psi(ji,jj,jk)-psi(ji,jj-1,jk) )*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
          vstar(ji,jj,jk) = - ( psi(ji,jj,jk)-psi(ji-1,jj,jk) )*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj)
       ENDDO
      ENDDO
     ENDDO
     mi = MAXVAL(ABS(ustar))
     ma = MAXVAL(ABS(vstar))
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ mi)
     IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
     IF(lwp) WRITE(numdiag,9304) lkt,'DEKE ABS CURRENT' ,mi,ma

     IF( ln_skeb_tune ) THEN
         DO jj=2,jpj
           DO ji=2,jpi
              mi = MAXVAL(ABS(ustar(ji,jj,:)))
              ma = MAXVAL(ABS(vstar(ji,jj,:)))
              stun(ji,jj,3) = stun(ji,jj,3) + max(mi,ma)
           ENDDO
         ENDDO
     ENDIF

     rn_mmax ( jk_skeb_deke ) =  MAX ( MAX( mi, ma), rn_mmax ( jk_skeb_deke ) )


9304 FORMAT(' it :', i8, ' ', A16, ' Min: ',d10.3 , ' Max: ',d10.3)

   ENDIF

   psi = 0._wp
   IF(ln_skeb_own_gauss) THEN
     DO jk=1,jpkm1
       psi(:,:,jk) = rn_skeb * sqrt( rn_beta_num * dnum(:,:,jk)+ rn_beta_con * dcon(:,:,jk) + rn_beta_eke * abs( deke(:,:,jk) )) * gauss_n_2d_k(:,:) 
     ENDDO
   ELSE
     DO jk=1,jpkm1
       psi(:,:,jk) = rn_skeb_stdev * rn_skeb * sqrt( rn_beta_num * dnum(:,:,jk)+ rn_beta_con * dcon(:,:,jk) + rn_beta_eke * abs( deke(:,:,jk) ) ) &
       & * gauss_n_2d(:,:) / rn_sppt_stdev
     ENDDO
   ENDIF

   ustar = 0._wp
   vstar = 0._wp
   ustar2= 0._wp
   vstar2= 0._wp
   zinthu= 0._wp
   zinthv= 0._wp

   IF ( ln_stopack_debug .AND. lwp ) THEN
      WRITE(numout,*)
      WRITE(numout,*) ' MIN/MAXVAL of deke : ',MINVAL( rn_beta_eke*deke ), MAXVAL( rn_beta_eke*deke )
      WRITE(numout,*) ' MIN/MAXVAL of psi  : ',MINVAL( psi ), MAXVAL( psi )
      WRITE(numout,*)
   ENDIF

   DO jk=1,jpkm1
    DO jj=2,jpj
     DO ji=2,jpi
        ustar(ji,jj,jk) =   ( psi(ji,jj,jk)-psi(ji,jj-1,jk) )*umask(ji,jj,jk)*umask(ji,jj-1,jk)/ e2u(ji,jj)
        vstar(ji,jj,jk) = - ( psi(ji,jj,jk)-psi(ji-1,jj,jk) )*vmask(ji,jj,jk)*vmask(ji-1,jj,jk)/ e1v(ji,jj)
     ENDDO
    ENDDO
#ifdef NEMO_V42
    ustar2(:,:) = ustar2(:,:) + ustar(:,:,jk)*umask(:,:,jk)*vmask(:,:,jk)*e3u(:,:,jk,Kmm)
    vstar2(:,:) = vstar2(:,:) + vstar(:,:,jk)*umask(:,:,jk)*vmask(:,:,jk)*e3v(:,:,jk,Kmm)
    zinthu(:,:) = zinthu(:,:) +               umask(:,:,jk)*vmask(:,:,jk)*e3u(:,:,jk,Kmm)
    zinthv(:,:) = zinthu(:,:) +               umask(:,:,jk)*vmask(:,:,jk)*e3v(:,:,jk,Kmm)
#else
    ustar2(:,:) = ustar2(:,:) + ustar(:,:,jk)*umask(:,:,jk)*vmask(:,:,jk)*fse3u(:,:,jk)
    vstar2(:,:) = vstar2(:,:) + vstar(:,:,jk)*umask(:,:,jk)*vmask(:,:,jk)*fse3v(:,:,jk)
    zinthu(:,:) = zinthu(:,:) +               umask(:,:,jk)*vmask(:,:,jk)*fse3u(:,:,jk)
    zinthv(:,:) = zinthu(:,:) +               umask(:,:,jk)*vmask(:,:,jk)*fse3v(:,:,jk)
#endif
   ENDDO
   WHERE ( zinthu > 0._wp ) ustar2 = ustar2 / zinthu
   WHERE ( zinthv > 0._wp ) vstar2 = vstar2 / zinthv

   CALL lbc_lnk(_LBCNAME_ ustar,'U',-1._wp)
   CALL lbc_lnk(_LBCNAME_ vstar,'V',-1._wp)
   CALL lbc_lnk(_LBCNAME_ ustar2,'U',-1._wp)
   CALL lbc_lnk(_LBCNAME_ vstar2,'V',-1._wp)

   IF( ln_stopack_diags ) THEN

      mi = MAXVAL(ABS(dnum))
      ma = MAXVAL(ABS(dcon))
      m2 = MAXVAL(ABS(deke))
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ mi)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ m2)
      IF(lwp) WRITE(numdiag,9303) lkt,'DNUM DCON DEKE MX' ,mi,ma,m2

      mi = MINVAL(ustar)
      ma = MAXVAL(ustar)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
      IF(lwp) WRITE(numdiag,9302) lkt,'SKEB U-CURRENT  ' ,mi,ma

      rn_mmax ( jk_skeb_tot ) =  MAX ( MAX( abs(mi), ma), rn_mmax ( jk_skeb_tot ) )

      mi = MINVAL(vstar)
      ma = MAXVAL(vstar)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_min(_LBCNAME_ mi)
      IF(lk_mpp .and. .not. ln_ctl) CALL mpp_max(_LBCNAME_ ma)
      IF(lwp) WRITE(numdiag,9302) lkt,'SKEB V-CURRENT  ' ,mi,ma
9302  FORMAT(' it :', i8, ' ', A16, ' Min: ',d10.3 , ' Max: ',d10.3)
9303  FORMAT(' it :', i8, ' ', A16, ' Max: ',d10.3 , ' Max: ',d10.3,' Max: ',d10.3)

      rn_mmax ( jk_skeb_tot ) =  MAX ( MAX( abs(mi), ma), rn_mmax ( jk_skeb_tot ) )

   ENDIF

   IF ( ln_stopack_debug .AND. lwp ) THEN
      WRITE(numout,*)
      WRITE(numout,*) ' Applying SKEB perturbation at timestep: ', mt
      WRITE(numout,*) ' MAXVAL of u, v perturbations: ',MAXVAL( ABS(ustar) ), MAXVAL( ABS(vstar) )
      WRITE(numout,*)
   ENDIF

   IF ( ln_skeb_tune ) THEN
      IF ( mt .eq. nitend ) THEN
         stun = stun / REAL(ktun,wp)
         CALL skeb_dia_write('skeb_tuning',mt)
      ENDIF
   ELSE

     WHERE( .NOT. abs(ustar) .lt. 1.e+6 ) ustar = 0._wp
     WHERE( .NOT. abs(vstar) .lt. 1.e+6 ) vstar = 0._wp
     WHERE( .NOT. abs(ustar2) .lt. 1.e+6 ) ustar2 = 0._wp
     WHERE( .NOT. abs(vstar2) .lt. 1.e+6 ) vstar2 = 0._wp

     ustar(:,:,:) = ustar(:,:,:) * umask(:,:,:)
     vstar(:,:,:) = vstar(:,:,:) * vmask(:,:,:)

     IF( ln_stopack_debug ) THEN
       DO jk=1,jpkm1
         DO jj=1,jpj
           DO ji=1,jpi
             IF(lwp) THEN
                   IF( .NOT. ABS( ustar(ji,jj,jk) ) .LT. 1._wp ) &
                           WRITE(numout,'(A,4I6,E12.3)') 'uprob',mt,ji,jj,jk,ustar(ji,jj,jk)
                   IF( .NOT. ABS( vstar(ji,jj,jk) ) .LT. 1._wp ) &
                           WRITE(numout,'(A,4I6,E12.3)') 'vprob',mt,ji,jj,jk,vstar(ji,jj,jk)
             ENDIF
           ENDDO
         ENDDO
       ENDDO
     ENDIF

     SELECT CASE ( nn_skst )
#ifdef NEMO_V42

     CASE (-1,-2)

         uu(:,:,:,Kmm) =  uu(:,:,:,Kmm) + ustar(:,:,:)
         vv(:,:,:,Kmm) =  vv(:,:,:,Kmm) + vstar(:,:,:)
         uu(:,:,:,Kbb) =  uu(:,:,:,Kbb) + ustar(:,:,:)
         vv(:,:,:,Kbb) =  vv(:,:,:,Kbb) + vstar(:,:,:)

     CASE (1,2)

         uu(:,:,:,Kmm) =  uu(:,:,:,Kmm) + ustar(:,:,:)
         vv(:,:,:,Kmm) =  vv(:,:,:,Kmm) + vstar(:,:,:)

     CASE (3)

         uu(:,:,:,Krhs) =  uu(:,:,:,Krhs) + ustar(:,:,:)
         vv(:,:,:,Krhs) =  vv(:,:,:,Krhs) + vstar(:,:,:)

#else
     CASE (-1,-2)

         un(:,:,:) = un(:,:,:) + ustar(:,:,:)
         vn(:,:,:) = vn(:,:,:) + vstar(:,:,:)
         ub(:,:,:) = ub(:,:,:) + ustar(:,:,:)
         vb(:,:,:) = vb(:,:,:) + vstar(:,:,:)

     CASE (1,2)

         un(:,:,:) = un(:,:,:) + ustar(:,:,:)
         vn(:,:,:) = vn(:,:,:) + vstar(:,:,:)

     CASE (3)

         ua(:,:,:) = ua(:,:,:) + ustar(:,:,:) * (1._wp / rdt )
         va(:,:,:) = va(:,:,:) + vstar(:,:,:) * (1._wp / rdt )

     CASE (4)

         un_e (:,:) = un_e (:,:) + ustar2(:,:) 
         vn_e (:,:) = vn_e (:,:) + vstar2(:,:) 

     CASE (5)

         ua_b (:,:) = ua_b (:,:) + ustar2(:,:) 
         va_b (:,:) = va_b (:,:) + vstar2(:,:) 
#endif
     CASE DEFAULT
         CALL ctl_stop('nn_skst is not recognized' ) 
     END SELECT

   ENDIF

#ifdef key_iomput
   IF (iom_use('ustar_skeb') ) CALL iom_put( "ustar_skeb" , ustar)
   IF (iom_use('vstar_skeb') ) CALL iom_put( "vstar_skeb" , vstar)
#endif

   IF ( mt .eq. nitend ) THEN
     IF(ALLOCATED(dnum)) DEALLOCATE ( dnum )
     IF(ALLOCATED(dcon)) DEALLOCATE ( dcon )
     IF(ALLOCATED(deke)) DEALLOCATE ( deke )
     IF (ln_stopack_diags .AND. lwp) CALL stopack_report
   ENDIF

   END SUBROUTINE

   !!======================================================================
   !! Random number generator, used in NEMO stochastic parameterization
   !!
   !!=====================================================================
   !! History :  3.3  ! 2011-10 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! The module is based on (and includes) the
   !! 64-bit KISS (Keep It Simple Stupid) random number generator
   !! distributed by George Marsaglia :
   !! http://groups.google.com/group/comp.lang.fortran/
   !!        browse_thread/thread/a85bf5f2a97f5a55
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   kiss          : 64-bit KISS random number generator (period ~ 2^250)
   !!   kiss_seed     : Define seeds for KISS random number generator
   !!   kiss_state    : Get current state of KISS random number generator
   !!   kiss_save     : Save current state of KISS (for future restart)
   !!   kiss_load     : Load the saved state of KISS
   !!   kiss_reset    : Reset the default seeds
   !!   kiss_check    : Check the KISS pseudo-random sequence
   !!   kiss_uniform  : Real random numbers with uniform distribution in [0,1]
   !!   kiss_gaussian : Real random numbers with Gaussian distribution N(0,1)
   !!   kiss_gamma    : Real random numbers with Gamma distribution Gamma(k,1)
   !!   kiss_sample   : Select a random sample from a set of integers
   !!
   !!   ---CURRENTLY NOT USED IN NEMO :
   !!   kiss_save, kiss_load, kiss_check, kiss_gamma, kiss_sample
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynhpg.F90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   FUNCTION kiss()
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION kiss  ***
      !!
      !! ** Purpose :   64-bit KISS random number generator
      !!
      !! ** Method  :   combine several random number generators:
      !!                (1) Xorshift (XSH), period 2^64-1,
      !!                (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
      !!                (3) Congruential generator (CNG), period 2^64.
      !!
      !!                overall period:
      !!                (2^250+2^192+2^64-2^186-2^129)/6
      !!                            ~= 2^(247.42) or 10^(74.48)
      !!
      !!                set your own seeds with 'kiss_seed'
      ! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: kiss, t

      t = ISHFT(x,58) + w
      IF (s(x).eq.s(t)) THEN
         w = ISHFT(x,-6) + s(x)
      ELSE
         w = ISHFT(x,-6) + 1 - s(x+t)
      ENDIF
      x = t + x
      y = m( m( m(y,13_8), -17_8 ), 43_8 )
      z = 6906969069_8 * z + 1234567_8

      kiss = x + y + z

      CONTAINS

         FUNCTION s(k)
            INTEGER(KIND=i8) :: s, k
            s = ISHFT(k,-63)
         END FUNCTION s

         FUNCTION m(k, n)
            INTEGER(KIND=i8) :: m, k, n
            m =  IEOR(k, ISHFT(k, n) )
         END FUNCTION m

   END FUNCTION kiss


   SUBROUTINE kiss_seed(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_seed  ***
      !!
      !! ** Purpose :   Define seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      IF(lwp) WRITE(numout,*) ' *** kiss_seed called with args:'
      IF(lwp) WRITE(numout,*) ' *** ix = ',ix
      IF(lwp) WRITE(numout,*) ' *** iy = ',iy
      IF(lwp) WRITE(numout,*) ' *** iz = ',iz
      IF(lwp) WRITE(numout,*) ' *** iw = ',iw

      x = ix
      y = iy
      z = iz
      w = iw

   END SUBROUTINE kiss_seed


   SUBROUTINE kiss_state(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_state  ***
      !!
      !! ** Purpose :   Get current state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      ix = x
      iy = y
      iz = z
      iw = w

   END SUBROUTINE kiss_state


   SUBROUTINE kiss_reset()
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_reset  ***
      !!
      !! ** Purpose :   Reset the default seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE

      x=1234567890987654321_8
      y=362436362436362436_8
      z=1066149217761810_8
      w=123456123456123456_8

    END SUBROUTINE kiss_reset

   SUBROUTINE kiss_uniform(uran)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_uniform  ***
      !!
      !! ** Purpose :   Real random numbers with uniform distribution in [0,1]
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: uran

      uran = half * ( one + REAL(kiss(),wp) / huge64 )

   END SUBROUTINE kiss_uniform


   SUBROUTINE kiss_gaussian(gran)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_gaussian  ***
      !!
      !! ** Purpose :   Real random numbers with Gaussian distribution N(0,1)
      !!
      !! ** Method  :   Generate 2 new Gaussian draws (gran1 and gran2)
      !!                from 2 uniform draws on [-1,1] (u1 and u2),
      !!                using the Marsaglia polar method
      !!                (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: gran, u1, u2, rsq, fac

      IF (ig.EQ.1) THEN
         rsq = two
         DO WHILE ( (rsq.GE.one).OR. (rsq.EQ.zero) )
            u1 = REAL(kiss(),wp) / huge64
            u2 = REAL(kiss(),wp) / huge64
            rsq = u1*u1 + u2*u2
         ENDDO
         fac = SQRT(-two*LOG(rsq)/rsq)
         gran1 = u1 * fac
         gran2 = u2 * fac
      ENDIF

      ! Output one of the 2 draws
      IF (ig.EQ.1) THEN
         gran = gran1 ; ig = 2
      ELSE
         gran = gran2 ; ig = 1
      ENDIF

   END SUBROUTINE kiss_gaussian


   SUBROUTINE kiss_gamma(gamr,k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_gamma  ***
      !!
      !! ** Purpose :   Real random numbers with Gamma distribution Gamma(k,1)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp), PARAMETER :: p1 = 4.5_8
      REAL(KIND=wp), PARAMETER :: p2 = 2.50407739677627_8  ! 1+LOG(9/2)
      REAL(KIND=wp), PARAMETER :: p3 = 1.38629436111989_8  ! LOG(4)
      REAL(KIND=wp) :: gamr, k, u1, u2, b, c, d, xx, yy, zz, rr, ee
      LOGICAL :: accepted

      IF (k.GT.one) THEN
         ! Cheng's rejection algorithm
         ! (see Devroye, Non-Uniform Random Variate Generation, p. 413)
         b = k - p3 ; d = SQRT(two*k-one) ; c = k + d

         accepted=.FALSE.
         DO WHILE (.NOT.accepted)
            CALL kiss_uniform(u1)
            yy = LOG(u1/(one-u1)) / d  ! Mistake in Devroye: "* k" instead of "/ d"
            xx = k * EXP(yy)
            rr = b + c * yy - xx
            CALL kiss_uniform(u2)
            zz = u1 * u1 * u2

            accepted = rr .GE. (zz*p1-p2)
            IF (.NOT.accepted) accepted =  rr .GE. LOG(zz)
         ENDDO

         gamr = xx

      ELSEIF (k.LT.one) THEN
        ! Rejection from the Weibull density
        ! (see Devroye, Non-Uniform Random Variate Generation, p. 415)
        c = one/k ; d = (one-k) * EXP( (k/(one-k)) * LOG(k) )

        accepted=.FALSE.
        DO WHILE (.NOT.accepted)
           CALL kiss_uniform(u1)
           zz = -LOG(u1)
           xx = EXP( c * LOG(zz) )
           CALL kiss_uniform(u2)
           ee = -LOG(u2)

           accepted = (zz+ee) .GE. (d+xx)  ! Mistake in Devroye: "LE" instead of "GE"
        ENDDO

        gamr = xx

      ELSE
         ! Exponential distribution
         CALL kiss_uniform(u1)
         gamr = -LOG(u1)

      ENDIF

   END SUBROUTINE kiss_gamma


   SUBROUTINE kiss_sample(a,n,k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_sample  ***
      !!
      !! ** Purpose :   Select a random sample of size k from a set of n integers
      !!
      !! ** Method  :   The sample is output in the first k elements of a
      !!                Set k equal to n to obtain a random permutation
      !!                  of the whole set of integers
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8), DIMENSION(:) :: a
      INTEGER(KIND=i8) :: n, k, i, j, atmp
      REAL(KIND=wp) :: uran

      ! Select the sample using the swapping method
      ! (see Devroye, Non-Uniform Random Variate Generation, p. 612)
      DO i=1,k
         ! Randomly select the swapping element between i and n (inclusive)
         CALL kiss_uniform(uran)
         j = i - 1 + CEILING( REAL(n-i+1,8) * uran )
         ! Swap elements i and j
         atmp = a(i) ; a(i) = a(j) ; a(j) = atmp
      ENDDO

   END SUBROUTINE kiss_sample

#ifdef NEMO_V34
SUBROUTINE ctl_nam(ier,cstr,lout)
INTEGER, INTENT(IN)            :: ier
LOGICAL, INTENT(IN)            :: lout
CHARACTER(LEN=*),INTENT(IN)    :: cstr
IF(lout) WRITE(numout,'(A,I4,X,A)') 'Error ',ier,TRIM(cstr)
END SUBROUTINE
#endif

   SUBROUTINE skeb_dia_write( cdfile_name, kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      INTEGER           , INTENT( in ) ::   kt               ! ocean time-step index
      !! 
      CHARACTER (len=32) :: clname
      CHARACTER (len=40) :: clop
      INTEGER  ::   id_i , nz_i, nh_i       
      INTEGER, DIMENSION(1) ::   idex             ! local workspace
      REAL(wp) ::   zsto, zout, zmax, zjulian, zdt
      INTEGER  ::   jstream
      !!----------------------------------------------------------------------
      ! 
      ! 0. Initialisation
      ! -----------------

      ! Define name, frequency of output and means
      clname = cdfile_name
      IF( .NOT. Agrif_Root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      zdt  = rdt
      zsto = rdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', clname, '.nc'


      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      ! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )         ! time axis 
      zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
      CALL histbeg( clname, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, nit000-1, zjulian, zdt, nh_i, id_i, domain_id=nidom, snc4chunks=snc4set ) ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
          "m", jpk, gdept_0, nz_i, "down")

      ! Declare all the output fields as NetCDF variables

      CALL histdef( id_i, "khtuning", "SKEB kh tuning"    , "m/s"    ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "kctuning", "SKEB kc tuning"    , "m/s"    ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "ketuning", "SKEB ke tuning"    , "m/s"    ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
      !
      CALL histend( id_i, snc4chunks=snc4set )
      ! 2. Start writing data
      ! ---------------------
      ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et idex la liste des indices a sortir
      idex(1) = 1   ! init to avoid compil warning

      ! Write all fields on T grid
      CALL histwrite( id_i, "khtuning", kt, stun(:,:,1)      , jpi*jpj    , idex )    ! SKEB tuning
      CALL histwrite( id_i, "kctuning", kt, stun(:,:,2)      , jpi*jpj    , idex )    ! SKEB tuning
      CALL histwrite( id_i, "ketuning", kt, stun(:,:,3)      , jpi*jpj    , idex )    ! SKEB tuning

      ! 3. Close the file
      ! -----------------
      CALL histclo( id_i )
      ! 
   END SUBROUTINE skeb_dia_write

END MODULE stopack
