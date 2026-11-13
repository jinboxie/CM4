module oro_drag_mod
!=======================================================================
!
!  Orographic drag package for GFDL AM4/ESM4 (FMS physics).
!
!  Four separately switchable components:
!    do_od_ls : large-scale orographic gravity wave drag  (Xie et al., 2020, JAMES)
!    do_od_bl : flow-blocking drag                        (Xie et al., 2020, JAMES)
!    do_od_ss : small-scale orographic gravity wave drag  (Tsiringakis et al., 2017, QJRMS)
!    do_od_fd : turbulent orographic form drag            (Beljaars et al., 2004, QJRMS)
!
!  Origin: E3SM EAM, components/eam/src/physics/cam/od_common.F90
!          (E3SM-Project/E3SM, commit 1bd1ee77ff5af2bcc7feded06faaf90ea400c8a5).
!  The E3SM routines od_gsd, od2d and dxygrid are kept in their original
!  form and with their original names (copied verbatim at the end of this
!  file).  Only the CAM-specific parts of od_common (oro_drag_readnl,
!  oro_drag_register, oro_drag_init, oro_drag_interface, grid_size,
!  pblh_get_level_idx) are replaced by the FMS routines oro_drag_init and
!  oro_drag, which call od_gsd exactly as oro_drag_interface does in E3SM.
!
!  How the E3SM inputs are supplied in ESM4:
!    OC, OA, OL        -> restart file INPUT/oro_drag.res.nc (see below)
!    sgh               -> sub-grid height std. dev. from INPUT/mg_drag.res.nc
!                         ('ghprime', the field mg_drag_mod uses) or from
!                         INPUT/oro_drag.res.nc ('SGH'), chosen by sgh_restart;
!                         used by the large-scale GWD, blocking and small-scale GWD
!    sgh30             -> std. dev. of ~1 km (30") sub-grid height, 'SGH30' from
!                         INPUT/oro_drag.res.nc; used by the turbulent form drag
!                         only (read when do_od_fd = .true.)
!    nm (Brunt-Vaisala)-> computed here, same algorithm as E3SM gw_prof
!    pblh              -> z_pbl from damping_driver (pbltop of previous step)
!    landfrac          -> frac_land from physics_driver_down
!    bulk Richardson   -> diagnosed here from lowest model level and t_surf
!                         (only its sign is used by the kernel)
!    dx, dy            -> sqrt(cell area)  (E3SM also uses a square cell)
!    etamid (znu)      -> pref(k)/pref(nlev+1) from the reference profile
!
!  od_gsd uses a single std. dev. (var2d) for all components.  To give the
!  form drag sgh30 while keeping od_gsd/od2d unchanged, od_gsd is called
!  twice: (1) var2d=sgh with ls/bl/ss, (2) var2d=sgh30 with fd only.  This
!  is exact: in od2d the form-drag block depends only on the winds, layer
!  heights, land fraction and var, and all outputs are reset on entry.
!
!  As in E3SM gw_drag.F90, the total tendency is multiplied by the land
!  fraction (scale_by_land_frac) and the lost kinetic energy is returned as
!  heating (do_dissipative_heating).
!
!  Input is read through restart files, as in mg_drag_mod (fms2_io domain
!  restart, so on the cubed sphere the files are *.res.tile[1-6].nc).
!  Dimensions as shown by ncdump (C order; Fortran order is reversed):
!
!    INPUT/oro_drag.res.nc
!      OC   (Time, yaxis_1, xaxis_1)       orographic convexity             [1]
!      OA   (Time, nOA, yaxis_1, xaxis_1)  orographic asymmetry; the first 2
!                                          of nOA (>= 2) are the (east, north)
!                                          components (E3SM files have nOA=3)[1]
!      OL   (Time, nOL, yaxis_1, xaxis_1)  effective orographic length,
!                                          nOL = 180 one-degree bins,
!                                          counter-clockwise from east      [1]
!      SGH  (Time, yaxis_1, xaxis_1)       sub-grid height std. dev. [m]
!                                          (only used if sgh_restart='oro_drag')
!      SGH30(Time, yaxis_1, xaxis_1)       std. dev. of ~1 km sub-grid height [m]
!                                          (required if do_od_fd = .true.)
!    INPUT/mg_drag.res.nc
!      ghprime (Time, yaxis_1, xaxis_1)    sub-grid height std. dev. [m]
!                                          (used if sgh_restart='mg_drag')
!
!  Directions must be relative to geographic east/north at each cell centre,
!  because the physics u,v in ESM4 are earth-relative.
!
!  Sub-grid height standard deviation, same options as mg_drag_nml:
!    source_of_sgsmtn = 'input'          from the restart file set by sgh_restart
!                       'computed'       get_topog_stdev (topography_mod, needs
!                                        INPUT/navy_topography.data.nc)
!                       'input/computed' restart if present, else computed
!    sgh_restart      = 'mg_drag'        INPUT/mg_drag.res.nc,  'ghprime' (default)
!                       'oro_drag'       INPUT/oro_drag.res.nc, 'SGH'
!
!  Restart output (oro_drag_restart, called from damping_driver at the end of
!  the run and at intermediate restarts), same layout as the input:
!    RESTART/oro_drag.res.nc   OC, OA, OL (+ SGH if sgh_restart='oro_drag',
!                              + SGH30 if it was read)
!    RESTART/mg_drag.res.nc    ghprime, only if sgh_restart='mg_drag' and
!                              mg_drag_mod is not active (otherwise mg_drag
!                              writes it itself)
!
!=======================================================================

use          mpp_mod, only: input_nml_file
use  mpp_domains_mod, only: domain2D
use          fms_mod, only: error_mesg, FATAL, NOTE, WARNING, mpp_pe, mpp_root_pe, &
                            stdlog, write_version_number, check_nml_error
use      fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file,         &
                            register_axis, register_field, register_restart_field, &
                            read_restart, write_restart, write_data, unlimited,   &
                            variable_exists, get_global_io_domain_indices,       &
                            dimension_exists, get_dimension_size
use   topography_mod, only: get_topog_stdev
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use time_manager_mod, only: time_type
use    constants_mod, only: GRAV, RDGAS, RVGAS, CP_AIR, PI, KAPPA

implicit none
private

public :: oro_drag_init, oro_drag, oro_drag_end, oro_drag_restart

! r8 is the default real kind, so the E3SM kernel below compiles at whatever
! precision ESM4 is built with (normally 8-byte via -r8).
integer, parameter :: r8 = kind(1.0)

integer, parameter :: ndir_asymmetry = 2    ! OA components used by od_gsd (first 2 of nOA)
integer, parameter :: ndir_efflength = 180  ! OL direction bins (nOL)

!---------------------------- namelist ---------------------------------
logical :: do_od_ls = .true.      ! large-scale orographic GWD
logical :: do_od_bl = .true.      ! flow-blocking drag
logical :: do_od_ss = .true.      ! small-scale orographic GWD
logical :: do_od_fd = .false.     ! turbulent orographic form drag (off in E3SM too)
real    :: od_ls_ncleff  = 3.0    ! tuning factor, oGWD   (E3SM default)
real    :: od_bl_ncd     = 3.0    ! tuning factor, FBD    (E3SM default)
real    :: od_ss_sncleff = 1.0    ! tuning factor, sGWD   (E3SM default)
logical :: scale_by_land_frac     = .true.  ! E3SM multiplies tendencies by landfrac
logical :: do_dissipative_heating = .true.  ! return lost KE as heating
character(len=128) :: source_of_sgsmtn = 'input'    ! as in mg_drag_nml
character(len=32)  :: sgh_restart      = 'mg_drag'  ! 'mg_drag' or 'oro_drag'

namelist /oro_drag_nml/ do_od_ls, do_od_bl, do_od_ss, do_od_fd,          &
                        od_ls_ncleff, od_bl_ncd, od_ss_sncleff,          &
                        scale_by_land_frac, do_dissipative_heating,      &
                        source_of_sgsmtn, sgh_restart

!---------------------------- module data ------------------------------
real, allocatable, dimension(:,:)   :: sgh         ! sub-grid height std. dev. [m]
real, allocatable, dimension(:,:)   :: sgh30       ! ~1 km sub-grid height std. dev. [m] (form drag)
real, allocatable, dimension(:,:)   :: oc          ! (nlon,nlat) compute domain
real, allocatable, dimension(:,:,:) :: oa          ! (nlon,nlat,nOA)
real, allocatable, dimension(:,:,:) :: ol          ! (nlon,nlat,nOL)
real, allocatable, dimension(:)     :: znu         ! reference eta, bottom-up (E3SM etamid(pver:1:-1))
real, parameter :: p00 = 1.0e5                     ! Exner reference pressure [Pa]

integer :: id_udt_oro, id_vdt_oro, id_tdt_diss_oro
integer :: id_udt_ls, id_vdt_ls, id_udt_bl, id_vdt_bl,  &
           id_udt_ss, id_vdt_ss, id_udt_fd, id_vdt_fd
integer :: id_taux_ls, id_tauy_ls, id_taux_bl, id_tauy_bl, &
           id_taux_ss, id_tauy_ss, id_taux_fd, id_tauy_fd
integer :: id_sgh, id_sgh30, id_oc, id_oa1, id_oa2
real    :: missing_value = -999.
character(len=8), parameter :: mod_name = 'oro_drag'

type(domain2D), pointer :: od_domain => NULL()
logical :: mg_on = .false.         ! mg_drag_mod active: it owns mg_drag.res.nc
logical :: have_sgh30 = .false.    ! SGH30 read from oro_drag.res.nc

logical :: module_is_initialized = .false.

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

contains

!#######################################################################

subroutine oro_drag_init (domain, lonb, latb, pref, axes, Time, hprime, mg_drag_active)

type(domain2D), target, intent(in) :: domain !< Atmosphere domain
real,            intent(in) :: lonb(:,:), latb(:,:)
real,            intent(in) :: pref(:)       ! (nlev+1): full-level ref. p, then ref. ps
integer,         intent(in) :: axes(4)
type(time_type), intent(in) :: Time
real, intent(out), dimension(:,:), optional :: hprime   ! returns sgh, like mg_drag_init
logical, intent(in), optional :: mg_drag_active          ! mg_drag_mod is also on

type(FmsNetcdfDomainFile_t) :: Od_restart, Mg_restart !< Fms2io domain decomposed fileobjs
integer :: io, ierr, logunit, nlon, nlat, nlev, k, noa, nol
character(len=256) :: msg
logical :: used, answer, sgh_found

  if (module_is_initialized) return

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------
  read (input_nml_file, nml=oro_drag_nml, iostat=io)
  ierr = check_nml_error(io, 'oro_drag_nml')
  call write_version_number (version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe()) write (logunit, nml=oro_drag_nml)

  if (.not.(do_od_ls .or. do_od_bl .or. do_od_ss .or. do_od_fd)) &
     call error_mesg ('oro_drag_init', 'all four components are switched off', WARNING)

  if (trim(sgh_restart) /= 'mg_drag' .and. trim(sgh_restart) /= 'oro_drag') &
     call error_mesg ('oro_drag_init', 'sgh_restart="'//trim(sgh_restart)//  &
                      '" must be "mg_drag" or "oro_drag"', FATAL)

  od_domain => domain
  if (present(mg_drag_active)) mg_on = mg_drag_active

!---------------------------------------------------------------------
! --- Allocate storage
!---------------------------------------------------------------------
  nlon = size(lonb,1) - 1
  nlat = size(latb,2) - 1
  nlev = size(pref)   - 1

  allocate (sgh(nlon,nlat)) ; sgh = 0.0
  allocate (sgh30(nlon,nlat)) ; sgh30 = 0.0
  allocate (oc (nlon,nlat)) ; oc  = 0.0

!--- reference eta at full levels, bottom-up, used by od_gsd for kpblmax
  allocate (znu(nlev))
  do k = 1, nlev
     znu(k) = pref(nlev+1-k)/pref(nlev+1)
  enddo

  module_is_initialized = .true.

!---------------------------------------------------------------------
! --- Input OC, OA, OL (and SGH) from INPUT/oro_drag.res.nc
!---------------------------------------------------------------------
  sgh_found = .false.
  if (open_file(Od_restart, "INPUT/oro_drag.res.nc", "read", od_domain, is_restart=.true.)) then
     if (mpp_pe() == mpp_root_pe()) call error_mesg ('oro_drag_init', &
          'Reading NetCDF formatted restart file: INPUT/oro_drag.res.nc', NOTE)
!--- direction dimensions of OA and OL, from the file
     if (.not.(dimension_exists(Od_restart, "nOA") .and. dimension_exists(Od_restart, "nOL"))) &
          call error_mesg ('oro_drag_init', &
          'INPUT/oro_drag.res.nc must have dimensions nOA and nOL', FATAL)
     call get_dimension_size (Od_restart, "nOA", noa)
     call get_dimension_size (Od_restart, "nOL", nol)
     if (noa < ndir_asymmetry .or. nol /= ndir_efflength) then
        write (msg,'(a,i0,a,i0,a,i0)') 'nOA = ', noa, ' (needs >= 2), nOL = ', nol, &
                                       ' (needs ', ndir_efflength
        call error_mesg ('oro_drag_init', trim(msg)//') in INPUT/oro_drag.res.nc', FATAL)
     endif
     allocate (oa (nlon,nlat,noa)) ; oa = 0.0
     allocate (ol (nlon,nlat,nol)) ; ol = 0.0
     if (trim(sgh_restart) == 'oro_drag' .and. trim(source_of_sgsmtn) /= 'computed') &
          sgh_found = variable_exists(Od_restart, 'SGH')
     if (do_od_fd) then
        have_sgh30 = variable_exists(Od_restart, 'SGH30')
        if (.not. have_sgh30) call error_mesg ('oro_drag_init', &
             'do_od_fd = .true. needs SGH30 in INPUT/oro_drag.res.nc', FATAL)
     endif
     call oro_register_restart (Od_restart, with_sgh=sgh_found, with_sgh30=have_sgh30)
     call read_restart (Od_restart)
     call close_file (Od_restart)
  else
     call error_mesg ('oro_drag_init', &
          'INPUT/oro_drag.res.nc (OC, OA, OL) does not exist', FATAL)
  endif

!---------------------------------------------------------------------
! --- Input sub-grid height std. dev. (same logic as mg_drag_init)
!---------------------------------------------------------------------
  if ( trim(source_of_sgsmtn) == 'computed' ) then
    answer = get_topog_stdev ( lonb, latb, sgh )
    if ( .not.answer ) then
      call error_mesg('oro_drag_init','source_of_sgsmtn="'//trim(source_of_sgsmtn)//'"'// &
                      ', but topography data file does not exist', FATAL)
    endif
  else if ( trim(source_of_sgsmtn) == 'input' .or. trim(source_of_sgsmtn) == 'input/computed' ) then
    if ( trim(sgh_restart) == 'mg_drag' ) then
      if (open_file(Mg_restart, "INPUT/mg_drag.res.nc", "read", od_domain, is_restart=.true.)) then
        if (mpp_pe() == mpp_root_pe()) call error_mesg ('oro_drag_init', &
             'Reading NetCDF formatted restart file: INPUT/mg_drag.res.nc', NOTE)
        call mg_register_restart (Mg_restart)
        call read_restart (Mg_restart)
        call close_file (Mg_restart)
        sgh_found = .true.
      endif
    endif
    if (.not. sgh_found) then
      if (trim(source_of_sgsmtn) == 'input') then
        call error_mesg ('oro_drag_init','source_of_sgsmtn="'//trim(source_of_sgsmtn)//'"'// &
             ', but no sub-grid height std. dev. in INPUT/'//trim(sgh_restart)//'.res.nc', FATAL)
      else
        answer = get_topog_stdev ( lonb, latb, sgh )
        if ( .not.answer ) then
          call error_mesg('oro_drag_init','source_of_sgsmtn="'//trim(source_of_sgsmtn)//'"'// &
                          ', but topography data file does not exist', FATAL)
        endif
      endif
    endif
  else
    call error_mesg ('oro_drag_init','"'//trim(source_of_sgsmtn)//'"'// &
          ' is not a valid value for source_of_sgsmtn', FATAL)
  endif

!--- E3SM fills missing points with NaN; here bad or fill values become 0,
!    which switches the drag off at that point.
  call clean_field2 (sgh) ; call clean_field2 (oc) ; call clean_field2 (sgh30)
  call clean_field3 (oa)  ; call clean_field3 (ol)
  sgh   = max(sgh,   0.)
  sgh30 = max(sgh30, 0.)

! return sub-grid scale topography?
  if (present(hprime)) hprime = sgh

!--- diagnostics
  id_udt_oro = register_diag_field (mod_name, 'udt_oro', axes(1:3), Time,         &
       'u tendency from orographic drag (all components)', 'm/s2', missing_value=missing_value)
  id_vdt_oro = register_diag_field (mod_name, 'vdt_oro', axes(1:3), Time,         &
       'v tendency from orographic drag (all components)', 'm/s2', missing_value=missing_value)
  id_tdt_diss_oro = register_diag_field (mod_name, 'tdt_diss_oro', axes(1:3), Time, &
       'dissipative heating from orographic drag', 'K/s', missing_value=missing_value)

  id_udt_ls = register_diag_field (mod_name, 'udt_od_ls', axes(1:3), Time, &
       'u tendency, large-scale orographic GWD', 'm/s2', missing_value=missing_value)
  id_vdt_ls = register_diag_field (mod_name, 'vdt_od_ls', axes(1:3), Time, &
       'v tendency, large-scale orographic GWD', 'm/s2', missing_value=missing_value)
  id_udt_bl = register_diag_field (mod_name, 'udt_od_bl', axes(1:3), Time, &
       'u tendency, flow-blocking drag', 'm/s2', missing_value=missing_value)
  id_vdt_bl = register_diag_field (mod_name, 'vdt_od_bl', axes(1:3), Time, &
       'v tendency, flow-blocking drag', 'm/s2', missing_value=missing_value)
  id_udt_ss = register_diag_field (mod_name, 'udt_od_ss', axes(1:3), Time, &
       'u tendency, small-scale orographic GWD', 'm/s2', missing_value=missing_value)
  id_vdt_ss = register_diag_field (mod_name, 'vdt_od_ss', axes(1:3), Time, &
       'v tendency, small-scale orographic GWD', 'm/s2', missing_value=missing_value)
  id_udt_fd = register_diag_field (mod_name, 'udt_od_fd', axes(1:3), Time, &
       'u tendency, turbulent orographic form drag', 'm/s2', missing_value=missing_value)
  id_vdt_fd = register_diag_field (mod_name, 'vdt_od_fd', axes(1:3), Time, &
       'v tendency, turbulent orographic form drag', 'm/s2', missing_value=missing_value)

  id_taux_ls = register_diag_field (mod_name, 'taux_od_ls', axes(1:2), Time, &
       'column zonal stress, large-scale orographic GWD', 'N/m2', missing_value=missing_value)
  id_tauy_ls = register_diag_field (mod_name, 'tauy_od_ls', axes(1:2), Time, &
       'column meridional stress, large-scale orographic GWD', 'N/m2', missing_value=missing_value)
  id_taux_bl = register_diag_field (mod_name, 'taux_od_bl', axes(1:2), Time, &
       'column zonal stress, flow-blocking drag', 'N/m2', missing_value=missing_value)
  id_tauy_bl = register_diag_field (mod_name, 'tauy_od_bl', axes(1:2), Time, &
       'column meridional stress, flow-blocking drag', 'N/m2', missing_value=missing_value)
  id_taux_ss = register_diag_field (mod_name, 'taux_od_ss', axes(1:2), Time, &
       'column zonal stress, small-scale orographic GWD', 'N/m2', missing_value=missing_value)
  id_tauy_ss = register_diag_field (mod_name, 'tauy_od_ss', axes(1:2), Time, &
       'column meridional stress, small-scale orographic GWD', 'N/m2', missing_value=missing_value)
  id_taux_fd = register_diag_field (mod_name, 'taux_od_fd', axes(1:2), Time, &
       'column zonal stress, turbulent orographic form drag', 'N/m2', missing_value=missing_value)
  id_tauy_fd = register_diag_field (mod_name, 'tauy_od_fd', axes(1:2), Time, &
       'column meridional stress, turbulent orographic form drag', 'N/m2', missing_value=missing_value)

  id_sgh = register_static_field (mod_name, 'od_sgh', axes(1:2), &
       'std. dev. of sub-grid orography (oro_drag input)', 'm')
  id_sgh30 = register_static_field (mod_name, 'od_sgh30', axes(1:2), &
       'std. dev. of ~1 km sub-grid orography (form drag input)', 'm')
  id_oc  = register_static_field (mod_name, 'od_oc',  axes(1:2), &
       'orographic convexity (oro_drag input)', '1')
  id_oa1 = register_static_field (mod_name, 'od_oa1', axes(1:2), &
       'orographic asymmetry, east component (oro_drag input)', '1')
  id_oa2 = register_static_field (mod_name, 'od_oa2', axes(1:2), &
       'orographic asymmetry, north component (oro_drag input)', '1')

  if (id_sgh > 0) used = send_data (id_sgh, sgh, Time)
  if (id_sgh30 > 0 .and. have_sgh30) used = send_data (id_sgh30, sgh30, Time)
  if (id_oc  > 0) used = send_data (id_oc,  oc,  Time)
  if (id_oa1 > 0) used = send_data (id_oa1, oa(:,:,1), Time)
  if (id_oa2 > 0) used = send_data (id_oa2, oa(:,:,2), Time)

end subroutine oro_drag_init

!#######################################################################
! Compute orographic drag tendencies for one physics window.
! All 3-D arrays are top-down (k=1 at model top), as everywhere in ESM4
! physics; columns are flipped to bottom-up for the E3SM/WRF kernel.

subroutine oro_drag (is, js, Time, delt, area, pfull, phalf, zfull, zhalf, &
                     u, v, t, q, z_pbl, frac_land, t_surf,                 &
                     utnd, vtnd, ttnd)

integer,         intent(in) :: is, js
type(time_type), intent(in) :: Time
real,            intent(in) :: delt
real, intent(in),  dimension(:,:)   :: area, z_pbl, frac_land, t_surf
real, intent(in),  dimension(:,:,:) :: pfull, phalf, zfull, zhalf, u, v, t, q
real, intent(out), dimension(:,:,:) :: utnd, vtnd, ttnd

integer :: idim, jdim, nlev, i, j, k, kr, jd, ie
real    :: fv, thv1, thvs, wsp2, z1, zb_top, zb_bot
real, dimension(size(u,1),size(u,3))   :: u_bu, v_bu, t_bu, q_bu, p_bu, pk_bu, &
                                          z_bu, dz_bu, nm_bu, nm_td,           &
                                          ut, vt, tht,                         &
                                          dx_ls, dy_ls, dx_bl, dy_bl,          &
                                          dx_ss, dy_ss, dx_fd, dy_fd
real, dimension(size(u,1),size(u,3)+1) :: pi_bu
real, dimension(size(u,1),size(u,3))   :: pk_e     ! Exner, passed as pi3d
real, dimension(size(u,1),size(u,3))   :: ut2, vt2, tht2,          & ! 2nd (form drag) call
                                          dx_x1, dy_x1, dx_x2, dy_x2, dx_x3, dy_x3
real, dimension(size(u,1)) :: su_x1, sv_x1, su_x2, sv_x2, su_x3, sv_x3, var30_c
real, dimension(size(u,1)) :: dxm, dym, hpbl, br, xland,                       &
                              su_ls, sv_ls, su_bl, sv_bl, su_ss, sv_ss, su_fd, sv_fd, &
                              var_c, oc_c
real, dimension(size(u,1),ndir_asymmetry) :: oa_c
real, dimension(size(u,1),ndir_efflength) :: ol_c
integer, dimension(size(u,1)) :: kpbl
real, dimension(size(u,1),size(u,2),size(u,3)) :: c_ls_u, c_ls_v, c_bl_u, c_bl_v, &
                                                  c_ss_u, c_ss_v, c_fd_u, c_fd_v
real, dimension(size(u,1),size(u,2)) :: s_ls_u, s_ls_v, s_bl_u, s_bl_v, &
                                        s_ss_u, s_ss_v, s_fd_u, s_fd_v
logical :: used

  if (.not. module_is_initialized) call error_mesg ('oro_drag', &
       'oro_drag_init must be called first', FATAL)

  idim = size(u,1) ; jdim = size(u,2) ; nlev = size(u,3)
  ie   = is + idim - 1
  fv   = RVGAS/RDGAS - 1.

  do j = 1, jdim
     jd = js + j - 1

!--- flip to bottom-up (k=1 lowest layer)
     do k = 1, nlev
        kr = nlev + 1 - k
        u_bu (:,k) = u(:,j,kr)
        v_bu (:,k) = v(:,j,kr)
        t_bu (:,k) = t(:,j,kr)
        q_bu (:,k) = q(:,j,kr)
        p_bu (:,k) = pfull(:,j,kr)
        pk_bu(:,k) = (pfull(:,j,kr)/p00)**KAPPA                 ! Exner function
        z_bu (:,k) = zhalf(:,j,kr+1)                            ! layer bottom height
        dz_bu(:,k) = zhalf(:,j,kr) - zhalf(:,j,kr+1)            ! layer thickness
     enddo
     do k = 1, nlev+1
        pi_bu(:,k) = phalf(:,j,nlev+2-k)                        ! interface pressure
     enddo

!--- mid-level Brunt-Vaisala frequency (E3SM gw_prof algorithm)
     call bv_freq (t(:,j,:), pfull(:,j,:), phalf(:,j,:), nm_td)
     do k = 1, nlev
        nm_bu(:,k) = nm_td(:,nlev+1-k)
     enddo

     do i = 1, idim
!--- PBL top: bottom-up index of the layer containing z_pbl
!    (E3SM pblh_get_level_idx; z_pbl is -999 before the first turbulence call)
        hpbl(i) = max(z_pbl(i,j), 0.)
        kpbl(i) = 1
        do k = nlev, 2, -1
           zb_bot = zhalf(i,j,k+1) - zhalf(i,j,nlev+1)
           zb_top = zhalf(i,j,k)   - zhalf(i,j,nlev+1)
           if (hpbl(i) >= zb_bot .and. hpbl(i) < zb_top) then
              kpbl(i) = nlev + 1 - k
              exit
           endif
        enddo

!--- bulk Richardson number of the surface layer (kernel only tests br > 0)
        thv1 = t(i,j,nlev)*(1. + fv*q(i,j,nlev))*(p00/pfull(i,j,nlev))**KAPPA
        thvs = t_surf(i,j)*(1. + fv*q(i,j,nlev))*(p00/phalf(i,j,nlev+1))**KAPPA
        z1   = max(zfull(i,j,nlev) - zhalf(i,j,nlev+1), 1.)
        wsp2 = max(u(i,j,nlev)**2 + v(i,j,nlev)**2, 1.e-2)
        br(i) = GRAV*z1*(thv1 - thvs)/(0.5*(thv1 + thvs)*wsp2)

        xland(i) = frac_land(i,j)
        dxm(i)   = sqrt(area(i,j))
        dym(i)   = dxm(i)
     enddo

     var_c = sgh(is:ie,jd)
     oc_c  = oc (is:ie,jd)
     oa_c  = oa (is:ie,jd,1:ndir_asymmetry)
     ol_c  = ol (is:ie,jd,:)

     ut = 0. ; vt = 0. ; tht = 0.
     pk_e = pk_bu
!--- (1) large-scale GWD, blocking, small-scale GWD with sgh
!    (same argument list as the od_gsd call in E3SM oro_drag_interface)
     call od_gsd (u3d=u_bu, v3d=v_bu, t3d=t_bu, qv3d=q_bu,                  &
                  p3d=p_bu, p3di=pi_bu, pi3d=pk_e, z=z_bu,                  &
                  od_ls_ncleff=od_ls_ncleff, od_bl_ncd=od_bl_ncd,           &
                  od_ss_sncleff=od_ss_sncleff,                              &
                  rublten=ut, rvblten=vt, rthblten=tht,                     &
                  dtaux3d_ls=dx_ls, dtauy3d_ls=dy_ls,                       &
                  dtaux3d_bl=dx_bl, dtauy3d_bl=dy_bl,                       &
                  dtaux3d_ss=dx_ss, dtauy3d_ss=dy_ss,                       &
                  dtaux3d_fd=dx_fd, dtauy3d_fd=dy_fd,                       &
                  dusfcg_ls=su_ls, dvsfcg_ls=sv_ls,                         &
                  dusfcg_bl=su_bl, dvsfcg_bl=sv_bl,                         &
                  dusfcg_ss=su_ss, dvsfcg_ss=sv_ss,                         &
                  dusfcg_fd=su_fd, dvsfcg_fd=sv_fd,                         &
                  xland=xland, br=br,                                       &
                  var2d=var_c, oc12d=oc_c, oa2d=oa_c, ol2d=ol_c,            &
                  znu=znu, dz=dz_bu, pblh=hpbl,                             &
                  cp=CP_AIR, g=GRAV, rd=RDGAS, rv=RVGAS, ep1=fv, pi=PI,     &
                  bnvbg=nm_bu, dt=delt, dx=dxm, dy=dym,                     &
                  kpbl2d=kpbl, gwd_opt=0,                                   &
                  ids=1, ide=idim, jds=0, jde=0, kds=1, kde=nlev,           &
                  ims=1, ime=idim, jms=0, jme=0, kms=1, kme=nlev,           &
                  its=1, ite=idim, jts=0, jte=0, kts=1, kte=nlev,           &
                  gwd_ls=do_od_ls, gwd_bl=do_od_bl,                         &
                  gwd_ss=do_od_ss, gwd_fd=.false.)

!--- (2) turbulent form drag with the ~1 km std. dev. sgh30
     if (do_od_fd) then
        var30_c = sgh30(is:ie,jd)
        ut2 = 0. ; vt2 = 0. ; tht2 = 0.
        call od_gsd (u3d=u_bu, v3d=v_bu, t3d=t_bu, qv3d=q_bu,                  &
                     p3d=p_bu, p3di=pi_bu, pi3d=pk_e, z=z_bu,                  &
                     od_ls_ncleff=od_ls_ncleff, od_bl_ncd=od_bl_ncd,           &
                     od_ss_sncleff=od_ss_sncleff,                              &
                     rublten=ut2, rvblten=vt2, rthblten=tht2,                  &
                     dtaux3d_ls=dx_x1, dtauy3d_ls=dy_x1,                       &
                     dtaux3d_bl=dx_x2, dtauy3d_bl=dy_x2,                       &
                     dtaux3d_ss=dx_x3, dtauy3d_ss=dy_x3,                       &
                     dtaux3d_fd=dx_fd, dtauy3d_fd=dy_fd,                       &
                     dusfcg_ls=su_x1, dvsfcg_ls=sv_x1,                         &
                     dusfcg_bl=su_x2, dvsfcg_bl=sv_x2,                         &
                     dusfcg_ss=su_x3, dvsfcg_ss=sv_x3,                         &
                     dusfcg_fd=su_fd, dvsfcg_fd=sv_fd,                         &
                     xland=xland, br=br,                                       &
                     var2d=var30_c, oc12d=oc_c, oa2d=oa_c, ol2d=ol_c,          &
                     znu=znu, dz=dz_bu, pblh=hpbl,                             &
                     cp=CP_AIR, g=GRAV, rd=RDGAS, rv=RVGAS, ep1=fv, pi=PI,     &
                     bnvbg=nm_bu, dt=delt, dx=dxm, dy=dym,                     &
                     kpbl2d=kpbl, gwd_opt=0,                                   &
                     ids=1, ide=idim, jds=0, jde=0, kds=1, kde=nlev,           &
                     ims=1, ime=idim, jms=0, jme=0, kms=1, kme=nlev,           &
                     its=1, ite=idim, jts=0, jte=0, kts=1, kte=nlev,           &
                     gwd_ls=.false., gwd_bl=.false.,                           &
                     gwd_ss=.false., gwd_fd=.true.)
        ut = ut + ut2
        vt = vt + vt2
     endif

     if (.not. scale_by_land_frac) xland = 1.

!--- flip back to top-down, apply land-fraction weighting
     do k = 1, nlev
        kr = nlev + 1 - k
        utnd  (:,j,kr) = ut   (:,k)*xland
        vtnd  (:,j,kr) = vt   (:,k)*xland
        c_ls_u(:,j,kr) = dx_ls(:,k)*xland ; c_ls_v(:,j,kr) = dy_ls(:,k)*xland
        c_bl_u(:,j,kr) = dx_bl(:,k)*xland ; c_bl_v(:,j,kr) = dy_bl(:,k)*xland
        c_ss_u(:,j,kr) = dx_ss(:,k)*xland ; c_ss_v(:,j,kr) = dy_ss(:,k)*xland
        c_fd_u(:,j,kr) = dx_fd(:,k)*xland ; c_fd_v(:,j,kr) = dy_fd(:,k)*xland
     enddo
     s_ls_u(:,j) = su_ls*xland ; s_ls_v(:,j) = sv_ls*xland
     s_bl_u(:,j) = su_bl*xland ; s_bl_v(:,j) = sv_bl*xland
     s_ss_u(:,j) = su_ss*xland ; s_ss_v(:,j) = sv_ss*xland
     s_fd_u(:,j) = su_fd*xland ; s_fd_v(:,j) = sv_fd*xland
  enddo

!--- kinetic energy removed by the drag is returned as heat (as in E3SM)
  if (do_dissipative_heating) then
     ttnd = -( utnd*(u + 0.5*delt*utnd) + vtnd*(v + 0.5*delt*vtnd) )/CP_AIR
  else
     ttnd = 0.
  endif

!--- diagnostics
  if (id_udt_oro      > 0) used = send_data (id_udt_oro,      utnd,   Time, is, js, 1)
  if (id_vdt_oro      > 0) used = send_data (id_vdt_oro,      vtnd,   Time, is, js, 1)
  if (id_tdt_diss_oro > 0) used = send_data (id_tdt_diss_oro, ttnd,   Time, is, js, 1)
  if (id_udt_ls > 0) used = send_data (id_udt_ls, c_ls_u, Time, is, js, 1)
  if (id_vdt_ls > 0) used = send_data (id_vdt_ls, c_ls_v, Time, is, js, 1)
  if (id_udt_bl > 0) used = send_data (id_udt_bl, c_bl_u, Time, is, js, 1)
  if (id_vdt_bl > 0) used = send_data (id_vdt_bl, c_bl_v, Time, is, js, 1)
  if (id_udt_ss > 0) used = send_data (id_udt_ss, c_ss_u, Time, is, js, 1)
  if (id_vdt_ss > 0) used = send_data (id_vdt_ss, c_ss_v, Time, is, js, 1)
  if (id_udt_fd > 0) used = send_data (id_udt_fd, c_fd_u, Time, is, js, 1)
  if (id_vdt_fd > 0) used = send_data (id_vdt_fd, c_fd_v, Time, is, js, 1)
  if (id_taux_ls > 0) used = send_data (id_taux_ls, s_ls_u, Time, is, js)
  if (id_tauy_ls > 0) used = send_data (id_tauy_ls, s_ls_v, Time, is, js)
  if (id_taux_bl > 0) used = send_data (id_taux_bl, s_bl_u, Time, is, js)
  if (id_tauy_bl > 0) used = send_data (id_tauy_bl, s_bl_v, Time, is, js)
  if (id_taux_ss > 0) used = send_data (id_taux_ss, s_ss_u, Time, is, js)
  if (id_tauy_ss > 0) used = send_data (id_tauy_ss, s_ss_v, Time, is, js)
  if (id_taux_fd > 0) used = send_data (id_taux_fd, s_fd_u, Time, is, js)
  if (id_tauy_fd > 0) used = send_data (id_tauy_fd, s_fd_v, Time, is, js)

end subroutine oro_drag

!#######################################################################

subroutine oro_drag_end
  if (.not. module_is_initialized) return
  if (allocated(sgh)) deallocate (sgh, sgh30, oc, oa, ol, znu)
  module_is_initialized = .false.
end subroutine oro_drag_end

!#######################################################################
!< oro_register_restart: registers the oro_drag.res.nc restart variables

subroutine oro_register_restart (Od_restart, with_sgh, with_sgh30)
type(FmsNetcdfDomainFile_t), intent(inout) :: Od_restart !< Fms2io domain decomposed fileobj
logical,                     intent(in)    :: with_sgh   !< also register SGH
logical,                     intent(in)    :: with_sgh30 !< also register SGH30
character(len=8), dimension(3) :: dim_names2 !< String array of dimension names
character(len=8), dimension(4) :: dim_names3 !< String array of dimension names (OA)
character(len=8), dimension(4) :: dim_names4 !< String array of dimension names (OL)

  ! Fortran order; ncdump shows (Time, [nOA|nOL,] yaxis_1, xaxis_1)
  dim_names2 = (/"xaxis_1", "yaxis_1", "Time   "/)
  dim_names3 = (/"xaxis_1", "yaxis_1", "nOA    ", "Time   "/)
  dim_names4 = (/"xaxis_1", "yaxis_1", "nOL    ", "Time   "/)

  call register_axis (Od_restart, "xaxis_1", "x")
  call register_axis (Od_restart, "yaxis_1", "y")
  call register_axis (Od_restart, "nOA", size(oa,3))
  call register_axis (Od_restart, "nOL", size(ol,3))
  call register_axis (Od_restart, "Time", unlimited)

  !< Register the domain decomposed dimensions as variables so that the combiner can work
  !! correctly
  call register_field (Od_restart, "xaxis_1", "double", (/"xaxis_1"/))
  call register_field (Od_restart, "yaxis_1", "double", (/"yaxis_1"/))

  call register_restart_field (Od_restart, 'OC', oc, dim_names2)
  call register_restart_field (Od_restart, 'OA', oa, dim_names3)
  call register_restart_field (Od_restart, 'OL', ol, dim_names4)
  if (with_sgh)   call register_restart_field (Od_restart, 'SGH',   sgh,   dim_names2)
  if (with_sgh30) call register_restart_field (Od_restart, 'SGH30', sgh30, dim_names2)

end subroutine oro_register_restart

!#######################################################################
!< mg_register_restart: registers ghprime, identical to mg_drag_mod

subroutine mg_register_restart (Mg_restart)
type(FmsNetcdfDomainFile_t), intent(inout) :: Mg_restart !< Fms2io domain decomposed fileobj
character(len=8), dimension(3) :: dim_names !< String array of dimension names

  dim_names(1) = "xaxis_1"
  dim_names(2) = "yaxis_1"
  dim_names(3) = "Time"

  call register_axis (Mg_restart, dim_names(1), "x")
  call register_axis (Mg_restart, dim_names(2), "y")
  call register_axis (Mg_restart, dim_names(3), unlimited)

  !< Register the domain decomposed dimensions as variables so that the combiner can work
  !! correctly
  call register_field (Mg_restart, dim_names(1), "double", (/dim_names(1)/))
  call register_field (Mg_restart, dim_names(2), "double", (/dim_names(2)/))

  call register_restart_field (Mg_restart, 'ghprime', sgh, dim_names)

end subroutine mg_register_restart

!#######################################################################
! <SUBROUTINE NAME="oro_drag_restart">
!
! write out restart file(s), same as mg_drag_restart.
! Arguments:
!   timestamp (optional, intent(in)) : A character string that represents the model time,
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix.

subroutine oro_drag_restart (timestamp)
character(len=*), intent(in), optional :: timestamp

type(FmsNetcdfDomainFile_t) :: Od_restart, Mg_restart !< Fms2io domain decomposed fileobjs
character(len=128) :: filename !< Restart filename

  if (.not. module_is_initialized) return

!--- oro_drag.res.nc
  if (present(timestamp)) then
    filename = "RESTART/"//trim(timestamp)//".oro_drag.res.nc"
  else
    filename = "RESTART/oro_drag.res.nc"
  endif
  if (mpp_pe() == mpp_root_pe()) call error_mesg ('oro_drag_mod', &
       'Writing NetCDF formatted restart file: '//trim(filename), NOTE)
  if (open_file(Od_restart, trim(filename), "overwrite", od_domain, is_restart=.true.)) then
    call oro_register_restart (Od_restart, with_sgh=(trim(sgh_restart) == 'oro_drag'), &
                               with_sgh30=have_sgh30)
    call write_restart (Od_restart)
    call add_domain_dimension_data (Od_restart)
    call close_file (Od_restart)
  endif

!--- mg_drag.res.nc, only if mg_drag_mod does not write it itself
  if (trim(sgh_restart) == 'mg_drag' .and. .not. mg_on) then
    if (present(timestamp)) then
      filename = "RESTART/"//trim(timestamp)//".mg_drag.res.nc"
    else
      filename = "RESTART/mg_drag.res.nc"
    endif
    if (mpp_pe() == mpp_root_pe()) call error_mesg ('oro_drag_mod', &
         'Writing NetCDF formatted restart file: '//trim(filename), NOTE)
    if (open_file(Mg_restart, trim(filename), "overwrite", od_domain, is_restart=.true.)) then
      call mg_register_restart (Mg_restart)
      call write_restart (Mg_restart)
      call add_domain_dimension_data (Mg_restart)
      call close_file (Mg_restart)
    endif
  endif

end subroutine oro_drag_restart
! </SUBROUTINE> NAME="oro_drag_restart"

!#######################################################################
!< Add_dimension_data: Adds dummy data for the domain decomposed axis

subroutine add_domain_dimension_data (fileobj)
type(FmsNetcdfDomainFile_t) :: fileobj !< Fms2io domain decomposed fileobj
integer, dimension(:), allocatable :: buffer !< Buffer with axis data
integer :: is, ie !< Starting and Ending indices for data

  call get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
  call write_data(fileobj, "xaxis_1", buffer)
  deallocate(buffer)

  call get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
  call write_data(fileobj, "yaxis_1", buffer)
  deallocate(buffer)

end subroutine add_domain_dimension_data

!#######################################################################
! Mid-level Brunt-Vaisala frequency, top-down column(s).  Same algorithm as
! E3SM gw_common::gw_prof (dry, n2min = 1e-8): interface values from
! dT/dp, top interface uses the isothermal value, bottom copies the one above,
! mid-levels are the average of the bounding interfaces.

subroutine bv_freq (t, pfull, phalf, nm)
real, intent(in)  :: t(:,:), pfull(:,:), phalf(:,:)
real, intent(out) :: nm(:,:)
real, parameter :: n2min = 1.e-8
real :: ni(size(t,1), 0:size(t,2)), ti, rhoi, dtdp, n2
integer :: i, k, nlev

  nlev = size(t,2)
  do i = 1, size(t,1)
     ni(i,0) = sqrt(GRAV*GRAV/(CP_AIR*t(i,1)))
     do k = 1, nlev-1                          ! interface between layers k, k+1
        ti   = 0.5*(t(i,k) + t(i,k+1))
        rhoi = phalf(i,k+1)/(RDGAS*ti)
        dtdp = (t(i,k+1) - t(i,k))/(pfull(i,k+1) - pfull(i,k))
        n2   = GRAV*GRAV/ti*(1./CP_AIR - rhoi*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     enddo
     ni(i,nlev) = ni(i,nlev-1)
     do k = 1, nlev
        nm(i,k) = 0.5*(ni(i,k-1) + ni(i,k))
     enddo
  enddo
end subroutine bv_freq

!#######################################################################

subroutine clean_field2 (x)
real, intent(inout) :: x(:,:)
  where (x /= x .or. abs(x) > 1.e20) x = 0.
end subroutine clean_field2

subroutine clean_field3 (x)
real, intent(inout) :: x(:,:,:)
  where (x /= x .or. abs(x) > 1.e20) x = 0.
end subroutine clean_field3

!#######################################################################
!=======================================================================
!  E3SM od_common.F90 routines, original form and names:
!  dxygrid, od_gsd, od2d
!=======================================================================

subroutine dxygrid(dx,dy,theta_in,dxy)

  IMPLICIT NONE
  real(r8),intent(in) :: dx,dy,theta_in
  real(r8),intent(out):: dxy
  !local variables
  real(r8) :: rad,theta,theta1

  rad=4.0_r8*atan(1.0_r8)/180.0_r8
  theta1=MOD(theta_in,360._r8)
  !set negative axis into 0~360
  if (theta1.ge.-360._r8.and.theta1.lt.0._r8) then
    theta1=theta1+360._r8
  endif
  !in case the angle is not into the judgement
  theta=theta1
  !transform of angle into first quadrant
  if      (theta1.ge.  0._r8.and.theta1.lt. 90._r8) then
    theta=theta1
  else if (theta1.gt. 90._r8.and.theta1.lt.180._r8) then
    theta=(180._r8-theta1)
  else if (theta1.gt.180._r8.and.theta1.lt.270._r8) then
    theta=(theta1-180._r8)
  else if (theta1.gt.270._r8.and.theta1.lt.360._r8) then
    theta=(360._r8-theta1)
  else if (theta1.eq.90._r8.or.theta1.eq.270._r8) then
    theta=90._r8
  else if (theta1.eq.0._r8.or.theta1.eq.180._r8) then
    theta=0._r8
  endif

  !get dxy
  if   (theta.ge. 0._r8.and.theta.lt.atan2(dy,dx)/rad) then
    dxy=dx/cos(theta*rad)
  else if (theta.ge.atan2(dy,dx)/rad.and.theta.le.90._r8)then
    dxy=dy/sin(theta*rad)
  endif

end subroutine dxygrid

!=======================================================================

subroutine od_gsd(u3d,v3d,t3d,qv3d,p3d,p3di,pi3d,z,                            &
                  od_ls_ncleff,od_bl_ncd,od_ss_sncleff,                        &
                  rublten,rvblten,rthblten,                                    &
                  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,                 &
                  dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,                 &
                  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,dusfcg_ss,dvsfcg_ss, &
                  dusfcg_fd,dvsfcg_fd,xland,br,                                &
                  var2d,oc12d,oa2d,ol2d,znu,znw,p_top,dz,pblh,                 &
                  cp,g,rd,rv,ep1,pi,bnvbg,                                     &
                  dt,dx,dy,kpbl2d,gwd_opt,                                     &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                  gwd_ls,gwd_bl,gwd_ss,gwd_fd)
  !-------------------------------------------------------------------------------
   implicit none
  !-------------------------------------------------------------------------------
  !                                                                       
  !-- u3d         3d u-velocity interpolated to theta points (m/s)
  !-- v3d         3d v-velocity interpolated to theta points (m/s)
  !-- t3d         temperature (k)
  !-- qv3d        3d water vapor mixing ratio (kg/kg)
  !-- p3d         3d pressure (pa)
  !-- p3di        3d pressure (pa) at interface level
  !-- pi3d        3d exner function (dimensionless)
  !-- rublten     u tendency due to pbl parameterization (m/s/s) 
  !-- rvblten     v tendency due to pbl parameterization (m/s/s)
  !-- rthblten    theta tendency due to pbl parameterization (K/s)
  !-- znu         eta values (sigma values)
  !-- cp          heat capacity at constant pressure for dry air (j/kg/k)
  !-- g           acceleration due to gravity (m/s^2)
  !-- rd          gas constant for dry air (j/kg/k)
  !-- z           height above sea level (m)
  !-- rv          gas constant for water vapor (j/kg/k)
  !-- dt          time step (s)
  !-- dx          model grid interval (m)
  !-- dz          height of model layers (m)
  !-- xland       land mask (1 for land, 2 for water)
  !-- br          bulk richardson number in surface layer
  !-- pblh        planetary boundary layer height (m)
  !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
  !-- ids         start index for i in domain
  !-- ide         end index for i in domain
  !-- jds         start index for j in domain
  !-- jde         end index for j in domain
  !-- kds         start index for k in domain
  !-- kde         end index for k in domain
  !-- ims         start index for i in memory
  !-- ime         end index for i in memory
  !-- jms         start index for j in memory
  !-- jme         end index for j in memory
  !-- kms         start index for k in memory
  !-- kme         end index for k in memory
  !-- its         start index for i in tile
  !-- ite         end index for i in tile
  !-- jts         start index for j in tile
  !-- jte         end index for j in tile
  !-- kts         start index for k in tile
  !-- kte         end index for k in tile
  !-------------------------------------------------------------------------------
  integer,  intent(in)   ::      ids, ide, jds, jde, kds, kde
  integer,  intent(in)   ::      ims, ime, jms, jme, kms, kme
  integer,  intent(in)   ::      its, ite, jts, jte, kts, kte
  integer,  intent(in)   ::      gwd_opt
  real(r8), intent(in)   ::      cp,g,rd,rv,ep1,pi
  !input model grid length for the grid dx,dy
  real(r8), dimension(:),                   intent(in)  :: dx
  real(r8), dimension(:),                   intent(in)  :: dy
  !input atmospheric variables
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  :: qv3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  p3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  :: pi3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  t3d
  real(r8), dimension( ims:ime, kms:kme+1), intent(in)  :: p3di
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  u3d
  real(r8), dimension( ims:ime, kms:kme ),  intent(in)  ::  v3d
  !input tunable parameters
  real(r8),                                 intent(in)  :: od_ls_ncleff
  real(r8),                                 intent(in)  :: od_bl_ncd
  real(r8),                                 intent(in)  :: od_ss_sncleff
  !logical variables for selection of 4 schemes
  logical,                                  intent(in)  :: gwd_ls
  logical,                                  intent(in)  :: gwd_bl
  logical,                                  intent(in)  :: gwd_ss
  logical,                                  intent(in)  :: gwd_fd
  !input variables
  integer,  dimension( ims:ime ),           intent(in)  :: kpbl2d
  real(r8), dimension( ims:ime ),           intent(in)  :: pblh, br, xland
  real(r8), dimension( ims:ime,kms:kme ),   intent(in)  ::          z
  real(r8), dimension( ims:ime,kms:kme ),   intent(in)  ::         dz
  !input topographic parameters
  real(r8), dimension( ims:ime ),           intent(in),   optional  :: var2d
  real(r8), dimension( ims:ime ),           intent(in),   optional  :: oc12d
  real(r8), dimension( ims:ime,ndir_efflength ),intent(in),   optional  :: ol2d
  real(r8), dimension( ims:ime,ndir_asymmetry ),intent(in),   optional  :: oa2d
  !input model parameters
  real(r8),                                 intent(in),   optional  :: dt
  real(r8),                                 intent(in),   optional  :: p_top
  real(r8), dimension(:),                   intent(in),   optional  :: znu
  real(r8), dimension(:),                   intent(in),   optional  :: znw
  !input bnv from gw_drag
  real(r8), dimension( ims:ime,kms:kme ),   intent(in),   optional  ::  bnvbg
  !output of vertical profile for momentum terms
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rublten
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rvblten
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   rthblten
  !output for 4 drag terms
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_ls
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_ls
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_bl
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_bl
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_ss
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_ss
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtaux3d_fd
  real(r8), dimension( ims:ime,kms:kme ),  intent(inout), optional  ::   dtauy3d_fd
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_ls
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_ls
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_bl
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_bl    
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_ss
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_ss
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dusfcg_fd
  real(r8), dimension( ims:ime ),          intent(inout), optional  ::   dvsfcg_fd
  !local drag terms
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_ls
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_ls
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_bl
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_bl
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_ss
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_ss
  real(r8), dimension( ims:ime, kms:kme )  ::  dtaux2d_fd
  real(r8), dimension( ims:ime, kms:kme )  ::  dtauy2d_fd
  real(r8), dimension( ims:ime ) ::  dusfc_ls
  real(r8), dimension( ims:ime ) ::  dvsfc_ls
  real(r8), dimension( ims:ime ) ::  dusfc_bl
  real(r8), dimension( ims:ime ) ::  dvsfc_bl
  real(r8), dimension( ims:ime ) ::  dusfc_ss
  real(r8), dimension( ims:ime ) ::  dvsfc_ss
  real(r8), dimension( ims:ime ) ::  dusfc_fd
  real(r8), dimension( ims:ime ) ::  dvsfc_fd
  !local variables
  real(r8),   dimension( its:ite, kts:kte )     ::  delprsi
  real(r8),   dimension( its:ite, kts:kte )     ::  pdh
  real(r8),   dimension( its:ite, kts:kte+1 )   ::  pdhi
  real(r8),   dimension( its:ite, ndir_asymmetry )  ::  oa4
  real(r8),   dimension( its:ite, ndir_efflength )  ::  ol4
  integer ::  i,j,k,kpblmax
  !determine the lowest level for planet boundary layer
  do k = kts,kte
    if( znu(k).gt.0.6_r8 ) kpblmax = k + 1
  enddo
  !get the interface pressure (hPa) for each level
  !for the last level over kte use mid-layer pressure instead
  do k = kts,kte+1
    do i = its,ite
      if(k.le.kte) pdh(i,k) = p3d(i,k)
        pdhi(i,k) = p3di(i,k)
    enddo
  enddo
  !get the layer pressure for each level
  do k = kts,kte
    do i = its,ite
      delprsi(i,k) = pdhi(i,k)-pdhi(i,k+1)
    enddo
  enddo
  !no need when there is no large drag
  if ( gwd_ls .or. gwd_bl ) then
    do i = its,ite
      oa4(i,:) = oa2d(i,:)
      ol4(i,:) = ol2d(i,:)
    enddo
  endif
  !call the od2d for calculatino of each grid
  call od2d(dudt=rublten(ims,kms),dvdt=rvblten(ims,kms)                   &
             ,dthdt=rthblten(ims,kms)                                       &
             ,ncleff=od_ls_ncleff,ncd=od_bl_ncd,sncleff=od_ss_sncleff                &
             ,dtaux2d_ls=dtaux2d_ls,dtauy2d_ls=dtauy2d_ls                   &
             ,dtaux2d_bl=dtaux2d_bl,dtauy2d_bl=dtauy2d_bl                   &
             ,dtaux2d_ss=dtaux2d_ss,dtauy2d_ss=dtauy2d_ss                   &
             ,dtaux2d_fd=dtaux2d_fd,dtauy2d_fd=dtauy2d_fd                   &
             ,u1=u3d(ims,kms),v1=v3d(ims,kms)                               &
             ,t1=t3d(ims,kms)                                               &
             ,q1=qv3d(ims,kms)                                              &
             ,del=delprsi(its,kts)                                          &
             ,prsi=pdhi(its,kts)                                            &
             ,prsl=pdh(its,kts),prslk=pi3d(ims,kms)                         &
             ,zl=z(ims,kms),rcl=1.0_r8                                      &
             ,xland1=xland(ims),br1=br(ims),hpbl=pblh(ims)                  &
             ,bnv_in=bnvbg(ims,kms)                                         &
             ,dz2=dz(ims,kms)                                               &
             ,kpblmax=kpblmax                                               &
             ,dusfc_ls=dusfc_ls,dvsfc_ls=dvsfc_ls                           &
             ,dusfc_bl=dusfc_bl,dvsfc_bl=dvsfc_bl                           &
             ,dusfc_ss=dusfc_ss,dvsfc_ss=dvsfc_ss                           &
             ,dusfc_fd=dusfc_fd,dvsfc_fd=dvsfc_fd                           &
             ,var=var2d(ims),oc1=oc12d(ims)                                 &
             ,oa4=oa4,ol4=ol4                                               &
             ,g=g,cp=cp,rd=rd,rv=rv,fv=ep1,pi=pi                            &
             ,dxmeter=dx,dymeter=dy,deltim=dt                               &
             ,kpbl=kpbl2d(ims)                                              &
             ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde             &
             ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme             &
             ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte             &
             ,gsd_gwd_ls=gwd_ls,gsd_gwd_bl=gwd_bl,gsd_gwd_ss=gwd_ss,gsd_gwd_fd=gwd_fd)
    !set the total stress output to each terms for the 4 drag schemes
    do i = its,ite
      dusfcg_ls(i) = dusfc_ls(i)
      dvsfcg_ls(i) = dvsfc_ls(i)
      dusfcg_bl(i) = dusfc_bl(i)
      dvsfcg_bl(i) = dvsfc_bl(i)
      dusfcg_ss(i) = dusfc_ss(i)
      dvsfcg_ss(i) = dvsfc_ss(i)
      dusfcg_fd(i) = dusfc_fd(i)
      dvsfcg_fd(i) = dvsfc_fd(i)
    enddo
    !set the 3D output tendencies to each terms for the 4 drag schemes
    dtaux3d_ls = dtaux2d_ls
    dtaux3d_bl = dtaux2d_bl
    dtauy3d_ls = dtauy2d_ls
    dtauy3d_bl = dtauy2d_bl
    dtaux3d_ss = dtaux2d_ss
    dtaux3d_fd = dtaux2d_fd
    dtauy3d_ss = dtauy2d_ss
    dtauy3d_fd = dtauy2d_fd

end subroutine od_gsd

!=======================================================================

subroutine od2d(dudt,dvdt,dthdt,ncleff,ncd,sncleff,                        &
                dtaux2d_ls,dtauy2d_ls,                                     &
                dtaux2d_bl,dtauy2d_bl,                                     &
                dtaux2d_ss,dtauy2d_ss,                                     &
                dtaux2d_fd,dtauy2d_fd,                                     &
                u1,v1,t1,q1,                                               &
                del,                                                       &
                prsi,prsl,prslk,zl,rcl,                                    &
                xland1,br1,hpbl,bnv_in,dz2,                                &
                kpblmax,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,               &
                dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,var,oc1,oa4,ol4,       &
                g,cp,rd,rv,fv,pi,dxmeter,dymeter,deltim,kpbl,              &
                ids,ide, jds,jde, kds,kde,                                 &
                ims,ime, jms,jme, kms,kme,                                 &
                its,ite, jts,jte, kts,kte,                                 &
                gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd)
  !  This code handles the time tendencies of u v due to the effect of mountain 
  !  induced gravity wave drag from sub-grid scale orography. It includes 4 parts:
  !  orographic gravity wave drag and flow-blocking drag (Xie et al.,2020),small-scale  
  !  orographic gravity wave drag (Tsiringakis et al. 2017), and turbulent orographic
  !  form drag (Beljaars et al.,2004).
  !
  !           Activation of each component is done by specifying the integer-parameters
  !           (defined below) to .true. (active) or .false. (inactive)
  !                    gsd_gwd_ls : large-scale
  !                    gsd_gwd_bl : blocking drag 
  !                    gsd_gwd_ss : small-scale gravity wave drag
  !                    gsd_gwd_fd : topographic form drag
  !
  !
  !        References:
  !        Xie et al. (2020), JAMES
  !        Tsiringakis et al. (2017), QJRMS
  !        Beljaars et al., (2004), QJRMS
  !-------------------------------------------------------------------------------
  !
  !  input                                                                
  !        dudt (ims:ime,kms:kme)  non-lin tendency for u wind component
  !        dvdt (ims:ime,kms:kme)  non-lin tendency for v wind component
  !        u1(ims:ime,kms:kme) zonal wind / sqrt(rcl)  m/sec  at t0-dt
  !        
  !        v1(ims:ime,kms:kme) meridional wind / sqrt(rcl) m/sec at t0-dt
  !        t1(ims:ime,kms:kme) temperature deg k at t0-dt
  !        q1(ims:ime,kms:kme) specific humidity at t0-dt
  !
  !        rcl     a scaling factor = reciprocal of square of cos(lat)
  !                for gmp.  rcl=1 if u1 and v1 are wind components.
  !        deltim  time step    secs                                       
  !        del(kts:kte)  positive increment of pressure across layer (pa)
  !                                                                       
  !  output
  !        dudt, dvdt    wind tendency due to gwdo
  !
  !-------------------------------------------------------------------------------
  implicit none
  !-------------------------------------------------------------------------------
  !intput parameters
  integer,  intent(in)         ::  ids,ide, jds,jde, kds,kde
  integer,  intent(in)         ::  ims,ime, jms,jme, kms,kme
  integer,  intent(in)         ::  its,ite, jts,jte, kts,kte
  !tunable parameter in oro_drag_nl, od_ls_ncleff,od_bl_ncd,od_ss_sncleff
  real(r8), intent(in) :: ncleff
  real(r8), intent(in) :: ncd
  real(r8), intent(in) :: sncleff
  !model timestep and other parameters
  real(r8), intent(in) ::  g,rd,rv,fv,cp,pi,deltim,rcl
  !input model grid length
  real(r8), dimension(:), intent(in)   ::  dxmeter
  real(r8), dimension(:), intent(in)   ::  dymeter
  !input topo variables 
  real(r8), dimension( ims:ime,ndir_asymmetry ), intent(in) ::  oa4
  real(r8), dimension( ims:ime,ndir_efflength ), intent(in) ::  ol4
  real(r8), dimension( ims:ime )           , intent(in) ::  var
  real(r8), dimension( ims:ime )           , intent(in) ::  oc1
  !input atmospheric variables
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  u1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  v1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  t1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  q1
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  zl
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  prsl
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  prslk
  real(r8), dimension( ims:ime,kms:kme+1),   intent(in) ::  prsi
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  del
  real(r8), dimension( ims:ime,kms:kme ),    intent(in) ::  bnv_in
  !input pbl level index
  integer,dimension( ims:ime ),              intent(in) ::  kpbl
  integer,                                   intent(in) ::  kpblmax
  !input for small-scale ogwd
  real(r8), dimension( ims:ime ),            intent(in) :: br1
  real(r8), dimension( ims:ime ),            intent(in) :: hpbl
  real(r8), dimension( ims:ime ),            intent(in) :: xland1
  real(r8), dimension( ims:ime, kms:kme ),   intent(in) :: dz2
  !variables for open/close process
  logical,                                   intent(in) ::  gsd_gwd_ls
  logical,                                   intent(in) ::  gsd_gwd_bl
  logical,                                   intent(in) ::  gsd_gwd_ss
  logical,                                   intent(in) ::  gsd_gwd_fd
  !input/output tendencies
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dudt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dvdt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dthdt
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_ls
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_ls
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_bl
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_bl
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_ss
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_ss
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtaux2d_fd
  real(r8), dimension( ims:ime,kms:kme ),    intent(inout) ::  dtauy2d_fd
  !input/output total column stress
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_ls
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_ls              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_bl
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_bl              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_ss
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_ss              
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dusfc_fd
  real(r8), dimension( ims:ime )        ,    intent(inout) ::  dvsfc_fd
  !
  ! added for small-scale orographic wave drag
  !
  integer                                  :: kpbl2,kvar
  real(r8), dimension(its:ite,kts:kte)     :: utendwave,vtendwave,thx,thvx,za
  real(r8), dimension(its:ite)             :: govrth
  real(r8), dimension(its:ite,kts:kte+1)   :: zq
  real(r8)                                 :: tauwavex0,tauwavey0,bnrf,density,tvcon,hpbl2
  real(r8), parameter                      :: varmax = 200._r8
  ! Variables for scale-awareness:
  ! Small-scale GWD + turbulent form drag
  real(r8), parameter   :: dxmin_ss = 1000._r8, dxmax_ss = 12000._r8  ! min,max range of tapering (m)
  ! Large-scale GWD
  real(r8), parameter   :: dxmin_ls = 3000._r8, dxmax_ls = 13000._r8  ! min,max range of tapering (m)
  !Add y axis for taper consider
  real(r8), parameter   :: dymin_ls = 3000._r8, dymax_ls = 13000._r8  ! min,maxrange of tapering (m)
  real(r8), parameter   :: dymin_ss = 3000._r8, dymax_ss = 13000._r8  ! min,maxrange of tapering (m)
  real(r8)              :: ss_taper, ls_taper  ! small- and large-scale tapering factors (-)
  !
  ! added Beljaars orographic form drag
  real(r8), dimension( its:ite,kts:kte )   :: utendform
  real(r8), dimension( its:ite,kts:kte )   :: vtendform
  real(r8)                                 :: a1,a2,wsp
  ! critical richardson number for wave breaking : ! larger drag with larger value
  real(r8),parameter       ::  ric     = 1._r8!original 0.25, but 1. seems better at drag profile
  real(r8),parameter       ::  ric_rig  = 0.25_r8
  real(r8),parameter       ::  dw2min  = 1._r8
  real(r8),parameter       ::  rimin   = -100._r8
  real(r8),parameter       ::  bnv2min = 1.0e-5_r8
  real(r8),parameter       ::  efmin   = 0.0_r8
  real(r8),parameter       ::  efmax   = 10.0_r8
  real(r8),parameter       ::  xl      = 4.0e4_r8
  real(r8),parameter       ::  critac  = 1.0e-5_r8
  real(r8),parameter       ::  gmax    = 1._r8
  real(r8),parameter       ::  veleps  = 1.0_r8
  real(r8),parameter       ::  factop  = 0.5_r8
  real(r8),parameter       ::  frc     = 1.0_r8
  real(r8),parameter       ::  ce      = 0.8_r8
  real(r8),parameter       ::  cg      = 0.5_r8
  real(r8),parameter       ::  tndmax  = 400._r8 / 86400._r8 ! convert 400 m/s/day to m/s/s
  integer,parameter        ::  kpblmin = 2
  !number of direction for ogwd
  integer,parameter        ::  mdir=2*ndir_efflength
  !  variables for flow-blocking drag
  real(r8),parameter       ::  frmax  = 10._r8
  real(r8),parameter       ::  olmin  = 1.0e-5_r8
  real(r8),parameter       ::  odmin  = 0.1_r8
  real(r8),parameter       ::  odmax  = 10._r8
  real(r8),parameter       ::  erad   = 6371.315e+3_r8

  !local variables
  integer                  ::  i,j,k,lcap,lcapp1,nwd,idir
  integer                  ::  klcap,kp1,ikount,kk,nwd1!added nwd1
  real(r8)                 ::  rcs,rclcs,csg,fdir,cleff,cs,rcsks
  real(r8)                 ::  wdir,ti,rdz,temp,tem2,dw2,shr2,bvf2,rdelks
  real(r8)                 ::  wtkbj,tem,gfobnv,hd,fro,rim,temc,tem1,efact
  real(r8)                 ::  temv,dtaux,dtauy,eng0,eng1,theta,rad,wdir1
  !logical variables
  logical,dimension( its:ite )            ::  ldrag
  logical,dimension( its:ite )            ::  icrilv
  logical,dimension( its:ite )            ::  flag
  logical,dimension( its:ite )            ::  kloop1
  !parameters for taking over the input
  real(r8),dimension( its:ite )           :: taub
  real(r8),dimension( its:ite )           :: xn
  real(r8),dimension( its:ite )           :: yn
  real(r8),dimension( its:ite )           :: ubar
  real(r8),dimension( its:ite )           :: vbar
  real(r8),dimension( its:ite )           :: fr
  real(r8),dimension( its:ite )           :: ulow
  real(r8),dimension( its:ite )           :: rulow
  real(r8),dimension( its:ite )           :: bnv
  real(r8),dimension( its:ite )           :: oa1
  real(r8),dimension( its:ite )           :: ol
  real(r8),dimension( its:ite )           :: roll
  real(r8),dimension( its:ite )           :: dtfac
  real(r8),dimension( its:ite )           :: brvf
  real(r8),dimension( its:ite )           :: xlinv
  real(r8),dimension( its:ite )           :: delks
  real(r8),dimension( its:ite )           :: delks1
  real(r8),dimension( its:ite,kts:kte )   :: bnv2
  real(r8),dimension( its:ite,kts:kte )   :: usqj
  real(r8),dimension( its:ite,kts:kte )   :: taud_ls
  real(r8),dimension( its:ite,kts:kte )   :: taud_bl
  real(r8),dimension( its:ite,kts:kte )   :: ro
  real(r8),dimension( its:ite,kts:kte )   :: vtk
  real(r8),dimension( its:ite,kts:kte )   :: vtj
  real(r8),dimension( its:ite )           :: zlowtop
  real(r8),dimension( its:ite )           :: coefm
  real(r8),dimension( its:ite,kts:kte+1 ) :: taup
  real(r8),dimension( its:ite,kts:kte-1 ) :: velco
  !pbl related variables
  integer ,dimension( its:ite )           :: kbl
  integer ,dimension( its:ite )           :: klowtop
  integer ,dimension( its:ite )           :: komax
  !flow-blocking related variables
  integer                                 :: kblk
  real(r8)                                :: cd
  real(r8)                                :: zblk,tautem
  real(r8)                                :: pe,ke
  !grid related variables
  real(r8),dimension( its:ite )           :: delx
  real(r8),dimension( its:ite )           :: dely
  real(r8),dimension( its:ite )           :: dxy
  real(r8),dimension( its:ite )           :: dxyp
  real(r8),dimension( its:ite,ndir_efflength ):: dxy4
  real(r8),dimension( its:ite,ndir_efflength ):: dxy4p
  !topo parameters
  real(r8),dimension( its:ite )           :: olp
  real(r8),dimension( its:ite )           :: od
  real(r8),dimension( its:ite,kts:kte+1 ) :: taufb
  !readdata for low-level determination of ogwd
  real(r8), dimension( its:ite )          :: zl_hint
  real(r8)                                :: l1,l2,S
  logical                                 :: iint
  !open/close low-level momentum adjustment according to scorer parameter
  logical  :: scorer_on=.false.
  !
  !---- constants                                                         
  !                                                                       
  rcs    = sqrt(rcl)
  cs     = 1._r8 / sqrt(rcl)
  csg    = cs * g              
  lcap   = kte
  lcapp1 = lcap + 1
  fdir   = mdir / (2.0_r8*pi)
  !
  !--- calculate scale-aware tapering factors, currently set as no taper
  !
  ls_taper=1._r8
  ss_taper=1._r8
  !
  !--- calculate length of grid for flow-blocking drag
  !
  delx(its:ite)=dxmeter(its:ite)
  dely(its:ite)=dymeter(its:ite)
  !
  !
  !-----initialize arrays
  dtaux = 0.0_r8
  dtauy = 0.0_r8
  do i = its,ite
    klowtop(i)    = 0
    kbl(i)        = 0
  enddo

  do i = its,ite
    xn(i)         = 0.0_r8
    yn(i)         = 0.0_r8
    ubar (i)      = 0.0_r8
    vbar (i)      = 0.0_r8
    roll (i)      = 0.0_r8
    taub (i)      = 0.0_r8
    oa1(i)        = 0.0_r8
    ol(i)         = 0.0_r8
    fr(i)         = 0.0_r8
    ulow (i)      = 0.0_r8
    dtfac(i)      = 1.0_r8
    ldrag(i)      = .false.
    icrilv(i)     = .false.
    flag(i)       = .true.
    zl_hint(i)    = 0.0_r8
  enddo
   
  do k = kts,kte
    do i = its,ite
      usqj(i,k) = 0.0_r8
      bnv2(i,k) = 0.0_r8
      vtj(i,k)  = 0.0_r8
      vtk(i,k)  = 0.0_r8
      taup(i,k) = 0.0_r8
      taud_ls(i,k) = 0.0_r8
      taud_bl(i,k) = 0.0_r8
      dtaux2d_ls(i,k)= 0.0_r8
      dtauy2d_ls(i,k)= 0.0_r8
      dtaux2d_bl(i,k)= 0.0_r8
      dtauy2d_bl(i,k)= 0.0_r8
      dtaux2d_ss(i,k)= 0.0_r8
      dtauy2d_ss(i,k)= 0.0_r8
      dtaux2d_fd(i,k)= 0.0_r8
      dtauy2d_fd(i,k)= 0.0_r8
    enddo
  enddo

  do i = its,ite
    dusfc_ls(i) = 0.0_r8
    dvsfc_ls(i) = 0.0_r8
    dusfc_bl(i) = 0.0_r8
    dvsfc_bl(i) = 0.0_r8
    dusfc_ss(i) = 0.0_r8
    dvsfc_ss(i) = 0.0_r8
    dusfc_fd(i) = 0.0_r8
    dvsfc_fd(i) = 0.0_r8
  enddo

  do i = its,ite
    taup(i,kte+1) = 0.0_r8
    xlinv(i)     = 1.0_r8/xl
  enddo
  !
  !  initialize array for flow-blocking drag
  !
  taufb(its:ite,kts:kte+1) = 0.0_r8
  komax(its:ite) = 0

  do k = kts,kte
    do i = its,ite
      vtj(i,k)  = t1(i,k)  * (1._r8+fv*q1(i,k))
      vtk(i,k)  = vtj(i,k) / prslk(i,k)
      ro(i,k)   = 1._r8/rd * prsl(i,k) / vtj(i,k) ! density kg/m**3
    enddo
  enddo
  !
  !  determine reference level: maximum of 2*var and pbl heights
  !
  do i = its,ite
    zlowtop(i) = 2._r8 * var(i)
  enddo

  do i = its,ite
    kloop1(i) = .true.
  enddo

  do k = kts+1,kte
    do i = its,ite
      if(kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
        klowtop(i) = k+1
        kloop1(i)  = .false.
      endif
    enddo
  enddo

  do i = its,ite
    kbl(i)   = max(kpbl(i), klowtop(i))
    kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
  enddo
  !
  !  determine the level of maximum orographic height
  !
  komax(:) = kbl(:)
  ! 
  do i = its,ite
    delks(i)  = 1.0_r8 / (prsi(i,1) - prsi(i,kbl(i)))
    delks1(i) = 1.0_r8 / (prsl(i,1) - prsl(i,kbl(i)))
  enddo
  !
  !  compute low level averages within pbl
  !
  do k = kts,kpblmax
    do i = its,ite
      if (k.lt.kbl(i)) then
        rcsks   = rcs     * del(i,k) * delks(i)
        rdelks  = del(i,k)  * delks(i)
        ubar(i) = ubar(i) + rcsks  * u1(i,k)      ! pbl u  mean
        vbar(i) = vbar(i) + rcsks  * v1(i,k)      ! pbl v  mean
        roll(i) = roll(i) + rdelks * ro(i,k)      ! ro mean
      endif
    enddo
  enddo
  !
  ! For ls and bl only
  IF (gsd_gwd_ls.or.gsd_gwd_bl) then
    !     figure out low-level horizontal wind direction 
    !     order into a counterclockwise index instead
    !
    do i = its,ite
      wdir  = atan2(vbar(i),ubar(i)) + pi!changed into y/x
      wdir1 = wdir-pi
      if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd  = MOD(nint(fdir*wdir1),mdir) + 1
      else!(-pi,0)
        nwd  = MOD(nint(fdir*(wdir1+2._r8*pi)),mdir) + 1
      endif
      !turn backwords because start is pi
      !need turning
      rad    = 4.0_r8*atan(1.0_r8)/180.0_r8
      theta  = (real(nwd,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
      !select OA
      oa1(i) = oa4(i,1)*cos(theta*rad)+oa4(i,2)*sin(theta*rad)
      !select OL
      ol(i)  = ol4(i,MOD(nwd-1,int(mdir/2))+1)
      !calculate dxygrid, not so slow
      call dxygrid(dxmeter(i),dymeter(i),theta,dxy(i))
      !
      !----- compute orographic width along (ol) and perpendicular (olp)
      !----- the direction of wind
      !put wdir inside the (0,2*pi) section
      !changing pi/2 either way is perpendicular
      !wdir1=wdir-pi
      if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd1  = MOD(nint(fdir*(wdir1+pi/2._r8)),mdir) + 1
        olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
      else!(-pi,0)
        nwd1  = MOD(nint(fdir*(wdir1-pi/2._r8+2._r8*pi)),mdir) + 1
        olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
      endif
      theta=(real(nwd1,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
      call dxygrid(dxmeter(i),dymeter(i),theta,dxyp(i))
      !
      !
      !----- compute orographic direction (horizontal orographic aspect ratio)
      !
      od(i) = olp(i)/max(ol(i),olmin)
      od(i) = min(od(i),odmax)
      od(i) = max(od(i),odmin)
      !
      !----- compute length of grid in the along(dxy) and cross(dxyp) wind directions
      !
    enddo
  ENDIF
  !============================================
  ! END INITIALIZATION; BEGIN GWD CALCULATIONS:
  !============================================
  IF (gsd_gwd_ls.or.gsd_gwd_bl.and.(ls_taper .GT. 1.E-02) ) THEN 
    !
    !---  saving richardson number in usqj for migwdi                       
    !
    do k = kts,kte-1
      do i = its,ite
        ti        = 2.0_r8 / (t1(i,k)+t1(i,k+1))
        rdz       = 1._r8/(zl(i,k+1) - zl(i,k))
        tem1      = u1(i,k) - u1(i,k+1)
        tem2      = v1(i,k) - v1(i,k+1)
        dw2       = rcl*(tem1*tem1 + tem2*tem2)
        shr2      = max(dw2,dw2min) * rdz * rdz
        bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti
        usqj(i,k) = max(bvf2/shr2,rimin)
        bnv2(i,k) = max(bnv_in(i,k)**2,bnv2min )
      enddo
    enddo
    !
    !----compute the "low level" or 1/3 wind magnitude (m/s)                
    !                                                                       
    do i = its,ite
      ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0_r8)
      rulow(i) = 1._r8/ulow(i)
    enddo

    do k = kts,kte-1
      do i = its,ite
        velco(i,k)  = (0.5_r8*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)                &
                                    + (v1(i,k)+v1(i,k+1)) * vbar(i))
        velco(i,k)  = velco(i,k) * rulow(i)
        if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0._r8)) then
          velco(i,k) = veleps
        endif
      enddo
    enddo
    !                                                                       
    !  no drag when critical level in the base layer                        
    !                                                                       
    do i = its,ite
      ldrag(i) = velco(i,1).le.0._r8
    enddo
    !
    !  no drag when velco.lt.0                                               
    ! 
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0._r8
      enddo
    enddo
    !                                                                       
    !  no drag when bnv2.lt.0                                               
    !                                                                       
    do k = kts,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0._r8
      enddo
    enddo
    !                                                                       
    !-----the low level weighted average ri is stored in usqj(1,1; im)      
    !-----the low level weighted average n**2 is stored in bnv2(1,1; im)    
    !---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2                      
    !---- rdelks (del(k)/delks) vert ave factor so we can * instead of /    
    !                                                                       
    do i = its,ite
      wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
      bnv2(i,1) = wtkbj * bnv2(i,1)
      usqj(i,1) = wtkbj * usqj(i,1)
    enddo
    
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) then
          rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
          bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
          usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
        endif
      enddo
    enddo
                                                                         
    do i = its,ite
      ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0_r8
      ldrag(i) = ldrag(i) .or. ulow(i)  .eq.1.0_r8
      ldrag(i) = ldrag(i) .or. var(i)   .le.0.0_r8
    enddo
    !                                                                       
    !  set all ri low level values to the low level value          
    !                                                                       
    do k = kpblmin,kpblmax
      do i = its,ite
        if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
      enddo
    enddo

    do i = its,ite
      if (.not.ldrag(i))   then
        bnv(i) = sqrt( bnv2(i,1) )
        fr(i) = bnv(i)  * rulow(i) * 2._r8 * var(i) * od(i)
        fr(i) = min(fr(i),frmax)
        xn(i)  = ubar(i) * rulow(i)
        yn(i)  = vbar(i) * rulow(i)
      endif
    enddo
    !
    !  compute the base level stress and store it in taub
    !  calculate enhancement factor, number of mountains & aspect        
    !  ratio const. use simplified relationship between standard            
    !  deviation & critical hgt                                          
    !
    do i = its,ite
      if (.not. ldrag(i))   then
        !maintain (oa+2) greater than or equal to 0
        efact    = max(oa1(i)+2._r8,0._r8) ** (ce*fr(i)/frc)
        efact    = min(max(efact,efmin),efmax)
        ! cleff (effective grid length) is highly tunable parameter
        ! the bigger (smaller) value produce weaker (stronger) wave drag
        cleff    = sqrt(dxy(i)**2._r8 + dxyp(i)**2._r8)
        !tune the times of drag
        cleff    = (3._r8/ncleff) * max(dxmax_ls,cleff)
        coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
        xlinv(i) = coefm(i) / cleff
        tem      = fr(i) * fr(i) * 1.!oc1(i)
        gfobnv   = gmax * tem / ((tem + cg)*bnv(i))
        !
        if (gsd_gwd_ls) then
           taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)                       &
                    * ulow(i) * gfobnv * efact
        else     ! We've gotten what we need for the blocking scheme
           taub(i) = 0.0_r8
        end if
      else
        taub(i) = 0.0_r8
        xn(i)   = 0.0_r8
        yn(i)   = 0.0_r8
      endif
    enddo

  ENDIF   ! (gsd_gwd_ls .eq. .true.).or.(gsd_gwd_bl .eq..true.)
  !=========================================================
  ! add small-scale wavedrag for stable boundary layer
  !=========================================================
  bnrf=0._r8
  tauwavex0=0._r8
  tauwavey0=0._r8
  utendwave=0._r8
  vtendwave=0._r8
  zq=0._r8

  IF (gsd_gwd_ss.and.(ss_taper.GT.1.E-02)) THEN
    !
    ! declaring potential temperature
    !
    do k = kts,kte
      do i = its,ite
        thx(i,k) = t1(i,k)/prslk(i,k)
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        tvcon = (1._r8+fv*q1(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      enddo
    enddo
    !
    ! Defining layer height
    !
    do k = kts,kte
      do i = its,ite
        zq(i,k+1) = dz2(i,k)+zq(i,k)
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
      enddo
    enddo

    do i=its,ite
      hpbl2 = hpbl(i)+10._r8
      kpbl2 = kpbl(i)
      kvar = 1
      do k=kts+1,MAX(kpbl(i),kts+1)
        IF (za(i,k)>300._r8) then
          kpbl2 = k
          IF (k == kpbl(i)) then
            hpbl2 = hpbl(i)+10._r8
          ELSE
            hpbl2 = za(i,k)+10._r8
          ENDIF
          exit
        ENDIF
      enddo

      if(xland1(i).gt.0._r8 .and. 2._r8*var(i).le.hpbl(i))then
        if(br1(i).gt.0._r8 .and. thvx(i,kpbl2)-thvx(i,kts) > 0._r8)then
          cleff    = sqrt(dxy(i)**2_r8 + dxyp(i)**2_r8)
          cleff    = (2.0_r8/sncleff) * max(dxmax_ss,cleff)
          coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
          xlinv(i) = coefm(i) / cleff
          govrth(i)=g/(0.5_r8*(thvx(i,kpbl2)+thvx(i,kts)))
          bnrf=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)

          if(abs(bnrf/u1(i,kpbl2)).gt.xlinv(i))then
            tauwavex0=0.5_r8*bnrf*xlinv(i)*(2._r8*MIN(var(i),varmax))**2_r8*ro(i,kvar)*u1(i,kvar)
            tauwavex0=tauwavex0*ss_taper   ! "Scale-awareness"
          else
            tauwavex0=0._r8
          endif

          if(abs(bnrf/v1(i,kpbl2)).gt.xlinv(i))then
            tauwavey0=0.5_r8*bnrf*xlinv(i)*(2._r8*MIN(var(i),varmax))**2._r8*ro(i,kvar)*v1(i,kvar)
            tauwavey0=tauwavey0*ss_taper   ! "Scale-awareness"
          else
            tauwavey0=0._r8
          endif

          do k=kts,kpbl(i) !MIN(kpbl2+1,kte-1)
            utendwave(i,k)=-1._r8*tauwavex0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
            vtendwave(i,k)=-1._r8*tauwavey0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
          enddo
        endif
      endif
    enddo ! end i loop

    do k = kts,kte
      do i = its,ite
        dudt(i,k)  = dudt(i,k) + utendwave(i,k)
        dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
        dtaux2d_ss(i,k) = utendwave(i,k)
        dtauy2d_ss(i,k) = vtendwave(i,k)
        dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
        dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
      enddo
    enddo

  ENDIF  ! end if gsd_gwd_ss == .true.
  !================================================================
  !add Beljaars et al. (2004, QJRMS, equ. 16) form drag:
  !================================================================
  IF (gsd_gwd_fd.and.(ss_taper.GT.1.E-02) ) THEN

    utendform=0._r8
    vtendform=0._r8
    zq=0._r8

    if (.not.gsd_gwd_ss.and.(ss_taper.GT.1.E-02) ) THEN
      ! Defining layer height. This is already done above is small-scale GWD is used
      do k = kts,kte
        do i = its,ite
          zq(i,k+1) = dz2(i,k)+zq(i,k)
        enddo
      enddo

      do k = kts,kte
        do i = its,ite
          za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
        enddo
      enddo
    endif

    do i=its,ite
      if (xland1(i) .gt. 0..and.2._r8*var(i).gt.0) then
        ! refer to Beljaars (2004) eq.16.
        a1=0.00026615161_r8*var(i)**2_r8
        a2=a1*0.005363_r8
        do k=kts,kte
          wsp=SQRT(u1(i,k)**2_r8 + v1(i,k)**2_r8)
          ! refer to Beljaars (2004) eq.16.
          ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759 
          utendform(i,k)=-0.0759_r8*wsp*u1(i,k)* &
                         EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
          vtendform(i,k)=-0.0759_r8*wsp*v1(i,k)* &
                         EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
        enddo
      endif
    enddo

    do k = kts,kte
      do i = its,ite
        dudt(i,k)  = dudt(i,k) + utendform(i,k)
        dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
        !limit drag tendency 
        !some tendency is likely to even overturn the wind,
        !making wind reverse in 1 timestep and reverse again in next,
        !this limitation may help to make model stable,
        !and no more wind reversal due to drag,
        !which is suppose to decelerate, not accelerate
        utendform(i,k)  = sign(min(abs(utendform(i,k)),abs(u1(i,k))/deltim),utendform(i,k))
        vtendform(i,k)  = sign(min(abs(vtendform(i,k)),abs(v1(i,k))/deltim),vtendform(i,k))
        dtaux2d_fd(i,k) = utendform(i,k)
        dtauy2d_fd(i,k) = vtendform(i,k)
        dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
        dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
      enddo
    enddo
  ENDIF  ! end if gsd_gwd_fd == .true.
  !=======================================================
  ! More for the large-scale gwd component
  !=======================================================
  IF (gsd_gwd_ls.and.(ls_taper.GT.1.E-02) ) THEN
    !                                                                       
    !   now compute vertical structure of the stress.
    !
    do k = kts,kpblmax
      do i = its,ite
        if (k .le. kbl(i)) taup(i,k) = taub(i)
      enddo
    enddo

    if (scorer_on) then
      !
      !determination of the interface height for scorer adjustment
      !
      do i=its,ite
        iint=.false.
        do k=kpblmin,kte-1
          if (k.gt.kbl(i).and.usqj(i,k)-usqj(i,k-1).lt.0.and.(.not.iint)) then
            iint=.true.
            zl_hint(i)=zl(i,k+1)
          endif
        enddo
      enddo
    endif

    do k = kpblmin, kte-1                   ! vertical level k loop!
      kp1 = k + 1
      do i = its,ite
        !
        !   unstablelayer if ri < ric
        !   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
        !   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
        !
        if (k .ge. kbl(i)) then
          !we modify the criteria for unstable layer
          !that the lv is critical under 0.25
          !while we keep wave breaking ric for
          !other larger lv
          icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric_rig)&
                                .or. (velco(i,k) .le. 0.0_r8)
          brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
          brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
        endif
      enddo

      do i = its,ite
        if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then
          if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0_r8 ) then
            temv = 1.0_r8 / velco(i,k)
            tem1 = coefm(i)/(dxy(i)/ncleff)*(ro(i,kp1)+ro(i,k))*brvf(i)*velco(i,k)*0.5_r8
            hd   = sqrt(taup(i,k) / tem1)
            fro  = brvf(i) * hd * temv
            !
            !  rim is the minimum-richardson number by shutts (1985)
            !
            tem2   = sqrt(usqj(i,k))
            tem    = 1._r8 + tem2 * fro
            rim    = usqj(i,k) * (1._r8-fro) / (tem * tem)

            !
            !  check stability to employ the 'saturation hypothesis'
            !  of lindzen (1981) except at tropospheric downstream regions
            !
            if (rim .le. ric) then  ! saturation hypothesis!
              if ((oa1(i) .le. 0._r8).or.(kp1 .ge. kpblmin )) then
                temc = 2.0_r8 + 1.0_r8 / tem2
                hd   = velco(i,k) * (2.0_r8*sqrt(temc)-temc) / brvf(i)
                taup(i,kp1) = tem1 * hd * hd
                !
                !  taup is restricted to monotoncally decrease
                !  to avoid unexpected high taup in calculation
                !
                taup(i,kp1)=min(tem1*hd*hd,taup(i,k))
                !
                !  add vertical decrease at low level below hint (Kim and Doyle 2005)
                !  where Ri first decreases
                !
                if (scorer_on.and.k.gt.klowtop(i).and.zl(i,k).le.zl_hint(i).and.k.lt.kte-1) then
                        l1=(9.81_r8*bnv2(i,kp1)/velco(i,kp1)**2)
                        l2=(9.81_r8*bnv2(i,k)/velco(i,k)**2)
                        taup(i,kp1)=min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd)
                endif
              endif
            else                    ! no wavebreaking!
              taup(i,kp1) = taup(i,k)
            endif
          endif
        endif
      enddo
    enddo

    if(lcap.lt.kte) then
      do klcap = lcapp1,kte
        do i = its,ite
          taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)
        enddo
      enddo
    endif

  ENDIF !END LARGE-SCALE TAU CALCULATION
  !===============================================================
  !COMPUTE BLOCKING COMPONENT 
  !===============================================================
  IF (gsd_gwd_bl.and.(ls_taper .GT. 1.E-02)) THEN

    do i = its,ite
      if(.not.ldrag(i)) then
        !
        !------- determine the height of flow-blocking layer
        !
        kblk = 0
        pe = 0.0_r8

        do k = kte, kpblmin, -1
          if(kblk.eq.0 .and. k.le.komax(i)) then
            !flow block appears within the reference level
            !compare potential energy and kinetic energy
            !divided by g*ro is to turn del(pa) into height
            pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g/ro(i,k)
            ke = 0.5_r8*((rcs*u1(i,k))**2._r8+(rcs*v1(i,k))**2._r8)
            !
            !---------- apply flow-blocking drag when pe >= ke 
            !
            if(pe.ge.ke) then
              kblk = k
              kblk = min(kblk,kbl(i))
              zblk = zl(i,kblk)-zl(i,kts)
            endif
          endif
        enddo

        if(kblk.ne.0) then
          !
          !--------- compute flow-blocking stress
          !

          !dxmax_ls is different than the usual one
          !because the taper is very different
          !dxy is a length scale mostly in the direction of the flow to the ridge
          !so it is good and not needed for an uneven grid area
          !ref Lott and Miller (1997) original scheme
          cd = max(2.0_r8-1.0_r8/od(i),0.0_r8)
          !
          !tuning of the drag magnitude
          cd=ncd*cd
          !
          taufb(i,kts) = 0.5_r8 * roll(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i)   &
                         * olp(i) * zblk * ulow(i)**2
          !changed grid box area into dy*dy
          tautem = taufb(i,kts)/float(kblk-kts)
          do k = kts+1, kblk
            taufb(i,k) = taufb(i,k-1) - tautem
          enddo

          !
          !----------sum orographic GW stress and flow-blocking stress
          !
          !taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
        endif
      endif
    enddo

  ENDIF   ! end blocking drag
!===========================================================
  IF (gsd_gwd_ls.OR.gsd_gwd_bl.and.(ls_taper .GT. 1.E-02)) THEN
    !                                                                       
    !  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
    !
    do k = kts,kte
      do i = its,ite
        taud_ls(i,k) = 1._r8 * (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
        taud_bl(i,k) = 1._r8 * (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
      enddo
    enddo
    !                                                                       
    !  limit de-acceleration (momentum deposition ) at top to 1/2 value 
    !  the idea is some stuff must go out the 'top'                     
    !                                                                       
    do klcap = lcap,kte
      do i = its,ite
        taud_ls(i,klcap) = taud_ls(i,klcap) * factop
        taud_bl(i,klcap) = taud_bl(i,klcap) * factop
      enddo
    enddo
    !                                                                       
    !  if the gravity wave drag would force a critical line             
    !  in the lower ksmm1 layers during the next deltim timestep,     
    !  then only apply drag until that critical line is reached.        
    !                                                                       
    do k = kts,kpblmax-1
      do i = its,ite
        if (k .le. kbl(i)) then
          if((taud_ls(i,k)+taud_bl(i,k)).ne.0._r8)                      &
            dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                 /(deltim*rcs*(taud_ls(i,k)+taud_bl(i,k)))))
        endif
      enddo
    enddo

    do k = kts,kte
      do i = its,ite
        taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper
        !apply limiter for ogwd
        !1.dudt < |c-u|/dt, so u-c cannot change sign(u^n+1 = u^n + du/dt * dt)
        !2.dudt<tndmax, eliminate ridiculous large tendency
        if (k.ne.kte) then!velco does not have top level value
          taud_ls(i,k)  = sign(min(abs(taud_ls(i,k)), 0.5*abs(velco(i,k))/deltim),taud_ls(i,k))
        endif
        taud_ls(i,k)  = sign(min(abs(taud_ls(i,k)),tndmax),taud_ls(i,k))
        taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper
        dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
        dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
        dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
        dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
        dudt(i,k)  = dtaux2d_ls(i,k) + dtaux2d_bl(i,k) + dudt(i,k)
        dvdt(i,k)  = dtauy2d_ls(i,k) + dtauy2d_bl(i,k) + dvdt(i,k)
        dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
        dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
        dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
        dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k)
      enddo
    enddo

  ENDIF

  !  Finalize dusfc and dvsfc diagnoses
  do i = its,ite
    dusfc_ls(i) = (-1._r8/g*rcs) * dusfc_ls(i)
    dvsfc_ls(i) = (-1._r8/g*rcs) * dvsfc_ls(i)
    dusfc_bl(i) = (-1._r8/g*rcs) * dusfc_bl(i)
    dvsfc_bl(i) = (-1._r8/g*rcs) * dvsfc_bl(i)
    dusfc_ss(i) = (-1._r8/g*rcs) * dusfc_ss(i)
    dvsfc_ss(i) = (-1._r8/g*rcs) * dvsfc_ss(i)
    dusfc_fd(i) = (-1._r8/g*rcs) * dusfc_fd(i)
    dvsfc_fd(i) = (-1._r8/g*rcs) * dvsfc_fd(i)
  enddo

  return

end subroutine od2d

end module oro_drag_mod
