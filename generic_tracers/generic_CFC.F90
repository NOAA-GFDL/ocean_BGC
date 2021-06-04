!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: CFC module
! This module contains the generic version of CFC Tracers and their chemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
! The chemistry calculations in this module are ported from MOM ocmip2_cfc.F90
! released in omsk_2008_03 
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 CFC
!       simulations as outlined in the CFC-HOWTO documentation,
!       revision 1.6, 1999/04/29.  Updated by JPD for revised Schmidt numbers
!       of Wanninkhof (2014).
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/CFC/HOWTO-CFC.html
! Wanninkhof, R. (2014) Relationship between wind speed and gas exchange over the ocean revisited.
! Limnology and Oceanography: Methods, 12, 351-362, DOI:10.4319/lom.2014.12.351.
! </REFERENCE>
! <DEVELOPER_NOTES>
! nnz: 
! A. Reproducing GOLD results
! The tracers in this module reproduce the corresponding ones
! in non-generic module GOLD_OMIP2_CFC.F90 with branch tag perth_gas_fluxes_nnz.
! The reproducing of non-generic tracers in this case is sufficient evidence for the consistency of vertical diffusion
! (tracer_vertdiff) routine used in the two cases. 
! 
! B. Reproducing MOM results
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_CFC

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, stdout
  use time_manager_mod, only: time_type
  use fm_util_mod,      only: fm_util_start_namelist, fm_util_end_namelist  
  use fms_mod,          only: stdout, stdlog, mpp_pe, mpp_root_pe

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  

  use g_tracer_utils, only : g_diag_type, g_diag_field_add
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use g_tracer_utils, only : g_send_data, is_root_pe

  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_CFC'
  character(len=fm_string_len), parameter :: package_name   = 'generic_cfc'

  public do_generic_CFC
  public generic_CFC_register
  public generic_CFC_init
  public generic_CFC_register_diag
  public generic_CFC_update_from_coupler
  public generic_CFC_update_from_source
  public generic_CFC_set_boundary_values
  public generic_CFC_end
  public as_param_cfc

  !The following variables for using this module
  ! are overwritten by generic_tracer_nml namelist
  logical, save :: do_generic_CFC = .false.
  character(len=10), save :: as_param_cfc   = 'gfdl_cmip6'

  real, parameter :: epsln=1.0e-30
  real, parameter :: missing_value1=-1.0e+10

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_CFC_params
     real :: a1_11, a2_11, a3_11, a4_11   ! Coefficients in the calculation of the
     real :: a1_12, a2_12, a3_12, a4_12   ! CFC11 and CFC12 Schmidt numbers, in
     
     real :: b1_11, b2_11, b3_11          ! More coefficients in the calculation of
     real :: b1_12, b2_12, b3_12          ! the CFC11 and CFC12 solubilities, in
                                          ! units of PSU-1, PSU-1 K-1, PSU-1 K-2.

     real :: sA_11, sB_11, sC_11, sD_11, sE_11   ! Coefficients in the calculation of the
     real :: sA_12, sB_12, sC_12, sD_12, sE_12   ! CFC11 and CFC12 Schmidt numbers, in
                                                 ! units of ND, degC-1, degC-2, degC-3.
     real :: Rho_0

     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_CFC_params


  type(generic_CFC_params) :: param

  !An auxiliary type for storing varible names
  type, public :: vardesc
     character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
     character(len=fm_string_len) :: longname ! The long name of that variable.
     character(len=1)  :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid   ! The vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid   ! The time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    ! The dimensions of the variable.
     character(len=1)  :: mem_size ! The size in memory: d or f.
  end type vardesc

  type generic_CFC_type
    integer :: &
      id_fgcfc11      = -1,    &
      id_fgcfc12      = -1,    &
      id_wc_vert_int_cfc11 = -1, &
      id_wc_vert_int_cfc12 = -1

    real, dimension(:,:), ALLOCATABLE :: &
      wc_vert_int_cfc11, &
      wc_vert_int_cfc12

    real, dimension (:,:), pointer :: &
      stf_gas_cfc11, &
      stf_gas_cfc12

    real, dimension(:,:,:,:), pointer :: &
      p_cfc11, &
      p_cfc12

  end type generic_CFC_type

  type(generic_CFC_type) :: cfc

contains

  subroutine generic_CFC_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_CFC_register

  ! <SUBROUTINE NAME="generic_CFC_init">
  !  <OVERVIEW>
  !   Initialize the generic CFC module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the CFC Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_CFC_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_CFC_init

  !   Register diagnostic fields to be used in this module.
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_CFC_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3)
    type(time_type):: init_time
    character(len=fm_string_len)          :: cmor_field_name
    character(len=fm_string_len)          :: cmor_long_name
    character(len=fm_string_len)          :: cmor_units
    character(len=fm_string_len)          :: cmor_standard_name
!    real                                  :: conversion

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes=axes,init_time=init_time)

    !   The following vardesc types contain a package of metadata about each tracer,
    ! including, in order, the following elements: name; longname; horizontal
    ! staggering ('h') for collocation with thickness points ; vertical staggering
    ! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
    ! and precision in non-restart output files ('f' for 32-bit float or 'd' for
    ! 64-bit doubles). For most tracers, only the name, longname and units should
    ! be changed.

    vardesc_temp = vardesc("fgcfc11","Surface Downward CFC11 flux",'h','1','s','mol sec-1 m-2','f')
    cfc%id_fgcfc11 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_downward_mole_flux_of_cfc11")

    vardesc_temp = vardesc("fgcfc12","Surface Downward CFC12 flux",'h','1','s','mol sec-1 m-2','f')
    cfc%id_fgcfc12 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_downward_mole_flux_of_cfc12")

    vardesc_temp = vardesc("wc_vert_int_cfc11","Total CFC11 vertical integral",'h','1','s','mol m-2','f')
    cfc%id_wc_vert_int_cfc11 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_cfc12","Total CFC12 vertical integral",'h','1','s','mol m-2','f')
    cfc%id_wc_vert_int_cfc12 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

  end subroutine generic_CFC_register_diag


  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,n
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    allocate(cfc%wc_vert_int_cfc11(isd:ied, jsd:jed))  ; cfc%wc_vert_int_cfc11=0.0
    allocate(cfc%wc_vert_int_cfc12(isd:ied, jsd:jed))  ; cfc%wc_vert_int_cfc12=0.0

  end subroutine user_allocate_arrays

  subroutine user_deallocate_arrays

    deallocate(cfc%wc_vert_int_cfc11)
    deallocate(cfc%wc_vert_int_cfc12)

  end subroutine user_deallocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_CFC_params type
    !==============================================================    

    !=============
    !Block Starts: g_tracer_add_param
    !=============
    !Add the known experimental parameters used for calculations
    !in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and g_tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.

    call g_tracer_start_param_list(package_name)

    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !      Use coefficients given by Zheng et al (1998), JGR vol 103, C1
    !         for CFC11 and CFC12
    !-----------------------------------------------------------------------
    !    g_tracer_add_param(name   , variable   ,  default_value)
    if (trim(as_param_cfc) == 'W92') then
        call g_tracer_add_param('sA_11', param%sA_11,  3501.8)
        call g_tracer_add_param('sB_11', param%sB_11, -210.31)
        call g_tracer_add_param('sC_11', param%sC_11,  6.1851)
        call g_tracer_add_param('sD_11', param%sD_11, -0.07513)
        call g_tracer_add_param('sE_11', param%sE_11,  0.0)      ! Not used for W92
        call g_tracer_add_param('sA_12', param%sA_12,  3845.4)
        call g_tracer_add_param('sB_12', param%sB_12, -228.95)
        call g_tracer_add_param('sC_12', param%sC_12,  6.1908)
        call g_tracer_add_param('sD_12', param%sD_12, -0.067430)
        call g_tracer_add_param('sE_12', param%sE_12,  0.0)      ! Not used for W92
        if (is_root_pe()) call mpp_error(NOTE,'generic_cfc: using Schmidt number coefficients for W92')
    else if ((trim(as_param_cfc) == 'W14') .or. (trim(as_param_cfc) == 'gfdl_cmip6')) then 
        call g_tracer_add_param('sA_11', param%sA_11,  3579.2)
        call g_tracer_add_param('sB_11', param%sB_11, -222.63)
        call g_tracer_add_param('sC_11', param%sC_11,  7.5749)
        call g_tracer_add_param('sD_11', param%sD_11, -0.14595)
        call g_tracer_add_param('sE_11', param%sE_11,  0.0011874)
        call g_tracer_add_param('sA_12', param%sA_12,  3828.1)
        call g_tracer_add_param('sB_12', param%sB_12, -249.86)
        call g_tracer_add_param('sC_12', param%sC_12,  8.7603)
        call g_tracer_add_param('sD_12', param%sD_12, -0.1716)
        call g_tracer_add_param('sE_12', param%sE_12,  0.001408)
        if (is_root_pe()) call mpp_error(NOTE,'generic_cfc: using Schmidt number coefficients for W14')
    else
        call mpp_error(FATAL,'Unable to set Schmidt number coefficients for as_param '//trim(as_param_cfc))
    endif

    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in mol/l/atm
    !      (1) for CFC11, (2) for CFC12
    !     after Warner and Weiss (1985) DSR, vol 32 for CFC11 and CFC12
    !-----------------------------------------------------------------------
    if ((trim(as_param_cfc) == 'W92') .or. (trim(as_param_cfc) == 'gfdl_cmip6')) then
        call g_tracer_add_param('A1_11', param%A1_11, -229.9261)
        call g_tracer_add_param('A2_11', param%A2_11,  319.6552)
        call g_tracer_add_param('A3_11', param%A3_11,  119.4471)
        call g_tracer_add_param('A4_11', param%A4_11, -1.39165)
        call g_tracer_add_param('B1_11', param%B1_11, -0.142382)
        call g_tracer_add_param('B2_11', param%B2_11,  0.091459)
        call g_tracer_add_param('B3_11', param%B3_11, -0.0157274)
        call g_tracer_add_param('A1_12', param%A1_12, -218.0971)
        call g_tracer_add_param('A2_12', param%A2_12,  298.9702)
        call g_tracer_add_param('A3_12', param%A3_12,  113.8049)
        call g_tracer_add_param('A4_12', param%A4_12, -1.39165)
        call g_tracer_add_param('B1_12', param%B1_12, -0.143566)
        call g_tracer_add_param('B2_12', param%B2_12,  0.091015)
        call g_tracer_add_param('B3_12', param%B3_12, -0.0153924)
        if (is_root_pe()) call mpp_error(NOTE,'generic_cfc: using solubility coefficients for W92')
    else if (trim(as_param_cfc) == 'W14') then
        call g_tracer_add_param('A1_11', param%A1_11, -134.1536)
        call g_tracer_add_param('A2_11', param%A2_11,  203.2156)
        call g_tracer_add_param('A3_11', param%A3_11,  56.2320)
        call g_tracer_add_param('A4_11', param%A4_11,  0.0)        ! Not used in W14
        call g_tracer_add_param('B1_11', param%B1_11, -0.144449)
        call g_tracer_add_param('B2_11', param%B2_11,  0.092952)
        call g_tracer_add_param('B3_11', param%B3_11, -0.0159977)
        call g_tracer_add_param('A1_12', param%A1_12, -122.3246)
        call g_tracer_add_param('A2_12', param%A2_12,  182.5306)
        call g_tracer_add_param('A3_12', param%A3_12,  50.5898)
        call g_tracer_add_param('A4_12', param%A4_12,  0.0)        ! Not used in W14
        call g_tracer_add_param('B1_12', param%B1_12, -0.145633)
        call g_tracer_add_param('B2_12', param%B2_12,  0.092509)
        call g_tracer_add_param('B3_12', param%B3_12, -0.0156627)
        if (is_root_pe()) call mpp_error(NOTE,'generic_cfc: using solubility coefficients for W14')
    else
        call mpp_error(FATAL,'Unable to set solubility coefficients for as_param '//trim(as_param_cfc))
    endif

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', param%Rho_0, 1035.0)

    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========

  end subroutine user_add_params

  !
  !   This is an internal sub, not a public interface.
  !   Add all the tracers to be used in this module. 
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    real :: as_coeff_cfc

    if ((trim(as_param_cfc) == 'W92') .or. (trim(as_param_cfc) == 'gfdl_cmip6')) then
        ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
        ! Value is 0.337 cm/hr in units of m/s.
        as_coeff_cfc = 9.36e-7
    else if (trim(as_param_cfc) == 'W14') then
        ! Value is 0.251 cm/hr in units of m/s
        as_coeff_cfc = 6.972e-7
    else
        call mpp_error(FATAL,'Unable to set wind speed coefficient coefficients for as_param '//trim(as_param_cfc))
    endif

    call g_tracer_start_param_list(package_name)!nnz: Does this append?
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip2_cfc.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip2_cfc.res.nc' )
    call g_tracer_add_param('IC_file'       , param%IC_file       , '')
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file=param%ice_restart_file, ocean_restart_file=param%ocean_restart_file )

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !prog_tracers: cfc11,cfc12 
    !diag_tracers: none
    !
    !cfc12

    call g_tracer_add(tracer_list,package_name,                      &
         name       = 'cfc12',                                       &
         longname   = 'Moles Per Unit Mass of CFC-12 in sea water',  &
         units      = 'mol/kg',                                      &
         prog       = .true.,                                        &
         requires_src_info  = .false.,                               &
         flux_gas       = .true.,                                    &
         flux_gas_type  = 'air_sea_gas_flux_generic',                &
         flux_gas_param = (/ as_coeff_cfc, 9.7561e-06 /),            &
         flux_gas_restart_file  = 'ocmip2_cfc_airsea_flux.res.nc',   &
         standard_name = "mole_concentration_of_cfc12_in_sea_water", &
         diag_field_units = 'mol m-3',                               &
         diag_field_scaling_factor = 1035.0)   ! rho = 1035.0 kg/m3, converts mol/kg to mol/m3

    !cfc11
    call g_tracer_add(tracer_list,package_name,                      &
         name       = 'cfc11',                                       &
         longname   = 'Moles Per Unit Mass of CFC-11 in sea water',  &
         units      = 'mol/kg',                                      &
         prog       = .true.,                                        &
         requires_src_info  = .false.,                               &
         flux_gas       = .true.,                                    &
         flux_gas_type  = 'air_sea_gas_flux_generic',                &
         flux_gas_param = (/ as_coeff_cfc, 9.7561e-06 /),            &
         flux_gas_restart_file  = 'ocmip2_cfc_airsea_flux.res.nc',   &
         standard_name = "mole_concentration_of_cfc11_in_sea_water", &
         diag_field_units = 'mol m-3',                               &
         diag_field_scaling_factor = 1035.0)   ! rho = 1035.0 kg/m3, converts mol/kg to mol/m3

  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_CFC_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_CFC_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_update_from_copler'
    !
    !Nothing specific to be done for CFC's
    !
    return
  end subroutine generic_CFC_update_from_coupler

  ! <SUBROUTINE NAME="generic_CFC_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for CFCs.
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_CFC_update_from_source(tracer_list,rho_dzt,dzt,hblt_depth,&
                                            ilb,jlb,tau,dt,grid_dat,model_time)

    type(g_tracer_type), pointer :: tracer_list
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dt
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time

    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau 
    integer :: i, j, k
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast, grid_kmt

    logical :: used, first

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    call g_tracer_get_pointer(tracer_list,'cfc11','stf_gas',cfc%stf_gas_cfc11)
    call g_tracer_get_pointer(tracer_list,'cfc12','stf_gas',cfc%stf_gas_cfc12)

    call g_tracer_get_pointer(tracer_list,'cfc11','field',  cfc%p_cfc11)
    call g_tracer_get_pointer(tracer_list,'cfc12','field',  cfc%p_cfc12)

    !-- calculate water column vertical integrals
    do j = jsc, jec ; do i = isc, iec !{
       cfc%wc_vert_int_cfc11(i,j) = 0.0
       cfc%wc_vert_int_cfc12(i,j) = 0.0
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
       cfc%wc_vert_int_cfc11(i,j) = cfc%wc_vert_int_cfc11(i,j) + (cfc%p_cfc11(i,j,k,tau) * rho_dzt(i,j,k))
       cfc%wc_vert_int_cfc12(i,j) = cfc%wc_vert_int_cfc12(i,j) + (cfc%p_cfc12(i,j,k,tau) * rho_dzt(i,j,k))
    enddo; enddo; enddo  !} i,j,k

    if (cfc%id_fgcfc11 .gt. 0)            &
        used = g_send_data(cfc%id_fgcfc11,  cfc%stf_gas_cfc11,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cfc%id_fgcfc12 .gt. 0)            &
        used = g_send_data(cfc%id_fgcfc12,  cfc%stf_gas_cfc12,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cfc%id_wc_vert_int_cfc11 .gt. 0)       &
       used = g_send_data(cfc%id_wc_vert_int_cfc11, cfc%wc_vert_int_cfc11, &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cfc%id_wc_vert_int_cfc12 .gt. 0)       &
       used = g_send_data(cfc%id_wc_vert_int_cfc12, cfc%wc_vert_int_cfc12, &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    return
  end subroutine generic_CFC_update_from_source

  ! <SUBROUTINE NAME="generic_CFC_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="SST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature   
  !  </IN>
  !  <IN NAME="SSS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  ! </SUBROUTINE>

  !User must provide the calculations for these boundary values.
  subroutine generic_CFC_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha_11,alpha_12,sc_11,sc_12
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_cfc11_field,g_cfc12_field
    real, dimension(:,:), ALLOCATABLE :: g_cfc11_alpha,g_cfc11_csurf,g_cfc12_alpha,g_cfc12_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no_11,sc_no_12

    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_set_boundary_values'


    !nnz: Can we treat these as source and move block to user_update_from_source?
    !
    !=============
    !Block Starts: Calculate the boundary values
    !=============
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'cfc11','field',g_cfc11_field)
    call g_tracer_get_pointer(tracer_list,'cfc12','field',g_cfc12_field)

    allocate(g_cfc11_alpha(isd:ied, jsd:jed)); g_cfc11_alpha=0.0
    allocate(g_cfc11_csurf(isd:ied, jsd:jed)); g_cfc11_csurf=0.0
    allocate(g_cfc12_alpha(isd:ied, jsd:jed)); g_cfc12_alpha=0.0
    allocate(g_cfc12_csurf(isd:ied, jsd:jed)); g_cfc12_csurf=0.0
    allocate(sc_no_11(isd:ied, jsd:jed))
    allocate(sc_no_12(isd:ied, jsd:jed))

    !The atmospheric code needs soluabilities in units of mol/m3/atm
    !
    !MOM
    !       The factor 1.0e+03 is for the conversion  
    !       from mol/(l * atm) to mol/(m3 * atm) 
    conv_fac = 1.0e+03
    !
    !GOLD
    !       The factor 1.e-09 converts 
    !       from mol/(l * atm) to mol/(m3 * pptv).
    !conv_fac = 1.0e-09

    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       ta = (ST(i,j) + 273.15) * 0.01 ! Why is this in dekaKelvin?
       sal = SSS(i,j) ; SST = ST(i,j)

       !---------------------------------------------------------------------
       !     Calculate solubilities
       !       Use Warner and Weiss (1985) DSR, vol 32, final result
       !       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
       !
       !       Use Bullister and Wisegavger for CCl4
       !---------------------------------------------------------------------
       if ((trim(as_param_cfc) == 'W92') .or. (trim(as_param_cfc) == 'gfdl_cmip6')) then
           alpha_11 = conv_fac * grid_tmask(i,j,1) * &
                exp(param%A1_11 + param%A2_11/ta + param%A3_11*log(ta) + param%A4_11*ta*ta +&
                sal * ((param%B3_11 * ta + param%B2_11) * ta + param%B1_11)&
                )
           alpha_12 = conv_fac * grid_tmask(i,j,1) * &
                exp(param%A1_12 + param%A2_12/ta + param%A3_12*log(ta) + param%A4_12*ta*ta +&
                sal * ((param%B3_12 * ta + param%B2_12) * ta + param%B1_12)&
                )
       else if (trim(as_param_cfc) == 'W14') then
           alpha_11 = conv_fac * grid_tmask(i,j,1) * &
                exp(param%A1_11 + param%A2_11/ta + param%A3_11*log(ta) +&
                sal * (param%B1_11 + ta * (param%B2_11 + ta * param%B3_11))&
                )
           alpha_12 = conv_fac * grid_tmask(i,j,1) * &
                exp(param%A1_12 + param%A2_12/ta + param%A3_12*log(ta) +&
                sal * (param%B1_12 + ta * (param%B2_12 + ta * param%B3_12))&
                )
       endif

       !---------------------------------------------------------------------
       !     Calculate Schmidt numbers
       !      use coefficients given by Zheng et al (1998), JGR vol 103, C1
       !---------------------------------------------------------------------
       if (trim(as_param_cfc) == 'W92') then
           sc_no_11(i,j) = param%sA_11 + SST * (param%sB_11 + SST * (param%sC_11 + SST * param%sD_11)) * &
                grid_tmask(i,j,1)
           sc_no_12(i,j) = param%sA_12 + SST * (param%sB_12 + SST * (param%sC_12 + SST * param%sD_12)) * &
                grid_tmask(i,j,1)
       else if ((trim(as_param_cfc) == 'W14') .or. (trim(as_param_cfc) == 'gfdl_cmip6')) then
           sc_no_11(i,j) = param%sA_11 + SST * (param%sB_11 + SST * (param%sC_11 + SST * (param%sD_11 + &
                SST * param%sE_11))) * grid_tmask(i,j,1)
           sc_no_12(i,j) = param%sA_12 + SST * (param%sB_12 + SST * (param%sC_12 + SST * (param%sD_12 + &
                SST * param%sE_12))) * grid_tmask(i,j,1)
       endif

       !sc_no_term = sqrt(660.0 / (sc_11 + epsln))
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       g_cfc11_alpha(i,j) = alpha_11              

       g_cfc11_csurf(i,j) = g_cfc11_field(i,j,1,taum1) *  param%Rho_0

       g_cfc12_alpha(i,j) = alpha_12              

       g_cfc12_csurf(i,j) = g_cfc12_field(i,j,1,taum1) *  param%Rho_0
 
    enddo; enddo
    !=============
    !Block Ends: Calculate the boundary values
    !=============

    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'cfc11','alpha',g_cfc11_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc11','csurf',g_cfc11_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc11','sc_no',sc_no_11,isd,jsd)

    call g_tracer_set_values(tracer_list,'cfc12','alpha',g_cfc12_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc12','csurf',g_cfc12_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'cfc12','sc_no',sc_no_12,isd,jsd)

    deallocate(g_cfc11_alpha,g_cfc11_csurf,g_cfc12_alpha,g_cfc12_csurf,sc_no_11,sc_no_12)

  end subroutine generic_CFC_set_boundary_values

  ! <SUBROUTINE NAME="generic_CFC_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_CFC_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_CFC_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_CFC_end'

    call user_deallocate_arrays

  end subroutine generic_CFC_end


end module generic_CFC
