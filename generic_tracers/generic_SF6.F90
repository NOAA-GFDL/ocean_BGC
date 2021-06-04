!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study: SF6 module
! This module contains the generic version of SF6 Tracer and its chemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP SF6
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
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_SF6

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, stdout
  use time_manager_mod, only : time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  
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

  character(len=fm_string_len), parameter :: mod_name       = 'generic_SF6'
  character(len=fm_string_len), parameter :: package_name   = 'generic_sf6'

  public do_generic_SF6
  public generic_SF6_register
  public generic_SF6_init
  public generic_SF6_register_diag
  public generic_SF6_update_from_coupler
  public generic_SF6_update_from_source
  public generic_SF6_set_boundary_values
  public generic_SF6_end
  public as_param_sf6

  !The following variables for using this module
  ! are overwritten by generic_tracer_nml namelist
  logical, save :: do_generic_SF6 = .false.
  character(len=10), save :: as_param_SF6   = 'gfdl_cmip6'

  real, parameter :: epsln=1.0e-30
  real, parameter :: missing_value1=-1.0e+10

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_SF6_params
     real :: sA, sB, sC, sD, sE   ! Coefficients in the calculation of SF6
                                  ! Schmidt number, in units of ND, degC-1, degC-2, degC-3.
     real :: A1, A2, A3, A4       ! Coefficients in the calculation of SF6
                                  ! solubility, in units of ND, K-1, log(K)^-1.
     real :: B1, B2, B3           ! More coefficients in the calculation of SF6
                                  ! solubility, in units of PSU-1, PSU-1 K-1, PSU-1 K-2.
      real :: Rho_0
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_SF6_params

  type(generic_SF6_params) :: param

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

  type generic_SF6_type
    integer :: &
      id_fgsf6 = -1, &
      id_wc_vert_int_sf6 = -1

    real, dimension(:,:), ALLOCATABLE :: &
      wc_vert_int_sf6

    real, dimension (:,:), pointer :: &
      stf_gas_sf6

    real, dimension(:,:,:,:), pointer :: &
      p_sf6

  end type generic_SF6_type

  type(generic_SF6_type) :: sf6

contains

  subroutine generic_SF6_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
    
  end subroutine generic_SF6_register

  ! <SUBROUTINE NAME="generic_SF6_init">
  !  <OVERVIEW>
  !   Initialize the generic SF6 module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the SF6 Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_SF6_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_SF6_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    call user_allocate_arrays !None for SF6 module currently

  end subroutine generic_SF6_init

  !   Register diagnostic fields to be used in this module.
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_SF6_register_diag(diag_list)
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

    vardesc_temp = vardesc("fgsf6","Surface Downward SF6 Flux",'h','1','s','mol sec-1 m-2','f')
    sf6%id_fgsf6 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_downward_mole_flux_of_sf6")

    vardesc_temp = vardesc("wc_vert_int_sf6","Total SF6 vertical integral",'h','1','s','mol m-2','f')
    sf6%id_wc_vert_int_sf6 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

  end subroutine generic_SF6_register_diag


  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,n
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    allocate(sf6%wc_vert_int_sf6(isd:ied, jsd:jed))  ; sf6%wc_vert_int_sf6=0.0

  end subroutine user_allocate_arrays

  subroutine user_deallocate_arrays

    deallocate(sf6%wc_vert_int_sf6)

  end subroutine user_deallocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_SF6_params type
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
    !      Use coefficients given by Wanninkhof (2014), L&O: Methods, 12, 351-362
    !         for SF6
    !-----------------------------------------------------------------------
    !    g_tracer_add_param(name   , variable   ,  default_value)
    if ((trim(as_param_SF6) == 'W14') .or. (trim(as_param_SF6) == 'gfdl_cmip6')) then
        call g_tracer_add_param('sA', param%sA,  3177.5)
        call g_tracer_add_param('sB', param%sB, -200.57)
        call g_tracer_add_param('sC', param%sC,  6.8865)
        call g_tracer_add_param('sD', param%sD, -0.13335)
        call g_tracer_add_param('sE', param%sE,  0.0010877)
        if (is_root_pe()) call mpp_error(NOTE,'generic_sf6: Using Schmidt number coefficients for W14')
    else
        call mpp_error(FATAL,'generic_sf6: Unable to set Schmidt number coefficients for as_param '//trim(as_param_SF6))
    endif

    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in mol/l/atm
    !      for SF6
    !     after Wanninkhof (2014), L&O: Methods, 12, 351-362
    !
    ! NOTE: Constants below DO NOT match those in Orr et al. 2017 GMDD
    !-----------------------------------------------------------------------
    if ((trim(as_param_SF6) == 'W92') .or. (trim(as_param_SF6) == 'gfdl_cmip6')) then
        call g_tracer_add_param('A1', param%A1, -80.0343)
        call g_tracer_add_param('A2', param%A2,  117.232)
        call g_tracer_add_param('A3', param%A3,  29.5817)
        call g_tracer_add_param('A4', param%A4,  0.0)        ! Not used for W92
        call g_tracer_add_param('B1', param%B1,  0.0335183)
        call g_tracer_add_param('B2', param%B2, -0.0373942)
        call g_tracer_add_param('B3', param%B3,  0.00774862)
        if (is_root_pe()) call mpp_error(NOTE,'generic_sf6: using solubility coefficients for W92')
    else if (trim(as_param_SF6) == 'W14') then
        call g_tracer_add_param('A1', param%A1, -96.5975)
        call g_tracer_add_param('A2', param%A2,  139.883)
        call g_tracer_add_param('A3', param%A3,  37.8193)
        call g_tracer_add_param('A4', param%A4,  0.0)        ! Not used for W92
        call g_tracer_add_param('B1', param%B1,  0.0310693)
        call g_tracer_add_param('B2', param%B2, -0.0356385)
        call g_tracer_add_param('B3', param%B3,  0.00743254)
        if (is_root_pe()) call mpp_error(NOTE,'generic_sf6: using solubility coefficients for W14')
    else
        call mpp_error(FATAL,'Unable to set solubility coefficients for as_param '//trim(as_param_SF6))
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
    real :: as_coeff_SF6

    if ((trim(as_param_SF6) == 'W92') .or. (trim(as_param_SF6) == 'gfdl_cmip6')) then
        ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
        ! Value is 0.337 cm/hr in units of m/s.
        as_coeff_SF6 = 9.36e-7
    else if (trim(as_param_SF6) == 'W14') then
        ! Value is 0.251 cm/hr in units of m/s
        as_coeff_SF6 = 6.972e-7
    else
        call mpp_error(FATAL,'Unable to set wind speed coefficient coefficients for as_param '//trim(as_param_SF6))
    endif

    call g_tracer_start_param_list(package_name)!nnz: Does this append?
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip_sf6.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip_sf6.res.nc' )
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
    !prog_tracers: sf6
    !diag_tracers: none
    !
    !sf6
    call g_tracer_add(tracer_list,package_name,                   &
         name = 'sf6',                                            &
         longname = 'Moles Per Unit Mass of SF6 in sea water',    &
         units = 'mol/kg',                                        &
         prog = .true.,                                           &
         flux_gas = .true.,                                       &
         flux_gas_type  = 'air_sea_gas_flux_generic',             &
         flux_gas_param = (/ as_coeff_sf6, 9.7561e-06 /),         &
         flux_gas_restart_file  = 'ocmip_sf6_airsea_flux.res.nc', &
         standard_name = "mole_concentration_of_sulfur_hexafluoride_in_sea_water", &
         diag_field_units = 'mol m-3', &
         diag_field_scaling_factor = 1035.0)   ! rho = 1035.0 kg/m3, converts mol/kg to mol/m3

  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_SF6_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for SF6s.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_SF6_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_SF6_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_update_from_copler'
    !
    !Nothing specific to be done for SF6's
    !
    return
  end subroutine generic_SF6_update_from_coupler

  ! <SUBROUTINE NAME="generic_SF6_update_from_source">
  !  <OVERVIEW>
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for SF6s.
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_SF6_update_from_source(tracer_list,rho_dzt,dzt,hblt_depth,&
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
    integer :: i,j,k
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast, grid_kmt

    logical :: used, first

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    call g_tracer_get_pointer(tracer_list,'sf6', 'stf_gas', sf6%stf_gas_sf6)
    call g_tracer_get_pointer(tracer_list,'sf6', 'field',   sf6%p_sf6)

    !-- calculate water column vertical integrals
    do j = jsc, jec ; do i = isc, iec !{
       sf6%wc_vert_int_sf6(i,j) = 0.0
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
       sf6%wc_vert_int_sf6(i,j) = sf6%wc_vert_int_sf6(i,j) + (sf6%p_sf6(i,j,k,tau) * rho_dzt(i,j,k))
    enddo; enddo; enddo  !} i,j,k

    if (sf6%id_fgsf6 .gt. 0)            &
        used = g_send_data(sf6%id_fgsf6,  sf6%stf_gas_sf6,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (sf6%id_wc_vert_int_sf6 .gt. 0)       &
       used = g_send_data(sf6%id_wc_vert_int_sf6, sf6%wc_vert_int_sf6, &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    return
  end subroutine generic_SF6_update_from_source

  ! <SUBROUTINE NAME="generic_SF6_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_SF6_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_SF6_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha,sc
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_sf6_field
    real, dimension(:,:), ALLOCATABLE :: g_sf6_alpha,g_sf6_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no

    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_set_boundary_values'


    !nnz: Can we treat these as source and move block to user_update_from_source?
    !
    !=============
    !Block Starts: Calculate the boundary values
    !=============
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'sf6','field',g_sf6_field)

    allocate(g_sf6_alpha(isd:ied, jsd:jed)); g_sf6_alpha=0.0
    allocate(g_sf6_csurf(isd:ied, jsd:jed)); g_sf6_csurf=0.0
    allocate(sc_no(isd:ied, jsd:jed))


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
       !---------------------------------------------------------------------

       if ((trim(as_param_sf6) == 'W92') .or. (trim(as_param_sf6) == 'gfdl_cmip6')) then
           alpha = conv_fac * grid_tmask(i,j,1) * &
               exp(param%A1 + param%A2/ta + param%A3*log(ta) +&
               sal * (param%B1 + ta * (param%B2 + ta * param%B3))&
               )

       else if (trim(as_param_sf6) == 'W14') then
           alpha = conv_fac * grid_tmask(i,j,1) * &
               exp(param%A1 + param%A2/ta + param%A3*log(ta) + param%A4*ta*ta +&
               sal * ((param%B3 * ta + param%B2) * ta + param%B1)&
               )
       endif

       !---------------------------------------------------------------------
       !     Calculate Schmidt numbers
       !---------------------------------------------------------------------
       if ((trim(as_param_sf6) == 'W14') .or. (trim(as_param_sf6) == 'gfdl_cmip6')) then
           sc_no(i,j) = param%sA + SST * (param%sB + SST * (param%sC + SST * (param%sD + &
               SST * param%sE))) * grid_tmask(i,j,1)
       endif

       !sc_no_term = sqrt(660.0 / (sc + epsln))
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       g_sf6_alpha(i,j) = alpha              

       g_sf6_csurf(i,j) = g_sf6_field(i,j,1,taum1) *  param%Rho_0
 
    enddo; enddo
    !=============
    !Block Ends: Calculate the boundary values
    !=============

    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'sf6','alpha',g_sf6_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'sf6','csurf',g_sf6_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'sf6','sc_no',sc_no,isd,jsd)

    deallocate(g_sf6_alpha,g_sf6_csurf,sc_no)

  end subroutine generic_SF6_set_boundary_values

  ! <SUBROUTINE NAME="generic_SF6_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_SF6_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_SF6_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_SF6_end'

    call user_deallocate_arrays

  end subroutine generic_SF6_end


end module generic_SF6
