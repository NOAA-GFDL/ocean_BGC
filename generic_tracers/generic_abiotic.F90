!----------------------------------------------------------------
! <CONTACT EMAIL="John.Krasting@noaa.gov"> John Krasting 
! </CONTACT>
! 
! <REVIEWER EMAIL="First.Last@noaa.gov"> NULL
! </REVIEWER>
!
! <OVERVIEW>
! Ocean Carbon Model Intercomparison Study: abiotic carbon module
! This module contains the generic version of the abiotic DIC and
! DI14C tracers along with associated chemistry.
! </OVERVIEW>
!
! <DESCRIPTION>
! Implementation of routines to solve the OCMIP abiotic
! simulations as outlined in the CMIP6 Ocen Model Intercomparison
! Project (OMIP) documentation
! </DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! Orr et al. Biogeochemical protocols and diagnostics for the CMIP6 
!   Ocean Model Intercomparison Project (OMIP).  Geosci. Model Dev
!   Discuss., doi:10.5194/gmd-2016-155  (2016).
! </REFERENCE>
! <DEVELOPER_NOTES>
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_abiotic

  use constants_mod,     only: WTMCO2
  use data_override_mod, only: data_override
  use field_manager_mod, only: fm_string_len
  use time_manager_mod,  only: time_type

  use g_tracer_utils,    only: g_diag_type,g_send_data
  use g_tracer_utils,    only: g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils,    only: g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils,    only: g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils,    only: g_tracer_get_values  
  use g_tracer_utils,    only: register_diag_field=>g_register_diag_field

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc,CO2_dope_vector

  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_abiotic'
  character(len=fm_string_len), parameter :: package_name   = 'generic_abiotic'

  public do_generic_abiotic
  public generic_abiotic_register
  public generic_abiotic_init
  public generic_abiotic_register_diag
  public generic_abiotic_update_from_source
  public generic_abiotic_set_boundary_values
  public generic_abiotic_end

  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_abiotic = .false.

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30
  real, parameter :: missing_value1=-1.0e+10
  real, parameter :: missing_value_diag=-1.0e+10

  !This type contains all the parameters and arrays used in this module.
  type generic_abiotic_params
     logical :: init              ! Initialize tracers
     logical :: force_update_fluxes  ! If OCMIP2 tracers fluxes should be updated 
                                     ! every coupling time step when update_from_source 
                                     ! is not called every coupling time step as is the 
                                     ! case with MOM6  THERMO_SPANS_COUPLING option
     real :: atm_delta_14c        ! Ratio of atmospheric 14C to 12C
     real :: half_life            ! Decay time scale of 14C (years)
     real :: lambda_14c           ! Radioactive decay constant for 14C (s-1)
     real :: htotal_in            ! Initial "first guess" for H+ concentration (mol/kg)
     real :: htotal_scale_hi      ! Initial upper limit of htotal range 
     real :: htotal_scale_lo      ! Initial lower limit of htotal range 
     real :: alkbar               ! Mean global alkalinity (eq/kg)
     real :: sio4_const           ! Silicate (SiO4) concentration (mol/kg)
     real :: po4_const            ! Phosphate (PO4) concentration (mol/kg)
     real :: a1_co2               ! Wanninkhof Coeffs.
     real :: a2_co2
     real :: a3_co2
     real :: a4_co2
     real :: Rho_0                ! Reference density (kg/m^3)

     ! Restart file names
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     ! Diagnostic Output IDs
     integer :: id_sfc_dissicabio=-1, id_sfc_dissi14cabio=-1, id_sfc_ab_htotal=-1
     integer :: id_ab_alk=-1, id_ab_po4=-1, id_ab_sio4=-1
     integer :: id_sfc_ab_alk=-1, id_sfc_ab_po4=-1, id_sfc_ab_sio4=-1
     integer :: id_ab_pco2surf=-1, id_ab_p14co2surf=-1
     integer :: id_jdecay_di14c=-1, id_delta_14catm=-1

     ! 2-dimensional fields
     real, dimension(:,:), ALLOCATABLE ::   &
          htotalhi, htotallo,               &! Upper and lower limits of htotal range 
          abco2_csurf,abco2_alpha,          &! Abiotic DIC csurf & alpha
          ab14co2_csurf,ab14co2_alpha,      &! Abiotic DI14C csurf & alpha
          abpco2_csurf,abp14co2_csurf,      &! Oceanic pCO2 (ppmv)
          delta_14catm                       ! Atm. Delta of 14C 
 
     ! 3-dimensional fields
     real, dimension(:,:,:), ALLOCATABLE :: &
          f_dissicabio, f_dissi14cabio,     &! Abiotic DIC & DI14C fields
          f_alk, f_po4, f_sio4,             &! Abiotic alkalinity, phosphate, silicate
          f_htotal,                         &! Abiotic H+ concentration
          jdecay_di14c                       ! DI14C decay rate

     ! 4-dimensional pointers
     real, dimension(:,:,:,:), pointer :: &
          p_dissicabio, p_dissi14cabio, p_htotal

  end type generic_abiotic_params

  ! An auxiliary type for storing varible names
  type, public :: vardesc
     character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
     character(len=fm_string_len) :: longname ! The long name of that variable.
     character(len=1)  :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid   ! The vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid   ! The time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    ! The dimensions of the variable.
     character(len=1)  :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(CO2_dope_vector)        :: CO2_dope_vec
  type(generic_abiotic_params) :: abiotic

contains

  subroutine generic_abiotic_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_abiotic_register



  ! <SUBROUTINE NAME="generic_abiotic_init">
  !  <OVERVIEW>
  !   Initialize the generic abiotic module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds the dissicabio, dissi14cabio, ab_htotal tracers to the list of generic tracers passed 
  !       to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_abiotic_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_abiotic_init(tracer_list, force_update_fluxes)
    type(g_tracer_type), pointer :: tracer_list
    logical          ,intent(in) :: force_update_fluxes

    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_init'

    !There are situations where the column_physics (update_from_source) is not called every timestep 
    ! such as when MOM6 THERMO_SPANS_COUPLING=True , yet we want the fluxes to be updated every timestep
    ! In that case we can force an update by setting the namelist generic_tracer_nml:force_update_fluxes=.true.
    abiotic%force_update_fluxes = force_update_fluxes

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    !    call user_allocate_arrays !None for abiotic module currently
    call user_allocate_arrays

  end subroutine generic_abiotic_init





  subroutine generic_abiotic_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3)
    type(time_type):: init_time 

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes=axes,init_time=init_time) 

    !   The following vardesc types contain a package of metadata about each tracer,
    ! including, in order, the following elements: name; longname; horizontal
    ! staggering ('h') for collocation with thickness points ; vertical staggering
    ! ('L') for a layer variable ; temporal staggering ('s' for snapshot) ; units ;
    ! and precision in non-restart output files ('f' for 32-bit float or 'd' for
    ! 64-bit doubles). For most tracers, only the name, longname and units should
    ! be changed.  

    vardesc_temp = vardesc("ab_alk","Abiotic Alkalinity",'h','1','s','eq kg-1','f')
    abiotic%id_ab_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ab_po4","Abiotic PO4",'h','1','s','mol kg-1','f')
    abiotic%id_ab_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ab_sio4","Abiotic SiO4",'h','1','s','mol kg-1','f')
    abiotic%id_ab_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("delta_14catm","Atmospheric Delta of 14C",'h','1','s','mil -1','f')
    abiotic%id_delta_14catm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ab_alk","Surface Abiotic Alkalinity",'h','1','s','eq kg-1','f')
    abiotic%id_sfc_ab_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ab_po4","Surface Abiotic PO4",'h','1','s','eq kg-1','f')
    abiotic%id_sfc_ab_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ab_sio4","Surface Abiotic SiO4",'h','1','s','eq kg-1','f')
    abiotic%id_sfc_ab_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_dissicabio","Surface Abiotic Dissolved Inorganic Carbon",&
                           'h','1','s','mol kg-1','f')
    abiotic%id_sfc_dissicabio = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_dissi14cabio","Surface Abiotic Dissolved Inorganic Radiocarbon",&
                           'h','1','s','mol kg-1','f')
    abiotic%id_sfc_dissi14cabio = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ab_htotal","Surface Abiotic Htotal",'h','1','s','mol kg-1','f')
    abiotic%id_sfc_ab_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ab_pco2surf","Oceanic Abiotic pCO2",'h','1','s','uatm','f')
    abiotic%id_ab_pco2surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ab_p14co2surf","Oceanic Abiotic p14CO2",'h','1','s','uatm','f')
    abiotic%id_ab_p14co2surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdecay_di14c","DI14C radioactive decay layer integral",&
                           'h','L','s','mol m-2 s-1','f')
    abiotic%id_jdecay_di14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

  end subroutine generic_abiotic_register_diag





  subroutine user_allocate_arrays
    !Allocate all the private arrays.
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    !Allocate all the private arrays.

    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc ; CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc ; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd ; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd ; CO2_dope_vec%jed = jed

    allocate(abiotic%htotallo(isd:ied,jsd:jed))
    allocate(abiotic%htotalhi(isd:ied,jsd:jed))
    allocate(abiotic%delta_14catm(isd:ied,jsd:jed))
    allocate(abiotic%f_alk(isd:ied,jsd:jed,1:nk))          ; abiotic%f_alk=0.0
    allocate(abiotic%f_htotal(isd:ied,jsd:jed,1:nk))       ; abiotic%f_htotal =0.0
    allocate(abiotic%f_dissicabio(isd:ied,jsd:jed,1:nk))   ; abiotic%f_dissicabio=0.0
    allocate(abiotic%f_dissi14cabio(isd:ied,jsd:jed,1:nk)) ; abiotic%f_dissi14cabio= 0.0
    allocate(abiotic%f_po4(isd:ied,jsd:jed,1:nk))          ; abiotic%f_po4=0.0
    allocate(abiotic%f_sio4(isd:ied,jsd:jed,1:nk))         ; abiotic%f_sio4=0.0
    allocate(abiotic%jdecay_di14c(isd:ied,jsd:jed, 1:nk))  ; abiotic%jdecay_di14c=0.0
    allocate(abiotic%abco2_csurf(isd:ied,jsd:jed))         ; abiotic%abco2_csurf=0.0
    allocate(abiotic%abco2_alpha(isd:ied,jsd:jed))         ; abiotic%abco2_alpha= 0.0
    allocate(abiotic%ab14co2_csurf(isd:ied,jsd:jed))       ; abiotic%ab14co2_csurf=0.0
    allocate(abiotic%ab14co2_alpha(isd:ied,jsd:jed))       ; abiotic%ab14co2_alpha=0.0
    allocate(abiotic%abpco2_csurf(isd:ied,jsd:jed))        ; abiotic%abpco2_csurf=0.0
    allocate(abiotic%abp14co2_csurf(isd:ied,jsd:jed))      ; abiotic%abp14co2_csurf=0.0

  end subroutine user_allocate_arrays





  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_abiotic_params type
    !==============================================================    

    !=============
    !Block Starts: g_tracer_add_param
    !=============
    !Add the known experimental parameters used for calculations
    !in this module.
    !All the g_tracer_add_param calls must happen between 
    !g_tracer_start_param_list and _tracer_end_param_list  calls.
    !This implementation enables runtime overwrite via field_table.

    call g_tracer_start_param_list(package_name)
    !-----------------------------------------------------------------------
    ! Initialization
    !-----------------------------------------------------------------------
    call g_tracer_add_param('init', abiotic%init, .false. )
    !-----------------------------------------------------------------------
    ! Constants for silicate and phosphate
    !-----------------------------------------------------------------------
    !    g_tracer_add_param(name   , variable   ,  default_value)
    call g_tracer_add_param('sio4_const', abiotic%sio4_const, 7.5e-6) !mol/kg
    call g_tracer_add_param('po4_const',  abiotic%po4_const,  5.0e-7) !mol/kg
    !-----------------------------------------------------------------------
    ! CO2 Solubility coefficients (Wanninkhof numbers)
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a1_co2', abiotic%a1_co2, 2068.9)
    call g_tracer_add_param('a2_co2', abiotic%a2_co2, -118.63)
    call g_tracer_add_param('a3_co2', abiotic%a3_co2, 2.9311)
    call g_tracer_add_param('a4_co2', abiotic%a4_co2, -0.027)
    !-----------------------------------------------------------------------
    ! H+ Concentration Parameters
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_in',       abiotic%htotal_in, 1.0e-8)
    call g_tracer_add_param('htotal_scale_lo', abiotic%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi', abiotic%htotal_scale_hi, 100.)
    !-----------------------------------------------------------------------
    ! Global Concentrations
    ! jpk: how-to doc specifies Alkbar as 2310 microeq/kg and value below 
    !      is from MOM4p1 OCMIP2 abiotic routine.  Is unit conversion correct?
    !-----------------------------------------------------------------------
    call g_tracer_add_param('half_life',    abiotic%half_life,    5700.)
    call g_tracer_add_param('lambda_14c',   abiotic%lambda_14c,   log(2.0) / &
                                                                  (abiotic%half_life * spery))
    call g_tracer_add_param('alkbar',       abiotic%alkbar,       2.31e-3)
    call g_tracer_add_param('atm_delta_14c',abiotic%atm_delta_14c,0.0)

    ! Rho_0 is used in the Boussinesq approximation to calculations of 
    ! pressure and pressure gradients, in units of kg m-3.
    call g_tracer_add_param('Rho_0', abiotic%Rho_0, 1035.0)

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


    call g_tracer_start_param_list(package_name)!nnz: Does this append?
    call g_tracer_add_param('ice_restart_file'   , abiotic%ice_restart_file   , 'ice_ocmip_abiotic.res.nc')
    call g_tracer_add_param('ocean_restart_file' , abiotic%ocean_restart_file , 'ocmip_abiotic.res.nc' )
    call g_tracer_add_param('IC_file'            , abiotic%IC_file            , '')
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file=abiotic%ice_restart_file, ocean_restart_file=abiotic%ocean_restart_file )

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !prog_tracers: abiotic
    !diag_tracers: none
    !

    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'dissicabio',                                    &
         longname   = 'Abiotic Dissolved Inorganic Carbon',            &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'abco2_flux',                                &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocmip_abiotic_airsea_flux.res.nc',  &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001 )          

    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'dissi14cabio',                                  &
         longname   = 'Abiotic Dissolved Inorganic Radioarbon',        &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'ab14co2_flux',                              &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocmip_abiotic_airsea_flux.res.nc',  &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001 )        

    call g_tracer_add(tracer_list,package_name,       &
         name       = 'ab_htotal',                    &
         longname   = 'abiotic H+ ion concentration', &
         units      = 'mol/kg',                       &
         prog       = .false.,                        &
         init_value = abiotic%htotal_in )

  end subroutine user_add_tracers


  ! <SUBROUTINE NAME="generic_abiotic_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_abiotic_update_from_source(tracer_list,Temp,Salt,sosga,rho_dzt,dzt,hblt_depth,&
       ilb,jlb,tau,dt,grid_dat,model_time)

    type(g_tracer_type),            pointer    :: tracer_list
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp,Salt,rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dt
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time
    real,                           intent(in) :: sosga

    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau, i, j, k , kblt, m, n, k_100, k_200 
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt
    !
    !------------------------------------------------------------------------
    ! Local Variables
    !------------------------------------------------------------------------
    !
    logical :: used
    real :: r_dt

    r_dt = 1.0 / dt

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    call g_tracer_get_values(tracer_list,'ab_htotal','field', abiotic%f_htotal,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'dissicabio'   ,'field', abiotic%f_dissicabio   ,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dissi14cabio' ,'field', abiotic%f_dissi14cabio ,isd,jsd,ntau=tau)

    !---------------------------------------------------------------------
    !Also calculate co2 fluxes csurf and alpha for the next round of exchnage
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       abiotic%htotallo(i,j) = abiotic%htotal_scale_lo * abiotic%f_htotal(i,j,k)
       abiotic%htotalhi(i,j) = abiotic%htotal_scale_hi * abiotic%f_htotal(i,j,k)
       abiotic%delta_14catm(i,j) = abiotic%atm_delta_14c ! default is 0.0
    enddo; enddo ; !} i, j
    
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      abiotic%f_alk(i,j,k)  = Salt(i,j,k) * (abiotic%alkbar/sosga)
      abiotic%f_po4(i,j,k)  = abiotic%po4_const
      abiotic%f_sio4(i,j,k) = abiotic%sio4_const
    enddo; enddo ; enddo ; !} i, j, k

    k=1
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         abiotic%f_dissicabio(:,:,k),                          &
         abiotic%f_po4(:,:,k),                          &  
         abiotic%f_sio4(:,:,k),                         &
         abiotic%f_alk(:,:,k),                          &
         abiotic%htotallo, abiotic%htotalhi,&
                                !InOut
         abiotic%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=abiotic%abco2_csurf(:,:), alpha=abiotic%abco2_alpha(:,:), &
         pCO2surf=abiotic%abpco2_csurf(:,:))

    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         abiotic%f_dissi14cabio(:,:,k),                          &
         abiotic%f_po4(:,:,k),                          &  
         abiotic%f_sio4(:,:,k),                         &
         abiotic%f_alk(:,:,k),                          &
         abiotic%htotallo, abiotic%htotalhi,&
                                !InOut
         abiotic%f_htotal(:,:,k),                       & 
                                !OUT
         co2star=abiotic%ab14co2_csurf(:,:), alpha=abiotic%ab14co2_alpha(:,:), &
         pCO2surf=abiotic%abp14co2_csurf(:,:))

    ! Update alpha based on the atmospheric 14C/12C ratio

    call data_override('OCN', 'delta_14catm', abiotic%delta_14catm(isc:iec,jsc:jec), model_time)
    do j = jsc, jec ; do i = isc, iec  !{
       abiotic%ab14co2_alpha(i,j) = abiotic%ab14co2_alpha(i,j) * &
                                    (1.0 + abiotic%delta_14catm(i,j) * 1.0e-03)
    enddo; enddo ; !} i, j

    call g_tracer_set_values(tracer_list,'ab_htotal','field',abiotic%f_htotal  ,isd,jsd,ntau=1)

    call g_tracer_set_values(tracer_list,'dissicabio','alpha',abiotic%abco2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissicabio','csurf',abiotic%abco2_csurf    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio','alpha',abiotic%ab14co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio','csurf',abiotic%ab14co2_csurf    ,isd,jsd)

    call g_tracer_get_pointer(tracer_list,'dissicabio',  'field',abiotic%p_dissicabio)
    call g_tracer_get_pointer(tracer_list,'dissi14cabio','field',abiotic%p_dissi14cabio)
    !call g_tracer_get_pointer(tracer_list,'ab_htotal','field',abiotic%p_htotal)

    ! Decay the radiocarbon
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

      abiotic%jdecay_di14c(i,j,k) = abiotic%f_dissi14cabio(i,j,k) * abiotic%lambda_14c 

      abiotic%p_dissi14cabio(i,j,k,tau) = abiotic%p_dissi14cabio(i,j,k,tau) - &
         (abiotic%jdecay_di14c(i,j,k) * (dt * grid_tmask(i,j,k)))

    enddo; enddo ; enddo !} i,j,k

    ! Send Diagnostics

   if (abiotic%id_delta_14catm .gt. 0)             &
       used = g_send_data(abiotic%id_delta_14catm, abiotic%delta_14catm(:,:),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   if (abiotic%id_ab_alk .gt. 0)             &
       used = g_send_data(abiotic%id_ab_alk, abiotic%f_alk(:,:,:),         &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ke_in=nk)

   if (abiotic%id_sfc_ab_alk .gt. 0)             &
       used = g_send_data(abiotic%id_sfc_ab_alk, abiotic%f_alk(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   if (abiotic%id_ab_po4 .gt. 0)             &
       used = g_send_data(abiotic%id_ab_po4, abiotic%f_po4(:,:,:),         &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ke_in=nk)

   if (abiotic%id_sfc_ab_po4 .gt. 0)             &
       used = g_send_data(abiotic%id_sfc_ab_po4, abiotic%f_po4(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   if (abiotic%id_ab_sio4 .gt. 0)             &
       used = g_send_data(abiotic%id_ab_sio4, abiotic%f_sio4(:,:,:),         &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

   if (abiotic%id_sfc_ab_sio4 .gt. 0)             &
       used = g_send_data(abiotic%id_sfc_ab_sio4, abiotic%f_sio4(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   if (abiotic%id_ab_pco2surf .gt. 0)              &
         used = g_send_data(abiotic%id_ab_pco2surf, abiotic%abpco2_csurf(:,:),           &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

   if (abiotic%id_ab_p14co2surf .gt. 0)              &
         used = g_send_data(abiotic%id_ab_p14co2surf, abiotic%abp14co2_csurf(:,:),           &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_sfc_dissicabio .gt. 0)  &
       used = g_send_data(abiotic%id_sfc_dissicabio, abiotic%p_dissicabio(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_sfc_dissi14cabio .gt. 0)  &
       used = g_send_data(abiotic%id_sfc_dissi14cabio, abiotic%p_dissi14cabio(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_sfc_ab_htotal .gt. 0)  &
       used = g_send_data(abiotic%id_sfc_ab_htotal, abiotic%f_htotal(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_jdecay_di14c .gt. 0)                                          &
         used = g_send_data(abiotic%id_jdecay_di14c, abiotic%jdecay_di14c(:,:,:)*rho_dzt, &
         model_time, rmask = grid_tmask(:,:,:),                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)


  end subroutine generic_abiotic_update_from_source





  ! <SUBROUTINE NAME="generic_abiotic_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_abiotic_set_boundary_values(tracer_list,SST,SSS,sosga,rho,ilb,jlb,tau)
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
  subroutine generic_abiotic_set_boundary_values(tracer_list,SST,SSS,sosga,rho,ilb,jlb,tau)
    type(g_tracer_type),            pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: SST, SSS
    real, intent(in)                           :: sosga
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,tau

    real    :: sal,ST
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: dissicabio_field, dissi14cabio_field
    real, dimension(:,:,:,:), pointer :: po4_field,sio4_field,alk_field
    real, dimension(:,:,:), ALLOCATABLE :: htotal_field
    real, dimension(:,:), ALLOCATABLE :: abco2_alpha,abco2_csurf,abco2_sc_no
    real, dimension(:,:), ALLOCATABLE :: ab14co2_alpha,ab14co2_csurf,ab14co2_sc_no
    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_set_boundary_values'

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    allocate(abco2_alpha(isd:ied, jsd:jed)); abco2_alpha=0.0
    allocate(abco2_csurf(isd:ied, jsd:jed)); abco2_csurf=0.0
    allocate(abco2_sc_no(isd:ied, jsd:jed)); abco2_sc_no=0.0
    allocate(ab14co2_alpha(isd:ied, jsd:jed)); ab14co2_alpha=0.0
    allocate(ab14co2_csurf(isd:ied, jsd:jed)); ab14co2_csurf=0.0
    allocate(ab14co2_sc_no(isd:ied, jsd:jed)); ab14co2_sc_no=0.0
    allocate(htotal_field(isd:ied,jsd:jed,nk)); htotal_field=0.0

    if (abiotic%init .OR. abiotic%force_update_fluxes) then
       !Get necessary fields
       call g_tracer_get_pointer(tracer_list,'dissicabio', 'field', dissicabio_field)
       call g_tracer_get_pointer(tracer_list,'dissi14cabio', 'field', dissi14cabio_field)
       call g_tracer_get_values(tracer_list,'ab_htotal' ,'field', htotal_field,isd,jsd,ntau=1)

       do j = jsc, jec ; do i = isc, iec  !{
          abiotic%htotallo(i,j) = abiotic%htotal_scale_lo * htotal_field(i,j,1)
          abiotic%htotalhi(i,j) = abiotic%htotal_scale_hi * htotal_field(i,j,1)
          abiotic%delta_14catm(i,j) = abiotic%atm_delta_14c ! default is 0.0
       enddo; enddo ; !} i, j

       do j = jsc, jec ; do i = isc, iec  !{
         abiotic%f_alk(i,j,1)  = SSS(i,j) * (abiotic%alkbar/sosga)
         abiotic%f_po4(i,j,1)  = abiotic%po4_const
         abiotic%f_sio4(i,j,1) = abiotic%sio4_const
       enddo; enddo ; !} i, j

       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dissicabio_field(:,:,1,tau),                          &
            abiotic%f_po4(:,:,1),                          &  
            abiotic%f_sio4(:,:,1),                         &
            abiotic%f_alk(:,:,1),                          &
            abiotic%htotallo, abiotic%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                           &
                                !OUT
            co2star=abco2_csurf(:,:), alpha=abco2_alpha(:,:),  &
            pCO2surf=abiotic%abpco2_csurf(:,:))

       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dissi14cabio_field(:,:,1,tau),                          &
            abiotic%f_po4(:,:,1),                          &  
            abiotic%f_sio4(:,:,1),                         &
            abiotic%f_alk(:,:,1),                          &
            abiotic%htotallo, abiotic%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                           &
                                !OUT
            co2star=ab14co2_csurf(:,:), alpha=ab14co2_alpha(:,:),  &
            pCO2surf=abiotic%abp14co2_csurf(:,:))

    !! Update alpha based on the atmospheric 14C/12C ratio
    !call data_override('OCN', 'delta_14catm', abiotic%delta_14catm(isc:iec,jsc:jec), model_time)
    !do j = jsc, jec ; do i = isc, iec  !{
    !   abiotic%ab14co2_alpha(i,j) = abiotic%ab14co2_alpha(i,j) * &
    !                                (1.0 + abiotic%delta_14catm(i,j) * 1.0e-03)
    !enddo; enddo ; !} i, j

       call g_tracer_set_values(tracer_list,'ab_htotal' ,'field',htotal_field,isd,jsd,ntau=1)
       call g_tracer_set_values(tracer_list,'dissicabio','alpha',abco2_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dissicabio','csurf',abco2_csurf    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dissi14cabio','alpha',ab14co2_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dissi14cabio','csurf',ab14co2_csurf    ,isd,jsd)
      
       abiotic%init = .false.

    endif

    call g_tracer_get_values(tracer_list,'dissicabio',  'alpha', abco2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dissicabio',  'csurf', abco2_csurf ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dissi14cabio',  'alpha', ab14co2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dissi14cabio',  'csurf', ab14co2_csurf ,isd,jsd)

    do j=jsc,jec ; do i=isc,iec
       sal = SSS(i,j) ; ST = SST(i,j)

       !nnz:
       !Note: In the following calculations in order to get results for co2 and o2
       !      identical with abiotic code in MOM abiotic%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !           abiotic%Rho_0 = rho(i,j,1,tau)
       !
       !      But since %Rho_0 plays the role of a unit conversion factor in this module
       !      it may be safer to keep it as a constant (1035.0) rather than the actual variable
       !      surface density rho(i,j,1,tau)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of CO2 in seawater using the
       !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
       !  7373-7382).
       !
       !     14CO2 - schmidt number is calculated same as before, as is alpha.
       !   Need to multiply alpha by frac_14catm to get the real alpha.
       !   NOTE: FOR NOW, D14C fixed at 0 permil!! Need to fix this later.
       !---------------------------------------------------------------------

       abco2_sc_no(i,j) = abiotic%a1_co2 + ST * (abiotic%a2_co2 + ST * (abiotic%a3_co2 + ST * abiotic%a4_co2)) * &
            grid_tmask(i,j,1)

       ab14co2_sc_no(i,j) = abiotic%a1_co2 + ST * (abiotic%a2_co2 + ST * (abiotic%a3_co2 + ST * abiotic%a4_co2)) * &
            grid_tmask(i,j,1)

       ! sc_no_term = sqrt(660.0 / (sc_co2 + epsln))
       !
       ! abco2_alpha(i,j) = abco2_alpha(i,j)* sc_no_term * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       ! abco2_csurf(i,j) = abco2_csurf(i,j)* sc_no_term * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       !
       ! in 'ocmip2_new' atmos_ocean_fluxes.F90 coupler formulation, the schmidt number is 
       ! carried in explicitly

       abco2_alpha(i,j) = abco2_alpha(i,j) * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       abco2_csurf(i,j) = abco2_csurf(i,j) * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)

       ab14co2_alpha(i,j) = ab14co2_alpha(i,j) * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       ab14co2_csurf(i,j) = ab14co2_csurf(i,j) * abiotic%Rho_0 !nnz: MOM has rho(i,j,1,tau)

    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. 
    !This will mark them for sending fluxes to coupler
    
    call g_tracer_set_values(tracer_list,'dissicabio',  'alpha',abco2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissicabio',  'csurf',abco2_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissicabio',  'sc_no',abco2_sc_no,isd,jsd)

    call g_tracer_set_values(tracer_list,'dissi14cabio',  'alpha',ab14co2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio',  'csurf',ab14co2_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio',  'sc_no',ab14co2_sc_no,isd,jsd)

    deallocate(abco2_alpha,abco2_csurf,abco2_sc_no,&
               ab14co2_alpha,ab14co2_csurf,ab14co2_sc_no)

  end subroutine generic_abiotic_set_boundary_values



  subroutine generic_abiotic_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_end'

  end subroutine generic_abiotic_end


end module generic_abiotic
