!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study: SF6 module
! This module contains the generic version of age Tracer and its chemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the ocean ideal age tracer.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
!  	Thiele, G., and J. L. Sarmiento, 1990: Tracer dating and ocean ventilation.
!       J. Geophys. Res., 95 (C6), 9377–9391. 
! </REFERENCE>
! <DEVELOPER_NOTES>
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_age

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, stdout
  use time_manager_mod, only : time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  


  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_age'
  character(len=fm_string_len), parameter :: package_name   = 'generic_age'

  public do_generic_age
  public generic_age_register
  public generic_age_init
  public generic_age_update_from_coupler
  public generic_age_update_from_source
  public generic_age_set_boundary_values
  public generic_age_end

  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_age = .false.

  real, parameter :: epsln=1.0e-30

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_age_params
      real :: Rho_0
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_age_params


  type(generic_age_params) :: param

contains

  subroutine generic_age_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_age_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
    
  end subroutine generic_age_register

  ! <SUBROUTINE NAME="generic_age_init">
  !  <OVERVIEW>
  !   Initialize the generic age module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the age Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_age_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_age_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_age_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    !    call user_allocate_arrays !None for age module currently

  end subroutine generic_age_init

  subroutine user_allocate_arrays
    !Allocate all the private arrays.
    !None for age module currently
  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_age_params type
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


    call g_tracer_start_param_list(package_name)!nnz: Does this append?
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip_age.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip_age.res.nc' )
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
    !prog_tracers: age
    !diag_tracers: none
    !
    !age
    call g_tracer_add(tracer_list,package_name,&
         name       = 'age',               &
         longname   = 'age Concentration', &
         units      = 'mol/kg',            &
         prog       = .true.)


  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_age_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for ages.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_age_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_age_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_age_update_from_copler'
    !
    !Nothing specific to be done for age's
    !
    return
  end subroutine generic_age_update_from_coupler

  ! <SUBROUTINE NAME="generic_age_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Sets age to zero in uppermost level and increments it by the time step elsewhere
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_age_update_from_source(tracer_list,tau,dt)
    type(g_tracer_type),            pointer    :: tracer_list
    integer,                        intent(in) :: tau
    real,                           intent(in) :: dt


    character(len=fm_string_len), parameter :: sub_name = 'generic_age_update_from_source'
    integer                                 :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,i,j,k
    real, parameter                         :: secs_in_year_r = 1.0 / (86400.0 * 365.25)
    real, dimension(:,:,:) ,pointer         :: grid_tmask
    real, dimension(:,:,:,:), pointer       :: g_age_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list,'age'    ,'field',g_age_field    )

   k = 1
    do j = jsc, jec ; do i = isc, iec
       g_age_field(i,j,k,tau) = 0.0
    enddo;  enddo !} i,j

    do k = 2, nk  ; do j = jsc, jec ; do i = isc, iec  
           g_age_field(i,j,k,tau) = g_age_field(i,j,k,tau) + dt * secs_in_year_r * grid_tmask(i,j,k)
    enddo;  enddo ;  enddo !} i,j,k

    return
  end subroutine generic_age_update_from_source

  ! <SUBROUTINE NAME="generic_age_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_age_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_age_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha,sc
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_age_field
    real, dimension(:,:), ALLOCATABLE :: g_age_alpha,g_age_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no

    character(len=fm_string_len), parameter :: sub_name = 'generic_age_set_boundary_values'


    return
  end subroutine generic_age_set_boundary_values

  ! <SUBROUTINE NAME="generic_age_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_age_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_age_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_age_end'

  end subroutine generic_age_end


end module generic_age
