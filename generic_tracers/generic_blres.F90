!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the amount of time that a water parcel 
!       at any give location has last resided in the mixed layer
!</DESCRIPTION>
!
! </DEVELOPER_NOTES>
! </INFO>
!
!----------------------------------------------------------------

module generic_blres

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  !use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, stdout
  use mpp_mod,           only: input_nml_file, mpp_error, stdlog, NOTE, WARNING, FATAL, stdout, mpp_chksum
  use fms_mod,           only: write_version_number, check_nml_error
  use time_manager_mod, only : time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_send_diag, g_tracer_get_values  


  implicit none ; private

  character(len=fm_string_len), parameter :: mod_name       = 'generic_blres'
  character(len=fm_string_len), parameter :: package_name   = 'generic_blres'

  public do_generic_blres
  public generic_blres_register
  public generic_blres_init
  public generic_blres_update_from_coupler
  public generic_blres_update_from_source
  public generic_blres_set_boundary_values
  public generic_blres_end

  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_blres = .false.

  real, parameter :: epsln=1.0e-30

  real :: reset_time=1.0

  namelist /generic_blres_nml/ reset_time
  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_blres_params
      real :: Rho_0
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_blres_params


  type(generic_blres_params) :: param

contains

  subroutine generic_blres_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

integer                                                 :: ioun
integer                                                 :: ierr
integer                                                 :: io_status
character(len=fm_string_len)                            :: name
integer                                                 :: stdoutunit,stdlogunit

    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_register'

! provide for namelist over-ride
! This needs to go before the add_tracers in order to allow the namelist 
! settings to switch tracers on and off.
!
stdoutunit=stdout();stdlogunit=stdlog()

read (input_nml_file, nml=generic_blres_nml, iostat=io_status)
ierr = check_nml_error(io_status,'generic_blres_nml')

write (stdoutunit,'(/)')
write (stdoutunit, generic_blres_nml)
write (stdlogunit, generic_blres_nml)

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)
    
    
  end subroutine generic_blres_register

  ! <SUBROUTINE NAME="generic_blres_init">
  !  <OVERVIEW>
  !   Initialize the generic mixed layer residence tracer
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the residence tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_blres_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_blres_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    !call user_allocate_arrays 

  end subroutine generic_blres_init

  subroutine user_allocate_arrays

    ! Nothing to do for mixed layer tracer
    ! call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

  end subroutine user_allocate_arrays

  subroutine user_deallocate_arrays
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  end subroutine user_deallocate_arrays
  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_blres_params type
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
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip_blres.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip_blres.res.nc' )
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
    !prog_tracers: blres
    !diag_tracers: none
    !
    !age
    call g_tracer_add(tracer_list,package_name,               &
         name          = 'blres',                             &
         longname      = 'residence time inside mixed layer', &
         units         = 'years',                             &
         init_value    = 0.0,                                 &
         prog          = .true.)

    call g_tracer_add(tracer_list,package_name,                &
         name          = 'blres_inv',                          &
         longname      = 'residence time outside mixed layer', &
         units         = 'years',                              &
         init_value    = 0.0,                                  &
         prog          = .true.)


  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_blres_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for ages.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_blres_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_blres_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_update_from_coupler'
    !
    !Nothing specific to be done for mixed layer residence tracers
    !
    return
  end subroutine generic_blres_update_from_coupler

  ! <SUBROUTINE NAME="generic_blres_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Sets age to zero in uppermost level and increments it by the time step elsewhere
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_blres_update_from_source( tracer_list,tau,dt,       &
                                               hblt_depth,dzt,ilb,jlb  )

    type(g_tracer_type),            pointer    :: tracer_list
    integer,                        intent(in) :: tau, ilb, jlb
    real,                           intent(in) :: dt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    real, dimension(ilb:,jlb:,:),   intent(in) :: dzt

    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_update_from_source'
    integer                                 :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,i,j,k
    real, parameter                         :: secs_in_year_r = 1.0 / (86400.0 * 365.25)
    real, dimension(:,:,:)  , pointer       :: grid_tmask
    real, dimension(:,:,:,:), pointer       :: p_blres_field, p_blres_inv_field

    real, dimension(:,:,:), allocatable :: zt

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    allocate( zt(isd:ied,jsd:jed,nk) ); zt=0.0

    call g_tracer_get_pointer(tracer_list, 'blres'    , 'field', p_blres_field    )
    call g_tracer_get_pointer(tracer_list, 'blres_inv', 'field', p_blres_inv_field)

    ! Compute depth of the bottom of the cell
    zt = 0.0
    do j = jsc, jec ; do i = isc, iec   !{
       zt(i,j,1) = dzt(i,j,1)
    enddo; enddo !} i,j
    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       zt(i,j,k) = zt(i,j,k-1) + dzt(i,j,k)
    enddo; enddo; enddo !} i,j

    do k = 1, nk; do j = jsc, jec ; do i = isc, iec

      if (zt(i,j,k) .le. hblt_depth(i,j) ) then
        !if ( p_blres_inv_field(i,j,k,tau) .gt. 1 ) then
        if ( p_blres_inv_field(i,j,k,tau) .gt. reset_time ) then
          ! Reset tracers
          p_blres_field    (i,j,k,tau) = 0.0
          p_blres_inv_field(i,j,k,tau) = 0.0
        else
          p_blres_field(i,j,k,tau) = p_blres_field(i,j,k,tau) + dt*secs_in_year_r*grid_tmask(i,j,k)
        endif
      else
        p_blres_inv_field(i,j,k,tau) = p_blres_inv_field(i,j,k,tau) + dt*secs_in_year_r*grid_tmask(i,j,k)
      endif

    enddo; enddo; enddo !} i,j,k

    deallocate (zt)

    return
  end subroutine generic_blres_update_from_source

  ! <SUBROUTINE NAME="generic_blres_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_blres_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_blres_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha,sc
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_blres_field
    real, dimension(:,:), ALLOCATABLE :: g_blres_alpha,g_blres_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no

    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_set_boundary_values'

    return
  end subroutine generic_blres_set_boundary_values

  ! <SUBROUTINE NAME="generic_blres_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_blres_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_blres_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_blres_end'
    !call user_deallocate_arrays

  end subroutine generic_blres_end

end module generic_blres
