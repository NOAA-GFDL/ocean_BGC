!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh 
! </CONTACT>
! 
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study: argon module
! This module contains the generic version of argon Tracer and its chemistry.
! It is designed so that both GFDL Ocean models, GOLD and MOM, can use it.
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP argon
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

module generic_argon

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

  character(len=fm_string_len), parameter :: mod_name       = 'generic_argon'
  character(len=fm_string_len), parameter :: package_name   = 'generic_argon'

  public do_generic_argon
  public generic_argon_register
  public generic_argon_init
  public generic_argon_update_from_coupler
  public generic_argon_update_from_source
  public generic_argon_set_boundary_values
  public generic_argon_end

  !The following logical for using this module is overwritten 
  ! by generic_tracer_nml namelist
  logical, save :: do_generic_argon = .false.

  real, parameter :: epsln=1.0e-30

  !
  !This type contains all the parameters and arrays used in this module.
  !
  !Note that there is no programatic reason for treating
  !the following as a type. These are the parameters used only in this module. 
  !It suffices for varables to be a declared at the top of the module. 
  !nnz: Find out about the timing overhead for using type%x rather than x

  type generic_argon_params
     real :: sA, sB, sC, sD, sE   ! Coefficients in the calculation of argon
                                  ! Schmidt number, in units of ND, degC-1, degC-2, degC-3.
     real :: A1, A2, A3           ! Coefficients in the calculation of argon
                                  ! solubility, in units of ND, K-1, log(K)^-1.
     real :: B1, B2, B3           ! More coefficients in the calculation of argon
                                  ! solubility, in units of PSU-1, PSU-1 K-1, PSU-1 K-2.
      real :: Rho_0
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file
  end type generic_argon_params


  type(generic_argon_params) :: param

contains

  subroutine generic_argon_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_argon_register'

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

    
    
  end subroutine generic_argon_register

  ! <SUBROUTINE NAME="generic_argon_init">
  !  <OVERVIEW>
  !   Initialize the generic argon module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the argon Tracers to the list of generic Tracers passed to it via utility subroutine g_tracer_add().
  !       Adds all the parameters used by this module via utility subroutine g_tracer_add_param().
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_argon_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>

  subroutine generic_argon_init(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_argon_init'

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate and initiate all the private work arrays used by this module.
    !    call user_allocate_arrays !None for argon module currently

  end subroutine generic_argon_init

  subroutine user_allocate_arrays
    !Allocate all the private arrays.
    !None for argon module currently
  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !==============================================================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_argon_params type
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
    !         for argon
    !-----------------------------------------------------------------------
    !    g_tracer_add_param(name   , variable   ,  default_value)
    call g_tracer_add_param('sA', param%sA,  2078.1)
    call g_tracer_add_param('sB', param%sB, -146.74)
    call g_tracer_add_param('sC', param%sC,  5.6403)
    call g_tracer_add_param('sD', param%sD, -0.11838)
    call g_tracer_add_param('sE', param%sE,  0.0010148)
    !-----------------------------------------------------------------------
    !     Solubility coefficients for alpha in mol/l/atm
    !      for argon
    !     after Wanninkhof (2014), L&O: Methods, 12, 351-362
    !-----------------------------------------------------------------------
    call g_tracer_add_param('A1', param%A1, -55.6578)
    call g_tracer_add_param('A2', param%A2,  82.0262)
    call g_tracer_add_param('A3', param%A3,  22.5929)
    call g_tracer_add_param('B1', param%B1,  0.036267)
    call g_tracer_add_param('B2', param%B2, -0.016241)
    call g_tracer_add_param('B3', param%B3, -0.0020114)

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
    call g_tracer_add_param('ice_restart_file'   , param%ice_restart_file   , 'ice_ocmip_argon.res.nc')
    call g_tracer_add_param('ocean_restart_file' , param%ocean_restart_file , 'ocmip_argon.res.nc' )
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
    !prog_tracers: argon
    !diag_tracers: none
    !
    !argon
    call g_tracer_add(tracer_list,package_name,&
         name       = 'argon',               &
         longname   = 'argon Concentration', &
         units      = 'mol/kg',            &
         prog       = .true.,              &
         flux_gas       = .true.,                      &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_param = (/ 9.36e-07, 9.7561e-06 /), &
         flux_gas_restart_file  = 'ocmip_argon_airsea_flux.res.nc' )


  end subroutine user_add_tracers

  ! <SUBROUTINE NAME="generic_argon_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for argons.
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_argon_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_argon_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'generic_argon_update_from_copler'
    !
    !Nothing specific to be done for argon's
    !
    return
  end subroutine generic_argon_update_from_coupler

  ! <SUBROUTINE NAME="generic_argon_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Currently an empty stub for argons.
  !  </DESCRIPTION>
  ! </SUBROUTINE>
  subroutine generic_argon_update_from_source(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    !
    !No source update for argon's currently exit in code.
    !
    return
  end subroutine generic_argon_update_from_source

  ! <SUBROUTINE NAME="generic_argon_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_argon_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau)
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
  subroutine generic_argon_set_boundary_values(tracer_list,ST,SSS,rho,ilb,jlb,taum1)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: ST, SSS 
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,taum1

    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: conv_fac,sal,ta,SST,alpha,sc
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: g_argon_field
    real, dimension(:,:), ALLOCATABLE :: g_argon_alpha,g_argon_csurf
    real, dimension(:,:), ALLOCATABLE :: sc_no

    character(len=fm_string_len), parameter :: sub_name = 'generic_argon_set_boundary_values'


    !nnz: Can we treat these as source and move block to user_update_from_source?
    !
    !=============
    !Block Starts: Calculate the boundary values
    !=============
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    call g_tracer_get_pointer(tracer_list,'argon','field',g_argon_field)

    allocate(g_argon_alpha(isd:ied, jsd:jed)); g_argon_alpha=0.0
    allocate(g_argon_csurf(isd:ied, jsd:jed)); g_argon_csurf=0.0
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
       !       Use Warner and Weiss (1985) DSR, vol 32, final result
       !       in mol/l/atm (note, atmospheric data may be in 1 part per trillion 1e-12, pptv)
       !
       !       Use Bullister and Wisegavger for CCl4
       !---------------------------------------------------------------------

       !nnz: MOM hmask=grid_tmask(i,j,1), GOLD hmask=G%hmask 
       alpha = conv_fac * grid_tmask(i,j,1) * &
            exp(param%A1 + param%A2/ta + param%A3*log(ta) +&
            sal * (param%B1 + ta * (param%B2 + ta * param%B3))&
            )

       !---------------------------------------------------------------------
       !     Calculate Schmidt numbers
       !      use coefficients given by Zheng et al (1998), JGR vol 103, C1
       !---------------------------------------------------------------------
       sc_no(i,j) = param%sA + SST * (param%sB + SST * (param%sC + SST * (param%sD + &
            SST * param%sE))) * grid_tmask(i,j,1)

       !sc_no_term = sqrt(660.0 / (sc + epsln))
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       g_argon_alpha(i,j) = alpha              

       g_argon_csurf(i,j) = g_argon_field(i,j,1,taum1) *  param%Rho_0
 
    enddo; enddo
    !=============
    !Block Ends: Calculate the boundary values
    !=============

    !
    !Set %csurf and %alpha for these tracers. This will mark them for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'argon','alpha',g_argon_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'argon','csurf',g_argon_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'argon','sc_no',sc_no,isd,jsd)

    deallocate(g_argon_alpha,g_argon_csurf,sc_no)

  end subroutine generic_argon_set_boundary_values

  ! <SUBROUTINE NAME="generic_argon_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_argon_end
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine generic_argon_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_argon_end'

  end subroutine generic_argon_end


end module generic_argon
