! DOXYGEN-compliant Documentation 
!
!-----------------------------------------------------------------------
!                   GNU General Public License
!
! This program is free software; you can redistribute it and/or modify it and
! are expected to follow the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! MOM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
! or see:   http://www.gnu.org/licenses/gpl.html
!-----------------------------------------------------------------------
!
!> \mainpage
!!
!! \brief  This module contains a generic tracer implementation of Ocean Model 
!! Intercomparison Project (OMIP) abiotic routines. The routine simulates dissolved
!! inorgantic carbon (dissicabio) and dissolved inorganic radiocarbon 
!! (dissi14cabio) tracers along with their associated flux exchange and chemistry.
!!
!! \author John.Krasting <John.Krasting@noaa.gov>
!! \author John.Dunne <John.Dunne@noaa.gov>
!!
!! The implementation of the abiotic routine follows the two-tracer approach outlined by
!! Orr et al. <a href="http://www.geosci-model-dev-discuss.net/gmd-2016-155/">
!! Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison 
!! Project (OMIP)</a>.
!!
!! The routine uses globally uniform, time-invariant concentrations of phosphate 
!! (PO<sub>4</sub>) and slicate (SiO<sub>4</sub>). Alkalintiy is based on a global
!! value that is scaled by the ratio of local sea surface salinity (sos) to the
!! global mean sea surface salinity (tosga). While the OMIP specification is that
!! tosga be computed as an annual mean, the GFDL implementation uses tosga at each
!! time step.
!!
!! In order to simulate a historical source of atmospheric 14C to the surface ocean,
!! an atmospheric delta (delta_14catm) is provided, defined by: delta = ((r14/r12)-1.)*1000.
!! This field is 0. by default, such as during a preindustrial crontrol run. A
!! time-varying field can be given to the routine by using a FMS data_override
!! entry for the 'OCN' 'delta_14catm' field.
!!
!! The tracers have independent CO<sub>2</sub> source fluxes from the atmosphere, by design, to
!! provide flexibility in the future.  The CO<sub>2</sub> fluxes for OMIP, however, should all
!! be identical.  The preindustrial value for OMIP is 284.65 ppmv.
!!
!! \section field_table_entries  Field Table Entries
!! A sample field_table is provided below.  In this example, the abiotic htotal field
!! is initialized to the same field from COBALT. The dissicabio and dissicabio fields
!! are both initialized from the preindustrial GLODAP fields.
!!
!! It may be preferrable to initialize the dissicabio field as the preindustrail GLODAP
!! DIC minus AOU*(C/O<sub>2</sub>) where C/O<sub>2</sub> = 106/160.  There is also a GLODAP 14C field that
!! can be used to initialize the dissi14cabio field.
!!
!! <pre>
!! "namelists","ocean_mod","generic_abiotic"
!! init = f
!! # min = 10^-8.4=3.981E-09,  max = 10^-7.9 = 1.259E-08 from input dataset
!! ab_htotal_src_file = INPUT/init_ocean_cobalt.res.nc
!! ab_htotal_src_var_name = htotal
!! ab_htotal_src_var_unit = none
!! ab_htotal_dest_var_name =  ab_htotal
!! ab_htotal_dest_var_unit = mol kg-1
!! ab_htotal_src_var_record = 1
!! ab_htotal_src_var_gridspec = NONE
!! ab_htotal_valid_min = 3.981E-09
!! ab_htotal_valid_max = 1.259E-08
!! #
!! ## divide by 1e6 to convert from umol/kg to mol/kg
!! dissicabio_src_file = INPUT/Preind_DIC.nc
!! dissicabio_src_var_name = Preind_DIC
!! dissicabio_src_var_unit = micromoles_per_kg
!! dissicabio_dest_var_name =  dissicabio
!! dissicabio_dest_var_unit = mol kg-1
!! dissicabio_src_var_record = 1
!! dissicabio_src_var_gridspec = NONE
!! dissicabio_valid_min = 0.0
!! #
!! ## divide by 1e6 to convert from umol/kg to mol/kg
!! dissi14cabio_src_file = INPUT/Preind_DIC.nc
!! dissi14cabio_src_var_name = Preind_DIC
!! dissi14cabio_src_var_unit = micromoles_per_kg
!! dissi14cabio_dest_var_name =  dissi14cabio
!! dissi14cabio_dest_var_unit = mol kg-1
!! dissi14cabio_src_var_record = 1
!! dissi14cabio_src_var_gridspec = NONE
!! dissi14cabio_valid_min = 0.0
!! /
!! </pre>
!!
!! \section data_table_entries  Data Table Entries
!! A sample data table is provided below.  A version time varying delta 14C file can be
!! found on Gaea at: /lustre/f1/unswept/John.Krasting/cwg/input/abiotic/atm_delta_13C_14C.nc 
!!
!! <pre>
!! "OCN", "delta_14catm", "D14C", "./INPUT/atm_delta_13C_14C.nc",  "bilinear", 1.0
!! "ATM", "abco2_flux_pcair_atm", "" , "", "none" ,284.262e-6
!! "ATM", "ab14co2_flux_pcair_atm", "", "", "none" ,284.262e-6
!! </pre>
!!
!! \section diag_table_entries Diag Table Entries
!! A sample diag table for the abiotic diagnostics is provided below.  Since some of the
!! nutrient fields are uniform in space and time, some fields may be commented out to save
!! space.  The full list of diagnostics are provided here to aid with debugging, if needed
!!
!! <pre>
!! "ocean_abiotic",           1, "months", 1, "days", "time",
!! "ocean_abiotic_z"          1, "months", 1, "days", "time",
!! "ocean_abiotic_annual",   12, "months", 1, "days", "time",
!! "ocean_abiotic_annual_z"  12, "months", 1, "days", "time",
!!
!! #-- CMIP Tracer Fields
!! "generic_abiotic",   "dissicabio",           "dissicabio",           "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "dissi14cabio",         "dissi14cabio",         "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "phabio",               "phabio",               "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic_z", "dissicabio",           "dissicabio",           "ocean_abiotic_z","all",.true.,"none", 2
!! "generic_abiotic_z", "dissi14cabio",         "dissi14cabio",         "ocean_abiotic_z","all",.true.,"none", 2
!! "generic_abiotic_z", "phabio",               "phabio",               "ocean_abiotic_z","all",.true.,"none", 2
!! "generic_abiotic",   "dissicabio",           "dissicabio",           "ocean_abiotic_annual","all",.true.,"none", 2
!! "generic_abiotic",   "dissi14cabio",         "dissi14cabio",         "ocean_abiotic_annual","all",.true.,"none", 2
!! "generic_abiotic",   "phabio",               "phabio",               "ocean_abiotic_annual","all",.true.,"none", 2
!! "generic_abiotic_z", "dissicabio",           "dissicabio",           "ocean_abiotic_annual_z","all",.true.,"none", 2
!! "generic_abiotic_z", "dissi14cabio",         "dissi14cabio",         "ocean_abiotic_annual_z","all",.true.,"none", 2
!! "generic_abiotic_z", "phabio",               "phabio",               "ocean_abiotic_annual_z","all",.true.,"none", 2
!!
!! #-- CMIP Surface Fields
!! "generic_abiotic",   "dissicabioos",         "dissicabioos",         "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissi14cabioos",       "dissi14cabioos",       "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dpco2abio",            "dpco2abio",            "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "fgco2abio",            "fgco2abio",            "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "fg14co2abio",          "fg14co2abio",          "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "phabioos",             "phabioos",             "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "spco2abio",            "spco2abio",            "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "fgco2abio",            "fgco2abio",            "ocean_abiotic_annual","all",.true.,"none",2
!! "generic_abiotic",   "fg14co2abio",          "fg14co2abio",          "ocean_abiotic_annual","all",.true.,"none",2
!!
!! #-- Non-CMIP Tracer Fields
!! "generic_abiotic",   "ab_alk",               "ab_alk",               "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "ab_htotal",            "ab_htotal",            "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "ab_po4",               "ab_po4",               "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "ab_sio4",              "ab_sio4",              "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic",   "jdecay_di14c",         "jdecay_di14c",         "ocean_abiotic","all",.true.,"none", 2
!! "generic_abiotic_z", "ab_alk",               "ab_alk",               "ocean_abiotic_z","all",.true.,"none",2
!! "generic_abiotic_z", "ab_htotal",            "ab_htotal",            "ocean_abiotic_z","all",.true.,"none",2
!! "generic_abiotic_z", "ab_po4",               "ab_po4",               "ocean_abiotic_z","all",.true.,"none",2
!! "generic_abiotic_z", "ab_sio4",              "ab_sio4",              "ocean_abiotic_z","all",.true.,"none",2

!! #-- Non-CMIP Surface Fields
!! "generic_abiotic",   "sfc_ab_alk",           "sfc_ab_alk",           "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "sfc_ab_htotal",        "sfc_ab_htotal",        "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "sfc_ab_po4",           "sfc_ab_po4",           "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "sfc_ab_sio4",          "sfc_ab_sio4",          "ocean_abiotic","all",.true.,"none",2
!! 
!! #-- Non-CMIP Gas exchange
!! "generic_abiotic",   "ab_pco2surf",          "ab_pco2surf",          "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "ab_p14co2surf",        "ab_p14co2surf",        "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "delta_14catm",         "delta_14catm",         "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissicabio_alpha",     "dissicabio_alpha",     "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissi14cabio_alpha",   "dissi14cabio_alpha",   "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissicabio_csurf",     "dissicabio_csurf",     "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissi14cabio_csurf",   "dissi14cabio_csurf",   "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissicabio_sc_no",     "dissicabio_sc_no",     "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissi14cabio_sc_no",   "dissi14cabio_sc_no",   "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissicabio_stf_gas",   "dissicabio_stf_gas",   "ocean_abiotic","all",.true.,"none",2
!! "generic_abiotic",   "dissi14cabio_stf_gas", "dissi14cabio_stf_gas", "ocean_abiotic","all",.true.,"none",2
!!</pre>
!!

module generic_abiotic

  use constants_mod,     only: WTMCO2
  use data_override_mod, only: data_override
  use field_manager_mod, only: fm_string_len
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist
  use fms_mod,           only: check_nml_error
  use mpp_mod,           only: input_nml_file, mpp_error, stdlog, NOTE, WARNING, FATAL, stdout, mpp_chksum
  use time_manager_mod,  only: time_type

  use g_tracer_utils,    only: g_diag_type,g_send_data
  use g_tracer_utils,    only: g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils,    only: g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils,    only: g_tracer_set_values,g_tracer_get_pointer,g_tracer_get_common
  use g_tracer_utils,    only: g_tracer_get_values  
  use g_tracer_utils,    only: register_diag_field=>g_register_diag_field
  use g_tracer_utils,    only: is_root_pe

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
  public as_param_abiotic

  !The following variables for using this module
  ! are overwritten by generic_tracer_nml namelist
  logical, save :: do_generic_abiotic = .false.
  character(len=10), save :: as_param_abiotic   = 'gfdl_cmip6'

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30
  real, parameter :: missing_value1=-1.0e+10
  real, parameter :: missing_value_diag=-1.0e+10

  !Namelist Options
  character(len=10) ::  co2_calc = 'ocmip2'  ! other option is 'mocsy'

namelist /generic_abiotic_nml/ co2_calc

  !This type contains all the parameters and arrays used in this module.
  type generic_abiotic_params
     logical :: init              !< Logical to initialize tracers
     logical :: force_update_fluxes  !< If OCMIP2 tracers fluxes should be updated 
                                     !! every coupling time step when update_from_source 
                                     !! is not called every coupling time step as is the 
                                     !! case with MOM6  THERMO_SPANS_COUPLING option
     real :: atm_delta_14c        !< Ratio of atmospheric 14C to 12C
     real :: half_life            !< Decay time scale of 14C (years)
     real :: lambda_14c           !< Radioactive decay constant for 14C (s-1)
     real :: htotal_in            !< Initial "first guess" for H+ concentration (mol/kg)
     real :: htotal_scale_hi      !< Initial upper limit of htotal range 
     real :: htotal_scale_lo      !< Initial lower limit of htotal range 
     real :: htotal14c_in         !< Initial "first guess" for H+ concentration used for 14C (mol/kg)
     real :: htotal14c_scale_hi   !< Initial upper limit of htotal range used for 14C
     real :: htotal14c_scale_lo   !< Initial lower limit of htotal range used for 14C
     real :: alkbar               !< Mean global alkalinity (eq/kg)
     real :: sio4_const           !< Silicate (SiO4) concentration (mol/kg)
     real :: po4_const            !< Phosphate (PO4) concentration (mol/kg)
     real :: sA_co2               !< Schmidt number Coeff.
     real :: sB_co2               !< Schmidt number Coeff.
     real :: sC_co2               !< Schmidt number Coeff.
     real :: sD_co2               !< Schmidt number Coeff.
     real :: sE_co2               !< Schmidt number Coeff.
     real :: Rho_0                !< Reference density (kg/m^3)

     ! Restart file names
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     ! Diagnostic Output IDs
     integer :: id_dissicabioos=-1, id_dissi14cabioos=-1
     integer :: id_sfc_ab_htotal=-1
     integer :: id_ab_alk=-1, id_ab_po4=-1, id_ab_sio4=-1
     integer :: id_sfc_ab_alk=-1, id_sfc_ab_po4=-1, id_sfc_ab_sio4=-1
     integer :: id_ab_pco2surf=-1, id_ab_p14co2surf=-1
     integer :: id_jdecay_di14c=-1, id_delta_14catm=-1
     integer :: id_phabio=-1, id_phabioos=-1
     integer :: id_fgco2abio=-1, id_fg14co2abio=-1
     integer :: id_spco2abio=-1, id_dpco2abio=-1

     ! 2-dimensional fields
     real, dimension(:,:), ALLOCATABLE ::   &
          htotalhi, htotallo,               &!< Upper and lower limits of htotal range 
          htotal14chi, htotal14clo,         &!< Upper and lower limits of htotal range 
          abco2_csurf,abco2_alpha,          &!< Abiotic DIC csurf & alpha
          ab14co2_csurf,ab14co2_alpha,      &!< Abiotic DI14C csurf & alpha
          abpco2_csurf,abp14co2_csurf,      &!< Oceanic pCO2 (ppmv)
          delta_14catm                       !< Atm. Delta of 14C 
 
     ! 2-dimensional pointers
     real, dimension(:,:), pointer :: &
          stf_gas_dissicabio,               &!< Surface partial pressure of DIC
          stf_gas_dissi14cabio,             &!< Sufface partial pressure of DI14C
          deltap_dissicabio                  !< Delta Abiotic pCO2

     ! 3-dimensional fields
     real, dimension(:,:,:), ALLOCATABLE :: &
          f_dissicabio, f_dissi14cabio,     &!< Abiotic DIC & DI14C fields
          f_alk, f_po4, f_sio4,             &!< Abiotic alkalinity, phosphate, silicate
          f_htotal,                         &!< Abiotic H+ concentration
          f_htotal14c,                      &!< Abiotic H+ concentration for 14C
          jdecay_di14c,                     &!< DI14C decay rate
          zt                                 !< Cell depth passed from ocean model

     ! 4-dimensional pointers
     real, dimension(:,:,:,:), pointer :: &
          p_dissicabio, p_dissi14cabio, p_htotal, p_htotal14c !< Pointers to tracer fields

  end type generic_abiotic_params

  ! An auxiliary type for storing varible names
  type, public :: vardesc
     character(len=fm_string_len) :: name     !< The variable name in a NetCDF file.
     character(len=fm_string_len) :: longname !< The long name of that variable.
     character(len=1)  :: hor_grid            !< The horiz. grid:  u, v, h, q, or 1.
     character(len=1)  :: z_grid              !< The vert. grid:  L, i, or 1.
     character(len=1)  :: t_grid              !< The time description: s, a, m, or 1.
     character(len=fm_string_len) :: units    !< The dimensions of the variable.
     character(len=1)  :: mem_size            !< The size in memory: d or f.
  end type vardesc

  type(CO2_dope_vector)        :: CO2_dope_vec
  type(generic_abiotic_params) :: abiotic
 
contains

  subroutine generic_abiotic_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    integer :: ioun, ierr, io_status, stdoutunit, stdlogunit

    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_register'
    character(len=256), parameter :: error_header =                                &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: warn_header =                               &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter   :: note_header =                               &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    ! provide for namelist over-ride
    ! This needs to go before the add_tracers in order to allow the namelist 
    ! settings to switch tracers on and off.

    stdoutunit=stdout();stdlogunit=stdlog()

    read (input_nml_file, nml=generic_abiotic_nml, iostat=io_status)
    ierr = check_nml_error(io_status,'generic_abiotic_nml')
    
    write (stdoutunit,'(/)')
    write (stdoutunit, generic_abiotic_nml)
    write (stdlogunit, generic_abiotic_nml)

    if (trim(co2_calc) == 'ocmip2') then
      write (stdoutunit,*) trim(note_header), 'Using FMS OCMIP2 CO2 routine'
    else if (trim(co2_calc) == 'mocsy') then
      write (stdoutunit,*) trim(note_header), 'Using Mocsy CO2 routine'
    else
      call mpp_error(FATAL,"Unknown co2_calc option specified in generic_abiotic_nml")
    endif

    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_abiotic_register


!> \brief   Initialize the generic abiotic module
!!
!! This subroutine adds the dissicabio, dissi14cabio, and ab_htotal tracers to the list of 
!! generic tracers passed to it via utility subroutine g_tracer_add().  Adds all the parameters 
!! used by this module via utility subroutine g_tracer_add_param(). Allocates all work arrays used 
!! in the module. 
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


!> \brief Register diagnostic fields to be used in this module. 
!!
!! Note that the tracer fields are automatically registered in user_add_tracers. User adds only 
!! diagnostics for fields that are not a member of g_tracer_type
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

    vardesc_temp = vardesc("dissicabioos","Surface Abiotic Dissolved Inorganic Carbon Concentration",&
                           'h','1','s','mol m-3','f')
    abiotic%id_dissicabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name = "mole_concentration_of_dissolved_inorganic_carbon_abiotic_analogue_in_sea_water")

    vardesc_temp = vardesc("dissi14cabioos","Surface Abiotic Dissolved Inorganic 14Carbon Concentration",&
                           'h','1','s','mol m-3','f')
    abiotic%id_dissi14cabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name = "mole_concentration_of_dissolved_inorganic_carbon14_in_sea_water")

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

    vardesc_temp = vardesc("phabio","Abiotic pH",'h','L','s','1.0','f')
    abiotic%id_phabio = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         standard_name="sea_water_ph_abiotic_analogue_reported_on_total_scale")

    vardesc_temp = vardesc("phabioos","Surface Abiotic pH",'h','L','s','1.0','f')
    abiotic%id_phabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1,&
         standard_name="sea_water_ph_abiotic_analogue_reported_on_total_scale")

    vardesc_temp = vardesc("fgco2abio","Surface Downward Abiotic CO2 Flux",'h','1','s','kg m-2 s-1','f')
    abiotic%id_fgco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_downward_mass_flux_of_carbon_dioxide_abiotic_analogue_expressed_as_carbon")

    vardesc_temp = vardesc("fg14co2abio","Surface Downward Abiotic 14CO2 Flux",'h','1','s','kg m-2 s-1','f')
    abiotic%id_fg14co2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_downward_mass_flux_of_carbon14_dioxide_abiotic_analogue_expressed_as_carbon")

    vardesc_temp = vardesc("spco2abio","Abiotic Surface Aqueous Partial Pressure of CO2",'h','1','s','Pa','f')
    abiotic%id_spco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_partial_pressure_of_carbon_dioxide_abiotic_analogue_in_sea_water")

    vardesc_temp = vardesc("dpco2abio","Abiotic Delta PCO2",'h','1','s','Pa','f')
    abiotic%id_dpco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         standard_name="surface_carbon_dioxide_abiotic_analogue_partial_pressure_difference_between_sea_water_and_air")


  end subroutine generic_abiotic_register_diag



!> \brief Internal subroutine to allocate user arrays
!!
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
    allocate(abiotic%htotal14clo(isd:ied,jsd:jed))
    allocate(abiotic%htotal14chi(isd:ied,jsd:jed))
    allocate(abiotic%delta_14catm(isd:ied,jsd:jed))
    allocate(abiotic%zt(isd:ied,jsd:jed,1:nk))             ; abiotic%zt=0.0
    allocate(abiotic%f_alk(isd:ied,jsd:jed,1:nk))          ; abiotic%f_alk=0.0
    allocate(abiotic%f_htotal(isd:ied,jsd:jed,1:nk))       ; abiotic%f_htotal =0.0
    allocate(abiotic%f_htotal14c(isd:ied,jsd:jed,1:nk))    ; abiotic%f_htotal14c =0.0
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


!> \brief Internal subroutine to allocate user parameters
!!
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
    ! CO2 Schmidt number coefficients
    !-----------------------------------------------------------------------
    if ((trim(as_param_abiotic) == 'W92') .or. (trim(as_param_abiotic) == 'gfdl_cmip6')) then
        call g_tracer_add_param('sA_co2', abiotic%sA_co2, 2068.9)
        call g_tracer_add_param('sB_co2', abiotic%sB_co2, -118.63)
        call g_tracer_add_param('sC_co2', abiotic%sC_co2, 2.9311)
        call g_tracer_add_param('sD_co2', abiotic%sD_co2, -0.027)
        call g_tracer_add_param('sE_co2', abiotic%sE_co2, 0.0)      ! Not used in W92
        if (is_root_pe()) call mpp_error(NOTE,'generic_abiotic: Using Schmidt number coefficients for W92')
    else if (trim(as_param_abiotic) == 'W14') then
        call g_tracer_add_param('sA_co2', abiotic%sA_co2,  2116.8)
        call g_tracer_add_param('sB_co2', abiotic%sB_co2, -136.25)
        call g_tracer_add_param('sC_co2', abiotic%sC_co2,  4.7353)
        call g_tracer_add_param('sD_co2', abiotic%sD_co2, -0.092307)
        call g_tracer_add_param('sE_co2', abiotic%sE_co2, -0.0007555)
        if (is_root_pe()) call mpp_error(NOTE,'generic_abiotic: Using Schmidt number coefficients for W14')
    else
        call mpp_error(FATAL,'generic_abiotic: unable to set Schmidt number coefficients for CO2.')
    endif

    !-----------------------------------------------------------------------
    ! H+ Concentration Parameters
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_in',          abiotic%htotal_in, 1.0e-8)
    call g_tracer_add_param('htotal_scale_lo',    abiotic%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi',    abiotic%htotal_scale_hi, 100.)
    call g_tracer_add_param('htotal14c_in',       abiotic%htotal14c_in, 1.0e-8)
    call g_tracer_add_param('htotal14c_scale_lo', abiotic%htotal14c_scale_lo, 0.01)
    call g_tracer_add_param('htotal14c_scale_hi', abiotic%htotal14c_scale_hi, 100.)
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

!> \brief Internal subroutine to add tracers to tracer_list
!!
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    real :: as_coeff_abiotic

    if ((trim(as_param_abiotic) == 'W92') .or. (trim(as_param_abiotic) == 'gfdl_cmip6')) then
        ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
        ! Value is 0.337 cm/hr in units of m/s.
        as_coeff_abiotic = 9.36e-7
    else if (trim(as_param_abiotic) == 'W14') then
        ! Value is 0.251 cm/hr in units of m/s
        as_coeff_abiotic = 6.972e-7
    else
        call mpp_error(FATAL,'Unable to set wind speed coefficient coefficients for as_param '//trim(as_param_abiotic))
    endif

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
         longname   = 'Abiotic Dissolved Inorganic Carbon Concentration',&
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'abco2_flux',                                &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ as_coeff_abiotic, 9.7561e-06 /),          &
         flux_gas_restart_file  = 'ocmip_abiotic_airsea_flux.res.nc',  &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001,                                           &
         standard_name = "mole_concentration_of_dissolved_inorganic_carbon_abiotic_analogue_in_sea_water", &
         diag_field_units = 'mol m-3',                                 &
         diag_field_scaling_factor = 1035.0)  !rho = 1035.0 kg/m3, converts mol/kg to mol/m3

    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'dissi14cabio',                                  &
         longname   = 'Abiotic Dissolved Inorganic 14Carbon Concentration',&
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'ab14co2_flux',                              &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ as_coeff_abiotic, 9.7561e-06 /),                  &
         flux_gas_restart_file  = 'ocmip_abiotic_airsea_flux.res.nc',  &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001,                                           &
         standard_name = "mole_concentration_of_dissolved_inorganic_carbon14_in_sea_water", &
         diag_field_units = 'mol m-3',                                &
         diag_field_scaling_factor = 1035.0) ! rho = 1035.0 kg/m3, converts mol/kg to mol/m3

    call g_tracer_add(tracer_list,package_name,       &
         name       = 'ab_htotal',                    &
         longname   = 'abiotic H+ ion concentration', &
         units      = 'mol/kg',                       &
         prog       = .false.,                        &
         init_value = abiotic%htotal_in )

  end subroutine user_add_tracers


!> \brief Update tracer concentration fields from the ocean model source
!!
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

    call g_tracer_get_values(tracer_list,'ab_htotal',    'field', abiotic%f_htotal       ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'dissicabio'   ,'field', abiotic%f_dissicabio   ,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dissi14cabio' ,'field', abiotic%f_dissi14cabio ,isd,jsd,ntau=tau)

    abiotic%zt = 0.0
    do j = jsc, jec ; do i = isc, iec   !{
       abiotic%zt(i,j,1) = dzt(i,j,1)
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       abiotic%zt(i,j,k) = abiotic%zt(i,j,k-1) + dzt(i,j,k)
    enddo; enddo ; enddo !} i,j,k


    !---------------------------------------------------------------------
    !Also calculate co2 fluxes csurf and alpha for the next round of exchnage
    !---------------------------------------------------------------------
   
    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       abiotic%htotallo(i,j) = abiotic%htotal_scale_lo * abiotic%f_htotal(i,j,k)
       abiotic%htotalhi(i,j) = abiotic%htotal_scale_hi * abiotic%f_htotal(i,j,k)
       abiotic%htotal14clo(i,j) = abiotic%htotal14c_scale_lo * abiotic%f_htotal14c(i,j,k)
       abiotic%htotal14chi(i,j) = abiotic%htotal14c_scale_hi * abiotic%f_htotal14c(i,j,k)
       abiotic%delta_14catm(i,j) = abiotic%atm_delta_14c ! default is 0.0
    enddo; enddo ; !} i, j
    
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      abiotic%f_alk(i,j,k)  = Salt(i,j,k) * (abiotic%alkbar/sosga)
      abiotic%f_po4(i,j,k)  = abiotic%po4_const
      abiotic%f_sio4(i,j,k) = abiotic%sio4_const
    enddo; enddo ; enddo ; !} i, j, k

    k=1
    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                      &
         abiotic%f_dissicabio(:,:,k),                   &
         abiotic%f_po4(:,:,k),                          &  
         abiotic%f_sio4(:,:,k),                         &
         abiotic%f_alk(:,:,k),                          &
         abiotic%htotallo, abiotic%htotalhi,&
                                !InOut
         abiotic%f_htotal(:,:,k),                       & 
                                !Optional In
         co2_calc=trim(co2_calc),                       & 
         zt=abiotic%zt(:,:,k),                           & 
                                !OUT
         co2star=abiotic%abco2_csurf(:,:), alpha=abiotic%abco2_alpha(:,:), &
         pCO2surf=abiotic%abpco2_csurf(:,:))

    ! Update csurf and alpha based on the atmospheric 14C/12C ratio

    ! The 14C surface concentration (ab14co2_csurf, 14CO2*) is proportional to the 14C surface concentration
    ! times the ratio of the *ocean* surface carbon concentration (i.e. DI14C/DIC). The 14C solubility (ab14co2_alpha) 
    ! is proportional to abco2_alpha times the 14C to C ratio of the *atmos* carbon concentration as the 
    ! ocean carbon solubility mainly depends on temperature and the atmos partial pressure of the gas.

    call data_override('OCN', 'delta_14catm', abiotic%delta_14catm(isc:iec,jsc:jec), model_time)
    do j = jsc, jec ; do i = isc, iec  !{
       abiotic%ab14co2_csurf(i,j) = abiotic%abco2_csurf(i,j) * &
                                   (abiotic%f_dissi14cabio(i,j,1)/(abiotic%f_dissicabio(i,j,1) + epsln))
       abiotic%ab14co2_alpha(i,j) = abiotic%abco2_alpha(i,j) * &
                                    (1.0 + abiotic%delta_14catm(i,j) * 1.0e-03)
    enddo; enddo ; !} i, j

    call g_tracer_set_values(tracer_list,'ab_htotal',   'field',abiotic%f_htotal   ,isd,jsd,ntau=1)

    call g_tracer_set_values(tracer_list,'dissicabio','alpha',abiotic%abco2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissicabio','csurf',abiotic%abco2_csurf    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio','alpha',abiotic%ab14co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dissi14cabio','csurf',abiotic%ab14co2_csurf    ,isd,jsd)

    call g_tracer_get_pointer(tracer_list,'dissicabio',  'field',abiotic%p_dissicabio)
    call g_tracer_get_pointer(tracer_list,'dissi14cabio','field',abiotic%p_dissi14cabio)
    call g_tracer_get_pointer(tracer_list,'ab_htotal','field',abiotic%p_htotal)

    call g_tracer_get_pointer(tracer_list,'dissicabio','stf_gas',abiotic%stf_gas_dissicabio)
    call g_tracer_get_pointer(tracer_list,'dissi14cabio','stf_gas',abiotic%stf_gas_dissi14cabio)

    call g_tracer_get_pointer(tracer_list,'dissicabio','deltap',abiotic%deltap_dissicabio)

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

    if (abiotic%id_dissicabioos .gt. 0)  &
       used = g_send_data(abiotic%id_dissicabioos, (abiotic%p_dissicabio(:,:,1,tau) * abiotic%Rho_0),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_dissi14cabioos .gt. 0)  &
       used = g_send_data(abiotic%id_dissi14cabioos, (abiotic%p_dissi14cabio(:,:,1,tau) * abiotic%Rho_0),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_sfc_ab_htotal .gt. 0)  &
       used = g_send_data(abiotic%id_sfc_ab_htotal, abiotic%f_htotal(:,:,1),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_jdecay_di14c .gt. 0)                                          &
         used = g_send_data(abiotic%id_jdecay_di14c, abiotic%jdecay_di14c(:,:,:), &
         model_time, rmask = grid_tmask(:,:,:),                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (abiotic%id_phabio .gt. 0) &
         used = g_send_data(abiotic%id_phabio, log10(abiotic%f_htotal(:,:,:))*(-1.),      &
         model_time, rmask = grid_tmask(:,:,:),                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (abiotic%id_phabioos .gt. 0) &
         used = g_send_data(abiotic%id_phabioos, log10(abiotic%f_htotal(:,:,1))*(-1.),    &
         model_time, rmask = grid_tmask(:,:,1),                                         & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    if (abiotic%id_fgco2abio .gt. 0)  &
        used = g_send_data(abiotic%id_fgco2abio, abiotic%stf_gas_dissicabio * 12e-3, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (abiotic%id_fg14co2abio .gt. 0)            &
        used = g_send_data(abiotic%id_fg14co2abio, abiotic%stf_gas_dissi14cabio * 12e-3, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (abiotic%id_spco2abio .gt. 0)            &
        used = g_send_data(abiotic%id_spco2abio,  abiotic%abpco2_csurf(:,:) * 0.1013,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (abiotic%id_dpco2abio .gt. 0)            &
        used = g_send_data(abiotic%id_dpco2abio,  abiotic%deltap_dissicabio * 0.1013,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


 end subroutine generic_abiotic_update_from_source



!> \brief Calculate and set coupler values at the surface
!!
!! This subroutine takes the pointer to the head of generic tracer list, the lower bounds of x and y 
!! extents of input arrays on data domain, sea surface temperature, sea surface salinity, global average
!! sea surface salinity, ocean denisty, and the time step index of %field as input.
  subroutine generic_abiotic_set_boundary_values(tracer_list,SST,SSS,sosga,rho,ilb,jlb,tau,model_time,dzt)
    type(g_tracer_type),            pointer    :: tracer_list
    real, dimension(ilb:,jlb:),     intent(in) :: SST, SSS
    real, intent(in)                           :: sosga
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,tau
    type(time_type),                intent(in) :: model_time
    real, dimension(ilb:,jlb:,:),   optional, intent(in) :: dzt

    real    :: sal,ST
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: dissicabio_field, dissi14cabio_field
    real, dimension(:,:,:,:), pointer :: po4_field,sio4_field,alk_field
    real, dimension(:,:,:), ALLOCATABLE :: htotal_field, htotal14c_field
    real, dimension(:,:), ALLOCATABLE :: abco2_alpha,abco2_csurf,abco2_sc_no
    real, dimension(:,:), ALLOCATABLE :: ab14co2_alpha,ab14co2_csurf,ab14co2_sc_no, delta_14catm
    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_set_boundary_values'

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    allocate(abco2_alpha(isd:ied, jsd:jed)); abco2_alpha=0.0
    allocate(abco2_csurf(isd:ied, jsd:jed)); abco2_csurf=0.0
    allocate(abco2_sc_no(isd:ied, jsd:jed)); abco2_sc_no=0.0
    allocate(ab14co2_alpha(isd:ied, jsd:jed)); ab14co2_alpha=0.0
    allocate(ab14co2_csurf(isd:ied, jsd:jed)); ab14co2_csurf=0.0
    allocate(ab14co2_sc_no(isd:ied, jsd:jed)); ab14co2_sc_no=0.0
    allocate(delta_14catm(isd:ied, jsd:jed)); delta_14catm=0.0
    allocate(htotal_field(isd:ied,jsd:jed,nk)); htotal_field=0.0
    allocate(htotal14c_field(isd:ied,jsd:jed,nk)); htotal14c_field=0.0

    if (abiotic%init .OR. abiotic%force_update_fluxes) then
       !Get necessary fields
       call g_tracer_get_pointer(tracer_list,'dissicabio', 'field', dissicabio_field)
       call g_tracer_get_pointer(tracer_list,'dissi14cabio', 'field', dissi14cabio_field)

       call g_tracer_get_values(tracer_list, 'ab_htotal', 'field', htotal_field,isd,jsd,ntau=1)

       do j = jsc, jec ; do i = isc, iec  !{
          abiotic%htotallo(i,j) = abiotic%htotal_scale_lo * htotal_field(i,j,1)
          abiotic%htotalhi(i,j) = abiotic%htotal_scale_hi * htotal_field(i,j,1)
          abiotic%htotal14clo(i,j) = abiotic%htotal14c_scale_lo * htotal14c_field(i,j,1)
          abiotic%htotal14chi(i,j) = abiotic%htotal14c_scale_hi * htotal14c_field(i,j,1)
          abiotic%delta_14catm(i,j) = abiotic%atm_delta_14c ! default is 0.0
       enddo; enddo ; !} i, j

       if(present(dzt)) then
         do j = jsc, jec ; do i = isc, iec  !{
          abiotic%zt(i,j,1) = dzt(i,j,1)
         enddo; enddo ; !} i, j
       elseif (trim(co2_calc) == 'mocsy') then
         call mpp_error(FATAL,"mocsy method of co2_calc needs dzt to be passed to the FMS_ocmip2_co2calc subroutine.")
       endif

       do j = jsc, jec ; do i = isc, iec  !{
         abiotic%f_alk(i,j,1)  = SSS(i,j) * (abiotic%alkbar/sosga)
         abiotic%f_po4(i,j,1)  = abiotic%po4_const
         abiotic%f_sio4(i,j,1) = abiotic%sio4_const
       enddo; enddo ; !} i, j

       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                                &
            dissicabio_field(:,:,1,tau),                       &
            abiotic%f_po4(:,:,1),                              &  
            abiotic%f_sio4(:,:,1),                             &
            abiotic%f_alk(:,:,1),                              &
            abiotic%htotallo, abiotic%htotalhi,                &
                                !InOut
            htotal_field(:,:,1),                               &
                                !Optional In
            co2_calc=trim(co2_calc),                           & 
            zt=abiotic%zt(:,:,1),                               & 
                                !OUT
            co2star=abco2_csurf(:,:), alpha=abco2_alpha(:,:),  &
            pCO2surf=abiotic%abpco2_csurf(:,:))

       ! Update csurf and alpha based on the atmospheric 14C/12C ratio

       ! The 14C surface concentration (ab14co2_csurf, 14CO2*) is proportional to the 14C surface concentration
       ! times the ratio of the *ocean* surface carbon concentration (i.e. DI14C/DIC). The 14C solubility (ab14co2_alpha) 
       ! is proportional to abco2_alpha times the 14C to C ratio of the *atmos* carbon concentration as the 
       ! ocean carbon solubility mainly depends on temperature and the atmos partial pressure of the gas.

       call data_override('OCN', 'delta_14catm', delta_14catm(isc:iec,jsc:jec), model_time)
       do j = jsc, jec ; do i = isc, iec  !{
          ab14co2_csurf(i,j) = abco2_csurf(i,j) * &
                              (dissi14cabio_field(i,j,1,tau)/(dissicabio_field(i,j,1,tau) + epsln))
          ab14co2_alpha(i,j) = abco2_alpha(i,j) * (1.0 + delta_14catm(i,j) * 1.0e-03)
       enddo; enddo ; !} i, j

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

       if ((trim(as_param_abiotic) == 'W92') .or. (trim(as_param_abiotic) == 'gfdl_cmip6')) then
           abco2_sc_no(i,j)  =  abiotic%sA_co2 + ST * (abiotic%sB_co2 + ST * (abiotic%sC_co2 + ST * abiotic%sD_co2)) * &
                                grid_tmask(i,j,1)

           ab14co2_sc_no(i,j) = abiotic%sA_co2 + ST * (abiotic%sB_co2 + ST * (abiotic%sC_co2 + ST * abiotic%sD_co2)) * &
                                grid_tmask(i,j,1)

       else if (trim(as_param_abiotic) == 'W14') then
           abco2_sc_no(i,j) = abiotic%sA_co2 + ST*(abiotic%sB_co2 + ST*(abiotic%sC_co2 + & 
                                               ST*(abiotic%sD_co2 + ST*abiotic%sE_co2))) * &
                                               grid_tmask(i,j,1)

           ab14co2_sc_no(i,j) = abiotic%sA_co2 + ST*(abiotic%sB_co2 + ST*(abiotic%sC_co2 + & 
                                               ST*(abiotic%sD_co2 + ST*abiotic%sE_co2))) * &
                                               grid_tmask(i,j,1)

       endif

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

    deallocate(abco2_alpha, abco2_csurf, abco2_sc_no,       &
               ab14co2_alpha, ab14co2_csurf, ab14co2_sc_no, &
               delta_14catm, htotal_field, htotal14c_field)

  end subroutine generic_abiotic_set_boundary_values


!> \brief Internal subroutine to deallocate user arrays
!!
  subroutine user_deallocate_arrays
    deallocate(abiotic%htotallo)
    deallocate(abiotic%htotalhi)
    deallocate(abiotic%htotal14clo)
    deallocate(abiotic%htotal14chi)
    deallocate(abiotic%delta_14catm)
    deallocate(abiotic%zt)
    deallocate(abiotic%f_alk)
    deallocate(abiotic%f_htotal)
    deallocate(abiotic%f_htotal14c)
    deallocate(abiotic%f_dissicabio)
    deallocate(abiotic%f_dissi14cabio)
    deallocate(abiotic%f_po4)
    deallocate(abiotic%f_sio4)
    deallocate(abiotic%jdecay_di14c)
    deallocate(abiotic%abco2_csurf)
    deallocate(abiotic%abco2_alpha)
    deallocate(abiotic%ab14co2_csurf)
    deallocate(abiotic%ab14co2_alpha)
    deallocate(abiotic%abpco2_csurf)
    deallocate(abiotic%abp14co2_csurf)
  end subroutine user_deallocate_arrays 

!> \brief Subroutine to end the abiotic module and deallocate all work arrays
!!
  subroutine generic_abiotic_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_abiotic_end'
    call user_deallocate_arrays
  end subroutine generic_abiotic_end

end module generic_abiotic
