!-----------------------------------------------------------------------
!
! <CONTACT EMAIL="Richard.Matear@csiro.au"> Richard Matear
! </CONTACT>
!
! <CONTACT EMAIL="Matthew.Chamberlain@csiro.au"> Matt Chamberlain
! </CONTACT>
!
! <CONTACT EMAIL="Dougal.Squire@anu.edu.au"> Dougie Squire
! </CONTACT>
!
! <CONTACT EMAIL="Pearse.Buchanan@csiro.au"> Pearse Buchanan
! </CONTACT>
!
! <OVERVIEW>
!  This module contains the generic version of WOMBATlite.
!  It is designed so that both GFDL Ocean models, GOLD and MOM, can use
!  it.
! </OVERVIEW>
!
! <DESCRIPTION>
!  Whole Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT) is
!  based on a NPZD (nutrient–phytoplankton–zooplankton–detritus) model.
!  This is the "lite" version of WOMBAT which includes one class each of
!  phyto-plankton and zooplankton and the addition of bio-available iron
!  limitation (Fe), dissolved inorganic carbon (DIC), calcium carbonate
!  (CaCO3), alkalinity (ALK), and oxygen (O).
!  Gas exchange follows OCMIP2 protocols.
! </DESCRIPTION>
!
! <INFO>
!  <REFERENCE>
!   This model is available for public use. Note that different tracers
!   may exist in different versions of the module.
!  </REFERENCE>
!
!  <DEVELOPER_NOTES>
!   This code was originally ported from WOMBAT v3 here:
!   https://github.com/mom-ocean/MOM5/tree/d7ba13a3f364ce130b6ad0ba813f01832cada7a2/src/mom5/ocean_csiro_bgc
!   using generic_BLING.F90 as a template.
!  </DEVELOPER_NOTES>
! </INFO>
!
! <NAMELIST NAME="generic_wombatlite_nml">
!  <DATA NAME="co2_calc" TYPE="character">
!   Defines the carbon equiliabration method.  Default is 'ocmip2' which
!   uses the FMS_ocmip2_co2calc routine.  The other option is 'mocsy',
!   which uses the set of routines authored by J. Orr. See reference at:
!   http://ocmip5.ipsl.jussieu.fr/mocsy/index.html
!  </DATA>
! </NAMELIST>
!
!-----------------------------------------------------------------------

module generic_WOMBATlite

  use field_manager_mod, only: fm_string_len
  use mpp_mod,           only: input_nml_file, mpp_error, FATAL
  use fms_mod,           only: write_version_number, check_nml_error, stdout, stdlog
  use time_manager_mod,  only: time_type
  use constants_mod,     only: WTMCO2, WTMO2

  use g_tracer_utils, only : g_diag_type, g_tracer_type
  use g_tracer_utils, only : g_tracer_start_param_list, g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add, g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_get_common, g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_values, g_tracer_set_values
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use g_tracer_utils, only : g_send_data

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  character(len=fm_string_len), parameter :: mod_name     = 'generic_WOMBATlite'
  character(len=fm_string_len), parameter :: package_name = 'generic_wombatlite'

  public do_generic_WOMBATlite
  public generic_WOMBATlite_register
  public generic_WOMBATlite_init
  public generic_WOMBATlite_register_diag
  public generic_WOMBATlite_update_from_coupler
  public generic_WOMBATlite_update_from_source
  public generic_WOMBATlite_update_from_bottom
  public generic_WOMBATlite_set_boundary_values
  public generic_WOMBATlite_end

  ! The following variable for using this module is overwritten by
  ! generic_tracer_nml namelist
  logical, save :: do_generic_WOMBATlite = .false.

  real, parameter :: missing_value1 = -1.0e+10

  !=======================================================================
  ! Namelist Options
  !=======================================================================
  character(len=10) :: co2_calc = 'ocmip2' ! other option is 'mocsy'

  namelist /generic_wombatlite_nml/ co2_calc

  !=======================================================================
  ! This type contains all the parameters and arrays used in this module
  !=======================================================================
  type generic_WOMBATlite_type
    !-----------------------------------------------------------------------
    ! Configurable parameters
    !-----------------------------------------------------------------------
    ! See user_add_params for descriptions of each parameter
    logical :: &
        init, &
        force_update_fluxes ! Set in generic_tracer_nml

    real :: &
        alphabio, &
        abio, &
        bbio, &
        cbio, &
        k1bio, &
        muepbio, &
        muepsbio, &
        gam1bio, &
        gbio, &
        epsbio, &
        muezbio, &
        gam2bio, &
        muedbio, &
        muedbio_sed, &
        wdetbio, &
        muecaco3, &
        muecaco3_sed, &
        wcaco3, &
        f_inorg, &
        tscav_fe, &
        fe_bkgnd, &
        dt_npzd, &
        sal_global, &
        dic_global, &
        adic_global, &
        alk_global, &
        no3_global, &
        sio2_surf, &
        htotal_scale_lo, &
        htotal_scale_hi, &
        htotal_in, &
        Rho_0, &
        a_0, a_1, a_2, a_3, a_4, a_5, &
        b_0, b_1, b_2, b_3, c_0, &
        a1_co2, a2_co2, a3_co2, a4_co2, &
        a1_o2, a2_o2, a3_o2, a4_o2

    character(len=fm_string_len) :: ice_restart_file
    character(len=fm_string_len) :: ocean_restart_file
    character(len=fm_string_len) :: IC_file

    !-----------------------------------------------------------------------
    ! Arrays for surface gas fluxes
    !-----------------------------------------------------------------------
    real, dimension(:,:), allocatable :: &
        htotallo, htotalhi, htotal,  &
        ahtotallo, ahtotalhi, ahtotal,  &
        sio2, &
        co2_csurf, co2_alpha, co2_sc_no, pco2_csurf, &
        aco2_csurf, aco2_alpha, paco2_csurf, &
        o2_csurf, o2_alpha, o2_sc_no, &
        no3_vstf, dic_vstf, adic_vstf, alk_vstf

    !-----------------------------------------------------------------------
    ! Arrays for tracer fields and source terms
    !-----------------------------------------------------------------------
    ! The prefixes "f_" refers to a "field", "j" to a volumetric rate, "b_"
    ! to a bottom flux and "p_" to a "pointer".
    real, dimension(:,:), allocatable :: &
        b_dic, &
        b_adic, &
        b_alk, &
        b_no3, &
        b_o2, &
        b_fe, &
        light_limit, &
        pprod_gross_2d, &
        wdet100, &
        export_prod, &
        export_inorg, &
        npp2d, &
        det_btm, &
        caco3_btm, &
        det_sed_remin, &
        caco3_sed_remin, &
        dic_intmld, &
        adic_intmld, &
        o2_intmld, &
        no3_intmld, &
        fe_intmld, &
        phy_intmld, &
        det_intmld, &
        pprod_gross_intmld, &
        npp_intmld, &
        radbio_intmld, &
        dic_int100, &
        adic_int100, &
        o2_int100, &
        no3_int100, &
        fe_int100, &
        phy_int100, &
        det_int100, &
        pprod_gross_int100, &
        npp_int100, &
        radbio_int100
      
    real, dimension(:,:,:), allocatable :: &
        f_dic, &
        f_adic, &
        f_alk, &
        f_no3, &
        f_phy, &
        f_zoo, &
        f_det, &
        f_o2, &
        f_caco3, &
        f_fe, &
        pprod_gross, &
        zprod_gross, &
        radbio3d, &
        npp3d, &
        vpbio, &
        avej, &
        no3_prev, &
        caco3_prev, &
        zw, &
        zm

    real, dimension(:,:,:,:), pointer :: &
        p_dic, &
        p_adic, &
        p_alk, &
        p_o2

    real, dimension(:,:,:), pointer :: &
        p_det_sediment, &
        p_caco3_sediment

    real, dimension(:,:), pointer :: &
        p_no3_stf, &
        p_dic_stf, &
        p_adic_stf, &
        p_alk_stf

    !-----------------------------------------------------------------------
    ! IDs for diagnostics
    !-----------------------------------------------------------------------
    ! See register_diagnostics for descriptions of each diagnostic
    integer :: &
        id_pco2 = -1, &
        id_paco2 = -1, &
        id_htotal = -1, &
        id_ahtotal = -1, &
        id_no3_vstf = -1, &
        id_dic_vstf = -1, &
        id_adic_vstf = -1, &
        id_alk_vstf = -1, &
        id_light_limit = -1, &
        id_radbio3d = -1, &
        id_radbio1 = -1, &
        id_pprod_gross = -1, &
        id_pprod_gross_2d = -1, &
        id_wdet100 = -1, &
        id_export_prod = -1, &
        id_export_inorg = -1, &
        id_npp3d = -1, &
        id_npp2d = -1, &
        id_npp1 = -1, &
        id_zprod_gross = -1, &
        id_dic_intmld = -1, &
        id_adic_intmld = -1, &
        id_o2_intmld = -1, &
        id_no3_intmld = -1, &
        id_fe_intmld = -1, &
        id_phy_intmld = -1, &
        id_det_intmld = -1, &
        id_pprod_gross_intmld = -1, &
        id_npp_intmld = -1, &
        id_radbio_intmld = -1, &
        id_dic_int100 = -1, &
        id_adic_int100 = -1, &
        id_o2_int100 = -1, &
        id_no3_int100 = -1, &
        id_fe_int100 = -1, &
        id_phy_int100 = -1, &
        id_det_int100 = -1, &
        id_pprod_gross_int100 = -1, &
        id_npp_int100 = -1, &
        id_radbio_int100 = -1, &
        id_det_sed_remin = -1, &
        id_det_sed_depst = -1, &
        id_caco3_sed_remin = -1, &
        id_caco3_sed_depst = -1

  end type generic_WOMBATlite_type

  type(generic_WOMBATlite_type), save :: wombat

  ! An auxiliary type for storing varible names
  type, public :: vardesc
    character(len=fm_string_len) :: name     ! The variable name in a NetCDF file.
    character(len=fm_string_len) :: longname ! The long name of that variable.
    character(len=1)             :: hor_grid ! The hor. grid:  u, v, h, q, or 1.
    character(len=1)             :: z_grid   ! The vert. grid:  L, i, or 1.
    character(len=1)             :: t_grid   ! The time description: s, a, m, or 1.
    character(len=fm_string_len) :: units    ! The dimensions of the variable.
    character(len=1)             :: mem_size ! The size in memory: d or f.
  end type vardesc

  type(CO2_dope_vector) :: CO2_dope_vec

  contains

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_register">
  !  <OVERVIEW>
  !   Register the generic WOMBATlite module
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This subroutine reads and checks the WOMBATlite namelist and adds all
  !   WOMBATlite tracers via subroutine user_add_tracers()
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_register(tracer_list, force_update_fluxes)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    integer                                 :: ierr
    integer                                 :: io_status
    integer                                 :: stdoutunit, stdlogunit
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATlite_register'
    character(len=256), parameter           :: error_header = &
        '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter           :: warn_header =  &
        '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
    character(len=256), parameter           :: note_header =  &
        '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '

    ! Provide for namelist over-ride
    ! This needs to go before user_add_tracers in order to allow the namelist
    ! settings to switch tracers on and off.
    stdoutunit = stdout(); stdlogunit = stdlog()

    read (input_nml_file, nml=generic_wombatlite_nml, iostat=io_status)
    ierr = check_nml_error(io_status, 'generic_wombatlite_nml')

    write (stdoutunit,'(/)')
    write (stdoutunit, generic_wombatlite_nml)
    write (stdlogunit, generic_wombatlite_nml)

    if (trim(co2_calc) == 'ocmip2') then
      write (stdoutunit,*) trim(note_header), 'Using FMS OCMIP2 CO2 routine'
    else if (trim(co2_calc) == 'mocsy') then
      write (stdoutunit,*) trim(note_header), 'Using Mocsy CO2 routine'
    else
      call mpp_error(FATAL,"Unknown co2_calc option specified in generic_wombatlite_nml")
    endif

    ! Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_WOMBATlite_register

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_init">
  !  <OVERVIEW>
  !   Initialize the generic WOMBATlite module
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This subroutine: adds all the WOMBATlite tracers to the list of
  !   generic tracers passed to it via utility subroutine g_tracer_add();
  !   adds all the parameters used by this module via utility subroutine
  !   g_tracer_add_param(); and allocates all work arrays used in the
  !   module.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_init(tracer_list, force_update_fluxes)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="force_update_fluxes" TYPE="logical">
  !   Flag to force update the fluxes every timestep. This maybe be
  !   necessary in situations where the column_physics (update_from_source)
  !   is not called every timestep such as when MOM6
  !   THERMO_SPANS_COUPLING=True
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_init(tracer_list, force_update_fluxes)
    type(g_tracer_type), pointer :: tracer_list
    logical, intent(in)          :: force_update_fluxes

    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATlite_init'

    wombat%force_update_fluxes = force_update_fluxes

    call write_version_number( version, tagname )

    ! Specify and initialize all parameters used by this package
    call user_add_params

    ! Allocate all the private work arrays used by this module.
    call user_allocate_arrays

  end subroutine generic_WOMBATlite_init

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_register_diag">
  !  <OVERVIEW>
  !   Register diagnostic fields to be used in this module.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Register diagnostic fields to be used in this module. Note that the
  !   tracer fields are automatically registered in user_add_tracers. User
  !   adds only diagnostics for fields that are not a member of
  !   g_tracer_type
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_register_diag(diag_list)
  !  </TEMPLATE>
  !
  !  <IN NAME="g_diag_type" TYPE="type(g_diag_type), pointer">
  !   Pointer to the head of generic diag list. Currently, this is not
  !   actually used.
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list ! dts: this is not actually used

    type(vardesc)   :: vardesc_temp
    integer         :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, axes(3)
    type(time_type) :: init_time

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        axes=axes, init_time=init_time)

    !=======================================================================
    ! Register all diagnostics in this module
    !=======================================================================
    !
    ! The following vardesc types contain a package of metadata about each
    ! tracer, including, in order, the following elements: name; longname;
    ! horizontal staggering ('h') for collocation with thickness points;
    ! vertical staggering ('L') for a layer variable ; temporal staggering
    ! ('s' for snapshot) ; units; and precision in non-restart output files
    ! ('f' for 32-bit float or 'd' for 64-bit doubles). For most tracers, only
    ! the name, longname and units should be changed.
    !
    ! Niki: The register_diag_field interface needs to be extended to take the
    ! MOM6 axes_grp as argument instead of this integer array axes_grp%handle.
    ! Currently the actual MOM6 diag axes is chosen to be T or Tl based on the
    ! size of the axes argument, 2 or 3. The actual values of these axes
    ! argument are not used, only their size is checked to determine the diag
    ! axes! This is not correct since axesTi and axesTl are both of size 3,
    ! likewise there are many axes of size 2. To accomodate axesTi with the
    ! least amount of code modification we can set and check for an input
    ! array of size 1.

    !=======================================================================
    ! Surface flux diagnostics
    !=======================================================================
    !
    ! dts: other gas exchange diagnostics are available via the "ocean_flux",
    ! "ocean_sfc", "atmos_sfc" diagnostic module_names
    vardesc_temp = vardesc( &
        'pco2', 'Surface aqueous partial pressure of CO2', 'h', '1', 's', 'uatm', 'f')
    wombat%id_pco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
  
    vardesc_temp = vardesc( &
        'paco2', 'Surface aqueous partial pressure of CO2 inc. anthropogenic', &
        'h', '1', 's', 'uatm', 'f')
    wombat%id_paco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'htotal', 'H+ concentration', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ahtotal', 'H+ concentration inc. anthropogenic', 'h', '1', 's', 'mol/kg', 'f')
    wombat%id_ahtotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'no3_vstf', 'Virtual flux of nitrate into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_no3_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'dic_vstf', 'Virtual flux of natural dissolved inorganic carbon into ocean due to '// &
        'salinity restoring/correction', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_dic_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'adic_vstf', 'Virtual flux of natural + anthropogenic dissolved inorganic carbon '// &
        'into ocean due to salinity restoring/correction', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_adic_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'alk_vstf', 'Virtual flux of alkalinity into ocean due to salinity restoring/correction', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_alk_vstf = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    !=======================================================================
    ! Tracer and source term diagnostics
    !=======================================================================
    vardesc_temp = vardesc( &
        'light_limit', 'Integrated light limitation of phytoplankton growth', &
        'h', '1', 's', 'm', 'f')
    wombat%id_light_limit = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio3d', 'Photosynthetically active radiation for phytoplankton growth', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radbio3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio1', 'Photosynthetically active radiation for phytoplankton growth at surface', &
        'h', '1', 's', 'W m-2', 'f')
    wombat%id_radbio1 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'det_sed_remin', 'Rate of remineralisation of detritus in accumulated sediment', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_det_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'det_sed_depst', 'Rate of deposition of detritus to sediment at base of water column', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_det_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caco3_sed_remin', 'Rate of remineralisation of CaCO3 in accumulated sediment', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_caco3_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caco3_sed_depst', 'Rate of deposition of CaCO3 to sediment at base of water column', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_caco3_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'wdet100', 'Detritus export at 100 m (det*sinking rate)', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_wdet100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'export_prod', 'Organic export through 100m', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_export_prod = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'export_inorg', 'Inorganic export through 100m', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_export_inorg = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp3d', 'Net primary productivity', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_npp3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp2d', 'Vertically integrated net primary productivity', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
  
    vardesc_temp = vardesc( &
        'npp1', 'Net primary productivity in the first ocean layer', &
        'h', '1', 's', 'mol/kg/s', 'f')
    wombat%id_npp1 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
      
    vardesc_temp = vardesc( &
        'pprod_gross', 'Gross phytoplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_pprod_gross = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pprod_gross_2d', 'Vertically integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zprod_gross', 'Gross zooplankton production', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_zprod_gross = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    ! MLD-integrated diagnostics
    !-----------------------------------------------------------------------
    vardesc_temp = vardesc( &
        'dic_intmld', 'MLD-integrated natural dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_dic_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'adic_intmld', 'MLD-integrated natural + anthropogenic dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_adic_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'o2_intmld', 'MLD-integrated dissolved oxygen', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_o2_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'no3_intmld', 'MLD-integrated nitrate', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_no3_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'fe_intmld', 'MLD-integrated iron', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_fe_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'phy_intmld', 'MLD-integrated phytoplankton', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_phy_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'det_intmld', 'MLD-integrated detritus', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_det_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'pprod_gross_intmld', 'MLD-integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp_intmld', 'MLD-integrated net primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio_intmld', 'MLD-integrated photosynthetically active radiation for phytoplankton '// &
        'growth', 'h', '1', 's', 'W m-1', 'f')
    wombat%id_radbio_intmld = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    ! 100m-integrated diagnostics
    !-----------------------------------------------------------------------
    vardesc_temp = vardesc( &
        'dic_int100', '100m-integrated natural dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_dic_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'adic_int100', '100m-integrated natural + anthropogenic dissolved inorganic carbon', &
        'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_adic_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'o2_int100', '100m-integrated dissolved oxygen', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_o2_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'no3_int100', '100m-integrated nitrate', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_no3_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'fe_int100', '100m-integrated iron', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_fe_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'phy_int100', '100m-integrated phytoplankton', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_phy_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'det_int100', '100m-integrated detritus', 'h', '1', 's', 'mol/m^2', 'f')
    wombat%id_det_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'pprod_gross_int100', '100m-integrated gross phytoplankton production', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_pprod_gross_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'npp_int100', '100m-integrated net primary productivity', 'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_npp_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radbio_int100', '100m-integrated photosynthetically active radiation for phytoplankton '// &
        'growth', 'h', '1', 's', 'W m-1', 'f')
    wombat%id_radbio_int100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

  end subroutine generic_WOMBATlite_register_diag

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Add all the parameters to be used in this module.
  !
  subroutine user_add_params

    !=======================================================================
    ! Specify all parameters used in this modules.
    !=======================================================================
    !
    ! Add the known experimental parameters used for calculations in this
    ! module. All the g_tracer_add_param calls must happen between
    ! g_tracer_start_param_list and g_tracer_end_param_list calls. This
    ! implementation enables runtime overwrite via field_table.
    ! dts: Note, some parameters are required by the user_add_tracers routine
    ! which is run _before_ this one. Those parameters are added in
    ! user_add_tracers.

    ! User adds one call for each parameter below with the template
    ! g_tracer_add_param(name, variable,  default_value)
    call g_tracer_start_param_list(package_name)

    !=======================================================================
    ! General parameters
    !=======================================================================
    !
    ! dts: This was copied from BLING to enable calculation of surface flux
    ! terms when the update_from_source routine is commented out for
    ! debugging.
    call g_tracer_add_param('init', wombat%init, .false. )

    ! Average density of sea water [kg/m^3]
    !-----------------------------------------------------------------------
    ! Rho_0 is used in the Boussinesq approximation to calculations of
    ! pressure and pressure gradients, in units of kg m-3.
    call g_tracer_add_param('Rho_0', wombat%Rho_0, 1035.0)
    
    !=======================================================================
    ! Surface gas flux parameters
    !=======================================================================

    ! Coefficients for O2 saturation [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', wombat%a_0, 2.00907)
    call g_tracer_add_param('a_1', wombat%a_1, 3.22014)
    call g_tracer_add_param('a_2', wombat%a_2, 4.05010)
    call g_tracer_add_param('a_3', wombat%a_3, 4.94457)
    call g_tracer_add_param('a_4', wombat%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', wombat%a_5, 3.88767)
    call g_tracer_add_param('b_0', wombat%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', wombat%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', wombat%b_2, -1.03410e-02)
    call g_tracer_add_param('b_3', wombat%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', wombat%c_0, -4.88682e-07)

    ! Schmidt number coefficients [1]
    !-----------------------------------------------------------------------
    ! Compute the Schmidt number of CO2 in seawater using the
    ! formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
    ! 7373-7382).
    call g_tracer_add_param('a1_co2', wombat%a1_co2,  2073.1)
    call g_tracer_add_param('a2_co2', wombat%a2_co2, -125.62)
    call g_tracer_add_param('a3_co2', wombat%a3_co2,  3.6276)
    call g_tracer_add_param('a4_co2', wombat%a4_co2, -0.043219)

    ! Compute the Schmidt number of O2 in seawater using the
    ! formulation proposed by Keeling et al. (1998, Global Biogeochem.
    ! Cycles, 12, 141-163).
    call g_tracer_add_param('a1_o2', wombat%a1_o2, 1638.0)
    call g_tracer_add_param('a2_o2', wombat%a2_o2, -81.83)
    call g_tracer_add_param('a3_o2', wombat%a3_o2, 1.483)
    call g_tracer_add_param('a4_o2', wombat%a4_o2, -0.008004)

    ! Initial H+ concentration [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_in', wombat%htotal_in, 1.e-8) ! dts: default conc from COBALT

    ! Scale factor to set lower limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_lo', wombat%htotal_scale_lo, 0.1)

    ! Scale factor to set upper limit of htotal range [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('htotal_scale_hi', wombat%htotal_scale_hi, 100.0)

    ! Global average surface concentration of inorganic silicate [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sio2_surf', wombat%sio2_surf, 35.0e-3 / 1035.0)

    !=======================================================================
    ! NPZD parameters
    !=======================================================================
    ! dts: note the parameter units and default values are as used by WOMBAT
    ! v3 in ACCESS-OM2 and ACCESS-ESM1.5. Unit conversions are done
    ! internally to account for the different units carried in this generic
    ! version of WOMBATlite.

    ! Initial slope of P-I curve [m2/W/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio', wombat%alphabio, 2.96296e-6)

    ! Phytoplankton maximum growth rate parameter a [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abio', wombat%abio, 3.125e-6)

    ! Phytoplankton maximum growth rate parameter b [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbio', wombat%bbio, 1.066)

    ! Phytoplankton maximum growth rate parameter c [1/K]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('cbio', wombat%cbio, 1.)

    ! Phytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('k1bio', wombat%k1bio, 0.7)

    ! Phytoplankton linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('muepbio', wombat%muepbio, 4.62963e-7)

    ! Phytoplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('muepsbio', wombat%muepsbio, 2.89352e-6)

    ! Zooplankton assimilation efficiency [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('gam1bio', wombat%gam1bio, 0.925)

    ! Zooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('gbio', wombat%gbio, 1.82292e-5)

    ! Zooplankton prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('epsbio', wombat%epsbio, 1.85185e-5)

    ! Zooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('muezbio', wombat%muezbio, 3.93519e-6)

    ! Zooplankton excretion rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('gam2bio', wombat%gam2bio, 1.15741e-7)

    ! Detritus remineralisation rate constant for <180 m; value for >=180m is half of this [1/s]
    !-----------------------------------------------------------------------
    ! Detritus remineralisation rate for (>=180 m) is hard-coded to be
    ! muedbio/2
    call g_tracer_add_param('muedbio', wombat%muedbio, 5.55556e-7)
    
    ! Detritus remineralisation rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    ! This would normally equal muedbio, but we set the default value to be
    ! consistent with what is in ACCESS-ESM1.5 (undocumented in Ziehn et al
    ! 2020)
    call g_tracer_add_param('muedbio_sed', wombat%muedbio_sed, 2.31481e-7)

    ! CaCO3 remineralisation rate constant [1/s]
    !-----------------------------------------------------------------------
    ! Default value matches 0.001714 day-1 in Ziehn et al 2020; differs from
    ! Hayashida et al 2020
    call g_tracer_add_param('muecaco3', wombat%muecaco3, 1.9838e-8)

    ! CaCO3 remineralization rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('muecaco3_sed', wombat%muecaco3_sed, 4.05093e-8)

    ! CaCO3 inorganic fraction [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('f_inorg', wombat%f_inorg, 0.062)

    ! Iron scavenging rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('tscav_fe', wombat%tscav_fe, 3.17098e-8)

    ! Iron background concentration [umol Fe m^-3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('fe_bkgnd', wombat%fe_bkgnd, 0.6)

    ! Nested timestep for the ecosystem model [s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dt_npzd', wombat%dt_npzd, 900.)

    ! Global average surface salinity used for virtual flux correction [g/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('sal_global', wombat%sal_global, 34.6)

    ! Global average surface dic used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('dic_global', wombat%dic_global, 1.90e-3)

    ! Global average surface adic used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('adic_global', wombat%adic_global, 1.90e-3)

    ! Global average surface alk used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alk_global', wombat%alk_global, 2.225e-3)

    ! Global average surface no3 used for virtual flux correction [mol/kg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('no3_global', wombat%no3_global, 4.512e-06)

    call g_tracer_end_param_list(package_name)

  end subroutine user_add_params

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Add all the tracers to be used in this module.
  !
  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    real                                    :: as_coeff_wombatlite

    !=======================================================================
    ! Parameters
    !=======================================================================
    ! Add here only the parameters that are required at the time of registeration
    ! (to make flux exchanging ocean tracers known for all PE's)

    ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
    !-----------------------------------------------------------------------
    ! From Wanninkhof 1992 for steady wind speed (in m/s)
    as_coeff_wombatlite = 0.31 / 3.6e5

    call g_tracer_start_param_list(package_name)

    ! CaCO3 sinking velocity [m/s]
    !-----------------------------------------------------------------------
    ! Default value matches Ziehn et al 2020 but differs from Hayashida et
    ! al 2020
    call g_tracer_add_param('wcaco3', wombat%wcaco3, 6.94444e-5)

    ! Detritus sinking velocity [m/s]
    !-----------------------------------------------------------------------
    ! Default value matches Ziehn et al 2020 but differs from Hayashida et
    ! al 2020
    call g_tracer_add_param('wdetbio', wombat%wdetbio, 2.77778e-4)

    call g_tracer_add_param('ice_restart_file', wombat%ice_restart_file, 'ice_wombatlite.res.nc')
    call g_tracer_add_param('ocean_restart_file', wombat%ocean_restart_file, 'ocean_wombatlite.res.nc')
    call g_tracer_add_param('IC_file', wombat%IC_file, '')

    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file = wombat%ice_restart_file, &
        ocean_restart_file = wombat%ocean_restart_file )

    !=======================================================================
    ! Specify all tracers of this module
    !=======================================================================
    !
    ! User adds one call for each tracer below!
    ! User should specify if fluxes must be extracted from boundary by passing
    ! one or more of the following methods as .true. and provide the corresponding
    ! parameters array methods: flux_gas, flux_runoff, flux_wetdep, flux_drydep.
    ! Pass an init_value arg if the tracers should be initialized to a nonzero
    ! value everywhere otherwise they will be initialized to zero.
    !
    ! dts: diagnostic (prog = .false.) tracers added here are automatically
    ! registered for restart but not for horizontal advection and diffusion. All
    ! tracer fields are registered for diag output.

    !=======================================================================
    ! Prognostic Tracers
    !=======================================================================

    ! Nitrate
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Nitrate
    call g_tracer_add(tracer_list, package_name, &
        name = 'no3', &
        longname = 'Nitrate', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Phytoplankton
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'phy', &
        longname = 'Phytoplankton', &
        units = 'mol/kg', &
        prog = .true.)

    ! Oxygen
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'o2', &
        longname = 'Oxygen', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true.,  &
        flux_bottom = .true., &
        flux_gas_name = 'o2_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMO2, &
        flux_gas_param = (/ as_coeff_wombatlite, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatlite_airsea_flux.res.nc')
    
    ! Zooplankton
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'zoo', &
        longname = 'Zooplankton', &
        units = 'mol/kg', &
        prog = .true.)

    ! Detritus
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'det', &
        longname = 'Detritus', &
        units = 'mol/kg', &
        prog = .true., &
        sink_rate = wombat%wdetbio, &
        btm_reservoir = .true.)

    ! CaCO3
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'caco3', &
        longname = 'CaCO3', &
        units = 'mol/kg', &
        prog = .true., &
        sink_rate = wombat%wcaco3, &
        btm_reservoir = .true.)

    ! ADIC (Natural + anthropogenic dissolved inorganic carbon)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'adic', &
        longname = 'Natural + anthropogenic Dissolved Inorganic Carbon', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true., &
        flux_bottom = .true., &
        flux_gas_name = 'co2_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMCO2, &
        flux_gas_param = (/ as_coeff_wombatlite, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatlite_airsea_flux.res.nc', &
        flux_virtual = .true.)
      
    ! DIC (Natural dissolved inorganic carbon)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'dic', &
        longname = 'Natural Dissolved Inorganic Carbon', &
        units = 'mol/kg', &
        prog = .true., &
        flux_gas = .true., &
        flux_bottom = .true., &
        flux_gas_name = 'co2_nat_flux', &
        flux_gas_type = 'air_sea_gas_flux_generic', &
        flux_gas_molwt = WTMCO2, &
        flux_gas_param = (/ as_coeff_wombatlite, 9.7561e-06 /), & ! dts: param(2) converts Pa -> atm
        flux_gas_restart_file = 'ocean_wombatlite_airsea_flux.res.nc', &
        flux_virtual = .true.)

    ! Alk (Total carbonate alkalinity)
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'alk', &
        longname = 'Alkalinity', &
        units = 'mol/kg', &
        prog = .true., &
        flux_bottom = .true., &
        flux_virtual = .true.)

    ! Dissolved Iron
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'fe', &
        longname = 'Dissolved Iron', &
        units = 'mol/kg', &
        prog = .true., &
        flux_drydep = .true., &
        flux_param  = (/ 1.0 /), & ! dts: conversion to mol/m2/s done in data_table
        flux_bottom = .true.)

    !=======================================================================
    ! Diagnostic Tracers
    !=======================================================================

    ! Detritus sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'det_sediment', &
        longname = 'Detritus at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

    ! CaCO3 sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name  = 'caco3_sediment', &
        longname = 'CaCO3 at base of column as sediment', &
        units = 'mol m-2', &
        prog = .false.)

  end subroutine user_add_tracers

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_update_from_coupler">
  !  <OVERVIEW>
  !     Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !    Some tracer fields could be modified after values are obtained from
  !    the coupler. This subroutine is the place for specific tracer
  !    manipulations. In WOMBATlite we apply virtual flux corrections due
  !    to salt flux restoring/correction here.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_update_from_coupler(tracer_list)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="salt_flux_added" TYPE="real, dimension(ilb:,jlb:), optional">
  !   Surface salt flux into ocean from restoring or flux adjustment [g/m2/s]
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_update_from_coupler(tracer_list, ilb, jlb, salt_flux_added)
    type(g_tracer_type), pointer           :: tracer_list
    integer, intent(in)                    :: ilb, jlb
    real, dimension(ilb:,jlb:), intent(in) :: salt_flux_added

    ! Account for virtual fluxes due to salt flux restoring/correction
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'no3', 'stf', wombat%p_no3_stf)
    wombat%no3_vstf(:,:) = (wombat%no3_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_no3_stf(:,:) = wombat%p_no3_stf(:,:) + wombat%no3_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'dic', 'stf', wombat%p_dic_stf)
    wombat%dic_vstf(:,:) = (wombat%dic_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_dic_stf(:,:) = wombat%p_dic_stf(:,:) + wombat%dic_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'adic', 'stf', wombat%p_adic_stf)
    wombat%adic_vstf(:,:) = (wombat%adic_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_adic_stf(:,:) = wombat%p_adic_stf(:,:) + wombat%adic_vstf(:,:) ! [mol/m2/s]

    call g_tracer_get_pointer(tracer_list, 'alk', 'stf', wombat%p_alk_stf)
    wombat%alk_vstf(:,:) = (wombat%alk_global / wombat%sal_global) * salt_flux_added(:,:) ! [mol/m2/s]
    wombat%p_alk_stf(:,:) = wombat%p_alk_stf(:,:) + wombat%alk_vstf(:,:) ! [mol/m2/s]

  end subroutine generic_WOMBATlite_update_from_coupler

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_update_from_bottom">
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Some tracers could have bottom fluxes and reservoirs. This subroutine
  !   is the place for specific tracer manipulations.
  !   In WOMBATlite, remineralization from the sediment tracers (which
  !   requires temperature) is done in update_from_source. Deposition from
  !   sinking is handled here.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_update_from_bottom(tracer_list, dt, tau)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real, intent(in)             :: dt
    integer, intent(in)          :: tau
    type(time_type), intent(in)  :: model_time

    integer                         :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau
    real, dimension(:,:,:), pointer :: grid_tmask
    logical                         :: used

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    ! Move bottom reservoirs to sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_values(tracer_list, 'det', 'btm_reservoir', wombat%det_btm, isd, jsd)
    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment)
    wombat%p_det_sediment(:,:,1) = wombat%p_det_sediment(:,:,1) + wombat%det_btm(:,:) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'det', 'btm_reservoir', 0.0)

    call g_tracer_get_values(tracer_list, 'caco3', 'btm_reservoir', wombat%caco3_btm, isd, jsd)
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment)
    wombat%p_caco3_sediment(:,:,1) =  wombat%p_caco3_sediment(:,:,1) + wombat%caco3_btm(:,:) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'caco3', 'btm_reservoir', 0.0)

    ! Send diagnostics
    !-----------------------------------------------------------------------
    if (wombat%id_det_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_det_sed_depst, wombat%det_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_caco3_sed_depst, wombat%caco3_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  end subroutine generic_WOMBATlite_update_from_bottom

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for
  !   calculating the interaction of tracers with each other and with outside
  !   forcings.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_update_from_source(tracer_list, Temp, Salt, &
  !     dzt, hblt_depth, ilb, jlb, tau, dt, grid_dat, sw_pen, opacity)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="Temp" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean temperature
  !  </IN>
  !
  !  <IN NAME="Salt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean salinity
  !  </IN>
  !
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  !
  !  <IN NAME="opacity" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean opacity
  !  </IN>
  !
  !  <IN NAME="sw_pen" TYPE="real, dimension(ilb:,jlb:)">
  !   Shortwave peneteration
  !  </IN>
  !
  !  <IN NAME="hblt_depth" TYPE="real, dimension(ilb:,jlb:)">
  !   Depth of actively mixing layer
  !  </IN>
  !
  !  <IN NAME="grid_dat" TYPE="real, dimension(ilb:,jlb:)">
  !   Grid area
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine generic_WOMBATlite_update_from_source(tracer_list, Temp, Salt,  &
      rho_dzt, dzt, hblt_depth, ilb, jlb, tau, dt, grid_dat, model_time, nbands, &
      max_wavelength_band, sw_pen_band, opacity_band)
    type(g_tracer_type), pointer               :: tracer_list
    real, dimension(ilb:,jlb:,:), intent(in)   :: Temp, Salt, rho_dzt, dzt
    real, dimension(ilb:,jlb:), intent(in)     :: hblt_depth
    integer, intent(in)                        :: ilb, jlb, tau
    real, intent(in)                           :: dt
    real, dimension(ilb:,jlb:), intent(in)     :: grid_dat
    type(time_type), intent(in)                :: model_time
    integer, intent(in)                        :: nbands
    real, dimension(:), intent(in)             :: max_wavelength_band
    real, dimension(:,ilb:,jlb:), intent(in)   :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band

    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATlite_update_from_source'
    integer                                 :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, tn
    integer                                 :: i, j, k, n, nz
    real, dimension(:,:,:), pointer         :: grid_tmask
    integer, dimension(:,:), pointer        :: grid_kmt
    integer, dimension(:,:), allocatable    :: kmeuph ! deepest level of euphotic zone
    integer, dimension(:,:), allocatable    :: k100 ! deepest level less than 100 m
    real                                    :: mmol_m3_to_mol_kg, umol_m3_to_mol_kg, mmol2_m6_to_mol2_kg2
    integer                                 :: ts_npzd ! number of time steps within NPZD model
    real                                    :: dtsb ! number of seconds per NPZD timestep
    real                                    :: rdtts ! 1 / dt
    real, dimension(nbands)                 :: sw_pen
    real                                    :: sw_pen_top, sw_pen_bottom
    real                                    :: u_npz, g_npz
    real                                    :: fbc
    real                                    :: f11, f21, f22, f23, f31, f32, f41, f51
    real                                    :: no3_bgc_change, caco3_bgc_change
    real                                    :: epsi = 1.0e-30
    logical                                 :: used

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask, grid_kmt=grid_kmt)

    ! dts: Note, other generic_tracer modules call this zt. However here we
    ! call this zw to be consistent with WOMBAT v3, which uses MOM5 terminology.
    ! zm here is halfway between interfaces
    wombat%zw(:,:,:) = 0.0
    wombat%zm(:,:,:) = 0.0
    do j = jsc,jec; do i = isc,iec
      wombat%zw(i,j,1) = dzt(i,j,1)
      wombat%zm(i,j,1) = 0.5 * dzt(i,j,1)
    enddo; enddo
    do k = 2,nk; do j = jsc,jec ; do i = isc,iec
      wombat%zw(i,j,k) = wombat%zw(i,j,k-1) + dzt(i,j,k)
      wombat%zm(i,j,k) = wombat%zw(i,j,k-1) + 0.5 * dzt(i,j,k)
    enddo; enddo ; enddo

    !=======================================================================
    ! Surface gas fluxes
    !=======================================================================
    !
    ! Calculate the surface gas fluxes for the next round of exchange. This
    ! is done here to align with other generic_tracer modules (e.g. BLING).
    ! dts: I think this done here in other modules because they calculate
    ! 3D Carbonate ion concentration (co3_ion) here using the FMS_ocmip2_co2calc
    ! routine. The FMS_ocmip2_co2calc routine also calculates co2star, alpha
    ! and pco2surf, so it makes sense to set these values here rather than
    ! recalculating them in set_boundary_values.
    
    call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=tau, &
        positive=.true.)
    call g_tracer_get_values(tracer_list, 'adic', 'field', wombat%f_adic, isd, jsd, ntau=tau, &
        positive=.true.)
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau, &
        positive=.true.)
    call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=tau, &
        positive=.true.)
  
    do j = jsc,jec; do i = isc,iec
      wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j)
      wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j)
      wombat%ahtotallo(i,j) = wombat%htotal_scale_lo * wombat%ahtotal(i,j)
      wombat%ahtotalhi(i,j) = wombat%htotal_scale_hi * wombat%ahtotal(i,j)
    enddo; enddo 

    call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
        Temp(:,:,1), Salt(:,:,1), &
        wombat%f_dic(:,:,1), &
        wombat%f_no3(:,:,1) / 16., &
        wombat%sio2(:,:), &
        wombat%f_alk(:,:,1), &
        wombat%htotallo(:,:), wombat%htotalhi(:,:), &
        wombat%htotal(:,:), &
        co2_calc=trim(co2_calc), &
        zt=wombat%zw(:,:,1), &
        co2star=wombat%co2_csurf(:,:), alpha=wombat%co2_alpha(:,:), &
        pCO2surf=wombat%pco2_csurf(:,:))

    call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
        Temp(:,:,1), Salt(:,:,1), &
        wombat%f_adic(:,:,1), &
        wombat%f_no3(:,:,1) / 16., &
        wombat%sio2(:,:), &
        wombat%f_alk(:,:,1), &
        wombat%ahtotallo(:,:), wombat%ahtotalhi(:,:), &
        wombat%ahtotal(:,:), &
        co2_calc=trim(co2_calc), &
        zt=wombat%zw(:,:,1), &
        co2star=wombat%aco2_csurf(:,:), alpha=wombat%aco2_alpha(:,:), &
        pCO2surf=wombat%paco2_csurf(:,:))

    call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'alpha', wombat%aco2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'csurf', wombat%aco2_csurf, isd, jsd)

    !=======================================================================
    ! Calculate the source terms
    !=======================================================================

    wombat%pprod_gross(:,:,:) = 0.0
    wombat%zprod_gross(:,:,:) = 0.0
    wombat%light_limit(:,:) = 0.0
    wombat%radbio3d(:,:,:) = 0.0
    wombat%npp3d(:,:,:) = 0.0
    wombat%vpbio(:,:,:) = 0.0
    wombat%avej(:,:,:) = 0.0
    wombat%export_prod(:,:) = 0.0
    wombat%export_inorg(:,:) = 0.0
    wombat%adic_intmld(:,:) = 0.0
    wombat%dic_intmld(:,:) = 0.0
    wombat%o2_intmld(:,:) = 0.0
    wombat%no3_intmld(:,:) = 0.0
    wombat%fe_intmld(:,:) = 0.0
    wombat%phy_intmld(:,:) = 0.0
    wombat%det_intmld(:,:) = 0.0
    wombat%pprod_gross_intmld(:,:) = 0.0
    wombat%npp_intmld(:,:) = 0.0
    wombat%radbio_intmld(:,:) = 0.0
    wombat%adic_int100(:,:) = 0.0
    wombat%dic_int100(:,:) = 0.0
    wombat%o2_int100(:,:) = 0.0
    wombat%no3_int100(:,:) = 0.0
    wombat%fe_int100(:,:) = 0.0
    wombat%phy_int100(:,:) = 0.0
    wombat%det_int100(:,:) = 0.0
    wombat%pprod_gross_int100(:,:) = 0.0
    wombat%npp_int100(:,:) = 0.0
    wombat%radbio_int100(:,:) = 0.0

    ! Some unit conversion factors
    mmol_m3_to_mol_kg = 1.e-3 / wombat%Rho_0
    umol_m3_to_mol_kg = 1.e-3 * mmol_m3_to_mol_kg
    mmol2_m6_to_mol2_kg2 = mmol_m3_to_mol_kg * mmol_m3_to_mol_kg
    
    ! Set the maximum index for euphotic depth
    ! dts: in WOMBAT v3, kmeuph and k100 are integers but here they are arrays since zw
    ! may vary spatially
    allocate(kmeuph(isc:iec, jsc:jec)); kmeuph(:,:)=1
    allocate(k100(isc:iec, jsc:jec)); k100(:,:)=1
    do j = jsc,jec; do i = isc,iec;
      nz = grid_kmt(i,j)
      do k = 1,nz
        if (wombat%zw(i,j,k) .le. 400) kmeuph(i,j)=k
        if (wombat%zw(i,j,k) .le. 100) k100(i,j)=k
      enddo
    enddo; enddo

    ! Get the timestep for the ecosystem model
    ts_npzd = max(1, nint(dt / wombat%dt_npzd)) ! number of ecosystem timesteps per model timestep
    rdtts = 1 / dt
    dtsb = dt / float(ts_npzd) ! number of seconds per nested ecosystem timestep

    ! Get the prognostic tracer values
    ! dts attn: we should probably update prognostic tracers via pointers to avoid
    ! having to allocate all these field arrays
    ! dts attn: do we really want/need to force these to be positive?
    call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]

    !-----------------------------------------------------------------------
    ! Available light and export production
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec
      ! dts: keep only shortwave wavelengths < 710 nm (taken from COBALT/BLING)
      do n = 1,nbands;
        if (max_wavelength_band(n) .lt. 710) then
          sw_pen(n) = max(0.0, sw_pen_band(n,i,j))
        else
          sw_pen(n) = 0.0
        endif
      enddo

      do k = 1,kmeuph(i,j)
        ! Shortwave intensity summed over bands and averaged over layer
        sw_pen_top = 0.0
        sw_pen_bottom = 0.0
        do n = 1,nbands;
            sw_pen_top = sw_pen_top + sw_pen(n) ! sw at top interface
            sw_pen(n) = sw_pen(n) * exp(-1.0 * opacity_band(n,i,j,k) * dzt(i,j,k))
            sw_pen_bottom = sw_pen_bottom + sw_pen(n) ! sw at bottom interface
        enddo
        wombat%radbio3d(i,j,k) = 0.5 * (sw_pen_top + sw_pen_bottom) ! [W/m2]

        wombat%vpbio(i,j,k) = wombat%abio * wombat%bbio ** (wombat%cbio * Temp(i,j,k)) ! [1/s]

        ! Growth term from Brian Griffiths
        wombat%avej(i,j,k) = wombat%vpbio(i,j,k) * &
            (1.0 - exp(-1.0 * (wombat%alphabio * wombat%radbio3d(i,j,k)) / wombat%vpbio(i,j,k))) ! [1/s]

        ! Calculate the average light limitation over the mixed layer
        if (wombat%zw(i,j,k) .le. hblt_depth(i,j)) &
          wombat%light_limit(i,j) = wombat%light_limit(i,j) + dzt(i,j,k) * &
              (1.0 - exp(-1.0 * (wombat%alphabio * wombat%radbio3d(i,j,k)) / wombat%vpbio(i,j,k))) ! [m]
      enddo

      ! dts: in WOMBAT v3, export production was calucated before updating the source terms so we do
      ! the same here
      k = k100(i,j)
      wombat%export_prod(i,j) = (wombat%Rho_0 * wombat%wdetbio) * wombat%f_det(i,j,k) ! [mol/m2/s]
      wombat%export_inorg(i,j) = (wombat%Rho_0 * wombat%wcaco3) * wombat%f_caco3(i,j,k) ! [mol/m2/s]
    enddo; enddo

    !-----------------------------------------------------------------------
    ! Calculate source terms using Euler forward timestepping
    !-----------------------------------------------------------------------
    ! chd: This is the NPZD model:
    !  (P: phytoplankton, Z: Zooplankton, N: Nitrate and D: Detritus)
    !  dP/dt = u(N,Temp.,Light) P - p_P P - g(P) P Z
    !  dZ/dt = a g(P) P Z - d Z - p_Z Z^2
    !  dN/dt = r D + d Z - u(N,Temp.,Light) P  [ + r_d DOC ]
    !  dD/dt = (1-s)[ (1-a) g(P) P Z + p_P P + p_Z Z^2] -r D + w_D dD/dz

    wombat%no3_prev(:,:,:) = wombat%f_no3(:,:,:)
    wombat%caco3_prev(:,:,:) = wombat%f_caco3(:,:,:)

    do tn = 1,ts_npzd
      do k = 1,nk; do j = jsc,jec; do i = isc,iec;
        ! chd: Use Liebigs Law of the Minimum (Liebig, 1845) for growth rate
        ! (minimum of light-limited and nutrient limited growth rates); although
        ! chlorophyll is not explicitly considered, this will later allow for a
        ! diagnostic determination of a Chl:N ratio depending on light- or
        ! nutrient-limited growth.

        ! Growth rate
        ! dts: convert k1bio from mmol/m3 to mol/kg
        u_npz = min(wombat%avej(i,j,k), wombat%vpbio(i,j,k) * wombat%f_no3(i,j,k) / &
            (mmol_m3_to_mol_kg * wombat%k1bio + wombat%f_no3(i,j,k))) ! [1/s]

        ! Iron limitation
        ! dts: convert 0.1 from umol/m3 to mol/kg
        u_npz = min(u_npz, wombat%vpbio(i,j,k) * wombat%f_fe(i,j,k) / &
            (umol_m3_to_mol_kg * 0.1 + wombat%f_fe(i,j,k))) ! [1/s]

        ! Grazing function
        ! dts: convert epsbio from m6/mmol2/s to kg2/mol2/s
        g_npz = wombat%gbio * &
            ((wombat%epsbio / mmol2_m6_to_mol2_kg2) * (wombat%f_phy(i,j,k) * wombat%f_phy(i,j,k))) / &
            (wombat%gbio + &
            (wombat%epsbio / mmol2_m6_to_mol2_kg2) * (wombat%f_phy(i,j,k) * wombat%f_phy(i,j,k))) ! [1/s]

        ! Temperature dependance of growth rates
        fbc = wombat%bbio ** (wombat%cbio * Temp(i,j,k)) ! [1]

        if (wombat%f_no3(i,j,k) .gt. epsi) then
          f11 = u_npz * wombat%f_phy(i,j,k) ! [mol/kg/s]
        else
          f11 = 0.0
        endif

        ! dts: convert muepsbio from m3/mmol/s to kg/mol/s
        if (wombat%f_phy(i,j,k) .gt. epsi) then
          f21 = g_npz * wombat%f_zoo(i,j,k) ! [mol/kg/s]
          f22 = wombat%muepbio * fbc * wombat%f_phy(i,j,k) ! [mol/kg/s]
          f23 = (wombat%muepsbio / mmol_m3_to_mol_kg) * wombat%f_phy(i,j,k) * wombat%f_phy(i,j,k) ! [mol/kg/s]
        else
          f21 = 0.0
          f22 = 0.0
          f23 = 0.0
        endif
        
        ! dts: convert muepsbio from m3/mmol/s to kg/mol/s
        if (wombat%f_zoo(i,j,k) .gt. epsi) then
          f31 = wombat%gam2bio * fbc * wombat%f_zoo(i,j,k) ! [mol/kg/s]
          f32 = (wombat%muezbio / mmol_m3_to_mol_kg) * wombat%f_zoo(i,j,k) * wombat%f_zoo(i,j,k) ! [mol/kg/s]
        else
          f31 = 0.0
          f32 = 0.0
        endif
        
        if (wombat%f_det(i,j,k) .gt. epsi) then
          f41 = wombat%muedbio * fbc * wombat%f_det(i,j,k) ! [mol/kg/s]

          if (wombat%zw(i,j,k) .ge. 180) f41 = 0.5 * f41 ! reduce decay below 180m
        else
          f41 = 0.0
        endif
        
        if (wombat%f_caco3(i,j,k) .gt. epsi) then
          f51 = wombat%muecaco3 * wombat%f_caco3(i,j,k) ! [mol/kg/s]
        else
          f51 = 0.0
        endif

        ! Nutrient equation
        !-----------------------------------------------------------------------
        wombat%f_no3(i,j,k) = wombat%f_no3(i,j,k) + dtsb * ((f41 + f31) + (f22 - f11)) ! [mol/kg]

        ! Phytoplankton equation
        !-----------------------------------------------------------------------
        wombat%f_phy(i,j,k)  = wombat%f_phy(i,j,k) + dtsb * ((f11 - f21) - (f22 + f23)) ! [mol/kg]

        ! Estimate primary productivity from phytoplankton growth
        wombat%pprod_gross(i,j,k) = wombat%pprod_gross(i,j,k) + dtsb * f11 ! [mol/kg]

        ! Net primary productivity (gross PP minus linear mortality)
        wombat%npp3d(i,j,k) = wombat%npp3d(i,j,k) + dtsb * (f11 - f21) ! [mol/kg]

        ! Zooplankton equation
        !-----------------------------------------------------------------------
        wombat%f_zoo(i,j,k)  = wombat%f_zoo(i,j,k)  + dtsb * &
            (wombat%gam1bio * f21 - (f31 + f32)) ! [mol/kg]

        ! Estimate secondary productivity from zooplankton growth
        wombat%zprod_gross(i,j,k) = wombat%zprod_gross(i,j,k) + dtsb * f21 ! [mol/kg]

        ! Detritus equation
        !-----------------------------------------------------------------------
        wombat%f_det(i,j,k) = wombat%f_det(i,j,k) + dtsb * &
            (((1 - wombat%gam1bio) * f21 + f23) + (f32 - f41)) ! [mol/kg]

        ! Oxygen equation
        !-----------------------------------------------------------------------
        if (wombat%f_o2(i,j,k) .gt. epsi) &
          wombat%f_o2(i,j,k) = wombat%f_o2(i,j,k) - 172. / 16. * dtsb * &
              ((f41 + f31) + (f22 - f11)) ! [mol/kg]

        ! Extra equation for caco3 - alkalinity
        !-----------------------------------------------------------------------
        ! rjm: 6.2% of POC by default 106/16*.062
        wombat%f_caco3(i,j,k) = wombat%f_caco3(i,j,k) + dtsb * &
            (((1 - wombat%gam1bio) * f21 + (f23 + f32)) * wombat%f_inorg * 106. / 16. - f51) ! [mol/kg]

        ! Extra equation for iron
        !-----------------------------------------------------------------------
        ! mac: molar Fe:N = 1.98e-5:1.0 (Christian et al. 2002)
        wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) + dtsb * 2.0e-5 * &
            ((f41 + f31) + (f22 - f11)) ! [mol/kg]

      enddo; enddo; enddo
    enddo

    ! Add biotically induced tendency to biotracers
    !-----------------------------------------------------------------------
    ! dts: update prognostic tracers via pointers
    call g_tracer_get_pointer(tracer_list, 'dic', 'field', wombat%p_dic) ! [mol/kg]
    call g_tracer_get_pointer(tracer_list, 'adic', 'field', wombat%p_adic) ! [mol/kg]
    call g_tracer_get_pointer(tracer_list, 'alk', 'field', wombat%p_alk) ! [mol/kg]

    do k = 1,nk; do j = jsc,jec; do i = isc,iec;
      no3_bgc_change = grid_tmask(i,j,k) * (wombat%f_no3(i,j,k) - wombat%no3_prev(i,j,k)) ! [mol/kg]
      caco3_bgc_change = grid_tmask(i,j,k) * (wombat%f_caco3(i,j,k) - wombat%caco3_prev(i,j,k)) ! [mol/kg]

      ! dts: convert fe_bkgnd from umol/m3 to mol/kg
      wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) - dt * wombat%tscav_fe * &
          max(0.0, (wombat%f_fe(i,j,k) - (umol_m3_to_mol_kg * wombat%fe_bkgnd))) ! [mol/kg]

      wombat%p_dic(i,j,k,tau) = wombat%p_dic(i,j,k,tau) + &
          (106. / 16. * no3_bgc_change - caco3_bgc_change) ! [mol/kg]

      wombat%p_adic(i,j,k,tau) = wombat%p_adic(i,j,k,tau) + &
          (106. / 16. * no3_bgc_change - caco3_bgc_change) ! [mol/kg]

      wombat%p_alk(i,j,k,tau) = wombat%p_alk(i,j,k,tau) + &
          (-2.0 * caco3_bgc_change - no3_bgc_change) ! [mol/kg]

      wombat%pprod_gross(i,j,k) = rdtts * wombat%pprod_gross(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%zprod_gross(i,j,k) = rdtts * wombat%zprod_gross(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]
      wombat%npp3d(i,j,k) = rdtts * wombat%npp3d(i,j,k) * grid_tmask(i,j,k) ! [mol/kg/s]

      if (wombat%zw(i,j,k) .le. hblt_depth(i,j)) then
        wombat%adic_intmld(i,j) = wombat%adic_intmld(i,j) + &
            wombat%p_adic(i,j,k,tau) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%dic_intmld(i,j) = wombat%dic_intmld(i,j) + wombat%p_dic(i,j,k,tau) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%o2_intmld(i,j)  = wombat%o2_intmld(i,j)  + wombat%f_o2(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%no3_intmld(i,j) = wombat%no3_intmld(i,j) + wombat%f_no3(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%fe_intmld(i,j)  = wombat%fe_intmld(i,j)  + wombat%f_fe(i,j,k)  * rho_dzt(i,j,k) ! [mol/m2]
        wombat%phy_intmld(i,j) = wombat%phy_intmld(i,j) + wombat%f_phy(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%det_intmld(i,j) = wombat%det_intmld(i,j) + wombat%f_det(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%pprod_gross_intmld(i,j) = wombat%pprod_gross_intmld(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%npp_intmld(i,j) = wombat%npp_intmld(i,j) + wombat%npp3d(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%radbio_intmld(i,j) = wombat%radbio_intmld(i,j) + wombat%radbio3d(i,j,k) * dzt(i,j,k) ! [W/m]
      endif

      if (wombat%zw(i,j,k) .le. 100) then
        wombat%adic_int100(i,j) = wombat%adic_int100(i,j) + &
            wombat%p_adic(i,j,k,tau) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%dic_int100(i,j) = wombat%dic_int100(i,j) + wombat%p_dic(i,j,k,tau) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%o2_int100(i,j)  = wombat%o2_int100(i,j)  + wombat%f_o2(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%no3_int100(i,j) = wombat%no3_int100(i,j) + wombat%f_no3(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%fe_int100(i,j)  = wombat%fe_int100(i,j)  + wombat%f_fe(i,j,k)  * rho_dzt(i,j,k) ! [mol/m2]
        wombat%phy_int100(i,j) = wombat%phy_int100(i,j) + wombat%f_phy(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%det_int100(i,j) = wombat%det_int100(i,j) + wombat%f_det(i,j,k) * rho_dzt(i,j,k) ! [mol/m2]
        wombat%pprod_gross_int100(i,j) = wombat%pprod_gross_int100(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%npp_int100(i,j) = wombat%npp_int100(i,j) + wombat%npp3d(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
        wombat%radbio_int100(i,j) = wombat%radbio_int100(i,j) + wombat%radbio3d(i,j,k) * dzt(i,j,k) ! [W/m]
      endif
    enddo; enddo; enddo

    ! Bottom iron fix
    !-----------------------------------------------------------------------
    ! mac: only apply this fix when the water is <= 200 m deep.  
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j) .gt. 0) then
        k = grid_kmt(i,j)
        if (wombat%zw(i,j,k) .le. 200) &
          wombat%f_fe(i,j,k)= umol_m3_to_mol_kg * 0.999 ! [mol/kg]
      endif
    enddo; enddo

    ! Set tracers values
    call g_tracer_set_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'phy', 'field', wombat%f_phy, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau)

    !-----------------------------------------------------------------------
    ! Remineralisation of sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment) ! [mol/m2]

    do j = jsc,jec; do i = isc,iec;
      k = grid_kmt(i,j)
      if (k .gt. 0) then
        fbc = wombat%bbio ** (wombat%cbio * Temp(i,j,k)) ! [1]
        wombat%det_sed_remin(i,j) = wombat%muedbio_sed * fbc * wombat%p_det_sediment(i,j,1) ! [mol/m2/s]
        wombat%caco3_sed_remin(i,j) = wombat%muecaco3_sed * fbc * wombat%p_caco3_sediment(i,j,1) ! [mol/m2/s]
        
        ! Remineralisation of sediments to supply nutrient fields.
        ! btf values are positive from the water column into the sediment.
        wombat%b_no3(i,j) = -1.0 * wombat%det_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_o2(i,j) = -172. / 16. * wombat%b_no3(i,j) ! [mol/m2/s]
        wombat%b_dic(i,j) = 106. / 16. * wombat%b_no3(i,j) - wombat%caco3_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_adic(i,j) = wombat%b_dic(i,j) ! [mol/m2/s]
        wombat%b_fe(i,j) = 2.0e-5 * wombat%b_no3(i,j) ! [mol/m2/s]
        wombat%b_alk(i,j) = -2.0 * wombat%caco3_sed_remin(i,j) - wombat%b_no3(i,j) ! [mol/m2/s]
      endif
    enddo; enddo

    ! Apply remineralisation rates to sediment tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j) .gt. 0) then
        wombat%p_det_sediment(i,j,1) = wombat%p_det_sediment(i,j,1) - dt * wombat%det_sed_remin(i,j) ! [mol/m2]
        wombat%p_caco3_sediment(i,j,1) = wombat%p_caco3_sediment(i,j,1) - &
            dt * wombat%caco3_sed_remin(i,j) ! [mol/m2]
      endif
    enddo; enddo

    call g_tracer_set_values(tracer_list, 'no3', 'btf', wombat%b_no3, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'btf', wombat%b_o2, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'btf', wombat%b_dic, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'btf', wombat%b_adic, isd, jsd)
    call g_tracer_set_values(tracer_list, 'fe', 'btf', wombat%b_fe, isd, jsd)
    call g_tracer_set_values(tracer_list, 'alk', 'btf', wombat%b_alk, isd, jsd)

    !=======================================================================
    ! Send diagnostics
    !=======================================================================

    if (wombat%id_pco2 .gt. 0) &
      used = g_send_data(wombat%id_pco2, wombat%pco2_csurf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_paco2 .gt. 0) &
      used = g_send_data(wombat%id_paco2, wombat%paco2_csurf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_htotal .gt. 0) &
      used = g_send_data(wombat%id_htotal, wombat%htotal, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_ahtotal .gt. 0) &
      used = g_send_data(wombat%id_ahtotal, wombat%ahtotal, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_no3_vstf .gt. 0) &
      used = g_send_data(wombat%id_no3_vstf, wombat%no3_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_vstf .gt. 0) &
      used = g_send_data(wombat%id_dic_vstf, wombat%dic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_adic_vstf .gt. 0) &
      used = g_send_data(wombat%id_adic_vstf, wombat%adic_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_alk_vstf .gt. 0) &
      used = g_send_data(wombat%id_alk_vstf, wombat%alk_vstf, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_light_limit .gt. 0) &
      used = g_send_data(wombat%id_light_limit, wombat%light_limit, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_radbio3d .gt. 0) &
      used = g_send_data(wombat%id_radbio3d, wombat%radbio3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radbio1 .gt. 0) &
      used = g_send_data(wombat%id_radbio1, wombat%radbio3d(:,:,1), model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross, wombat%pprod_gross, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pprod_gross_2d .gt. 0) then
      wombat%pprod_gross_2d = 0.0
      do k = 1,nk; do j = jsc,jec; do i = isc,iec;
        wombat%pprod_gross_2d(i,j) = wombat%pprod_gross_2d(i,j) + &
            wombat%pprod_gross(i,j,k) * rho_dzt(i,j,k) ! [mol/m2/s]
      enddo; enddo; enddo
      used = g_send_data(wombat%id_pprod_gross_2d, wombat%pprod_gross_2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_wdet100 .gt. 0) then
      do j = jsc,jec; do i = isc,iec;
        wombat%wdet100(i,j) = (wombat%Rho_0 * wombat%wdetbio) * wombat%f_det(i,j, &
            minloc(abs(wombat%zm(i,j,:)-100), dim=1)) ! [mol/m2/s]
      enddo; enddo
      used = g_send_data(wombat%id_wdet100, wombat%wdet100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_export_prod .gt. 0) &
      used = g_send_data(wombat%id_export_prod, wombat%export_prod, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_export_inorg .gt. 0) &
      used = g_send_data(wombat%id_export_inorg, wombat%export_inorg, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp3d .gt. 0) &
      used = g_send_data(wombat%id_npp3d, wombat%npp3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_npp2d .gt. 0) then
      wombat%npp2d = 0.0
      do k = 1,nk
        wombat%npp2d(isc:iec,jsc:jec) = wombat%npp2d(isc:iec,jsc:jec) + &
            wombat%npp3d(isc:iec,jsc:jec,k) * rho_dzt(isc:iec,jsc:jec,k) ! [mol/m2/s]
      enddo
      used = g_send_data(wombat%id_npp2d, wombat%npp2d, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if (wombat%id_npp1 .gt. 0) &
      used = g_send_data(wombat%id_npp1, wombat%npp3d(:,:,1), model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_zprod_gross .gt. 0) &
      used = g_send_data(wombat%id_zprod_gross, wombat%zprod_gross, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_adic_intmld .gt. 0) &
      used = g_send_data(wombat%id_adic_intmld, wombat%adic_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_intmld .gt. 0) &
      used = g_send_data(wombat%id_dic_intmld, wombat%dic_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_o2_intmld .gt. 0) &
      used = g_send_data(wombat%id_o2_intmld, wombat%o2_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_no3_intmld .gt. 0) &
      used = g_send_data(wombat%id_no3_intmld, wombat%no3_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fe_intmld .gt. 0) &
      used = g_send_data(wombat%id_fe_intmld, wombat%fe_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_phy_intmld .gt. 0) &
      used = g_send_data(wombat%id_phy_intmld, wombat%phy_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_intmld .gt. 0) &
      used = g_send_data(wombat%id_det_intmld, wombat%det_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross_intmld .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross_intmld, wombat%pprod_gross_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp_intmld .gt. 0) &
      used = g_send_data(wombat%id_npp_intmld, wombat%npp_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_radbio_intmld .gt. 0) &
      used = g_send_data(wombat%id_radbio_intmld, wombat%radbio_intmld, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_adic_int100 .gt. 0) &
      used = g_send_data(wombat%id_adic_int100, wombat%adic_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_dic_int100 .gt. 0) &
      used = g_send_data(wombat%id_dic_int100, wombat%dic_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_o2_int100 .gt. 0) &
      used = g_send_data(wombat%id_o2_int100, wombat%o2_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_no3_int100 .gt. 0) &
      used = g_send_data(wombat%id_no3_int100, wombat%no3_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_fe_int100 .gt. 0) &
      used = g_send_data(wombat%id_fe_int100, wombat%fe_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_phy_int100 .gt. 0) &
      used = g_send_data(wombat%id_phy_int100, wombat%phy_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_int100 .gt. 0) &
      used = g_send_data(wombat%id_det_int100, wombat%det_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross_int100 .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross_int100, wombat%pprod_gross_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_npp_int100 .gt. 0) &
      used = g_send_data(wombat%id_npp_int100, wombat%npp_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_radbio_int100 .gt. 0) &
      used = g_send_data(wombat%id_radbio_int100, wombat%radbio_int100, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_det_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_det_sed_remin, wombat%det_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_caco3_sed_remin, wombat%caco3_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    deallocate(kmeuph, k100)

  end subroutine generic_WOMBATlite_update_from_source

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Calculate and set coupler values at the surface / bottom of the ocean.
  !   User must provide the calculations for these boundary values.
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_set_boundary_values(tracer_list, SST, SSS, rho, ilb, jlb, tau, dzt)
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !
  !  <IN NAME="SST" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Temperature
  !  </IN>
  !
  !  <IN NAME="SSS" TYPE="real, dimension(ilb:,jlb:)">
  !   Sea Surface Salinity
  !  </IN>
  !
  !  <IN NAME="rho" TYPE="real, dimension(ilb:,jlb:,:,:)">
  !   Ocean density
  !  </IN>
  !
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Layer thickness
  !  </IN>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_set_boundary_values(tracer_list, SST, SSS, rho, ilb, jlb, tau, dzt)
    type(g_tracer_type), pointer                       :: tracer_list
    real, dimension(ilb:,jlb:), intent(in)             :: SST, SSS
    real, dimension(ilb:,jlb:,:,:), intent(in)         :: rho
    integer, intent(in)                                :: ilb, jlb, tau
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt

    integer                                 :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, i, j
    real                                    :: sal, ST, o2_saturation
    real                                    :: tt, tk, ts, ts2, ts3, ts4, ts5
    real, dimension(:,:,:), pointer         :: grid_tmask
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATlite_set_boundary_values'

    ! Get the necessary properties
    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
        grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list, 'o2', 'field', wombat%p_o2)
    
    ! nnz: Since the generic_WOMBATlite_update_from_source() subroutine is called by this time
    ! the following if block is not really necessary (since this calculation is already done in
    ! source).
    ! It is only neccessary if source routine is commented out for debugging.
    ! Note: In order for this to work we should NOT zero out the coupler values for generic tracers
    ! This zeroing is done for non-generic TOPAZ by calling zero_ocean_sfc.
    ! Since the coupler values here are non-cumulative there is no need to zero them out anyway.
    if (wombat%init .OR. wombat%force_update_fluxes) then
      ! Get necessary fields
      call g_tracer_get_values(tracer_list, 'dic', 'field', wombat%f_dic, isd, jsd, ntau=1, &
          positive=.true.)
      call g_tracer_get_values(tracer_list, 'adic', 'field', wombat%f_adic, isd, jsd, ntau=1, &
          positive=.true.)
      call g_tracer_get_values(tracer_list, 'no3', 'field', wombat%f_no3, isd, jsd, ntau=1, &
          positive=.true.)
      call g_tracer_get_values(tracer_list, 'alk', 'field', wombat%f_alk, isd, jsd, ntau=1, &
          positive=.true.)

      do j = jsc, jec; do i = isc, iec
          wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j)
          wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j)
          wombat%ahtotallo(i,j) = wombat%htotal_scale_lo * wombat%ahtotal(i,j)
          wombat%ahtotalhi(i,j) = wombat%htotal_scale_hi * wombat%ahtotal(i,j)
      enddo; enddo 

      if ((trim(co2_calc) == 'mocsy') .and. (.not. present(dzt))) then
        call mpp_error(FATAL,"mocsy method of co2_calc needs dzt to be passed to the "// &
            "FMS_ocmip2_co2calc subroutine.")
      endif

      call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
          SST(:,:), SSS(:,:), &
          wombat%f_dic(:,:,1), &
          wombat%f_no3(:,:,1) / 16., &
          wombat%sio2(:,:), &
          wombat%f_alk(:,:,1), &
          wombat%htotallo(:,:), wombat%htotalhi(:,:), &
          wombat%htotal(:,:), &
          co2_calc=trim(co2_calc), &
          zt=dzt(:,:,1), &
          co2star=wombat%co2_csurf(:,:), alpha=wombat%co2_alpha(:,:), &
          pCO2surf=wombat%pco2_csurf(:,:))

      call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,1), &
          SST(:,:), SSS(:,:), &
          wombat%f_adic(:,:,1), &
          wombat%f_no3(:,:,1) / 16., &
          wombat%sio2(:,:), &
          wombat%f_alk(:,:,1), &
          wombat%ahtotallo(:,:), wombat%ahtotalhi(:,:), &
          wombat%ahtotal(:,:), &
          co2_calc=trim(co2_calc), &
          zt=dzt(:,:,1), &
          co2star=wombat%aco2_csurf(:,:), alpha=wombat%aco2_alpha(:,:), &
          pCO2surf=wombat%paco2_csurf(:,:))

      call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
      call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)
      call g_tracer_set_values(tracer_list, 'adic', 'alpha', wombat%aco2_alpha, isd, jsd)
      call g_tracer_set_values(tracer_list, 'adic', 'csurf', wombat%aco2_csurf, isd, jsd)

      ! nnz: If source is called uncomment the following
      wombat%init = .false. !nnz: This is necessary since the above calls appear in source subroutine too.
    endif

    call g_tracer_get_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha ,isd, jsd)
    call g_tracer_get_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf ,isd, jsd)
    call g_tracer_get_values(tracer_list, 'adic', 'alpha', wombat%aco2_alpha ,isd, jsd)
    call g_tracer_get_values(tracer_list, 'adic', 'csurf', wombat%aco2_csurf ,isd, jsd)

    do j=jsc,jec ; do i=isc,iec
      !-----------------------------------------------------------------------
      ! Compute the Schmidt number of CO2 in seawater using the formulation
      ! presented by Wanninkhof (1992, J. Geophys. Res., 97, 7373-7382).
      !-----------------------------------------------------------------------
      ST = SST(i,j)
      wombat%co2_sc_no(i,j) = wombat%a1_co2 + ST*(wombat%a2_co2 + ST*(wombat%a3_co2 + &
          ST*wombat%a4_co2)) * grid_tmask(i,j,1)
      
      wombat%co2_alpha(i,j) = wombat%co2_alpha(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
      wombat%co2_csurf(i,j) = wombat%co2_csurf(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
      wombat%aco2_alpha(i,j) = wombat%aco2_alpha(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
      wombat%aco2_csurf(i,j) = wombat%aco2_csurf(i,j) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'dic', 'sc_no', wombat%co2_sc_no, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'alpha', wombat%aco2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'csurf', wombat%aco2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'adic', 'sc_no', wombat%co2_sc_no, isd, jsd)

    call g_tracer_get_values(tracer_list, 'o2', 'alpha', wombat%o2_alpha, isd, jsd)
    call g_tracer_get_values(tracer_list, 'o2', 'csurf', wombat%o2_csurf ,isd, jsd)

    do j=jsc,jec ; do i=isc,iec
      !-----------------------------------------------------------------------
      ! Compute the oxygen saturation concentration at 1 atm total
      ! pressure in mol/kg given the temperature (t, in deg C) and
      ! the salinity (s, in permil)
      !
      ! From Garcia and Gordon (1992), Limnology and Oceonography.
      ! The formula used is from page 1310, eq (8).
      !
      ! *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
      ! *** It shouldn't be there.                                ***
      !
      ! o2_saturation is defined between T(freezing) <= T <= 40 deg C and
      !                                   0 permil <= S <= 42 permil
      ! We impose these bounds here.
      !
      ! check value: T = 10 deg C, S = 35 permil,
      !              o2_saturation = 0.282015 mol m-3
      !-----------------------------------------------------------------------
      sal = SSS(i,j) ; ST = SST(i,j)

      ! jgj 2015/05/14 impose temperature and salinity bounds for o2sat
      sal = min(42.0, max(0.0, sal))
      tt = 298.15 - min(40.0, max(0.0, ST))
      tk = 273.15 + min(40.0, max(0.0, ST))
      ts = log(tt / tk)
      ts2 = ts  * ts
      ts3 = ts2 * ts
      ts4 = ts3 * ts
      ts5 = ts4 * ts

      o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
          exp( wombat%a_0 + wombat%a_1*ts + wombat%a_2*ts2 + wombat%a_3*ts3 + wombat%a_4*ts4 + &
              wombat%a_5*ts5 + (wombat%b_0 + wombat%b_1*ts + wombat%b_2*ts2 + wombat%b_3*ts3 + &
              wombat%c_0*sal)*sal)

      !-----------------------------------------------------------------------
      !  Compute the Schmidt number of O2 in seawater using the
      !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
      !  Cycles, 12, 141-163).
      !-----------------------------------------------------------------------
      wombat%o2_sc_no(i,j)  = wombat%a1_o2 + ST * (wombat%a2_o2 + ST * (wombat%a3_o2 + ST * &
          wombat%a4_o2 )) * grid_tmask(i,j,1)

      ! renormalize the alpha value for atm o2
      ! data table override for o2_flux_pcair_atm is now set to 0.21
      wombat%o2_alpha(i,j) = (o2_saturation / 0.21)
      wombat%o2_csurf(i,j) = wombat%p_o2(i,j,1,tau) * wombat%Rho_0 !nnz: MOM has rho(i,j,1,tau)
    enddo; enddo

    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    call g_tracer_set_values(tracer_list, 'o2', 'alpha', wombat%o2_alpha, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'csurf', wombat%o2_csurf, isd, jsd)
    call g_tracer_set_values(tracer_list, 'o2', 'sc_no', wombat%o2_sc_no, isd, jsd)

  end subroutine generic_WOMBATlite_set_boundary_values

  !#######################################################################
  !
  ! <SUBROUTINE NAME="generic_WOMBATlite_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_WOMBATlite_end
  !  </TEMPLATE>
  ! </SUBROUTINE>
  !
  subroutine generic_WOMBATlite_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_WOMBATlite_end'

    call user_deallocate_arrays

  end subroutine generic_WOMBATlite_end

  !#######################################################################
  !
  ! This is an internal sub, not a public interface.
  ! Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau)

    ! Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc; CO2_dope_vec%iec = iec
    CO2_dope_vec%jsc = jsc; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd; CO2_dope_vec%jed = jed

    allocate(wombat%htotallo(isd:ied, jsd:jed)); wombat%htotallo(:,:)=0.0
    allocate(wombat%htotalhi(isd:ied, jsd:jed)); wombat%htotalhi(:,:)=0.0
    allocate(wombat%htotal(isd:ied, jsd:jed)); wombat%htotal(:,:)=wombat%htotal_in
    allocate(wombat%ahtotallo(isd:ied, jsd:jed)); wombat%ahtotallo(:,:)=0.0
    allocate(wombat%ahtotalhi(isd:ied, jsd:jed)); wombat%ahtotalhi(:,:)=0.0
    allocate(wombat%ahtotal(isd:ied, jsd:jed)); wombat%ahtotal(:,:)=wombat%htotal_in
    allocate(wombat%sio2(isd:ied, jsd:jed)); wombat%sio2(:,:)=wombat%sio2_surf
    allocate(wombat%co2_csurf(isd:ied, jsd:jed)); wombat%co2_csurf(:,:)=0.0
    allocate(wombat%co2_alpha(isd:ied, jsd:jed)); wombat%co2_alpha(:,:)=0.0
    allocate(wombat%co2_sc_no(isd:ied, jsd:jed)); wombat%co2_sc_no(:,:)=0.0
    allocate(wombat%pco2_csurf(isd:ied, jsd:jed)); wombat%pco2_csurf(:,:)=0.0
    allocate(wombat%aco2_csurf(isd:ied, jsd:jed)); wombat%aco2_csurf(:,:)=0.0
    allocate(wombat%aco2_alpha(isd:ied, jsd:jed)); wombat%aco2_alpha(:,:)=0.0
    allocate(wombat%paco2_csurf(isd:ied, jsd:jed)); wombat%paco2_csurf(:,:)=0.0
    allocate(wombat%o2_csurf(isd:ied, jsd:jed)); wombat%o2_csurf(:,:)=0.0
    allocate(wombat%o2_alpha(isd:ied, jsd:jed)); wombat%o2_alpha(:,:)=0.0
    allocate(wombat%o2_sc_no(isd:ied, jsd:jed)); wombat%o2_sc_no(:,:)=0.0
    allocate(wombat%no3_vstf(isd:ied, jsd:jed)); wombat%no3_vstf(:,:)=0.0
    allocate(wombat%dic_vstf(isd:ied, jsd:jed)); wombat%dic_vstf(:,:)=0.0
    allocate(wombat%adic_vstf(isd:ied, jsd:jed)); wombat%adic_vstf(:,:)=0.0
    allocate(wombat%alk_vstf(isd:ied, jsd:jed)); wombat%alk_vstf(:,:)=0.0

    allocate(wombat%f_dic(isd:ied, jsd:jed, 1:nk)); wombat%f_dic(:,:,:)=0.0
    allocate(wombat%f_adic(isd:ied, jsd:jed, 1:nk)); wombat%f_adic(:,:,:)=0.0
    allocate(wombat%f_alk(isd:ied, jsd:jed, 1:nk)); wombat%f_alk(:,:,:)=0.0
    allocate(wombat%f_no3(isd:ied, jsd:jed, 1:nk)); wombat%f_no3(:,:,:)=0.0
    allocate(wombat%f_phy(isd:ied, jsd:jed, 1:nk)); wombat%f_phy(:,:,:)=0.0
    allocate(wombat%f_zoo(isd:ied, jsd:jed, 1:nk)); wombat%f_zoo(:,:,:)=0.0
    allocate(wombat%f_det(isd:ied, jsd:jed, 1:nk)); wombat%f_det(:,:,:)=0.0
    allocate(wombat%f_o2(isd:ied, jsd:jed, 1:nk)); wombat%f_o2(:,:,:)=0.0
    allocate(wombat%f_caco3(isd:ied, jsd:jed, 1:nk)); wombat%f_caco3(:,:,:)=0.0
    allocate(wombat%f_fe(isd:ied, jsd:jed, 1:nk)); wombat%f_fe(:,:,:)=0.0

    allocate(wombat%b_no3(isd:ied, jsd:jed)); wombat%b_no3(:,:)=0.0
    allocate(wombat%b_o2(isd:ied, jsd:jed)); wombat%b_o2(:,:)=0.0
    allocate(wombat%b_dic(isd:ied, jsd:jed)); wombat%b_dic(:,:)=0.0
    allocate(wombat%b_adic(isd:ied, jsd:jed)); wombat%b_adic(:,:)=0.0
    allocate(wombat%b_fe(isd:ied, jsd:jed)); wombat%b_fe(:,:)=0.0
    allocate(wombat%b_alk(isd:ied, jsd:jed)); wombat%b_alk(:,:)=0.0

    allocate(wombat%light_limit(isd:ied, jsd:jed)); wombat%light_limit(:,:)=0.0
    allocate(wombat%radbio3d(isd:ied, jsd:jed, 1:nk)); wombat%radbio3d(:,:,:)=0.0
    allocate(wombat%pprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%pprod_gross(:,:,:)=0.0
    allocate(wombat%pprod_gross_2d(isd:ied, jsd:jed)); wombat%pprod_gross_2d(:,:)=0.0
    allocate(wombat%zprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%zprod_gross(:,:,:)=0.0
    allocate(wombat%wdet100(isd:ied, jsd:jed)); wombat%wdet100(:,:)=0.0
    allocate(wombat%export_prod(isd:ied, jsd:jed)); wombat%export_prod(:,:)=0.0
    allocate(wombat%export_inorg(isd:ied, jsd:jed)); wombat%export_inorg(:,:)=0.0
    allocate(wombat%npp2d(isd:ied, jsd:jed)); wombat%npp2d(:,:)=0.0
    allocate(wombat%npp3d(isd:ied, jsd:jed, 1:nk)); wombat%npp3d(:,:,:)=0.0
    allocate(wombat%vpbio(isd:ied, jsd:jed, 1:nk)); wombat%vpbio(:,:,:)=0.0
    allocate(wombat%avej(isd:ied, jsd:jed, 1:nk)); wombat%avej(:,:,:)=0.0
    allocate(wombat%no3_prev(isd:ied, jsd:jed, 1:nk)); wombat%no3_prev(:,:,:)=0.0
    allocate(wombat%caco3_prev(isd:ied, jsd:jed, 1:nk)); wombat%caco3_prev(:,:,:)=0.0
    allocate(wombat%det_sed_remin(isd:ied, jsd:jed)); wombat%det_sed_remin(:,:)=0.0
    allocate(wombat%det_btm(isd:ied, jsd:jed)); wombat%det_btm(:,:)=0.0
    allocate(wombat%caco3_sed_remin(isd:ied, jsd:jed)); wombat%caco3_sed_remin(:,:)=0.0
    allocate(wombat%caco3_btm(isd:ied, jsd:jed)); wombat%caco3_btm(:,:)=0.0
    allocate(wombat%zw(isd:ied, jsd:jed, 1:nk)); wombat%zw(:,:,:)=0.0
    allocate(wombat%zm(isd:ied, jsd:jed, 1:nk)); wombat%zm(:,:,:)=0.0

    allocate(wombat%dic_intmld(isd:ied, jsd:jed)); wombat%dic_intmld(:,:)=0.0
    allocate(wombat%adic_intmld(isd:ied, jsd:jed)); wombat%adic_intmld(:,:)=0.0
    allocate(wombat%o2_intmld(isd:ied, jsd:jed)); wombat%o2_intmld(:,:)=0.0
    allocate(wombat%no3_intmld(isd:ied, jsd:jed)); wombat%no3_intmld(:,:)=0.0
    allocate(wombat%fe_intmld(isd:ied, jsd:jed)); wombat%fe_intmld(:,:)=0.0
    allocate(wombat%phy_intmld(isd:ied, jsd:jed)); wombat%phy_intmld(:,:)=0.0
    allocate(wombat%det_intmld(isd:ied, jsd:jed)); wombat%det_intmld(:,:)=0.0
    allocate(wombat%pprod_gross_intmld(isd:ied, jsd:jed)); wombat%pprod_gross_intmld(:,:)=0.0
    allocate(wombat%npp_intmld(isd:ied, jsd:jed)); wombat%npp_intmld(:,:)=0.0
    allocate(wombat%radbio_intmld(isd:ied, jsd:jed)); wombat%radbio_intmld(:,:)=0.0

    allocate(wombat%dic_int100(isd:ied, jsd:jed)); wombat%dic_int100(:,:)=0.0
    allocate(wombat%adic_int100(isd:ied, jsd:jed)); wombat%adic_int100(:,:)=0.0
    allocate(wombat%o2_int100(isd:ied, jsd:jed)); wombat%o2_int100(:,:)=0.0
    allocate(wombat%no3_int100(isd:ied, jsd:jed)); wombat%no3_int100(:,:)=0.0
    allocate(wombat%fe_int100(isd:ied, jsd:jed)); wombat%fe_int100(:,:)=0.0
    allocate(wombat%phy_int100(isd:ied, jsd:jed)); wombat%phy_int100(:,:)=0.0
    allocate(wombat%det_int100(isd:ied, jsd:jed)); wombat%det_int100(:,:)=0.0
    allocate(wombat%pprod_gross_int100(isd:ied, jsd:jed)); wombat%pprod_gross_int100(:,:)=0.0
    allocate(wombat%npp_int100(isd:ied, jsd:jed)); wombat%npp_int100(:,:)=0.0
    allocate(wombat%radbio_int100(isd:ied, jsd:jed)); wombat%radbio_int100(:,:)=0.0

  end subroutine user_allocate_arrays

  !#######################################################################
  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays

    deallocate( &
        wombat%htotallo, &
        wombat%htotalhi, &
        wombat%htotal, &
        wombat%ahtotallo, &
        wombat%ahtotalhi, &
        wombat%ahtotal, &
        wombat%co2_csurf, &
        wombat%co2_alpha, &
        wombat%co2_sc_no, &
        wombat%pco2_csurf, &
        wombat%aco2_csurf, &
        wombat%aco2_alpha, &
        wombat%paco2_csurf, &
        wombat%o2_csurf, &
        wombat%o2_alpha, &
        wombat%o2_sc_no, &
        wombat%no3_vstf, &
        wombat%dic_vstf, &
        wombat%adic_vstf, &
        wombat%alk_vstf)

    deallocate( &
        wombat%f_dic, &
        wombat%f_adic, &
        wombat%f_alk, &
        wombat%f_no3, &
        wombat%f_phy, &
        wombat%f_zoo, &
        wombat%f_det, &
        wombat%f_o2, &
        wombat%f_caco3, &
        wombat%f_fe)

    deallocate( &
        wombat%b_no3, &
        wombat%b_o2, &
        wombat%b_dic, &
        wombat%b_adic, &
        wombat%b_fe, &
        wombat%b_alk)

    deallocate( &
        wombat%light_limit, &
        wombat%radbio3d, &
        wombat%pprod_gross, &
        wombat%pprod_gross_2d, &
        wombat%zprod_gross, &
        wombat%wdet100, &
        wombat%export_prod, &
        wombat%export_inorg, &
        wombat%npp2d, &
        wombat%npp3d, &
        wombat%vpbio, &
        wombat%avej, &
        wombat%no3_prev, &
        wombat%caco3_prev, &
        wombat%det_sed_remin, &
        wombat%det_btm, &
        wombat%caco3_sed_remin, &
        wombat%caco3_btm, &
        wombat%zw, &
        wombat%zm)

    deallocate( &
        wombat%dic_intmld, &
        wombat%adic_intmld, &
        wombat%o2_intmld, &
        wombat%no3_intmld, &
        wombat%fe_intmld, &
        wombat%phy_intmld, &
        wombat%det_intmld, &
        wombat%pprod_gross_intmld, &
        wombat%npp_intmld, &
        wombat%radbio_intmld, &
        wombat%dic_int100, &
        wombat%adic_int100, &
        wombat%o2_int100, &
        wombat%no3_int100, &
        wombat%fe_int100, &
        wombat%phy_int100, &
        wombat%det_int100, &
        wombat%pprod_gross_int100, &
        wombat%npp_int100, &
        wombat%radbio_int100)

  end subroutine user_deallocate_arrays

end module generic_WOMBATlite
