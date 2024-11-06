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
!
!
!     (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
!     / o o \  : :.-.: :: ,. :: `' :: .; :: .; :`-. .-'
!    (   "   ) : :: :: :: :: :: .. ::   .':    :  : :  
!     \__ __/  : `' `' ;: :; :: :; :: .; :: :: :  : :  
!               `.,`.,' `.__.':_;:_;:___.':_;:_;  :_;     
!
!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
!
!  World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT) is
!  based on a NPZD (nutrient–phytoplankton–zooplankton–detritus) model.
!  This is the "lite" version of WOMBAT which includes one class each of
!  phytoplankton, zooplankton and sinking detritus, as well as nitrate
!  (NO3),  bio-available iron (Fe), dissolved inorganic carbon (DIC), 
!  calcium carbonate (CaCO3), alkalinity (ALK), and oxygen (O2). Fe is 
!  carried through the zooplankton and detrital pools as well.
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
        abioa, &
        bbioa, &
        bbioh, &
        phykn, &
        phykf, &
        phyminqc, &
        phyoptqc, &
        phyoptqf, &
        phymaxqf, &
        phylmor, &
        phyqmor, &
        zooassi, &
        zookz, &
        zoogmax, &
        zooepsi, &
        zoolmor, &
        zooqmor, &
        detlrem, &
        detlrem_sed, &
        wdetbio, &
        caco3lrem, &
        caco3lrem_sed, &
        wcaco3, &
        f_inorg, &
        tscav_fe, &
        ligand, &
        knano_dfe, &
        kscav_dfe, &
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
        htotallo, htotalhi,  &
        ahtotallo, ahtotalhi,  &
        sio2, &
        co2_csurf, co2_alpha, co2_sc_no, pco2_csurf, &
        aco2_csurf, aco2_alpha, paco2_csurf, &
        o2_csurf, o2_alpha, o2_sc_no, &
        no3_vstf, dic_vstf, adic_vstf, alk_vstf

    real, dimension(:,:,:), allocatable :: &
        htotal, ahtotal

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
        detfe_btm, &
        caco3_btm, &
        det_sed_remin, &
        detfe_sed_remin, &
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
        radbio_int100, &
        zeuphot, &
        phy_leup
      
    real, dimension(:,:,:), allocatable :: &
        f_dic, &
        f_adic, &
        f_alk, &
        f_no3, &
        f_phy, &
        f_pchl, &
        f_phyfe, &
        f_zoo, &
        f_zoofe, &
        f_det, &
        f_detfe, &
        f_o2, &
        f_caco3, &
        f_fe, &
        pprod_gross, &
        zprod_gross, &
        radbio3d, &
        radmid3d, &
        npp3d, &
        phy_mumax, &
        phy_mu, &
        pchl_mu, &
        phy_lpar, &
        phy_lnit, &
        phy_lfer, &
        phy_dfeupt, &
        feIII, &
        felig, &
        fecol, &
        feprecip, &
        fescaven, &
        fescadet, &
        fecoag2det, &
        fesources, &
        fesinks, &
        phygrow, &
        phyresp, &
        phymort, &
        zoograz, &
        zooresp, &
        zoomort, &
        zooexcr, &
        zooslop, &
        detremi, &
        caldiss, &
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
        p_detfe_sediment, &
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
        id_radmid3d = -1, &
        id_radbio1 = -1, &
        id_pprod_gross = -1, &
        id_phy_lpar = -1, &
        id_phy_lnit = -1, &
        id_phy_lfer = -1, &
        id_phy_dfeupt = -1, &
        id_feIII = -1, &
        id_felig = -1, &
        id_fecol = -1, &
        id_feprecip = -1, &
        id_fescaven = -1, &
        id_fescadet = -1, &
        id_fecoag2det = -1, &
        id_fesources = -1, &
        id_fesinks = -1, &
        id_phygrow = -1, &
        id_phyresp = -1, &
        id_phymort = -1, &
        id_zoograz = -1, &
        id_zooresp = -1, &
        id_zoomort = -1, &
        id_zooexcr = -1, &
        id_zooslop = -1, &
        id_detremi = -1, &
        id_caldiss = -1, &
        id_phy_mumax = -1, &
        id_phy_mu = -1, &
        id_pchl_mu = -1, &
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
        id_detfe_sed_remin = -1, &
        id_detfe_sed_depst = -1, &
        id_caco3_sed_remin = -1, &
        id_caco3_sed_depst = -1, &
        id_zeuphot = -1, &
        id_phy_leup = -1

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
        'htotal', 'H+ concentration', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'ahtotal', 'H+ concentration inc. anthropogenic', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_ahtotal = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'radbio3d', 'Photosynthetically active radiation available for phytoplankton growth', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radbio3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)
    
    vardesc_temp = vardesc( &
        'radmid3d', 'Photosynthetically active radiation at centre point of grid cell', &
        'h', 'L', 's', 'W m-2', 'f')
    wombat%id_radmid3d = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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
        'detfe_sed_remin', 'Rate of remineralisation of detrital iron in accumulated sediment', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_detfe_sed_remin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detfe_sed_depst', 'Rate of deposition of detrital iron to sediment at base of water column', &
        'h', '1', 's', 'mol/m^2/s', 'f')
    wombat%id_detfe_sed_depst = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
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
        'phy_mumax', 'Maximum growth rate of phytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_phy_mumax = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_mu', 'Realised growth rate of phytoplankton', 'h', 'L', 's', '/s', 'f')
    wombat%id_phy_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'pchl_mu', 'Realised growth rate of phytoplankton chlorophyll', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_pchl_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lpar', 'Limitation of phytoplankton by light', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lpar = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lnit', 'Limitation of phytoplankton by nitrogen', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lnit = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_lfer', 'Limitation of phytoplankton by iron', 'h', 'L', 's', '[0-1]', 'f')
    wombat%id_phy_lfer = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_dfeupt', 'Uptake of dFe by phytoplankton', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_phy_dfeupt = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'feIII', 'free iron (Fe3+)', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_feIII = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'felig', 'ligand-bound dissolved iron', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_felig = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecol', 'Colloidal dissolved iron', 'h', 'L', 's', 'mol/kg', 'f')
    wombat%id_fecol = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'feprecip', 'Precipitation of free Fe onto nanoparticles', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_feprecip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescaven', 'Scavenging of free Fe onto detritus (organic + inorganic)', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescaven = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fescadet', 'Scavenging of free Fe onto organic detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fescadet = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fecoag2det', 'Coagulation of colloidal dFe onto detritus', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fecoag2det = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fesources', 'Total source of dFe in water column', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fesources = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'fesinks', 'Total sink of dFe in water column', 'h', 'L', 's', 'mol/kg/s', 'f')
    wombat%id_fesinks = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phygrow', 'Growth of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phygrow = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phyresp', 'Respiration of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phyresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phymort', 'Mortality of phytoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_phymort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoograz', 'Grazing rate of zootoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoograz = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooresp', 'Respiration of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooresp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zoomort', 'Mortality of zooplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zoomort = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooexcr', 'Excretion rate of zootoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooexcr = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'zooslop', 'Sloppy feeding of zootoplankton', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_zooslop = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'detremi', 'Remineralisation rate of detritus', 'h', 'L', 's', 'molC/kg/s', 'f')
    wombat%id_detremi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'caldiss', 'Dissolution rate of CaCO3', 'h', 'L', 's', 'molCaCO3/kg/s', 'f')
    wombat%id_caldiss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
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

    vardesc_temp = vardesc( &
        'zeuphot', 'Depth of the euphotic zone (defined as 1% of incident light)', &
        'h', '1', 's', 'm', 'f')
    wombat%id_zeuphot = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
        init_time, vardesc_temp%longname, vardesc_temp%units, missing_value=missing_value1)

    vardesc_temp = vardesc( &
        'phy_leup', 'Limitation factor due to depth of euphotic zone / MLD', &
        'h', '1', 's', 'm', 'f')
    wombat%id_phy_leup = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
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

    ! Initial slope of P-I curve [(mg Chl m-3)-1 (W m-2)-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('alphabio', wombat%alphabio, 2.25)

    ! Autotrophy maximum growth rate parameter a [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('abioa', wombat%abioa, 1.0/86400.0)

    ! Autotrophy maximum growth rate parameter b [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioa', wombat%bbioa, 1.050)

    ! Heterotrophy maximum growth rate parameter b [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('bbioh', wombat%bbioh, 1.066)

    ! Phytoplankton half saturation constant for nitrogen uptake [mmolN/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykn', wombat%phykn, 2.0)

    ! Phytoplankton half saturation constant for iron uptake [umolFe/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phykf', wombat%phykf, 2.5)

    ! Phytoplankton minimum quota of chlorophyll to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyminqc', wombat%phyminqc, 0.004)

    ! Phytoplankton optimal quota of chlorophyll to carbon [mg/mg]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyoptqc', wombat%phyoptqc, 0.036)

    ! Phytoplankton optimal quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyoptqf', wombat%phyoptqf, 10e-6)

    ! Phytoplankton maximum quota of iron to carbon [mol/mol]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phymaxqf', wombat%phymaxqf, 50e-6)

    ! Phytoplankton linear mortality rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phylmor', wombat%phylmor, 0.005/86400.0)

    ! Phytoplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('phyqmor', wombat%phyqmor, 0.025/86400.0)

    ! Zooplankton assimilation efficiency [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooassi', wombat%zooassi, 0.6)

    ! Zooplankton half saturation coefficient for linear mortality
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zookz', wombat%zookz, 0.25)

    ! Zooplankton maximum grazing rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoogmax', wombat%zoogmax, 3.0/86400.0)

    ! Zooplankton prey capture rate constant [m6/mmol2/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooepsi', wombat%zooepsi, 0.05/86400.0)

    ! Zooplankton respiration rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zoolmor', wombat%zoolmor, 0.01/86400.0)

    ! Zooplankton quadratic mortality rate constant [m3/mmolN/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('zooqmor', wombat%zooqmor, 0.5/86400.0)

    ! Detritus remineralisation rate constant for <180 m; value for >=180m is half of this [1/s]
    !-----------------------------------------------------------------------
    ! Detritus remineralisation rate for (>=180 m) is hard-coded to be
    ! detlrem/2
    call g_tracer_add_param('detlrem', wombat%detlrem, 0.09/86400.0)
    
    ! Detritus remineralisation rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    ! This would normally equal detlrem, but we set the default value to be
    ! consistent with what is in ACCESS-ESM1.5 (undocumented in Ziehn et al
    ! 2020)
    call g_tracer_add_param('detlrem_sed', wombat%detlrem_sed, 0.02/86400.0)

    ! CaCO3 remineralisation rate constant [1/s]
    !-----------------------------------------------------------------------
    ! Default value matches 0.001714 day-1 in Ziehn et al 2020; differs from
    ! Hayashida et al 2020
    call g_tracer_add_param('caco3lrem', wombat%caco3lrem, 0.01/86400.0)

    ! CaCO3 remineralization rate constant in sediments [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('caco3lrem_sed', wombat%caco3lrem_sed, 0.0035/86400.0)

    ! CaCO3 inorganic fraction [1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('f_inorg', wombat%f_inorg, 0.062)

    ! Iron scavenging rate constant [1/s]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('tscav_fe', wombat%tscav_fe, 3.17098e-8)
    
    ! Background concentration of iron-binding ligand [umol/m3]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('ligand', wombat%ligand, 0.7)

    ! Precipitation of Fe` as nanoparticles (in excess of solubility) [/d]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('knano_dfe', wombat%knano_dfe, 0.01)

    ! Scavenging of Fe` onto biogenic particles [(mmolC/m3)-1 d-1]
    !-----------------------------------------------------------------------
    call g_tracer_add_param('kscav_dfe', wombat%kscav_dfe, 0.005)

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
    call g_tracer_add_param('wcaco3', wombat%wcaco3, 5.0/86400.0)

    ! Detritus sinking velocity [m/s]
    !-----------------------------------------------------------------------
    ! Default value matches Ziehn et al 2020 but differs from Hayashida et
    ! al 2020
    call g_tracer_add_param('wdetbio', wombat%wdetbio, 18.0/86400.0)

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

    ! Phytoplankton Chlorophyll
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'pchl', &
        longname = 'Phytoplankton chlorophyll', &
        units = 'mol/kg', &
        prog = .true.)

    ! Phytoplankton Iron content
    !-----------------------------------------------------------------------
    ! dts: There is currently no sea-ice coupling of Phytoplankton
    call g_tracer_add(tracer_list, package_name, &
        name = 'phyfe', &
        longname = 'Phytoplankton iron content', &
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

    ! Zooplankton iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'zoofe', &
        longname = 'Zooplankton iron content', &
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

    ! Detrital iron content
    !-----------------------------------------------------------------------
    call g_tracer_add(tracer_list, package_name, &
        name = 'detfe', &
        longname = 'Detrital iron content', &
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

    ! Detrital irons sitting at base of column as sediment
    !-----------------------------------------------------------------------
    ! dts: included here so included in restart
    call g_tracer_add(tracer_list, package_name, &
        name = 'detfe_sediment', &
        longname = 'Detrital iron at base of column as sediment', &
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

    call g_tracer_get_values(tracer_list, 'detfe', 'btm_reservoir', wombat%detfe_btm, isd, jsd)
    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment)
    wombat%p_detfe_sediment(:,:,1) = wombat%p_detfe_sediment(:,:,1) + wombat%detfe_btm(:,:) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'detfe', 'btm_reservoir', 0.0)

    call g_tracer_get_values(tracer_list, 'caco3', 'btm_reservoir', wombat%caco3_btm, isd, jsd)
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment)
    wombat%p_caco3_sediment(:,:,1) =  wombat%p_caco3_sediment(:,:,1) + wombat%caco3_btm(:,:) ! [mol/m2]
    call g_tracer_set_values(tracer_list, 'caco3', 'btm_reservoir', 0.0)

    ! Send diagnostics
    !-----------------------------------------------------------------------
    if (wombat%id_det_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_det_sed_depst, wombat%det_btm / dt, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_detfe_sed_depst .gt. 0) &
      used = g_send_data(wombat%id_detfe_sed_depst, wombat%detfe_btm / dt, model_time, &
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
    real                                    :: mmol_m3_to_mol_kg, umol_m3_to_mol_kg
    integer                                 :: ts_npzd ! number of time steps within NPZD model
    real                                    :: dtsb ! number of seconds per NPZD timestep
    real                                    :: rdtts ! 1 / dt
    real, dimension(nbands)                 :: sw_pen
    real                                    :: swpar
    real                                    :: u_npz, g_npz
    real                                    :: biophy, biozoo, biodet, biono3, biofer, biocaco3
    real                                    :: biophyfe 
    real                                    :: fbc
    real                                    :: no3_bgc_change, caco3_bgc_change
    real                                    :: epsi = 1.0e-30
    real                                    :: pi = 3.14159265358979
    integer                                 :: ichl
    real                                    :: par_phy_mldsum, par_z_mldsum
    real                                    :: chl, zchl, zval, phy_chlc
    real                                    :: phy_k_nit, phy_k_fer
    real                                    :: phy_pisl, phy_pisl2 
    real                                    :: pchl_pisl, pchl_lpar, pchl_mumin, pchl_muopt
    real, dimension(:,:), allocatable       :: ek_bgr, par_bgr_mid, par_bgr_top
    real, dimension(:), allocatable         :: par_mid
    real, dimension(:,:,:), allocatable     :: par_phy, par_phymld
    real, dimension(4,61)                   :: zbgr
    real                                    :: ztemk, fe_keq, fe_par, fe_sfe, fe_tfe, partic
    real                                    :: fesol1, fesol2, fesol3, fesol4, fesol5, hp, fe3sol
    real                                    :: biof, biodoc
    real                                    :: phy_Fe2C, zoo_Fe2C, det_Fe2C
    real                                    :: phy_minqfe, phy_maxqfe, phy_dFeupt_upreg, phy_dFeupt_doreg
    real                                    :: zoo_slmor
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
    ! Attenuation coefficients for blue, green and red light
    !=======================================================================
    ! Chlorophyll      ! Blue attenuation    ! Green attenuation   ! Red attenuation
    zbgr(1, 1) =  0.010; zbgr(2, 1) = 0.01618; zbgr(3, 1) = 0.07464; zbgr(4, 1) = 0.3780
    zbgr(1, 2) =  0.011; zbgr(2, 2) = 0.01654; zbgr(3, 2) = 0.07480; zbgr(4, 2) = 0.37823
    zbgr(1, 3) =  0.013; zbgr(2, 3) = 0.01693; zbgr(3, 3) = 0.07499; zbgr(4, 3) = 0.37840
    zbgr(1, 4) =  0.014; zbgr(2, 4) = 0.01736; zbgr(3, 4) = 0.07518; zbgr(4, 4) = 0.37859
    zbgr(1, 5) =  0.016; zbgr(2, 5) = 0.01782; zbgr(3, 5) = 0.07539; zbgr(4, 5) = 0.37879
    zbgr(1, 6) =  0.018; zbgr(2, 6) = 0.01831; zbgr(3, 6) = 0.07562; zbgr(4, 6) = 0.37900
    zbgr(1, 7) =  0.020; zbgr(2, 7) = 0.01885; zbgr(3, 7) = 0.07586; zbgr(4, 7) = 0.37923
    zbgr(1, 8) =  0.022; zbgr(2, 8) = 0.01943; zbgr(3, 8) = 0.07613; zbgr(4, 8) = 0.37948
    zbgr(1, 9) =  0.025; zbgr(2, 9) = 0.02005; zbgr(3, 9) = 0.07641; zbgr(4, 9) = 0.37976
    zbgr(1,10) =  0.028; zbgr(2,10) = 0.02073; zbgr(3,10) = 0.07672; zbgr(4,10) = 0.38005
    zbgr(1,11) =  0.032; zbgr(2,11) = 0.02146; zbgr(3,11) = 0.07705; zbgr(4,11) = 0.38036
    zbgr(1,12) =  0.035; zbgr(2,12) = 0.02224; zbgr(3,12) = 0.07741; zbgr(4,12) = 0.38070
    zbgr(1,13) =  0.040; zbgr(2,13) = 0.02310; zbgr(3,13) = 0.07780; zbgr(4,13) = 0.38107
    zbgr(1,14) =  0.045; zbgr(2,14) = 0.02402; zbgr(3,14) = 0.07821; zbgr(4,14) = 0.38146
    zbgr(1,15) =  0.050; zbgr(2,15) = 0.02501; zbgr(3,15) = 0.07866; zbgr(4,15) = 0.38189
    zbgr(1,16) =  0.056; zbgr(2,16) = 0.02608; zbgr(3,16) = 0.07914; zbgr(4,16) = 0.38235
    zbgr(1,17) =  0.063; zbgr(2,17) = 0.02724; zbgr(3,17) = 0.07967; zbgr(4,17) = 0.38285
    zbgr(1,18) =  0.071; zbgr(2,18) = 0.02849; zbgr(3,18) = 0.08023; zbgr(4,18) = 0.38338
    zbgr(1,19) =  0.079; zbgr(2,19) = 0.02984; zbgr(3,19) = 0.08083; zbgr(4,19) = 0.38396
    zbgr(1,20) =  0.089; zbgr(2,20) = 0.03131; zbgr(3,20) = 0.08149; zbgr(4,20) = 0.38458
    zbgr(1,21) =  0.100; zbgr(2,21) = 0.03288; zbgr(3,21) = 0.08219; zbgr(4,21) = 0.38526
    zbgr(1,22) =  0.112; zbgr(2,22) = 0.03459; zbgr(3,22) = 0.08295; zbgr(4,22) = 0.38598
    zbgr(1,23) =  0.126; zbgr(2,23) = 0.03643; zbgr(3,23) = 0.08377; zbgr(4,23) = 0.38676
    zbgr(1,24) =  0.141; zbgr(2,24) = 0.03842; zbgr(3,24) = 0.08466; zbgr(4,24) = 0.38761
    zbgr(1,25) =  0.158; zbgr(2,25) = 0.04057; zbgr(3,25) = 0.08561; zbgr(4,25) = 0.38852
    zbgr(1,26) =  0.178; zbgr(2,26) = 0.04289; zbgr(3,26) = 0.08664; zbgr(4,26) = 0.38950
    zbgr(1,27) =  0.200; zbgr(2,27) = 0.04540; zbgr(3,27) = 0.08775; zbgr(4,27) = 0.39056
    zbgr(1,28) =  0.224; zbgr(2,28) = 0.04811; zbgr(3,28) = 0.08894; zbgr(4,28) = 0.39171
    zbgr(1,29) =  0.251; zbgr(2,29) = 0.05103; zbgr(3,29) = 0.09023; zbgr(4,29) = 0.39294
    zbgr(1,30) =  0.282; zbgr(2,30) = 0.05420; zbgr(3,30) = 0.09162; zbgr(4,30) = 0.39428
    zbgr(1,31) =  0.316; zbgr(2,31) = 0.05761; zbgr(3,31) = 0.09312; zbgr(4,31) = 0.39572
    zbgr(1,32) =  0.355; zbgr(2,32) = 0.06130; zbgr(3,32) = 0.09474; zbgr(4,32) = 0.39727
    zbgr(1,33) =  0.398; zbgr(2,33) = 0.06529; zbgr(3,33) = 0.09649; zbgr(4,33) = 0.39894
    zbgr(1,34) =  0.447; zbgr(2,34) = 0.06959; zbgr(3,34) = 0.09837; zbgr(4,34) = 0.40075
    zbgr(1,35) =  0.501; zbgr(2,35) = 0.07424; zbgr(3,35) = 0.10040; zbgr(4,35) = 0.40270
    zbgr(1,36) =  0.562; zbgr(2,36) = 0.07927; zbgr(3,36) = 0.10259; zbgr(4,36) = 0.40480
    zbgr(1,37) =  0.631; zbgr(2,37) = 0.08470; zbgr(3,37) = 0.10495; zbgr(4,37) = 0.40707
    zbgr(1,38) =  0.708; zbgr(2,38) = 0.09056; zbgr(3,38) = 0.10749; zbgr(4,38) = 0.40952
    zbgr(1,39) =  0.794; zbgr(2,39) = 0.09690; zbgr(3,39) = 0.11024; zbgr(4,39) = 0.41216
    zbgr(1,40) =  0.891; zbgr(2,40) = 0.10374; zbgr(3,40) = 0.11320; zbgr(4,40) = 0.41502
    zbgr(1,41) =  1.000; zbgr(2,41) = 0.11114; zbgr(3,41) = 0.11639; zbgr(4,41) = 0.41809
    zbgr(1,42) =  1.122; zbgr(2,42) = 0.11912; zbgr(3,42) = 0.11984; zbgr(4,42) = 0.42142
    zbgr(1,43) =  1.259; zbgr(2,43) = 0.12775; zbgr(3,43) = 0.12356; zbgr(4,43) = 0.42500
    zbgr(1,44) =  1.413; zbgr(2,44) = 0.13707; zbgr(3,44) = 0.12757; zbgr(4,44) = 0.42887
    zbgr(1,45) =  1.585; zbgr(2,45) = 0.14715; zbgr(3,45) = 0.13189; zbgr(4,45) = 0.43304
    zbgr(1,46) =  1.778; zbgr(2,46) = 0.15803; zbgr(3,46) = 0.13655; zbgr(4,46) = 0.43754
    zbgr(1,47) =  1.995; zbgr(2,47) = 0.16978; zbgr(3,47) = 0.14158; zbgr(4,47) = 0.44240
    zbgr(1,48) =  2.239; zbgr(2,48) = 0.18248; zbgr(3,48) = 0.14701; zbgr(4,48) = 0.44765
    zbgr(1,49) =  2.512; zbgr(2,49) = 0.19620; zbgr(3,49) = 0.15286; zbgr(4,49) = 0.45331
    zbgr(1,50) =  2.818; zbgr(2,50) = 0.21102; zbgr(3,50) = 0.15918; zbgr(4,50) = 0.45942
    zbgr(1,51) =  3.162; zbgr(2,51) = 0.22703; zbgr(3,51) = 0.16599; zbgr(4,51) = 0.46601
    zbgr(1,52) =  3.548; zbgr(2,52) = 0.24433; zbgr(3,52) = 0.17334; zbgr(4,52) = 0.47313
    zbgr(1,53) =  3.981; zbgr(2,53) = 0.26301; zbgr(3,53) = 0.18126; zbgr(4,53) = 0.48080
    zbgr(1,54) =  4.467; zbgr(2,54) = 0.28320; zbgr(3,54) = 0.18981; zbgr(4,54) = 0.48909
    zbgr(1,55) =  5.012; zbgr(2,55) = 0.30502; zbgr(3,55) = 0.19903; zbgr(4,55) = 0.49803
    zbgr(1,56) =  5.623; zbgr(2,56) = 0.32858; zbgr(3,56) = 0.20898; zbgr(4,56) = 0.50768
    zbgr(1,57) =  6.310; zbgr(2,57) = 0.35404; zbgr(3,57) = 0.21971; zbgr(4,57) = 0.51810
    zbgr(1,58) =  7.079; zbgr(2,58) = 0.38154; zbgr(3,58) = 0.23129; zbgr(4,58) = 0.52934
    zbgr(1,59) =  7.943; zbgr(2,59) = 0.41125; zbgr(3,59) = 0.24378; zbgr(4,59) = 0.54147
    zbgr(1,60) =  8.912; zbgr(2,60) = 0.44336; zbgr(3,60) = 0.25725; zbgr(4,60) = 0.55457
    zbgr(1,61) = 10.000; zbgr(2,61) = 0.47804; zbgr(3,61) = 0.27178; zbgr(4,61) = 0.56870

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
 
    do k = 1,nk !{ 
     do j = jsc,jec; do i = isc,iec
       wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,k)
       wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,k)
       wombat%ahtotallo(i,j) = wombat%htotal_scale_lo * wombat%ahtotal(i,j,k)
       wombat%ahtotalhi(i,j) = wombat%htotal_scale_hi * wombat%ahtotal(i,j,k)
     enddo; enddo 

     if (k.eq.1) then !{
       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), Salt(:,:,k), &
           wombat%f_dic(:,:,k), &
           wombat%f_no3(:,:,k) / 16., &
           wombat%sio2(:,:), &
           wombat%f_alk(:,:,k), &
           wombat%htotallo(:,:), wombat%htotalhi(:,:), &
           wombat%htotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k), &
           co2star=wombat%co2_csurf(:,:), alpha=wombat%co2_alpha(:,:), &
           pCO2surf=wombat%pco2_csurf(:,:))
 
       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), Salt(:,:,k), &
           wombat%f_adic(:,:,k), &
           wombat%f_no3(:,:,k) / 16., &
           wombat%sio2(:,:), &
           wombat%f_alk(:,:,k), &
           wombat%ahtotallo(:,:), wombat%ahtotalhi(:,:), &
           wombat%ahtotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k), &
           co2star=wombat%aco2_csurf(:,:), alpha=wombat%aco2_alpha(:,:), &
           pCO2surf=wombat%paco2_csurf(:,:))

       call g_tracer_set_values(tracer_list, 'dic', 'alpha', wombat%co2_alpha, isd, jsd)
       call g_tracer_set_values(tracer_list, 'dic', 'csurf', wombat%co2_csurf, isd, jsd)
       call g_tracer_set_values(tracer_list, 'adic', 'alpha', wombat%aco2_alpha, isd, jsd)
       call g_tracer_set_values(tracer_list, 'adic', 'csurf', wombat%aco2_csurf, isd, jsd)

     else

       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), Salt(:,:,k), &
           wombat%f_dic(:,:,k), &
           wombat%f_no3(:,:,k) / 16., &
           wombat%sio2(:,:), &
           wombat%f_alk(:,:,k), &
           wombat%htotallo(:,:), wombat%htotalhi(:,:), &
           wombat%htotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k))
 
       call FMS_ocmip2_co2calc(CO2_dope_vec, grid_tmask(:,:,k), &
           Temp(:,:,k), Salt(:,:,k), &
           wombat%f_adic(:,:,k), &
           wombat%f_no3(:,:,k) / 16., &
           wombat%sio2(:,:), &
           wombat%f_alk(:,:,k), &
           wombat%ahtotallo(:,:), wombat%ahtotalhi(:,:), &
           wombat%ahtotal(:,:,k), &
           co2_calc=trim(co2_calc), &
           zt=wombat%zw(:,:,k))

     endif !} if k.eq.1
   enddo !} do k = 1,nk

    !=======================================================================
    ! Calculate the source terms
    !=======================================================================

    wombat%pprod_gross(:,:,:) = 0.0
    wombat%zprod_gross(:,:,:) = 0.0
    wombat%light_limit(:,:) = 0.0
    wombat%radbio3d(:,:,:) = 0.0
    wombat%radmid3d(:,:,:) = 0.0
    wombat%npp3d(:,:,:) = 0.0
    wombat%phy_mumax(:,:,:) = 0.0
    wombat%phy_mu(:,:,:) = 0.0
    wombat%pchl_mu(:,:,:) = 0.0
    wombat%phy_lpar(:,:,:) = 1.0
    wombat%phy_lnit(:,:,:) = 1.0
    wombat%phy_lfer(:,:,:) = 1.0
    wombat%phy_dfeupt(:,:,:) = 0.0
    wombat%feIII(:,:,:) = 0.0
    wombat%felig(:,:,:) = 0.0
    wombat%fecol(:,:,:) = 0.0
    wombat%feprecip(:,:,:) = 0.0
    wombat%fescaven(:,:,:) = 0.0
    wombat%fescadet(:,:,:) = 0.0
    wombat%fesources(:,:,:) = 0.0
    wombat%fesinks(:,:,:) = 0.0
    wombat%fecoag2det(:,:,:) = 0.0
    wombat%phygrow(:,:,:) = 0.0
    wombat%phyresp(:,:,:) = 0.0
    wombat%phymort(:,:,:) = 0.0
    wombat%zoograz(:,:,:) = 0.0
    wombat%zooresp(:,:,:) = 0.0
    wombat%zoomort(:,:,:) = 0.0
    wombat%zooexcr(:,:,:) = 0.0
    wombat%zooslop(:,:,:) = 0.0
    wombat%detremi(:,:,:) = 0.0
    wombat%caldiss(:,:,:) = 0.0
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
    wombat%zeuphot(:,:) = 0.0
    wombat%phy_leup(:,:) = 1.0

    ! Some unit conversion factors
    mmol_m3_to_mol_kg = 1.e-3 / wombat%Rho_0
    umol_m3_to_mol_kg = 1.e-3 * mmol_m3_to_mol_kg
   
    ! Allocate and initialise some multi-dimensional variables
    allocate(ek_bgr(nk,3)); ek_bgr(:,:)=0.0
    allocate(par_bgr_mid(nk,3)); par_bgr_mid(:,:)=0.0
    allocate(par_bgr_top(nk,3)); par_bgr_top(:,:)=0.0
    allocate(par_mid(nk)); par_mid(:)=0.0
    allocate(par_phy(isc:iec, jsc:jec, nk)); par_phy(:,:,:)=0.0
    allocate(par_phymld(isc:iec, jsc:jec, nk)); par_phymld(:,:,:)=0.0

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
    call g_tracer_get_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'detfe', 'field', wombat%f_detfe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    call g_tracer_get_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau, &
        positive=.true.) ! [mol/kg]
    
    
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !     (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.                 !
    !     / o o \  : :.-.: :: ,. :: `' :: .; :: .; :`-. .-'                 !
    !    (   "   ) : :: :: :: :: :: .. ::   .':    :  : :                   !
    !     \__ __/  : `' `' ;: :; :: :; :: .; :: :: :  : :                   !
    !               `.,`.,' `.__.':_;:_;:___.':_;:_;  :_;                   !
    !                                                                       !
    ! World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)    !
    !                                                                       !
    !  Steps:                                                               !
    !    1.  Light attenuation through water column                         !
    !    2.  Nutrient limitation of phytoplankton                           !
    !    3.  Temperature-dependence of heterotrophy                         !
    !    4.  Light limitation of phytoplankton                              !
    !    5.  Growth of chlorophyll                                          !
    !    6.  Phytoplankton uptake of iron                                   !
    !    7.  Iron chemistry                                                 !
    !    8.  Mortality scalings and grazing                                 !
    !    9.  Sources and sinks                                              !
    !    10. Tracer tendencies                                              !
    !                                                                       !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!


    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !  [Step 1] Light attenutation through water column                     !
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    !-----------------------------------------------------------------------!
    do j = jsc,jec; do i = isc,iec
      !--- Retrieve incident light ---!
      ! dts: keep only shortwave wavelengths < 710 nm (taken from COBALT/BLING)
      swpar = 0.0
      do n = 1,nbands;
        if (max_wavelength_band(n) .lt. 710) then
          sw_pen(n) = max(0.0, sw_pen_band(n,i,j))
        else
          sw_pen(n) = 0.0
        endif
        swpar = swpar + sw_pen(n)
      enddo

      !--- Light field ---!
     
      par_bgr_top(:,:) = 0.0
      par_bgr_mid(:,:) = 0.0
      par_mid(:) = 0.0

      do k = 1,grid_kmt(i,j)  !{

        ! chlorophyll concentration conversion from mol/kg --> mg/m3
        chl = wombat%f_pchl(i,j,k) * 12.0 / mmol_m3_to_mol_kg 

        ! Attenuation coefficients given chlorophyll concentration
        zchl = max(0.05, min(10.0, chl) )
        ichl = nint( 41 + 20.0*log10(zchl) + epsi )
        ek_bgr(k,1) = zbgr(2,ichl) * dzt(i,j,k)
        ek_bgr(k,2) = zbgr(3,ichl) * dzt(i,j,k)
        ek_bgr(k,3) = zbgr(4,ichl) * dzt(i,j,k)

        ! BGR light available in the water column
        if (swpar.gt.0.0) then
          if (k.eq.1) then
            par_bgr_top(k,1) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is blue (coming through top of cell)
            par_bgr_top(k,2) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is green (coming through top of cell)
            par_bgr_top(k,3) = swpar * 1.0/3.0  ! Assumes 1/3 of PAR is red (coming through top of cell)
            par_bgr_mid(k,1) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,1)) ! Assumes 1/3 of PAR is blue (at centre point)
            par_bgr_mid(k,2) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,2)) ! Assumes 1/3 of PAR is green (at centre point)
            par_bgr_mid(k,3) = swpar * 1.0/3.0 * exp(-0.5 * ek_bgr(k,3)) ! Assumes 1/3 of PAR is red (at centre point)
          else
            par_bgr_top(k,1) = par_bgr_top(k-1,1) * exp(-ek_bgr(k-1,1))
            par_bgr_top(k,2) = par_bgr_top(k-1,2) * exp(-ek_bgr(k-1,2))
            par_bgr_top(k,3) = par_bgr_top(k-1,3) * exp(-ek_bgr(k-1,3))
            par_bgr_mid(k,1) = par_bgr_mid(k-1,1) * exp(-0.5 * (ek_bgr(k-1,1)+ek_bgr(k,1)))
            par_bgr_mid(k,2) = par_bgr_mid(k-1,2) * exp(-0.5 * (ek_bgr(k-1,2)+ek_bgr(k,2)))
            par_bgr_mid(k,3) = par_bgr_mid(k-1,3) * exp(-0.5 * (ek_bgr(k-1,3)+ek_bgr(k,3)))
          endif
        endif
        
        ! Light attenuation at mid point of the cells
        par_mid(k) = par_bgr_mid(k,1) + par_bgr_mid(k,2) + par_bgr_mid(k,3)

      enddo !} k

      ! initialise some variables
      par_phy_mldsum = 0.0
      par_z_mldsum = 0.0
      wombat%zeuphot(i,j) = 0.0  ! Euphotic zone depth set at 0.0 meters when there is no light

      !--- Collect euphotic zone depth, light at mid points, and mean light ---!
      do k = 1,grid_kmt(i,j)  !{

        if (swpar.gt.0.0) then
          ! Euphotic zone
          if (par_mid(k).gt.(swpar*0.01) .and. par_mid(k).gt.1.0) then
            wombat%zeuphot(i,j) = wombat%zw(i,j,k)
          endif
          ! Light attenuation mean over the grid cells (Eq. 19 of Baird et al., 2020 GMD)
          if (k.lt.grid_kmt(i,j)) then
            par_phy(i,j,k) = ((par_bgr_top(k,1) - par_bgr_top(k+1,1)) / ek_bgr(k,1)) + &
                             ((par_bgr_top(k,2) - par_bgr_top(k+1,2)) / ek_bgr(k,2)) + &
                             ((par_bgr_top(k,3) - par_bgr_top(k+1,3)) / ek_bgr(k,3))
          else
            par_phy(i,j,k) = ((par_bgr_top(k,1) - par_bgr_top(k,1)*exp(-ek_bgr(k,1))) / ek_bgr(k,1)) + &
                             ((par_bgr_top(k,2) - par_bgr_top(k,2)*exp(-ek_bgr(k,2))) / ek_bgr(k,2)) + &
                             ((par_bgr_top(k,3) - par_bgr_top(k,3)*exp(-ek_bgr(k,3))) / ek_bgr(k,3))
          endif
        else
          par_phy(i,j,k) = 0.0
        endif

        ! Integrated light field in the mixed layer
        if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
          par_phy_mldsum = par_phy_mldsum + par_phy(i,j,k) * dzt(i,j,k)
          par_z_mldsum = par_z_mldsum + dzt(i,j,k)
        endif

      enddo !}

      ! Calculate impact of euphotic zone depth on phytoplankton and chlorophyll production
      wombat%phy_leup(i,j) = MIN(1.0, wombat%zeuphot(i,j)/(hblt_depth(i,j) + epsi))

      !--- Aggregate light in mixed layer and calculate maximum growth rates ---!
      do k = 1,grid_kmt(i,j)  !{

        ! Calculate average light level in the mixed layer
        if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
          zval = 1.0/(par_z_mldsum + epsi)
          par_phymld(i,j,k) = par_phy_mldsum * zval
        else
          par_phymld(i,j,k) = par_phy(i,j,k)
        endif

        ! Save total PAR to radbio and radmid arrays for diagnostic output
        wombat%radbio3d(i,j,k) = par_phy(i,j,k)
        wombat%radmid3d(i,j,k) = par_mid(k)

        ! Temperature-dependent maximum growth rate (Eppley curve)
        wombat%phy_mumax(i,j,k) = wombat%abioa * wombat%bbioa ** (Temp(i,j,k)) ! [1/s]
         
      enddo  !} k

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

      ! Initialise some values and ratios (put into nicer units than mol/kg)
      biophy   = max(epsi, wombat%f_phy(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biophyfe = max(epsi, wombat%f_phyfe(i,j,k))/ mmol_m3_to_mol_kg  ![mmol/m3]
      biozoo   = max(epsi, wombat%f_zoo(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biodet   = max(epsi, wombat%f_det(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biono3   = max(epsi, wombat%f_no3(i,j,k) ) / mmol_m3_to_mol_kg  ![mmol/m3]
      biofer   = max(epsi, wombat%f_fe(i,j,k)  ) / umol_m3_to_mol_kg  ![umol/m3]
      biocaco3 = max(epsi, wombat%f_caco3(i,j,k))/ mmol_m3_to_mol_kg  ![mmol/m3]
      phy_chlc = max(epsi, wombat%f_pchl(i,j,k)) / max(epsi, wombat%f_phy(i,j,k))
      phy_Fe2C = max(epsi, wombat%f_phyfe(i,j,k))/ max(epsi, wombat%f_phy(i,j,k))
      zoo_Fe2C = max(epsi, wombat%f_zoofe(i,j,k))/ max(epsi, wombat%f_zoo(i,j,k))
      det_Fe2C = max(epsi, wombat%f_detfe(i,j,k))/ max(epsi, wombat%f_det(i,j,k))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 2] Nutrient limitation of phytoplankton                        !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Allometric scaling of 0.37 (Wickman et al., 2024; Science)
      ! 2. Apply variable K to set limitation term

      phy_k_nit = wombat%phykn * max(biophy, 0.1)**(0.37)
      phy_k_fer = wombat%phykf * max(biophy, 0.1)**(0.37)
      ! Nitrogen limitation (currently Michaelis-Menten term)
      wombat%phy_lnit(i,j,k) = biono3 / (biono3 + phy_k_nit + epsi)
      ! Iron limitation (Quota model, constants from Flynn & Hipkin 1999)
      if (biophy.gt.1e-3) then
        phy_minqfe = 0.00167 / 55.85 * phy_chlc + &
                     1.21e-5 * 14.0 / 55.85 / 7.625 * 0.5 * 1.5 * wombat%phy_lnit(i,j,k) + &
                     1.15e-4 * 14.0 / 55.85 / 7.625 * 0.5 * wombat%phy_lnit(i,j,k)
        wombat%phy_lfer(i,j,k) = min(1.0, max(0.0, (phy_Fe2C - phy_minqfe) / wombat%phyoptqf ))
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 3] Temperature-dependence of heterotrophy                      !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Temperature dependance of heterotrophy
      fbc = wombat%bbioh ** (Temp(i,j,k))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 4] Light limitation of phytoplankton                           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      
      ! 1. initial slope of Photosynthesis-Irradiance curve
      ! 2. Alter the slope to account for respiration and daylength limitation
      ! 3. Light limitation 
      ! 4. Apply light and nutrient limitations to maximum growth rate

      phy_pisl  = max(wombat%alphabio * phy_chlc, wombat%alphabio * wombat%phyminqc)
      phy_pisl2 = phy_pisl / ( (1. + wombat%phylmor*86400.0 * fbc) ) ! add daylength estimate here
      wombat%phy_lpar(i,j,k) = (1. - exp(-phy_pisl2 * par_phy(i,j,k))) * wombat%phy_leup(i,j)
      wombat%phy_mu(i,j,k) = wombat%phy_mumax(i,j,k) * wombat%phy_lpar(i,j,k) * & 
                             min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))

      ! Save the average light limitation over the euphotic zone
      if (wombat%zw(i,j,k) .le. wombat%zeuphot(i,j)) then 
        wombat%light_limit(i,j) = wombat%light_limit(i,j) + dzt(i,j,k) &
                                  * dtsb*rdtts * wombat%phy_lpar(i,j,k)
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 5] Growth of chlorophyll                                       !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Light limitation of chlorophyll production
      ! 2. minimum and optimal rates of chlorophyll growth
      ! 3. Calculate mg Chl m-3 s-1
      
      pchl_pisl = phy_pisl / ( wombat%phy_mumax(i,j,k) * 86400.0 * & 
                               (1. - min(wombat%phy_lnit(i,j,k), &
                                         wombat%phy_lfer(i,j,k))) + epsi )
      pchl_lpar = (1. - exp(-pchl_pisl * par_phymld(i,j,k))) * wombat%phy_leup(i,j)
      pchl_mumin = wombat%phyminqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s] 
      pchl_muopt = wombat%phyoptqc * wombat%phy_mu(i,j,k) * biophy * 12.0   ![mg/m3/s]
      wombat%pchl_mu(i,j,k) = (pchl_muopt - pchl_mumin) * pchl_lpar * &
                min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))
      if ( (phy_pisl * par_phymld(i,j,k)) .gt. 0.0 ) then
        wombat%pchl_mu(i,j,k) = pchl_mumin + wombat%pchl_mu(i,j,k) / &
                                (phy_pisl * par_phymld(i,j,k))
      endif
      wombat%pchl_mu(i,j,k) = wombat%pchl_mu(i,j,k) / 12.0 * mmol_m3_to_mol_kg  ![mol/kg/s]


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 6] Phytoplankton uptake of iron                                !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! 1. Maximum iron content of phytoplankton cell
      ! 2. Ensure that dFe uptake increases or decreases in response to cell quota
      ! 3. Iron uptake of phytoplankton

      phy_maxqfe = biophy * wombat%phymaxqf
      phy_dFeupt_upreg = (4.0 - 4.5 * wombat%phy_lfer(i,j,k) / &
                          (wombat%phy_lfer(i,j,k) + 0.5) )
      phy_dFeupt_doreg = max(0.0, (1.0 - biophyfe/phy_maxqfe) / &
                                  abs(1.05 - biophyfe/phy_maxqfe) )
      wombat%phy_dfeupt(i,j,k) = (biophy * wombat%phy_mumax(i,j,k) * wombat%phymaxqf * &
                                  biofer / (biofer + phy_k_fer) * &
                                  phy_dFeupt_doreg * phy_dFeupt_upreg) * mmol_m3_to_mol_kg


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 7] Iron chemistry (Aumont et al., 2015; GMD)                   !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Estimate solubility of Fe3+ (free Fe) in solution using temperature, pH and salinity
      ztemk = max(5.0, Temp(i,j,k)) + 273.15    ! temperature in kelvin
      zval = 19.924 * Salt(i,j,k) / ( 1000. - 1.005 * Salt(i,j,k))
      fesol1 = 10**(-13.486 - 0.1856*zval**0.5 + 0.3073*zval + 5254.0/ztemk)
      fesol2 = 10**(2.517 - 0.8885*zval**0.5 + 0.2139*zval - 1320.0/ztemk)
      fesol3 = 10**(0.4511 - 0.3305*zval**0.5 - 1996.0/ztemk)
      fesol4 = 10**(-0.2965 - 0.7881*zval**0.5 - 4086.0/ztemk)
      fesol5 = 10**(4.4466 - 0.8505*zval**0.5 - 7980.0/ztemk)
      hp = wombat%ahtotal(i,j,k)
      fe3sol = fesol1 * ( hp**3 + fesol2 * hp**2 + fesol3 * hp + fesol4 + fesol5 / hp ) *1e9

      ! Estimate total colloidal iron
      ! ... for now, we assume that 50% of all dFe is colloidal
      wombat%fecol(i,j,k) = 0.5 * biofer 

      ! Determine equilibriuim fractionation of the remain dFe (non-colloidal fraction) into Fe' and L-Fe
      fe_keq = 10**( 17.27 - 1565.7 / ztemk ) * 1e-9 ! Temperature reduces solubility
      fe_par = 4.77e-7 * wombat%radbio3d(i,j,k) * 0.5 / 10**(-6.3) ! Light increases solubility
      fe_sfe = max(0.0, biofer - wombat%fecol(i,j,k))
      fe_tfe = (1.0 + fe_par) * fe_sfe
      wombat%feIII(i,j,k) = ( -( 1. + wombat%ligand * fe_keq + fe_par - fe_sfe * fe_keq ) &
                             + SQRT( ( 1. + wombat%ligand * fe_keq + fe_par - fe_sfe * fe_keq )**2 &
                                      + 4. * fe_tfe * fe_keq) ) / ( 2. * fe_keq + epsi )
      wombat%felig(i,j,k) = max(0.0, fe_sfe - wombat%feIII(i,j,k))

      ! Precipitation of Fe' (creation of nanoparticles)
      wombat%feprecip(i,j,k) = max(0.0, ( wombat%feIII(i,j,k) - fe3sol ) ) * wombat%knano_dfe/86400.0

      ! Scavenging of Fe` onto biogenic particles 
      partic = (biodet + biocaco3)
      wombat%fescaven(i,j,k) = wombat%feIII(i,j,k) * (3e-5 + wombat%kscav_dfe * partic) / 86400.0
      wombat%fescadet(i,j,k) = wombat%fescaven(i,j,k) * biodet / (partic+epsi) 

      ! Coagulation of colloidal Fe (umol/m3) to form sinking particles (mmol/m3)
      ! Following Tagliabue et al. (2023), make coagulation rate dependent on DOC and Phytoplankton biomass
      biof = max(1/3., biophy / (biophy + 0.03))
      biodoc = 2.0 + (1.0 - min(wombat%phy_lnit(i,j,k), wombat%phy_lfer(i,j,k))) * 38.0 ! proxy of DOC (mmol/m3)
      if (wombat%zw(i,j,k).le.hblt_depth(i,j)) then
        zval = (      (12.*biof*biodoc + 9.*biodet) + 2.5*biodet + 128.*biof*biodoc + 725.*biodet )*1e-8
      else
        zval = ( 0.01*(12.*biof*biodoc + 9.*biodet) + 2.5*biodet + 128.*biof*biodoc + 725.*biodet )*1e-8
      endif
      wombat%fecoag2det(i,j,k) = wombat%fecol(i,j,k) * zval / 86400.0

      ! Convert the sink terms to mol/kg
      wombat%feprecip(i,j,k) = wombat%feprecip(i,j,k) * umol_m3_to_mol_kg
      wombat%fescaven(i,j,k) = wombat%fescaven(i,j,k) * umol_m3_to_mol_kg
      wombat%fescadet(i,j,k) = wombat%fescadet(i,j,k) * umol_m3_to_mol_kg
      wombat%fecoag2det(i,j,k) = wombat%fecoag2det(i,j,k) * umol_m3_to_mol_kg


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 8] Mortality scalings and grazing                              !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      zoo_slmor = biozoo / (biozoo + wombat%zookz)
      
      ! Grazing function ! [1/s]
      g_npz = wombat%zoogmax * fbc * (wombat%zooepsi * biophy * biophy) / &
              (wombat%zoogmax * fbc + (wombat%zooepsi * biophy * biophy))


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 9] Sources and sinks                                           !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      if (wombat%f_no3(i,j,k) .gt. epsi) then
        wombat%phygrow(i,j,k) = wombat%phy_mu(i,j,k) * wombat%f_phy(i,j,k) ! [molC/kg/s]
      else
        wombat%phygrow(i,j,k) = 0.0
      endif

      if (biophy .gt. 1e-3) then
        wombat%zoograz(i,j,k) = g_npz * wombat%f_zoo(i,j,k) ! [molC/kg/s]
        wombat%phyresp(i,j,k) = wombat%phylmor * fbc * wombat%f_phy(i,j,k) ! [molC/kg/s]
        wombat%phymort(i,j,k) = wombat%phyqmor / mmol_m3_to_mol_kg * wombat%f_phy(i,j,k) * wombat%f_phy(i,j,k) ! [molC/kg/s]
      else
        wombat%zoograz(i,j,k) = 0.0
        wombat%phyresp(i,j,k) = 0.0
        wombat%phymort(i,j,k) = 0.0
      endif
        wombat%zooexcr(i,j,k) = wombat%zoograz(i,j,k) * (1.0 - wombat%zooassi)*0.75
        wombat%zooslop(i,j,k) = wombat%zoograz(i,j,k) * (1.0 - wombat%zooassi)*0.25
      
      if (biozoo .gt. 1e-3) then
        wombat%zooresp(i,j,k) = wombat%zoolmor * fbc * wombat%f_zoo(i,j,k) * zoo_slmor ! [molC/kg/s]
        wombat%zoomort(i,j,k) = wombat%zooqmor / mmol_m3_to_mol_kg * wombat%f_zoo(i,j,k) * wombat%f_zoo(i,j,k) ! [molC/kg/s]
      else
        wombat%zooresp(i,j,k) = 0.0
        wombat%zoomort(i,j,k) = 0.0
      endif
      
      if (wombat%f_det(i,j,k) .gt. epsi) then
        wombat%detremi(i,j,k) = wombat%detlrem * fbc * wombat%f_det(i,j,k) ! [molC/kg/s]
        if (wombat%zw(i,j,k) .ge. 180.0) wombat%detremi(i,j,k) = 0.5 * wombat%detremi(i,j,k) ! reduce decay below 180m
      else
        wombat%detremi(i,j,k) = 0.0
      endif
      
      if (wombat%f_caco3(i,j,k) .gt. epsi) then
        wombat%caldiss(i,j,k) = wombat%caco3lrem * wombat%f_caco3(i,j,k) ! [mol/kg/s]
      else
        wombat%caldiss(i,j,k) = 0.0
      endif


      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !  [Step 10] Tracer tendencies                                          !
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!
      !-----------------------------------------------------------------------!

      ! Nitrate equation ! [molN/kg]
      !----------------------------------------------------------------------
      wombat%f_no3(i,j,k) = wombat%f_no3(i,j,k) + dtsb * 16./122. * ( &
                              wombat%detremi(i,j,k) + &
                              wombat%zooresp(i,j,k) + &
                              wombat%zooexcr(i,j,k) + &
                              wombat%phyresp(i,j,k) - &
                              wombat%phygrow(i,j,k) )

      ! Phytoplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_phy(i,j,k)  = wombat%f_phy(i,j,k) + dtsb * ( &
                               wombat%phygrow(i,j,k) - &
                               wombat%phyresp(i,j,k) - &
                               wombat%phymort(i,j,k) - &
                               wombat%zoograz(i,j,k) )
      
      ! Phytoplankton chlorophyll equation ! [molChl/kg]
      !-----------------------------------------------------------------------
      wombat%f_pchl(i,j,k)  = wombat%f_pchl(i,j,k) + dtsb * ( &
                                wombat%pchl_mu(i,j,k) - &
                                wombat%phyresp(i,j,k) * phy_chlc - &
                                wombat%phymort(i,j,k) * phy_chlc - &
                                wombat%zoograz(i,j,k) * phy_chlc )

      ! Phytoplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_phyfe(i,j,k)  = wombat%f_phyfe(i,j,k) + dtsb * ( & 
                                 wombat%phy_dfeupt(i,j,k) - &
                                 wombat%phyresp(i,j,k) * phy_Fe2C - &
                                 wombat%phymort(i,j,k) * phy_Fe2C - &
                                 wombat%zoograz(i,j,k) * phy_Fe2C )

      ! Estimate primary productivity from phytoplankton growth ! [molC/kg/s]
      wombat%pprod_gross(i,j,k) = wombat%pprod_gross(i,j,k) + dtsb * wombat%phygrow(i,j,k)

      ! Net primary productivity (gross PP minus linear mortality) ! [molC/kg/s]
      wombat%npp3d(i,j,k) = wombat%npp3d(i,j,k) + dtsb * ( &
                              wombat%phygrow(i,j,k) - wombat%phyresp(i,j,k) )

      ! Zooplankton equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoo(i,j,k)  = wombat%f_zoo(i,j,k) + dtsb * ( &
                               wombat%zooassi * wombat%zoograz(i,j,k) - &
                               wombat%zooresp(i,j,k) - &
                               wombat%zoomort(i,j,k) )

      ! Zooplankton iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_zoofe(i,j,k)  = wombat%f_zoofe(i,j,k) + dtsb * ( &
                                 wombat%zooassi * wombat%zoograz(i,j,k) * phy_Fe2C - &
                                 wombat%zooresp(i,j,k) * zoo_Fe2C - &
                                 wombat%zoomort(i,j,k) * zoo_Fe2C )

      ! Estimate secondary productivity from zooplankton growth ! [molC/kg/s]
      wombat%zprod_gross(i,j,k) = wombat%zprod_gross(i,j,k) + dtsb * &
                                    wombat%zooassi * wombat%zoograz(i,j,k)

      ! Detritus equation ! [molC/kg]
      !-----------------------------------------------------------------------
      wombat%f_det(i,j,k) = wombat%f_det(i,j,k) + dtsb * ( & 
                              wombat%zooslop(i,j,k) + &
                              wombat%phymort(i,j,k) + &
                              wombat%zoomort(i,j,k) - &
                              wombat%detremi(i,j,k) )

      ! Detrital iron equation ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_detfe(i,j,k) = wombat%f_detfe(i,j,k) + dtsb * ( & 
                                wombat%zooslop(i,j,k) * phy_Fe2C + &
                                wombat%phymort(i,j,k) * phy_Fe2C + &
                                wombat%zoomort(i,j,k) * zoo_Fe2C - &
                                wombat%detremi(i,j,k) * det_Fe2C + &
                                wombat%fescadet(i,j,k) + &
                                wombat%fecoag2det(i,j,k) )

      ! Oxygen equation ! [molO2/kg]
      !-----------------------------------------------------------------------
      if (wombat%f_o2(i,j,k) .gt. epsi) &
        wombat%f_o2(i,j,k) = wombat%f_o2(i,j,k) - 172./122. * dtsb * ( &
                               wombat%detremi(i,j,k) + &
                               wombat%zooresp(i,j,k) + &
                               wombat%phyresp(i,j,k) - &
                               wombat%phygrow(i,j,k) )

      ! Extra equation for caco3 - alkalinity ! [molCaCO3/kg]
      !-----------------------------------------------------------------------
      ! rjm: 6.2% of POC by default 
      wombat%f_caco3(i,j,k) = wombat%f_caco3(i,j,k) + dtsb * ( &
                                wombat%zooslop(i,j,k) * wombat%f_inorg + &
                                wombat%phymort(i,j,k) * wombat%f_inorg + &
                                wombat%zoomort(i,j,k) * wombat%f_inorg - &
                                wombat%caldiss(i,j,k) )

      ! Extra equation for iron ! [molFe/kg]
      !-----------------------------------------------------------------------
      wombat%f_fe(i,j,k) = wombat%f_fe(i,j,k) + dtsb * ( &
                             wombat%detremi(i,j,k) * det_Fe2C + &
                             wombat%zooresp(i,j,k) * zoo_Fe2C + &
                             wombat%zooexcr(i,j,k) * zoo_Fe2C + &
                             wombat%phyresp(i,j,k) * phy_Fe2C - &
                             wombat%phy_dfeupt(i,j,k) - &
                             wombat%feprecip(i,j,k) - &
                             wombat%fescaven(i,j,k) - &
                             wombat%fecoag2det(i,j,k) )

      ! Collect dFe sources and sinks for diagnostic output
      wombat%fesources(i,j,k) = wombat%fesources(i,j,k) + dtsb*rdtts * ( &
                                  wombat%detremi(i,j,k) * det_Fe2C + &
                                  wombat%zooresp(i,j,k) * zoo_Fe2C + &
                                  wombat%zooexcr(i,j,k) * zoo_Fe2C + &
                                  wombat%phyresp(i,j,k) * phy_Fe2C)
      wombat%fesinks(i,j,k) = wombat%fesinks(i,j,k) + dtsb*rdtts * ( & 
                                wombat%phy_dfeupt(i,j,k) + &
                                wombat%feprecip(i,j,k) + &
                                wombat%fescaven(i,j,k) + &
                                wombat%fecoag2det(i,j,k)) 

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

      wombat%p_dic(i,j,k,tau) = wombat%p_dic(i,j,k,tau) + &
          (122./16. * no3_bgc_change - caco3_bgc_change) ! [mol/kg]

      wombat%p_adic(i,j,k,tau) = wombat%p_adic(i,j,k,tau) + &
          (122./16. * no3_bgc_change - caco3_bgc_change) ! [mol/kg]

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
    call g_tracer_set_values(tracer_list, 'pchl', 'field', wombat%f_pchl, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'phyfe', 'field', wombat%f_phyfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoo', 'field', wombat%f_zoo, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'zoofe', 'field', wombat%f_zoofe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'det', 'field', wombat%f_det, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'detfe', 'field', wombat%f_detfe, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'o2', 'field', wombat%f_o2, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'caco3', 'field', wombat%f_caco3, isd, jsd, ntau=tau)
    call g_tracer_set_values(tracer_list, 'fe', 'field', wombat%f_fe, isd, jsd, ntau=tau)

    !-----------------------------------------------------------------------
    ! Remineralisation of sediment tracers
    !-----------------------------------------------------------------------
    call g_tracer_get_pointer(tracer_list, 'det_sediment', 'field', wombat%p_det_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'detfe_sediment', 'field', wombat%p_detfe_sediment) ! [mol/m2]
    call g_tracer_get_pointer(tracer_list, 'caco3_sediment', 'field', wombat%p_caco3_sediment) ! [mol/m2]

    do j = jsc,jec; do i = isc,iec;
      k = grid_kmt(i,j)
      if (k .gt. 0) then
        fbc = wombat%bbioh ** (Temp(i,j,k)) ! [1]
        wombat%det_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_det_sediment(i,j,1) ! [mol/m2/s]
        wombat%detfe_sed_remin(i,j) = wombat%detlrem_sed * fbc * wombat%p_detfe_sediment(i,j,1) ! [mol/m2/s]
        wombat%caco3_sed_remin(i,j) = wombat%caco3lrem_sed * fbc * wombat%p_caco3_sediment(i,j,1) ! [mol/m2/s]
        
        ! Remineralisation of sediments to supply nutrient fields.
        ! btf values are positive from the water column into the sediment.
        wombat%b_no3(i,j) = -16./122. * wombat%det_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_o2(i,j) = -172./16. * wombat%b_no3(i,j) ! [mol/m2/s]
        wombat%b_dic(i,j) = 122./16. * wombat%b_no3(i,j) - wombat%caco3_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_adic(i,j) = wombat%b_dic(i,j) ! [mol/m2/s]
        wombat%b_fe(i,j) = -1.0 * wombat%detfe_sed_remin(i,j) ! [mol/m2/s]
        wombat%b_alk(i,j) = -2.0 * wombat%caco3_sed_remin(i,j) - wombat%b_no3(i,j) ! [mol/m2/s]
      endif
    enddo; enddo

    ! Apply remineralisation rates to sediment tracers
    !-----------------------------------------------------------------------
    do j = jsc,jec; do i = isc,iec;
      if (grid_kmt(i,j) .gt. 0) then
        wombat%p_det_sediment(i,j,1) = wombat%p_det_sediment(i,j,1) - dt * wombat%det_sed_remin(i,j) ! [mol/m2]
        wombat%p_detfe_sediment(i,j,1) = wombat%p_detfe_sediment(i,j,1) - dt * wombat%detfe_sed_remin(i,j) ! [mol/m2]
        wombat%p_caco3_sediment(i,j,1) = wombat%p_caco3_sediment(i,j,1) - dt * wombat%caco3_sed_remin(i,j) ! [mol/m2]
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
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_ahtotal .gt. 0) &
      used = g_send_data(wombat%id_ahtotal, wombat%ahtotal, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

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

    if (wombat%id_radmid3d .gt. 0) &
      used = g_send_data(wombat%id_radmid3d, wombat%radmid3d, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_radbio1 .gt. 0) &
      used = g_send_data(wombat%id_radbio1, wombat%radbio3d(:,:,1), model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_pprod_gross .gt. 0) &
      used = g_send_data(wombat%id_pprod_gross, wombat%pprod_gross, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mumax .gt. 0) &
      used = g_send_data(wombat%id_phy_mumax, wombat%phy_mumax, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_mu .gt. 0) &
      used = g_send_data(wombat%id_phy_mu, wombat%phy_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_pchl_mu .gt. 0) &
      used = g_send_data(wombat%id_pchl_mu, wombat%pchl_mu, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lpar .gt. 0) &
      used = g_send_data(wombat%id_phy_lpar, wombat%phy_lpar, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lnit .gt. 0) &
      used = g_send_data(wombat%id_phy_lnit, wombat%phy_lnit, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_lfer .gt. 0) &
      used = g_send_data(wombat%id_phy_lfer, wombat%phy_lfer, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phy_dfeupt .gt. 0) &
      used = g_send_data(wombat%id_phy_dfeupt, wombat%phy_dfeupt, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_feIII .gt. 0) &
      used = g_send_data(wombat%id_feIII, wombat%feIII, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_felig .gt. 0) &
      used = g_send_data(wombat%id_felig, wombat%felig, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecol .gt. 0) &
      used = g_send_data(wombat%id_fecol, wombat%fecol, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_feprecip .gt. 0) &
      used = g_send_data(wombat%id_feprecip, wombat%feprecip, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescaven .gt. 0) &
      used = g_send_data(wombat%id_fescaven, wombat%fescaven, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fescadet .gt. 0) &
      used = g_send_data(wombat%id_fescadet, wombat%fescadet, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fecoag2det .gt. 0) &
      used = g_send_data(wombat%id_fecoag2det, wombat%fecoag2det, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesources .gt. 0) &
      used = g_send_data(wombat%id_fesources, wombat%fesources, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_fesinks .gt. 0) &
      used = g_send_data(wombat%id_fesinks, wombat%fesinks, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phygrow .gt. 0) &
      used = g_send_data(wombat%id_phygrow, wombat%phygrow, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phyresp .gt. 0) &
      used = g_send_data(wombat%id_phyresp, wombat%phyresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_phymort .gt. 0) &
      used = g_send_data(wombat%id_phymort, wombat%phymort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoograz .gt. 0) &
      used = g_send_data(wombat%id_zoograz, wombat%zoograz, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooresp .gt. 0) &
      used = g_send_data(wombat%id_zooresp, wombat%zooresp, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zoomort .gt. 0) &
      used = g_send_data(wombat%id_zoomort, wombat%zoomort, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooexcr .gt. 0) &
      used = g_send_data(wombat%id_zooexcr, wombat%zooexcr, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_zooslop .gt. 0) &
      used = g_send_data(wombat%id_zooslop, wombat%zooslop, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_detremi .gt. 0) &
      used = g_send_data(wombat%id_detremi, wombat%detremi, model_time, &
          rmask=grid_tmask, is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

    if (wombat%id_caldiss .gt. 0) &
      used = g_send_data(wombat%id_caldiss, wombat%caldiss, model_time, &
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

    if (wombat%id_detfe_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_detfe_sed_remin, wombat%detfe_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_caco3_sed_remin .gt. 0) &
      used = g_send_data(wombat%id_caco3_sed_remin, wombat%caco3_sed_remin, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_zeuphot .gt. 0) &
      used = g_send_data(wombat%id_zeuphot, wombat%zeuphot, model_time, &
          rmask=grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (wombat%id_phy_leup .gt. 0) &
      used = g_send_data(wombat%id_phy_leup, wombat%phy_leup, model_time, &
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
          wombat%htotallo(i,j) = wombat%htotal_scale_lo * wombat%htotal(i,j,1)
          wombat%htotalhi(i,j) = wombat%htotal_scale_hi * wombat%htotal(i,j,1)
          wombat%ahtotallo(i,j) = wombat%htotal_scale_lo * wombat%ahtotal(i,j,1)
          wombat%ahtotalhi(i,j) = wombat%htotal_scale_hi * wombat%ahtotal(i,j,1)
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
          wombat%htotal(:,:,1), &
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
          wombat%ahtotal(:,:,1), &
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
    allocate(wombat%htotal(isd:ied, jsd:jed, 1:nk)); wombat%htotal(:,:,:)=wombat%htotal_in
    allocate(wombat%ahtotallo(isd:ied, jsd:jed)); wombat%ahtotallo(:,:)=0.0
    allocate(wombat%ahtotalhi(isd:ied, jsd:jed)); wombat%ahtotalhi(:,:)=0.0
    allocate(wombat%ahtotal(isd:ied, jsd:jed, 1:nk)); wombat%ahtotal(:,:,:)=wombat%htotal_in
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
    allocate(wombat%f_pchl(isd:ied, jsd:jed, 1:nk)); wombat%f_pchl(:,:,:)=0.0
    allocate(wombat%f_phyfe(isd:ied, jsd:jed, 1:nk)); wombat%f_phyfe(:,:,:)=0.0
    allocate(wombat%f_zoo(isd:ied, jsd:jed, 1:nk)); wombat%f_zoo(:,:,:)=0.0
    allocate(wombat%f_zoofe(isd:ied, jsd:jed, 1:nk)); wombat%f_zoofe(:,:,:)=0.0
    allocate(wombat%f_det(isd:ied, jsd:jed, 1:nk)); wombat%f_det(:,:,:)=0.0
    allocate(wombat%f_detfe(isd:ied, jsd:jed, 1:nk)); wombat%f_detfe(:,:,:)=0.0
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
    allocate(wombat%radmid3d(isd:ied, jsd:jed, 1:nk)); wombat%radmid3d(:,:,:)=0.0
    allocate(wombat%pprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%pprod_gross(:,:,:)=0.0
    allocate(wombat%pprod_gross_2d(isd:ied, jsd:jed)); wombat%pprod_gross_2d(:,:)=0.0
    allocate(wombat%zprod_gross(isd:ied, jsd:jed, 1:nk)); wombat%zprod_gross(:,:,:)=0.0
    allocate(wombat%wdet100(isd:ied, jsd:jed)); wombat%wdet100(:,:)=0.0
    allocate(wombat%export_prod(isd:ied, jsd:jed)); wombat%export_prod(:,:)=0.0
    allocate(wombat%export_inorg(isd:ied, jsd:jed)); wombat%export_inorg(:,:)=0.0
    allocate(wombat%npp2d(isd:ied, jsd:jed)); wombat%npp2d(:,:)=0.0
    allocate(wombat%npp3d(isd:ied, jsd:jed, 1:nk)); wombat%npp3d(:,:,:)=0.0
    allocate(wombat%phy_mumax(isd:ied, jsd:jed, 1:nk)); wombat%phy_mumax(:,:,:)=0.0
    allocate(wombat%phy_mu(isd:ied, jsd:jed, 1:nk)); wombat%phy_mu(:,:,:)=0.0
    allocate(wombat%pchl_mu(isd:ied, jsd:jed, 1:nk)); wombat%pchl_mu(:,:,:)=0.0
    allocate(wombat%phy_lpar(isd:ied, jsd:jed, 1:nk)); wombat%phy_lpar(:,:,:)=0.0
    allocate(wombat%phy_lnit(isd:ied, jsd:jed, 1:nk)); wombat%phy_lnit(:,:,:)=0.0
    allocate(wombat%phy_lfer(isd:ied, jsd:jed, 1:nk)); wombat%phy_lfer(:,:,:)=0.0
    allocate(wombat%phy_dfeupt(isd:ied, jsd:jed, 1:nk)); wombat%phy_dfeupt(:,:,:)=0.0
    allocate(wombat%feIII(isd:ied, jsd:jed, 1:nk)); wombat%feIII(:,:,:)=0.0
    allocate(wombat%felig(isd:ied, jsd:jed, 1:nk)); wombat%felig(:,:,:)=0.0
    allocate(wombat%fecol(isd:ied, jsd:jed, 1:nk)); wombat%fecol(:,:,:)=0.0
    allocate(wombat%feprecip(isd:ied, jsd:jed, 1:nk)); wombat%feprecip(:,:,:)=0.0
    allocate(wombat%fescaven(isd:ied, jsd:jed, 1:nk)); wombat%fescaven(:,:,:)=0.0
    allocate(wombat%fescadet(isd:ied, jsd:jed, 1:nk)); wombat%fescadet(:,:,:)=0.0
    allocate(wombat%fecoag2det(isd:ied, jsd:jed, 1:nk)); wombat%fecoag2det(:,:,:)=0.0
    allocate(wombat%fesources(isd:ied, jsd:jed, 1:nk)); wombat%fesources(:,:,:)=0.0
    allocate(wombat%fesinks(isd:ied, jsd:jed, 1:nk)); wombat%fesinks(:,:,:)=0.0
    allocate(wombat%phygrow(isd:ied, jsd:jed, 1:nk)); wombat%phygrow(:,:,:)=0.0
    allocate(wombat%phyresp(isd:ied, jsd:jed, 1:nk)); wombat%phyresp(:,:,:)=0.0
    allocate(wombat%phymort(isd:ied, jsd:jed, 1:nk)); wombat%phymort(:,:,:)=0.0
    allocate(wombat%zoograz(isd:ied, jsd:jed, 1:nk)); wombat%zoograz(:,:,:)=0.0
    allocate(wombat%zooresp(isd:ied, jsd:jed, 1:nk)); wombat%zooresp(:,:,:)=0.0
    allocate(wombat%zoomort(isd:ied, jsd:jed, 1:nk)); wombat%zoomort(:,:,:)=0.0
    allocate(wombat%zooexcr(isd:ied, jsd:jed, 1:nk)); wombat%zooexcr(:,:,:)=0.0
    allocate(wombat%zooslop(isd:ied, jsd:jed, 1:nk)); wombat%zooslop(:,:,:)=0.0
    allocate(wombat%detremi(isd:ied, jsd:jed, 1:nk)); wombat%detremi(:,:,:)=0.0
    allocate(wombat%caldiss(isd:ied, jsd:jed, 1:nk)); wombat%caldiss(:,:,:)=0.0
    allocate(wombat%no3_prev(isd:ied, jsd:jed, 1:nk)); wombat%no3_prev(:,:,:)=0.0
    allocate(wombat%caco3_prev(isd:ied, jsd:jed, 1:nk)); wombat%caco3_prev(:,:,:)=0.0
    allocate(wombat%det_sed_remin(isd:ied, jsd:jed)); wombat%det_sed_remin(:,:)=0.0
    allocate(wombat%det_btm(isd:ied, jsd:jed)); wombat%det_btm(:,:)=0.0
    allocate(wombat%detfe_sed_remin(isd:ied, jsd:jed)); wombat%detfe_sed_remin(:,:)=0.0
    allocate(wombat%detfe_btm(isd:ied, jsd:jed)); wombat%detfe_btm(:,:)=0.0
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

!    allocate(wombat%yt(isd:ied, jsd:jed)); wombat%yt(:,:)=0.0
    allocate(wombat%zeuphot(isd:ied, jsd:jed)); wombat%zeuphot(:,:)=0.0
    allocate(wombat%phy_leup(isd:ied, jsd:jed)); wombat%phy_leup(:,:)=0.0

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
        wombat%f_pchl, &
        wombat%f_phyfe, &
        wombat%f_zoo, &
        wombat%f_zoofe, &
        wombat%f_det, &
        wombat%f_detfe, &
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
        wombat%radmid3d, &
        wombat%pprod_gross, &
        wombat%pprod_gross_2d, &
        wombat%zprod_gross, &
        wombat%wdet100, &
        wombat%export_prod, &
        wombat%export_inorg, &
        wombat%npp2d, &
        wombat%npp3d, &
        wombat%phy_mumax, &
        wombat%phy_mu, &
        wombat%pchl_mu, &
        wombat%phy_lpar, &
        wombat%phy_lnit, &
        wombat%phy_lfer, &
        wombat%phy_dfeupt, &
        wombat%feIII, &
        wombat%felig, &
        wombat%fecol, &
        wombat%feprecip, &
        wombat%fescaven, &
        wombat%fescadet, &
        wombat%fecoag2det, &
        wombat%fesources, &
        wombat%fesinks, &
        wombat%phygrow, &
        wombat%phyresp, &
        wombat%phymort, &
        wombat%zoograz, &
        wombat%zooresp, &
        wombat%zoomort, &
        wombat%zooexcr, &
        wombat%zooslop, &
        wombat%detremi, &
        wombat%caldiss, &
        wombat%no3_prev, &
        wombat%caco3_prev, &
        wombat%det_sed_remin, &
        wombat%det_btm, &
        wombat%detfe_sed_remin, &
        wombat%detfe_btm, &
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

    deallocate( &
!        wombat%yt, &
        wombat%phy_leup, &
        wombat%zeuphot)

  end subroutine user_deallocate_arrays

end module generic_WOMBATlite
