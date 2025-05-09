
! <CONTACT EMAIL="Charles.Stock@noaa.gov"> Charles Stock 
! </CONTACT>
!
! <OVERVIEW>
! This module contains the generic version of the COBALT 2.0 model: "Carbon Ocean
! Biogeochemistry and Lower Trophics".  COBALT augments the foodweb dynamics 
! in TOPAZ to enable anaylisis of the energy flow through the planktonic
! foodweb and improve the mechanistic resolution of foodweb dynamics that
! influence biogeochemical processes.
! </OVERVIEW>
!<DESCRIPTION>
!       COBALT simulates the biogeochemical cycling of carbon, nitrogen,
!       phosphorous, iron, silica, calcium carbonate, and lithogenic
!       material in the ocean.  The code is built upon the TOPAZ code
!       developed by John Dunne.  The primary changes to TOPAZ are:
!
!          1) the addition of three zooplankton groups
!          2) The addition of bacteria
!          3) The expansion of the dissolved organic nitrogen and 
!             phosphorous groups to include three types each: labile,
!             semi-labile, and refractory
!          4) Constant Stoichiometry by plankton functional type 
!       
!       The primary COBALT reference is:
!
!       Stock, CA, Dunne, JP, John, JG. 2014. Global-scale carbon and energy
!         flows through the marine planktonic food web: An analysis with a'
!         coupled physical-biological model.  Progress in Oceanography 120, 1-18.  
!
!       Version 2.0 has a number of refinements:
!
!           1) Ammonia uptake parameters are now based on the "high-affinity" 
!              settings from Paulot et al., 2015; GBC; 29(8)
!           2) Phytoplankton aggregation is initiated only when growth rates
!              fall below 1/4 maximum values
!           3) The default parameterization has elevated N:P ratios for both
!              diazotrophs and small phytoplankton
!	    4) Remineralization of sinking detritus is now based on the temperature
!              and oxygen dependences described in Laufkotter et al., 2017; O2
!              dependence of other aerobic processes have also been adjusted
!              for consistency.
!           5) The default carbon chemistry routine is now MOCSY
!           6) The iron scavenging has been re-tuned to new atmospheric (Ginoux-AM4),
!              sediment, river and hydrothermal vent (Tagliabue) sources.
!       
!       parameterizations used herein, as well as definitions for all variables
!       and parameters.  The 33 model state variables are:
!
!       alk: alkalinity
!       cadet_arag: calcium carbonate detritus (aragonite)                
!       cadet_calc: calcium carbonate detritus (calcite)                  
!       dic: dissolved inorganic carbon                                   
!       fed: dissolved iron                                               
!       fedi: diazotroph iron                                             
!       felg: large phytoplankton iron
!       fedet: iron detritus                                              
!       fesm: small phytoplankton iron
!       ldon: labile dissolved organic nitrogen                           
!       ldop: labile dissolved organic phosphorous
!       lith: lithogenic aluminosilicate particles                        
!       lithdet: lithogenic detritus                                      
!       nbact: bacteria
!       ndet: nitrogen detritus                                           
!       ndi: diazotroph nitrogen                                          
!       nlg: large phyto nitrogen
!       nsm: high-light adapted small phyto nitrogen
!       nh4: ammonia                                                      
!       no3: nitrate                                                      
!       o2: oxygen                                                        
!       pdet: phosphorous detritus                                        
!       po4: phosphate                                                    
!       srdon: semi-refractory dissolved organic nitrogen
!             (decays over years to decades)
!       srdop: semi-refractory dissolved organic phosphorous
!             (decays over years to decades)
!       sldon: semi-labile dissolved organic nitrogen 
!             (decays on monthly time scales)               
!       sldop: semi-labile dissolved organic phosphorous                
!             (decays on monthly time scales)
!       sidet: silica detritus                                            
!       silg: large phyto silica
!       sio4: silicate                                                    
!       nsmz: small zooplankton nitrogen
!       nmdz: medium zooplankton nitrogen
!       nlgz: large zooplankton nitrogen
!   
!<NAMELIST NAME="generic_COBALT_nml">
!
!  <DATA NAME="do_14c" TYPE="logical">
!  If true, then simulate radiocarbon. Includes 2 prognostic tracers, DI14C
! and DO14C. Requires that do_carbon = .true. Note that 14C is not taken up
! by CaCO3 at the current time, but cycles only through the soft tissue.
! This is a mistake that will be fixed later.
!  </DATA> 
!
!  <DATA NAME="co2_calc" TYPE="character">
!  Defines the carbon equilibration method.  Default is 'ocmip2' which uses
! the FMS_ocmip2_co2calc routine.  The other option is 'mocsy', which uses
! the set of routines authored by J. Orr. See reference at: 
! http://ocmip5.ipsl.jussieu.fr/mocsy/index.html
!
!</NAMELIST>
!
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! Stock, Charles A., John P Dunne, and Jasmin G John, 2014: Global-scale carbon and energy flows
! through the marine food web: an analysis with a coupled physical-biological mode. Progress in Oceanography,
! 120, DOI:10.1016/j.pocean.2013.07.001.
!
! Stock, Charles A., and John P Dunne, 2010: Controls on the ratio of mesozooplankton production to
! primary production in marine ecosystems. Deep-Sea Research, Part I, 57(1), DOI:10.1016/j.dsr.2009.10.006.
!
! Dunne, John P., Jasmin G John, Elena Shevliakova, Ronald J Stouffer, John P Krasting, Sergey Malyshev, 
! P C D Milly, Lori T Sentman, Alistair Adcroft, William F Cooke, Krista A Dunne, Stephen M Griffies,
! Robert Hallberg, Matthew J Harrison, Hiram Levy II, Andrew T Wittenberg, Peter Phillipps, and Niki Zadeh,
! 2013: GFDL's ESM2 global coupled climate-carbon Earth System Models Part II: Carbon system formulation and
! baseline simulation characteristics. Journal of Climate, 26(7), DOI:10.1175/JCLI-D-12-00150.1.
! </REFERENCE>
! <DEVELOPER_NOTES>
! </DEVELOPER_NOTES>
! </INFO>
!----------------------------------------------------------------

module generic_COBALT

  use coupler_types_mod, only: coupler_2d_bc_type
  use field_manager_mod, only: fm_string_len
  use mpp_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod,           only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE
  use mpp_mod,           only: input_nml_file, mpp_error, stdlog, NOTE, WARNING, FATAL, stdout, mpp_chksum
  use time_manager_mod,  only: time_type
  use fm_util_mod,       only: fm_util_start_namelist, fm_util_end_namelist  
  use constants_mod,     only: WTMCO2, WTMO2,WTMN,rdgas,wtmair
  use data_override_mod, only: data_override
  use fms_mod,           only: write_version_number, FATAL, WARNING, stdout, stdlog,mpp_pe,mpp_root_pe
  use fms_mod,           only: check_nml_error

  use g_tracer_utils, only : g_tracer_type,g_tracer_start_param_list,g_tracer_end_param_list
  use g_tracer_utils, only : g_tracer_add,g_tracer_add_param, g_tracer_set_files
  use g_tracer_utils, only : g_tracer_set_values,g_tracer_get_pointer
  use g_tracer_utils, only : g_tracer_get_common,g_tracer_set_common 
  use g_tracer_utils, only : g_tracer_coupler_set,g_tracer_coupler_get
  use g_tracer_utils, only : g_tracer_get_values  
  use g_tracer_utils, only : g_diag_type, g_diag_field_add
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use g_tracer_utils, only : g_send_data, is_root_pe

  use FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector

  implicit none ; private
!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: generic_COBALT.F90,v 20.0.2.1.2.1 2014/09/29 16:40:08 Niki.Zadeh Exp $'
  character(len=128) :: tag = '$Name: bugfix_nnz $'
!-----------------------------------------------------------------------

  character(len=fm_string_len), parameter :: mod_name       = 'generic_COBALT'
  character(len=fm_string_len), parameter :: package_name   = 'generic_cobalt'

  public do_generic_COBALT
  public generic_COBALT_register
  public generic_COBALT_init
  public generic_COBALT_register_diag
  public generic_COBALT_update_from_coupler
  public generic_COBALT_update_from_source
  public generic_COBALT_update_from_bottom
  public generic_COBALT_set_boundary_values
  public generic_COBALT_end
  public as_param_cobalt

  !The following variables for using this module 
  ! are overwritten by generic_tracer_nml namelist
  logical, save :: do_generic_COBALT = .false.
  character(len=10), save :: as_param_cobalt = 'gfdl_cmip6'

  real, parameter :: sperd = 24.0 * 3600.0
  real, parameter :: spery = 365.25 * sperd
  real, parameter :: epsln=1.0e-30
  real,parameter :: missing_value1=-1.0e+10
  real, parameter :: missing_value_diag=-1.0e+10

  real, parameter :: vb_nh3 = 25.

  logical :: do_nh3_diag

! Namelist Options

  character(len=10) ::  co2_calc = 'ocmip2'  ! other option is 'mocsy'
  logical :: do_14c             = .false.
  logical :: debug              = .false.
  logical :: do_nh3_atm_ocean_exchange = .false. 
  real    :: k_nh4_small = 1.e-8
  real    :: k_nh4_diazo = 1.e-7
  real    :: k_nh4_large = 5.e-8
  real    :: k_no3_small = 5.e-7
  real    :: k_no3_diazo = 5.0e-6
  real    :: k_no3_large = 2.5e-6
  real    :: o2_min_nit= 0.01e-6
  real    :: k_o2_nit  = 3.9e-6
  real    :: irr_inhibit = 10.
  real    :: gamma_nitrif= 3.5e6 !month(-1)
  real    :: k_nh3_nitrif= 3.1e-9 !mol/kg
  real    :: imbalance_tolerance=1.0e-10 !the tolerance for non-conservation in C,N,P,Sc,Fe

  integer :: scheme_no3_nh4_lim = 2 !1-Frost and Franzen (1992)
                                    !2-O'Neill

  integer :: scheme_nitrif = 3 !1-default COBALT
                               !2-update with no temperature dependence
                               !3-update with temperature dependence

namelist /generic_COBALT_nml/ do_14c, co2_calc, debug, do_nh3_atm_ocean_exchange, scheme_nitrif, &
     k_nh4_small,k_nh4_large,k_nh4_diazo,scheme_no3_nh4_lim,k_no3_small,k_no3_large,k_no3_diazo, &
     o2_min_nit,k_o2_nit,irr_inhibit,k_nh3_nitrif,gamma_nitrif,imbalance_tolerance

  ! Declare phytoplankton, zooplankton and cobalt variable types, which contain
  ! the vast majority of all variables used in this module. 

  type phytoplankton
     real :: alpha,   &			
          fe_2_n_max,    &
          p_2_n_static,  &
          k_fe_2_n,      &
          k_fed,         &
          k_nh4,         &
          k_no3,         &
          k_po4,         &
          k_sio4,        &
          P_C_max,       &
          si_2_n_max,    &
          si_2_n_static, &
          thetamax,      &     
          bresp,         &
          agg,           &
          frac_mu_agg,   &
          vir,           &            
          exu 
     real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      & 
          jprod_n_new_100,  & 
          jprod_n_n2_100,   &   
          jzloss_n_100,     &
          jaggloss_n_100,   &
          jvirloss_n_100,   &
          jexuloss_n_100,   &
          f_n_100,          &
          juptake_fe_100,   &
          juptake_po4_100,  &
          juptake_sio4_100, &
          nlim_bw_100,      &
          plim_bw_100,      &
          def_fe_bw_100,    &
          irrlim_bw_100        
     real, ALLOCATABLE, dimension(:,:,:)  :: &
          def_fe      , & 
          def_p       , & 
          f_fe        , & 
          f_n         , & 
          felim       , & 
          irrlim      , & 
          jzloss_fe   , & 
          jzloss_n    , & 
          jzloss_p    , & 
          jzloss_sio2 , & 
          jaggloss_fe , &  
          jaggloss_n  , & 
          jaggloss_p  , &
          jaggloss_sio2,&
          agg_lim      ,& 
          jvirloss_fe , & 
          jvirloss_n  , & 
          jvirloss_p  , & 
          jvirloss_sio2,&
          jexuloss_fe , &
          jexuloss_n  , &
          jexuloss_p  , &
          jhploss_fe  , & 
          jhploss_n   , &
          jhploss_p   , & 
          jhploss_sio2, & 
          juptake_n2  , & 
          juptake_fe  , & 
          juptake_nh4 , & 
          juptake_no3 , & 
          juptake_po4 , & 
          juptake_sio4, &
          jprod_n     , & 
          liebig_lim  , & 
          mu          , &
          f_mu_mem    , &
          mu_mix      , & 
          nh4lim      , & 
          no3lim      , & 
          po4lim      , &
          o2lim       , & 
          q_fe_2_n    , & 
          q_p_2_n     , & 
          silim       , & 
          q_si_2_n    , & 
          theta            
     integer ::            &
          id_def_fe       = -1, & 
          id_def_p        = -1, &
          id_felim        = -1, &
          id_irrlim       = -1, &
          id_jzloss_fe    = -1, &
          id_jzloss_n     = -1, & 
          id_jzloss_p     = -1, & 
          id_jzloss_sio2  = -1, &
          id_jaggloss_fe  = -1, &
          id_jaggloss_n   = -1, &
          id_jaggloss_p   = -1, &
          id_jaggloss_sio2= -1, &
          id_agg_lim      = -1, & 
          id_jvirloss_fe  = -1, & 
          id_jvirloss_n   = -1, &
          id_jvirloss_p   = -1, &
          id_jvirloss_sio2= -1, &
          id_jexuloss_n   = -1, &
          id_jexuloss_p   = -1, &
          id_jexuloss_fe  = -1, &
          id_jhploss_fe   = -1, & 
          id_jhploss_n    = -1, & 
          id_jhploss_p    = -1, &
          id_jhploss_sio2 = -1, &
          id_juptake_n2   = -1, &
          id_juptake_fe   = -1, &
          id_juptake_nh4  = -1, &
          id_juptake_no3  = -1, & 
          id_juptake_po4  = -1, &
          id_juptake_sio4 = -1, &
          id_jprod_n      = -1, & 
          id_liebig_lim   = -1, &
          id_mu           = -1, &
          id_f_mu_mem     = -1, &
          id_mu_mix       = -1, &
          id_nh4lim       = -1, &
          id_no3lim       = -1, &
          id_po4lim       = -1, &
          id_o2lim        = -1, &
          id_q_fe_2_n     = -1, &
          id_q_p_2_n      = -1, &
          id_silim        = -1, &
          id_q_si_2_n     = -1, & 
          id_theta        = -1, &
          id_jprod_n_100  = -1, &
          id_jprod_n_new_100  = -1, &     
          id_jprod_n_n2_100 = -1, &
          id_jzloss_n_100     = -1, &
          id_jaggloss_n_100   = -1, &
          id_jvirloss_n_100   = -1, &
          id_jexuloss_n_100   = -1, &
          id_f_n_100          = -1, &
          id_sfc_f_n          = -1, &
          id_sfc_chl          = -1, &
          id_sfc_def_fe       = -1, &
          id_sfc_felim        = -1, &
          id_sfc_q_fe_2_n     = -1, &
          id_sfc_nh4lim       = -1, &
          id_sfc_no3lim       = -1, &
          id_sfc_po4lim       = -1, &
          id_sfc_irrlim       = -1, &
          id_sfc_theta        = -1, &
          id_sfc_mu           = -1
  end type phytoplankton

  type zooplankton
    real ::  &
	  imax,             & ! maximum ingestion rate (sec-1)         
          ki,               & ! half-sat for ingestion (moles N m-3)
          gge_max,          & ! max gross growth efficiciency (approached as i >> bresp, dimensionless)
          nswitch,          & ! switching parameter (dimensionless)
          mswitch,          & ! switching parameter (dimensionless)
          bresp,            & ! basal respiration rate (sec-1)
          ktemp,            & ! temperature dependence of zooplankton rates (C-1)
          phi_det,          & ! fraction of ingested N to detritus
          phi_ldon,         & ! fraction of ingested N/P to labile don
          phi_sldon,        & ! fraction of ingested N/P to semi-labile don
          phi_srdon,        & ! fraction of ingested N/P to semi-refractory don
          phi_ldop,         & ! fraction of ingested N/P to labile dop
          phi_sldop,        & ! fraction of ingested N/P to semi-labile dop
          phi_srdop,        & ! fraction of ingested N/P to semi-refractory dop 
          phi_nh4,          & ! fraction of ingested N to nh4 due to ingestion-related metabolism
          phi_po4,	    & ! fraction of ingested N to po4 due to ingestion-related metabolism
          q_p_2_n,          & ! p:n ratio of zooplankton
          ipa_smp,          & ! innate prey availability of low-light adapt. small phytos 
          ipa_lgp,          & ! innate prey availability of large phytoplankton
          ipa_diaz,         & ! innate prey availability of diazotrophs 
          ipa_smz,          & ! innate prey availability of small zooplankton
          ipa_mdz,          & ! innate prey availability of large zooplankton
          ipa_lgz,          & ! innate prey availability of x-large zooplankton
          ipa_det,          & ! innate prey availability of detritus
          ipa_bact            ! innate prey availability for bacteria
    real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      &
          jingest_n_100,    &
          jzloss_n_100,     &
          jhploss_n_100,    &
          jprod_ndet_100,   &
          jprod_don_100,    &
          jremin_n_100,     &
          f_n_100          
    real, ALLOCATABLE, dimension(:,:,:) :: &
          f_n,              & ! zooplankton biomass
          jzloss_n,         & ! Losses of n due to consumption by other zooplankton groups
          jzloss_p,	    & ! Losses of p due to consumption by other zooplankton groups 
          jhploss_n,        & ! Losses of n due to consumption by unresolved higher preds
          jhploss_p,	    & ! Losses of p due to consumption by unresolved higher preds
          jingest_n,        & ! Total ingestion of n
          jingest_p,        & ! Total ingestion of p
          jingest_sio2,     & ! Total ingestion of silicate
          jingest_fe,	    & ! Total ingestion of iron
          jprod_ndet,       & ! production of nitrogen detritus by zooplankton group 
          jprod_pdet,       & ! production of phosphorous detritus by zooplankton group
          jprod_ldon,       & ! production of labile dissolved organic N by zooplankton group
          jprod_ldop,       & ! production of labile dissolved organic P by zooplankton group
          jprod_srdon,      & ! production of semi-refractory dissolved organic N by zooplankton group
          jprod_srdop,      & ! production of semi-refractory dissolved organic P by zooplankton group 
          jprod_sldon,      & ! production of semi-labile dissolved organic N by zooplankton group
          jprod_sldop,      & ! production of semi-labile dissolved organic P by zooplankton group
          jprod_fed,	    & ! production of dissolved iron
          jprod_fedet,      & ! production of iron detritus
          jprod_sidet,	    & ! production of silica detritus
          jprod_sio4,       & ! production of silicate via rapid dissolution at surface
          jprod_po4,        & ! phosphate production by zooplankton
          jprod_nh4,        & ! ammonia production by zooplankton
          jprod_n,          & ! zooplankton production
          o2lim,            & ! oxygen limitation of zooplankton activity
          temp_lim            ! Temperature limitation
    integer ::		    &
          id_jzloss_n       = -1, &
          id_jzloss_p       = -1, &
          id_jhploss_n      = -1, &
          id_jhploss_p      = -1, &
          id_jingest_n      = -1, &
          id_jingest_p      = -1, &
          id_jingest_sio2   = -1, &
          id_jingest_fe     = -1, &
          id_jprod_ndet     = -1, &
          id_jprod_pdet     = -1, &
          id_jprod_ldon     = -1, &
          id_jprod_ldop     = -1, &
          id_jprod_srdon    = -1, &
          id_jprod_srdop    = -1, &
          id_jprod_sldon    = -1, &
          id_jprod_sldop    = -1, &
          id_jprod_fed      = -1, &
          id_jprod_fedet    = -1, &
          id_jprod_sidet    = -1, &
          id_jprod_sio4     = -1, &
          id_jprod_po4      = -1, &
          id_jprod_nh4      = -1, &
          id_jprod_n        = -1, &
          id_o2lim          = -1, &
          id_temp_lim       = -1, &
          id_jprod_n_100    = -1, &
          id_jingest_n_100  = -1, &
          id_jzloss_n_100   = -1, &
          id_jhploss_n_100  = -1, &
          id_jprod_ndet_100 = -1, &
          id_jprod_don_100  = -1, &
          id_jremin_n_100   = -1, &
          id_f_n_100        = -1
  end type zooplankton

  type bacteria
    real ::  &
          mu_max,           & ! maximum bacterial growth rate (sec-1)
          k_ldon,           & ! half-sat for nitrogen-limited growth (mmoles N m-3)
          gge_max,          & ! max gross growth efficiciency (dimensionless)
          bresp,            & ! basal respiration rate (sec-1)
          ktemp,            & ! temperature dependence of bacterial rates (C-1)
          vir,              & ! virus-driven loss rate for bacteria (sec-1 mmole N m-3)
          q_p_2_n             ! p:n ratio for bacteria 
    real, ALLOCATABLE, dimension(:,:)  :: &
          jprod_n_100,      &
          jzloss_n_100,     &
          jvirloss_n_100,   &
          jremin_n_100,     &
          juptake_ldon_100, &
          f_n_100
    real, ALLOCATABLE, dimension(:,:,:) :: &
          f_n,              & ! bacteria biomass
          jzloss_n,         & ! Losses of n due to consumption by zooplankton 
          jzloss_p,         & ! Losses of p due to consumption by zooplankton
          jhploss_n,        & ! Losses of n due to consumption by unresolved higher preds
          jhploss_p,        & ! Losses of p due to consumption by unresolved higher preds
          jvirloss_n  ,     & ! nitrogen losses via viruses
          jvirloss_p  ,     & ! phosphorous losses via viruses
          juptake_ldon,     & ! Total uptake of ldon
          juptake_ldop,     & ! Total uptake of sldon
          jprod_nh4,        & ! production of ammonia bacteria  
          jprod_po4,        & ! production of phosphate by bacteria 
          jprod_n,          & ! bacterial production
          ldonlim,          & ! limitation due to organic substrate
          o2lim,            & ! limitation due to oxygen	 
          temp_lim            ! Temperature limitation
    integer ::              &
          id_jzloss_n       = -1, &
          id_jzloss_p       = -1, &
          id_jhploss_n      = -1, &
          id_jhploss_p      = -1, &
          id_jvirloss_n     = -1, &
          id_jvirloss_p     = -1, &
          id_juptake_ldon   = -1, &
          id_juptake_ldop   = -1, &
          id_jprod_nh4      = -1, &
          id_jprod_po4      = -1, &
          id_jprod_n        = -1, &
          id_temp_lim       = -1, &
          id_o2lim          = -1, &
          id_ldonlim        = -1, &
          id_jprod_n_100    = -1, &
          id_jzloss_n_100   = -1, &
          id_jvirloss_n_100 = -1, &
          id_jremin_n_100   = -1, &
          id_juptake_ldon_100 = -1, &
          id_f_n_100
  end type bacteria

  integer, parameter :: NUM_PHYTO  = 3
  !
  ! Array allocations and flux calculations assume that phyto(1) is the
  ! only phytoplankton group cabable of nitrogen uptake by N2 fixation while phyto(2:NUM_PHYTO) 
  ! are only cabable of nitrgen uptake by NH4 and NO3 uptake
  !
  integer, parameter :: DIAZO      = 1
  integer, parameter :: LARGE      = 2
  integer, parameter :: SMALL      = 3
  type(phytoplankton), dimension(NUM_PHYTO) :: phyto

  ! define three zooplankton types
  integer, parameter :: NUM_ZOO = 3
  type(zooplankton), dimension(NUM_ZOO) :: zoo

  type(bacteria), dimension(1) :: bact

  integer, parameter :: NUM_PREY = 8

  type generic_COBALT_type

     logical  ::       &
          init,             &                  ! If tracers should be initializated
          force_update_fluxes,&                ! If OCMIP2 tracers fluxes should be updated every coupling timesteps
                                               !    when update_from_source is not called every coupling timesteps
                                               !    as is the case with MOM6  THERMO_SPANS_COUPLING option
          p_2_n_static,     &                  ! If P:N is fixed in phytoplankton
          cased_steady,     &                  ! steady state approximation for cased
          tracer_debug

     real  ::          &
          atm_co2_flux,     &
          c_2_n,            & 
          ca_2_n_arag,      &
          ca_2_n_calc,      & 
          caco3_sat_max,    &
          doc_background,   &
          fe_2_n_upt_fac,   &
          fe_2_n_sed,       &
          ffe_sed_max,      &
          ffe_geotherm_ratio,&
          ffe_iceberg_ratio,&
          fe_coast,         &
          felig_2_don,      &
          felig_bkg ,       &
          gamma_cadet_arag, & 
          gamma_cadet_calc, & 
          gamma_irr_mem,    &
          gamma_mu_mem,     &
          gamma_ndet,       &
          gamma_nitrif,     &
          k_nh3_nitrif,     &
          gamma_sidet,      &
          gamma_srdon,      &
          gamma_srdop,      &
          gamma_sldon,      &
          gamma_sldop,      &
          kappa_sidet,      &
          irr_inhibit,      &
          k_n_inhib_di,     &
          k_o2,             &
          k_o2_nit,         &
          kappa_eppley,     &
          kappa_remin,      &
          remin_ramp_scale, &
          kfe_eq_lig_hl,    &
          kfe_eq_lig_ll,    &
          alpha_fescav,     &
          beta_fescav,      &
          io_fescav,        &
          remin_eff_fedet,  &
          half_life_14c,    &
          lambda_14c,       &
          k_lith,           &
          phi_lith,         &
          alk_2_n_denit,    &
          n_2_n_denit,      &
          k_no3_denit,      &
          z_burial,         &
          phi_surfresp_cased, & 
          phi_deepresp_cased, &
          alpha_cased,      &
          beta_cased,       &
          gamma_cased,      &
          Co_cased,         &
          o2_min,           &
          o2_min_nit,       &
          o2_2_nfix,        & 
          o2_2_nh4,         &
          o2_2_no3,         &
          o2_2_nitrif,      &
          o2_inhib_di_pow,  &
          o2_inhib_di_sat,  &
          rpcaco3,          &
          rplith,           &
          rpsio2,           &
          thetamin,         &
          thetamin_nolim,   &
          vir_ktemp,        &
          lysis_phi_ldon,   &
          lysis_phi_srdon,  &
          lysis_phi_sldon,  & 
          lysis_phi_ldop,   & 
          lysis_phi_srdop,  &
          lysis_phi_sldop,  &
          wsink,            &
          z_sed,            &
          zeta,             &
          refuge_conc,      &
          imax_hp,          & ! unresolved higher pred. max ingestion rate
          ki_hp,            & ! unresolved higher pred. half-sat
          ktemp_hp,         & ! temperature dependence for higher predators
          coef_hp,          & ! scaling between unresolved preds and available prey
          nswitch_hp,	    & ! higher predator switching behavior
          mswitch_hp,       & ! higher predator switching behavior
          hp_ipa_smp,       & ! innate prey availability of small phytos to hp's
          hp_ipa_lgp,       & ! "  "  "  "  "  "  "  "  "   large phytos to hp's
          hp_ipa_diaz,      & ! "  "  "  "  "  "  "  "  "   diazotrophs to hp's  
          hp_ipa_bact,      & ! "  "  "  "  "  "  "  "  "   bacteria to hp's
          hp_ipa_smz,       & ! "  "  "  "  "  "  "  "  "   small zooplankton to hp's
          hp_ipa_mdz,       & ! "  "  "  "  "  "  "  "  "   medium zooplankton to hp's
          hp_ipa_lgz,       & ! "  "  "  "  "  "  "  "  "   large zooplankton to hp's
          hp_ipa_det,       & ! "  "  "  "  "  "  "  "  "   detritus to hp's
          hp_phi_det          ! fraction of ingested N to detritus
          
     real, dimension(3)                    :: total_atm_co2

     real    :: htotal_scale_lo, htotal_scale_hi, htotal_in
     real    :: Rho_0, a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
     real    :: a1_co2, a2_co2, a3_co2, a4_co2, a5_co2
     real    :: a1_o2, a2_o2, a3_o2, a4_o2, a5_o2

     logical, dimension(:,:), ALLOCATABLE ::  &
          mask_z_sat_arag,&
          mask_z_sat_calc

     real, dimension(:,:,:), ALLOCATABLE ::  &
          f_alk,&				! Other prognostic variables
          f_cadet_arag,&
          f_cadet_calc,&
          f_dic,&
          f_fed,&
          f_fedet,&
          f_ldon,&
          f_ldop,&
          f_lith,&
          f_lithdet,&
          f_ndet,&
          f_nh4,&
          f_no3,&
          f_o2,&
          f_pdet,&
          f_po4,&
          f_srdon,&
          f_srdop,&
          f_sldon,&
          f_sldop,&
          f_sidet,&
          f_silg,&
          f_sio4,&
          co3_sol_arag,&
          co3_sol_calc,&
          f_chl,&
          f_nh3,&
          f_co3_ion,&
          f_htotal,&
          f_irr_mem,&
          f_cased,&
          f_cadet_arag_btf,&
          f_cadet_calc_btf,&
          f_fedet_btf, &
          f_lithdet_btf, &
          f_ndet_btf,&
          f_pdet_btf,&
          f_sidet_btf,&
          jnbact,&
          jndi,&
          jnsm,&
          jnlg,&
          jnsmz,&
          jnmdz,&
          jnlgz,&
          jalk,&
          jalk_plus_btm,&
          jcadet_arag,&
          jcadet_calc,&
          jdic,&
          jdic_plus_btm,&
          jdin_plus_btm,&
          jfed,&
          jfed_plus_btm,&
          jfedi,&
          jfelg,&
          jfesm,&
          jfedet,&
          jldon,&
          jldop,&
          jlith,&
          jlithdet,&
          jndet,&
          jnh4,&
          jnh4_plus_btm,&
          jno3,&
          jno3_plus_btm,&
          jo2,&
          jo2_plus_btm,&
          jpdet,&
          jpo4,&
          jpo4_plus_btm,&
          jsrdon,&
          jsrdop,&
          jsldon,&
          jsldop,&
          jsidet,&
          jsilg,&
          jsio4,&                 
          jsio4_plus_btm,&
          jprod_ndet,&
          jprod_pdet,&
          jprod_ldon,&
          jprod_ldop,&
          jprod_sldon,&
          jprod_sldop,&
          jprod_srdon,&
          jprod_srdop,&
          jprod_fed,&
          jprod_fedet,&
          jprod_sidet,&
          jprod_sio4, &
          jprod_lithdet,&
          jprod_cadet_arag,&
          jprod_cadet_calc,&
          jprod_nh4,&
          jprod_nh4_plus_btm,&
          jprod_po4,&
          det_jzloss_n,&
          det_jzloss_p,&
          det_jzloss_fe,&
          det_jzloss_si,&
          det_jhploss_n,&
          det_jhploss_p,&
          det_jhploss_fe,&
          det_jhploss_si,&
          jdiss_cadet_arag,&
          jdiss_cadet_arag_plus_btm,&
          jdiss_cadet_calc,&
          jdiss_cadet_calc_plus_btm,&
          jdiss_sidet,&
          jremin_ndet,&
          jremin_pdet,&
          jremin_fedet,&
          jfe_ads,&
          jfe_coast,&
          kfe_eq_lig,&
          feprime,&
          ligand,&
          fe_sol,&
          expkT,&
          expkreminT,&
          hp_temp_lim,&
          hp_o2lim,&
          hp_jingest_n,&
          hp_jingest_p,&
          hp_jingest_fe,&
          hp_jingest_sio2,&
          irr_inst,&
          irr_mix,&
          jno3denit_wc,&
          jo2resp_wc, &
          jnitrif,&
          omega_arag,&
          omega_calc,&                                                  
          omegaa,&                                                  
          omegac,&                                                  
          tot_layer_int_c,&
          tot_layer_int_fe,&
          tot_layer_int_n,&
          tot_layer_int_p,&
          tot_layer_int_si,&
          tot_layer_int_o2,&
          tot_layer_int_alk,&
          total_filter_feeding,&
          nlg_diatoms,&
          q_si_2_n_lg_diatoms,&
          zt, &
          zm, &
          c14_2_n,&
          f_di14c,&
          f_do14c,&
          fpo14c,&
          j14c_decay_dic,&
          j14c_decay_doc,&
          j14c_reminp,&
          jdi14c,&
          jdo14c, &
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 Ocnbgc 
! CAS: added tot_layer_int_dic 
          dissoc, &
          o2sat, &
          remoc, &
          tot_layer_int_doc,&
          tot_layer_int_poc,&
          tot_layer_int_dic
 
!==============================================================================================================

     real, dimension(:,:), ALLOCATABLE :: &
          b_alk,b_dic,b_fed,b_nh4,b_no3,b_o2,b_po4,b_sio4,b_di14c,&	! bottom flux terms
          co2_csurf,pco2_csurf,co2_alpha,c14o2_csurf,c14o2_alpha,&
          nh3_csurf,nh3_alpha,pnh3_csurf,&
          fcadet_arag_btm,&
          fcadet_calc_btm,&
          ffedet_btm,&
          flithdet_btm,&
          fpdet_btm,&
          fndet_btm,&
          fsidet_btm,&      
          fcased_burial,&
          fcased_redis,&
          fcased_redis_surfresp,&
          cased_redis_coef,&
          cased_redis_delz,&
          ffe_sed,&
          ffe_geotherm,&
          ffe_iceberg,&
          fnfeso4red_sed,&
          fno3denit_sed,&
          fnoxic_sed,&
          frac_burial,&
          fndet_burial,&
          fpdet_burial,&
          jprod_allphytos_100,&
          jprod_diat_100,&
          htotallo, htotalhi,&
          hp_jingest_n_100,&
          hp_jremin_n_100,&
          hp_jprod_ndet_100,&
          jprod_lithdet_100,&
          jprod_sidet_100,&
          jprod_cadet_calc_100,&
          jprod_cadet_arag_100,&
          jprod_mesozoo_200, &
          jremin_ndet_100, &
          f_ndet_100, &
          f_don_100, &
          f_silg_100, &
          f_mesozoo_200, &
          fndet_100, &
          fpdet_100, &
          fsidet_100, &
          fcadet_calc_100, &
          fcadet_arag_100, &
          ffedet_100, &
          flithdet_100, &
          btm_temp,     &
          btm_o2,       &
          btm_htotal,   &
          btm_co3_ion,  &
          btm_co3_sol_arag, &
          btm_co3_sol_calc, &
          cased_2d,     &
          o2min, & 
          z_o2min, & 
          z_sat_arag,&
          z_sat_calc,&
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 Ocnbgc 
          f_alk_int_100, &
          f_dic_int_100, &
          f_din_int_100, &
          f_fed_int_100, &
          f_po4_int_100, &
          f_sio4_int_100, &
          jalk_100, &
          jdic_100, &
          jdin_100, &
          jfed_100, &
          jpo4_100, &
          jsio4_100, &
          jprod_ptot_100, &
          wc_vert_int_c,&
          wc_vert_int_dic,&
          wc_vert_int_doc,&
          wc_vert_int_poc,&
          wc_vert_int_n,&
          wc_vert_int_p,&
          wc_vert_int_fe,&
          wc_vert_int_si,&
          wc_vert_int_o2,&
          wc_vert_int_alk,&
          wc_vert_int_jdiss_sidet,&
          wc_vert_int_jdiss_cadet,&
          wc_vert_int_jo2resp,&
          wc_vert_int_jprod_cadet,&
          wc_vert_int_jno3denit,&
          wc_vert_int_jnitrif,&
          wc_vert_int_juptake_nh4,&
          wc_vert_int_jprod_nh4,&
          wc_vert_int_juptake_no3,&
          wc_vert_int_nfix
!==============================================================================================================

     real, dimension(:,:,:,:), pointer :: &
          p_alk,&
          p_cadet_arag,&
          p_cadet_calc,&
          p_dic,&
          p_di14c,&
          p_do14c,&
          p_fed,&
          p_fedet,&
          p_fedi,&
          p_felg,&
          p_fesm,&
          p_ldon,&
          p_ldop,&
          p_lith,&
          p_lithdet,&          
          p_nbact,&
          p_ndet,&
          p_ndi,&
          p_nlg,&
          p_nsm,&
          p_nh4,&
          p_no3,&
          p_o2,&
          p_pdet,&
          p_po4,&
          p_srdon,&
          p_srdop,&
          p_sldon,&
          p_sldop,&
          p_sidet,&
          p_silg,&
          p_sio4,&
          p_nsmz,&
          p_nmdz,&
          p_nlgz

      real, dimension (:,:), allocatable :: &
          runoff_flux_alk,&
          runoff_flux_dic,&
          runoff_flux_di14c,&
          runoff_flux_lith,&
          runoff_flux_fed,&
          runoff_flux_no3,&
          runoff_flux_ldon,&
          runoff_flux_sldon,&
          runoff_flux_srdon,&
          runoff_flux_ndet,&
          runoff_flux_po4,&
          runoff_flux_ldop,&
          runoff_flux_sldop,&
          runoff_flux_srdop,&
          dry_fed, wet_fed,&
          dry_lith, wet_lith,&
          dry_no3, wet_no3,&
          dry_nh4, wet_nh4,&
          dry_po4, wet_po4, &
          stf_gas_dic,&
          stf_gas_o2,&
          deltap_dic,&
          deltap_o2

     integer :: nkml
     character(len=fm_string_len)          :: file
     character(len=fm_string_len) :: ice_restart_file
     character(len=fm_string_len) :: ocean_restart_file,IC_file

     integer               ::          &
          id_ndi           = -1,       &
          id_nlg           = -1,       &
          id_nsm           = -1,       &
          id_nsmz          = -1,       &
          id_nmdz          = -1,       &
          id_nlgz          = -1,       & 
          id_nbact         = -1,       &
          id_alk           = -1,       &
          id_cadet_arag    = -1,       &
          id_cadet_calc    = -1,       &
          id_dic           = -1,       &
          id_fed           = -1,       &
          id_fedet         = -1,       &
          id_fedi          = -1,       &
          id_felg          = -1,       &
          id_fesm          = -1,       &
          id_ldon          = -1,       &
          id_ldop          = -1,       &
          id_lith          = -1,       &
          id_lithdet       = -1,       &
          id_ndet          = -1,       &
          id_nh4           = -1,       &
          id_no3           = -1,       &
          id_o2            = -1,       &
          id_pdet          = -1,       & 
          id_po4           = -1,       &
          id_srdop         = -1,       &
          id_srdon         = -1,       &
          id_sldon         = -1,       &
          id_sldop         = -1,       &
          id_sidet         = -1,       &
          id_silg          = -1,       &
          id_sio4          = -1,       &
          id_co3_sol_arag  = -1,       &
          id_co3_sol_calc  = -1,       &
          id_dep_dry_fed   = -1,       &
          id_dep_dry_nh4   = -1,       & 
          id_dep_dry_no3   = -1,       &
          id_dep_dry_po4   = -1,       & 
          id_dep_wet_fed   = -1,       & 
          id_dep_wet_nh4   = -1,       &
          id_dep_wet_no3   = -1,       &
          id_dep_wet_po4   = -1,       &
          id_dep_wet_lith  = -1,       &
          id_dep_dry_lith  = -1,       &
          id_omega_arag    = -1,       &
          id_omega_calc    = -1,       &
          id_chl           = -1,       &
          id_co3_ion       = -1,       &
          id_htotal        = -1,       &
          id_irr_mem       = -1,       &
          id_cased         = -1,       &
	  id_cadet_arag_btf = -1,      & 
          id_cadet_calc_btf = -1,      &
          id_fedet_btf     = -1,       & 
          id_lithdet_btf   = -1,       & 
          id_ndet_btf      = -1,       & 
          id_pdet_btf      = -1,       & 
          id_sidet_btf     = -1,       &
          id_jfed          = -1,       &
          id_jprod_ndet    = -1,       &
          id_jprod_pdet    = -1,       &
          id_jprod_sldon   = -1,       &
          id_jprod_ldon    = -1,       &
          id_jprod_srdon   = -1,       &
          id_jprod_sldop   = -1,       &
          id_jprod_ldop    = -1,       &
          id_jprod_srdop   = -1,       &
          id_jprod_fed     = -1,       &
          id_jprod_fedet   = -1,       &
          id_jprod_sidet   = -1,       &
          id_jprod_sio4    = -1,       &
          id_jprod_lithdet = -1,       &
          id_jprod_cadet_arag = -1,    &
          id_jprod_cadet_calc = -1,    & 
          id_jprod_po4     = -1,       &
          id_jprod_nh4     = -1,       &
          id_jprod_nh4_plus_btm = -1,  &
          id_det_jzloss_n  = -1,       &
          id_det_jzloss_p  = -1,       &
          id_det_jzloss_fe = -1,       &
          id_det_jzloss_si = -1,       &
          id_det_jhploss_n = -1,       &
          id_det_jhploss_p = -1,       &
          id_det_jhploss_fe = -1,      &
          id_det_jhploss_si = -1,      &
          id_jdiss_sidet   = -1,       &
          id_jdiss_cadet_arag = -1,    &
          id_jdiss_cadet_arag_plus_btm = -1, &
          id_jdiss_cadet_calc = -1,    &
          id_jdiss_cadet_calc_plus_btm = -1, &
          id_jremin_ndet   = -1,       &
          id_jremin_pdet   = -1,       & 
          id_jremin_fedet  = -1,       &
          id_jfe_ads       = -1,       &
          id_jfe_coast     = -1,       &
          id_kfe_eq_lig    = -1,       &
          id_feprime       = -1,       &
          id_ligand        = -1,       &
          id_fe_sol        = -1,       &
          id_expkT         = -1,       &
          id_expkreminT    = -1,       &
          id_hp_temp_lim   = -1,       &
          id_hp_o2lim      = -1,       &
          id_hp_jingest_n  = -1,       &
          id_hp_jingest_p  = -1,       &
          id_hp_jingest_fe = -1,       &
          id_hp_jingest_sio2 = -1,     &                
          id_irr_inst      = -1,       &
          id_irr_mix       = -1,       &
          id_jalk          = -1,       & 
          id_jalk_plus_btm = -1,       & 
          id_jdic          = -1,       & 
          id_jdic_plus_btm = -1,       & 
          id_jnh4          = -1,       & 
          id_jndet         = -1,       & 
          id_jnh4_plus_btm = -1,       & 
          id_jno3denit_wc  = -1,       &
          id_jo2resp_wc    = -1,       &
          id_jnitrif       = -1,       &
          id_co2_csurf     = -1,       & 
          id_pco2_csurf    = -1,       &
          id_co2_alpha     = -1,       &
          id_nh3_csurf     = -1,       & 
          id_nh3_alpha     = -1,       &
          id_fcadet_arag   = -1,       &
          id_fcadet_calc   = -1,       &
          id_ffedet        = -1,       &
          id_fndet         = -1,       &
          id_fpdet         = -1,       &
          id_fsidet        = -1,       & 
          id_flithdet      = -1,       &
          id_fcadet_arag_btm = -1,     &
          id_fcadet_calc_btm = -1,     &
          id_ffedet_btm    = -1,       &
          id_flithdet_btm  = -1,       &
          id_fndet_btm     = -1,       &
          id_fpdet_btm     = -1,       &
          id_fsidet_btm    = -1,       &
          id_fcased_burial = -1,       &
          id_fcased_redis  = -1,       &
          id_fcased_redis_surfresp  = -1, &
          id_cased_redis_coef  = -1,   &
          id_cased_redis_delz  = -1,   &
          id_ffe_sed       = -1,       &
          id_ffe_geotherm  = -1,       &
          id_ffe_iceberg = -1,         &
          id_fnfeso4red_sed= -1,       &
          id_fno3denit_sed = -1,       &
          id_fnoxic_sed    = -1,       &
          id_frac_burial   = -1,       &
          id_fndet_burial  = -1,       &
          id_fpdet_burial  = -1,       &
          id_nphyto_tot    = -1,       &
          id_no3_in_source = -1,       &
          id_pco2surf      = -1,       &
          id_pnh3surf      = -1,       &
          id_sfc_alk       = -1,       &
          id_sfc_cadet_arag= -1,       & 
          id_sfc_cadet_calc= -1,       & 
          id_sfc_dic       = -1,       & 
          id_sfc_fed       = -1,       & 
          id_sfc_ldon      = -1,       &
          id_sfc_sldon     = -1,       &
          id_sfc_srdon     = -1,       &
          id_sfc_no3       = -1,       &
          id_sfc_nh4       = -1,       &
          id_sfc_po4       = -1,       &
          id_sfc_sio4      = -1,       &
          id_sfc_htotal    = -1,       &
          id_sfc_o2        = -1,       &
          id_sfc_chl       = -1,       &
          id_sfc_irr       = -1,       &
          id_sfc_irr_mem   = -1,       &
          id_sfc_temp      = -1,       &
          id_btm_temp      = -1,       &
          id_btm_o2        = -1,       &
          id_btm_htotal    = -1,       &
          id_btm_co3_sol_arag = -1,    &
          id_btm_co3_sol_calc = -1,    &
          id_btm_co3_ion      = -1,    &
          id_cased_2d      = -1,       &
          id_sfc_co3_ion   = -1,       &
          id_sfc_co3_sol_arag = -1,    &
          id_sfc_co3_sol_calc = -1,    &
          id_runoff_flux_alk = -1,     &
          id_runoff_flux_dic = -1,     &
          id_runoff_flux_di14c = -1,     &
          id_runoff_flux_fed = -1,     &
          id_runoff_flux_lith = -1,    &
          id_runoff_flux_no3 = -1,     &
          id_runoff_flux_ldon = -1,    &
          id_runoff_flux_sldon = -1,   &
          id_runoff_flux_srdon = -1,   &
          id_runoff_flux_ndet = -1,    &
          id_runoff_flux_po4 = -1,     &
          id_runoff_flux_ldop = -1,    &
          id_runoff_flux_sldop = -1,   &
          id_runoff_flux_srdop = -1,   &
          id_tot_layer_int_c = -1,     & 
          id_tot_layer_int_fe = -1,    & 
          id_tot_layer_int_n = -1,     & 
          id_tot_layer_int_p = -1,     & 
          id_tot_layer_int_si = -1,    &
          id_tot_layer_int_o2 = -1,    &
          id_tot_layer_int_alk = -1,   &
          id_wc_vert_int_c = -1,       &
          id_wc_vert_int_dic = -1,     &
          id_wc_vert_int_doc = -1,     &
          id_wc_vert_int_poc = -1,     &
          id_wc_vert_int_n = -1,       &
          id_wc_vert_int_p = -1,       &
          id_wc_vert_int_fe = -1,      &
          id_wc_vert_int_si = -1,      &
          id_wc_vert_int_o2 = -1,      &
          id_wc_vert_int_alk = -1,     &
          id_wc_vert_int_jdiss_sidet = -1, &
          id_wc_vert_int_jdiss_cadet = -1, &
          id_wc_vert_int_jo2resp = -1,     &
          id_wc_vert_int_jprod_cadet = -1, &
          id_wc_vert_int_jno3denit = -1,   &
          id_wc_vert_int_jnitrif = -1,     &
          id_wc_vert_int_juptake_nh4 = -1, &
          id_wc_vert_int_jprod_nh4 = -1, &
          id_wc_vert_int_juptake_no3 = -1, &
          id_wc_vert_int_nfix = -1,        &
          id_total_filter_feeding = -1,&
          id_nlg_diatoms = -1,         &
          id_jprod_allphytos_100 = -1, &
          id_jprod_diat_100 = -1,      &
          id_q_si_2_n_lg_diatoms = -1, &
          id_hp_jingest_n_100 = -1,    &
          id_hp_jremin_n_100 = -1,     &
          id_hp_jprod_ndet_100 = -1,   &
          id_jprod_lithdet_100 = -1,   &
          id_jprod_sidet_100 = -1,     &
          id_jprod_cadet_calc_100 = -1, &
          id_jprod_cadet_arag_100 = -1, &
          id_jprod_mesozoo_200 = -1,   &
          id_jremin_ndet_100 = -1,     &
          id_f_ndet_100 = -1,          &
          id_f_don_100 = -1,           &
          id_f_silg_100 = -1,          &
          id_f_mesozoo_200 = -1,       &
          id_fndet_100 = -1,           &
          id_fpdet_100 = -1,           &
          id_ffedet_100 = -1,          &
          id_fcadet_calc_100 = -1,     &
          id_fcadet_arag_100 = -1,     &
          id_flithdet_100 = -1,        &
          id_fsidet_100 = -1,          &
          id_o2min         = -1,       &
          id_z_o2min       = -1,       &
          id_z_sat_arag    = -1,       & ! Depth of Aragonite saturation
          id_z_sat_calc    = -1,       & ! Depth of Calcite saturation
          id_b_di14c       = -1,       & ! Bottom flux of DI14C
          id_c14_2_n       = -1,       & ! DI14C to PO4 uptake ratio
          id_c14o2_csurf   = -1,       & ! Surface water 14CO2*
          id_c14o2_alpha   = -1,       & ! Surface water 14CO2* solubility 
          id_fpo14c        = -1,       & ! PO14C sinking flux
          id_j14c_decay_dic= -1,       & ! Radioactive decay of DI14C
          id_j14c_decay_doc= -1,       & ! Radioactive decay of DO14C
          id_j14c_reminp   = -1,       & ! 14C particle remineralization layer integral
          id_jdi14c        = -1,       & ! DI14C source layer integral
          id_jdo14c        = -1,       & ! Semilabile DO14C source layer integral
          id_di14c         = -1,       & ! Dissolved inorganic radiocarbon Prognostic tracer
          id_do14c         = -1,       & ! Dissolved organic radiocarbon Prognostic tracer
          id_f_alk_int_100  = -1, &
          id_f_dic_int_100  = -1, &
          id_f_din_int_100  = -1, &
          id_f_fed_int_100  = -1, &
          id_f_po4_int_100  = -1, &
          id_f_sio4_int_100 = -1, &
          id_jo2_plus_btm   = -1, &
          id_jalk_100       = -1, &
          id_jdic_100       = -1, &
          id_jdin_100       = -1, &
          id_jfed_100       = -1, &
          id_jpo4_100       = -1, &
          id_jsio4_100      = -1, &
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem 
          id_thetao         = -1, &        ! for testing
          id_dissic         = -1, & 
          id_dissicnat      = -1, & 
          id_dissicabio     = -1, & 
          id_dissi14cabio   = -1, & 
          id_dissoc         = -1, &
          id_phyc           = -1, &
          id_zooc           = -1, &
          id_bacc           = -1, &
          id_detoc          = -1, &
          id_calc           = -1, &
          id_arag           = -1, &
          id_phydiat        = -1, &
          id_phydiaz        = -1, &
          id_phypico        = -1, &
          id_phymisc        = -1, &
          id_zmicro         = -1, &
          id_zmeso          = -1, &
          id_talk           = -1, &
          id_talknat        = -1, &
          id_ph             = -1, &
          id_phnat          = -1, &
          id_phabio         = -1, &
          id_o2_cmip        = -1, &
          id_o2sat          = -1, &
          id_no3_cmip       = -1, &
          id_pka_nh3        = -1, &
          id_nh4_cmip       = -1, &
          id_po4_cmip       = -1, &
          id_dfe            = -1, &
          id_si             = -1, &
          id_chl_cmip       = -1, &
          id_chldiat        = -1, &
          id_chldiaz        = -1, &
          id_chlpico        = -1, &
          id_chlmisc        = -1, &
          id_poc            = -1, &
          id_pon            = -1, &
          id_pop            = -1, &
          id_bfe            = -1, &
          id_bsi            = -1, &
          id_phyn           = -1, &
          id_phyp           = -1, &
          id_phyfe          = -1, &
          id_physi          = -1, &
          id_co3            = -1, &
          id_co3nat         = -1, &
          id_co3abio        = -1, &
          id_co3satcalc     = -1, &
          id_co3satarag     = -1, &
          id_pp             = -1, &
          id_pnitrate       = -1, &
          id_pphosphate     = -1, &
          id_pbfe           = -1, &
          id_pbsi           = -1, &
          id_parag          = -1, &
          id_pcalc          = -1, &
          id_expc           = -1, &
          id_expn           = -1, &
          id_expp           = -1, &
          id_expfe          = -1, &
          id_expsi          = -1, &
          id_expcalc        = -1, &
          id_exparag        = -1, &
          id_remoc          = -1, &
          id_dcalc          = -1, &
          id_darag          = -1, &
          id_ppdiat         = -1, &
          id_ppdiaz         = -1, &
          id_pppico         = -1, &
          id_ppmisc        = -1, &
          id_bddtdic        = -1, &
          id_bddtdin        = -1, &
          id_bddtdip        = -1, &
          id_bddtdife       = -1, &
          id_bddtdisi       = -1, &
          id_bddtalk        = -1, &
          id_fescav         = -1, &
          id_fediss         = -1, &
          id_graz           = -1, &
          id_dissicos           = -1, & 
          id_dissicnatos        = -1, & 
          id_dissicabioos       = -1, & 
          id_dissi14cabioos     = -1, & 
          id_dissocos           = -1, &
          id_phycos             = -1, &
          id_zoocos             = -1, &
          id_baccos             = -1, &
          id_detocos            = -1, &
          id_calcos             = -1, &
          id_aragos             = -1, &
          id_phydiatos          = -1, &
          id_phydiazos          = -1, &
          id_phypicoos          = -1, &
          id_phymiscos          = -1, &
          id_zmicroos           = -1, &
          id_zmesoos            = -1, &
          id_talkos             = -1, &
          id_talknatos          = -1, &
          id_phos               = -1, &
          id_phnatos            = -1, &
          id_phabioos           = -1, &
          id_o2os               = -1, &
          id_o2satos            = -1, &
          id_no3os              = -1, &
          id_nh4os              = -1, &
          id_po4os              = -1, &
          id_dfeos              = -1, &
          id_sios               = -1, &
          id_chlos              = -1, &
          id_chldiatos          = -1, &
          id_chldiazos          = -1, &
          id_chlpicoos          = -1, &
          id_chlmiscos          = -1, &
          id_ponos              = -1, &
          id_popos              = -1, &
          id_bfeos              = -1, &
          id_bsios              = -1, &
          id_phynos             = -1, &
          id_phypos             = -1, &
          id_phyfeos            = -1, &
          id_physios            = -1, &
          id_co3os              = -1, &
          id_co3natos           = -1, &
          id_co3abioos          = -1, &
          id_co3satcalcos       = -1, &
          id_co3sataragos       = -1, &
          id_limndiat           = -1, &
          id_limndiaz           = -1, &
          id_limnpico           = -1, &
          id_limnmisc           = -1, &
          id_limirrdiat         = -1, &
          id_limirrdiaz         = -1, &
          id_limirrpico         = -1, &
          id_limirrmisc         = -1, &
          id_limfediat          = -1, &
          id_limfediaz          = -1, &
          id_limfepico          = -1, &
          id_limfemisc          = -1, &
          id_limpdiat           = -1, &
          id_limpdiaz           = -1, &
          id_limppico           = -1, &
          id_limpmisc           = -1, &
          id_intpp              = -1, &
          id_intppnitrate       = -1, &
          id_intppdiat          = -1, &
          id_intppdiaz          = -1, &
          id_intpppico          = -1, &
          id_intppmisc          = -1, &
          id_intpbn             = -1, &
          id_intpbp             = -1, &
          id_intpbfe            = -1, &
          id_intpbsi            = -1, &
          id_intpcalcite        = -1, &
          id_intparag           = -1, &
          id_epc100             = -1, &
          id_epn100             = -1, &
          id_epp100             = -1, &
          id_epfe100            = -1, &
          id_epsi100            = -1, &
          id_epcalc100          = -1, &
          id_eparag100          = -1, &
          id_intdic             = -1, &
          id_intdoc             = -1, &
          id_intpoc             = -1, &
          id_spco2              = -1, &
          id_spco2nat           = -1, &
          id_spco2abio          = -1, &
          id_dpco2              = -1, &
          id_dpco2nat           = -1, &
          id_dpco2abio          = -1, &
          id_dpo2               = -1, &
          id_fgco2              = -1, &
          id_fgco2nat           = -1, &
          id_fgco2abio          = -1, &
          id_fg14co2abio        = -1, &
          id_fgo2               = -1, &
          id_icfriver           = -1, &
          id_fric               = -1, &
          id_ocfriver           = -1, &
          id_froc               = -1, &
          id_intpn2             = -1, &
          id_fsn                = -1, &
          id_frn                = -1, &
          id_fsfe               = -1, &
          id_frfe               = -1, &
!          id_o2min              = -1, &  ! previously defined
          id_zo2min             = -1, &
          id_zsatcalc           = -1, &
          id_zsatarag           = -1, &
          id_fddtdic            = -1, &
          id_fddtdin            = -1, &
          id_fddtdip            = -1, &
          id_fddtdife           = -1, &
          id_fddtdisi           = -1, &
          id_fddtalk            = -1, &
          id_fbddtdic           = -1, &
          id_fbddtdin           = -1, &
          id_fbddtdip           = -1, &
          id_fbddtdife          = -1, &
          id_fbddtdisi          = -1, &
          id_fbddtalk           = -1

!==============================================================================================================
  end type generic_COBALT_type

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

  type(generic_COBALT_type) :: cobalt 
  
  type(CO2_dope_vector) :: CO2_dope_vec

  ! identification numbers for mpp clocks
  integer :: id_clock_carbon_calculations
  integer :: id_clock_phyto_growth
  integer :: id_clock_bacteria_growth
  integer :: id_clock_zooplankton_calculations
  integer :: id_clock_other_losses
  integer :: id_clock_production_loop
  integer :: id_clock_ballast_loops
  integer :: id_clock_source_sink_loop1
  integer :: id_clock_source_sink_loop2
  integer :: id_clock_source_sink_loop3
  integer :: id_clock_source_sink_loop4
  integer :: id_clock_source_sink_loop5
  integer :: id_clock_source_sink_loop6
  integer :: id_clock_cobalt_send_diagnostics
  integer :: id_clock_cobalt_calc_diagnostics

contains

  subroutine generic_COBALT_register(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
integer                                                 :: ioun
integer                                                 :: ierr
integer                                                 :: io_status
integer                                                 :: stdoutunit,stdlogunit
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!
    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_register'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '): '
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '): '



! provide for namelist over-ride
! This needs to go before the add_tracers in order to allow the namelist 
! settings to switch tracers on and off.
!
stdoutunit=stdout();stdlogunit=stdlog()

read (input_nml_file, nml=generic_COBALT_nml, iostat=io_status)
ierr = check_nml_error(io_status,'generic_COBALT_nml')

write (stdoutunit,'(/)')
write (stdoutunit, generic_COBALT_nml)
write (stdlogunit, generic_COBALT_nml)

  if (do_14c) then
    write (stdoutunit,*) trim(note_header), 'Simulating radiocarbon'
  endif

  if (trim(co2_calc) == 'ocmip2') then
    write (stdoutunit,*) trim(note_header), 'Using FMS OCMIP2 CO2 routine'
  else if (trim(co2_calc) == 'mocsy') then
    write (stdoutunit,*) trim(note_header), 'Using Mocsy CO2 routine'
  else
    call mpp_error(FATAL,"Unknown co2_calc option specified in generic_COBALT_nml")
  endif
    !Specify all prognostic and diagnostic tracers of this modules.
    call user_add_tracers(tracer_list)

  end subroutine generic_COBALT_register

  !  <SUBROUTINE NAME="generic_COBALT_init">
  !  <OVERVIEW>
  !   Initialize the generic COBALT module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine: 
  !       Adds all the COBALT Tracers to the list of generic Tracers
  !       passed to it via utility subroutine g_tracer_add().
  !
  !       Adds all the parameters used by this module via utility
  !       subroutine g_tracer_add_param().
  !
  !       Allocates all work arrays used in the module. 
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_init(tracer_list)
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_COBALT_init(tracer_list, force_update_fluxes)
    type(g_tracer_type), pointer :: tracer_list
    logical          ,intent(in) :: force_update_fluxes

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_init'

    !There are situations where the column_physics (update_from_source) is not called every timestep 
    ! such as when MOM6 THERMO_SPANS_COUPLING=True , yet we want the fluxes to be updated every timestep
    ! In that case we can force an update by setting the namelist generic_tracer_nml:force_update_fluxes=.true.
    cobalt%force_update_fluxes = force_update_fluxes

    if (do_nh3_atm_ocean_exchange .or. scheme_nitrif.eq.2 .or. scheme_nitrif.eq.3) then
       do_nh3_diag=.true.
    else
       do_nh3_diag=.false.
    end if

    !Specify and initialize all parameters used by this package
    call user_add_params

    !Allocate all the private work arrays used by this module.
    call user_allocate_arrays

    id_clock_carbon_calculations = mpp_clock_id('(Cobalt: carbon calcs)' ,grain=CLOCK_MODULE)
    id_clock_phyto_growth = mpp_clock_id('(Cobalt: phytoplankton growth calcs)',grain=CLOCK_MODULE)
    id_clock_bacteria_growth = mpp_clock_id('(Cobalt: bacteria growth calcs)',grain=CLOCK_MODULE)
    id_clock_zooplankton_calculations = mpp_clock_id('(Cobalt: zooplankton calculations)',grain=CLOCK_MODULE)
    id_clock_other_losses = mpp_clock_id('(Cobalt: other losses)',grain=CLOCK_MODULE)
    id_clock_production_loop = mpp_clock_id('(Cobalt: production loop)',grain=CLOCK_MODULE)
    id_clock_ballast_loops = mpp_clock_id('(Cobalt: ballasting loops)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop1 = mpp_clock_id('(Cobalt: source/sink loop 1)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop2 = mpp_clock_id('(Cobalt: source/sink loop 2)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop3 = mpp_clock_id('(Cobalt: source/sink loop 3)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop4 = mpp_clock_id('(Cobalt: source/sink loop 4)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop5 = mpp_clock_id('(Cobalt: source/sink loop 5)',grain=CLOCK_MODULE)
    id_clock_source_sink_loop6 = mpp_clock_id('(Cobalt: source/sink loop 6)',grain=CLOCK_MODULE)
    id_clock_cobalt_send_diagnostics = mpp_clock_id('(Cobalt: send diagnostics)',grain=CLOCK_MODULE)
    id_clock_cobalt_calc_diagnostics = mpp_clock_id('(Cobalt: calculate diagnostics)',grain=CLOCK_MODULE)

  end subroutine generic_COBALT_init

  !   Register diagnostic fields to be used in this module. 
  !   Note that the tracer fields are automatically registered in user_add_tracers
  !   User adds only diagnostics for fields that are not a member of g_tracer_type
  !
  subroutine generic_COBALT_register_diag(diag_list)
    type(g_diag_type), pointer :: diag_list
    type(vardesc)  :: vardesc_temp
    integer        :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau, axes(3), axesTi(3)
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


    ! Register the diagnostics for the various phytoplankton 
    !
    ! Register Limitation Diagnostics
    !
    vardesc_temp = vardesc("def_fe_Di","Diaz. Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Lg","Large Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("def_fe_Sm","Small Phyto. Fe Deficiency",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Di","Diaz. Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Lg","Large Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("felim_Sm","Small Phyto. Fed uptake Limitation",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Di","Diaz. Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Lg","Large Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irrlim_Sm","Small Phyto. Light Limitation",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Di","Diaz. Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(DIAZO)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Lg","Large Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(LARGE)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("theta_Sm","Small Phyto. Chl:C",'h','L','s','g Chl (g C)-1','f')
    phyto(SMALL)%id_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Di","Diaz. Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(DIAZO)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Lg","Large Phyto. Overall Growth Rate",'h','L','s','s-1','f')
    phyto(LARGE)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_Sm","Small Phyto. Growth Rate",'h','L','s','s-1','f')
    phyto(SMALL)%id_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mem_Di","Diaz. Phyto. Growth Memory",'h','L','s','s-1','f')
    phyto(DIAZO)%id_f_mu_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mem_Lg","Large Phyto. Growth memory",'h','L','s','s-1','f')
    phyto(LARGE)%id_f_mu_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mem_Sm","Small Phyto. Growth Memory",'h','L','s','s-1','f')
    phyto(SMALL)%id_f_mu_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mix_Di","Diaz. Phyto. ML ave",'h','L','s','s-1','f')
    phyto(DIAZO)%id_mu_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mix_Lg","Large Phyto. ML ave",'h','L','s','s-1','f')
    phyto(LARGE)%id_mu_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mu_mix_Sm","Small Phyto. ML ave",'h','L','s','s-1','f')
    phyto(SMALL)%id_mu_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Lg","Ammonia Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Sm","Ammonia Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nh4lim_Di","Ammonia Limitation of Diazo",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Lg","Nitrate Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Sm","Nitrate Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("no3lim_Di","Ammonia Limitation of Diazo",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Di","Phosphate Limitation of Diaz. Phyto",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Lg","Phosphate Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("po4lim_Sm","Phosphate Limitation of Small Phyto",'h','L','s','dimensionless','f')
    phyto(SMALL)%id_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Di","Oxygen Limitation of Diaz. Phyto",'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Di","Fe:N ratio of Diaz. Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(DIAZO)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Lg","Fe:N ratio of Large Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(LARGE)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_fe_2_n_Sm","Fe:N ratio of Small Phyto",'h','L','s','mol Fe/mol N','f')
    phyto(SMALL)%id_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("silim_Lg","SiO4 Limitation of Large Phyto",'h','L','s','dimensionless','f')
    phyto(LARGE)%id_silim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_si_2_n_Lg","Si:N ratio of Large Phyto",'h','L','s','mol Si/mol N','f')
    phyto(LARGE)%id_q_si_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for phytoplankton loss terms: zooplankton 
    ! CAS: loss diagnostics simplified to just N

    vardesc_temp = vardesc("jzloss_n_Di","Diazotroph nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Lg","Large phyto nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Sm","Small phyto nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Register diagnostics for phytoplankton loss terms: aggregation 
    !

    vardesc_temp = vardesc("jaggloss_n_Di","Diazotroph nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_n_Lg","Large phyto nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_n_Sm","Small phyto nitrogen loss to aggregation layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jaggloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("agg_lim_Di","Diazotroph aggregation limitation",&
                           'h','L','s','dimensionless','f')
    phyto(DIAZO)%id_agg_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("agg_lim_Lg","Large phyto aggregation limitation",&
                           'h','L','s','dimensionless','f')
    phyto(LARGE)%id_agg_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("agg_lim_Sm","Small phyto aggregation limitation",&
                           'h','L','s','dimensionless','f')
    phyto(SMALL)%id_agg_lim= register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)


    !
    !  Register diagnostics for phytoplankton loss terms: viruses
    !

    vardesc_temp = vardesc("jvirloss_n_Di","Diazotroph nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_n_Lg","Large phyto nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_n_Sm","Small phyto nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for phytoplankton exudation
    !
    vardesc_temp = vardesc("jexuloss_n_Di","Diazotroph nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_n_Lg","Large phyto nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(LARGE)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_n_Sm","Small phyto nitrogen loss via exudation",&
                           'h','L','s','mol N m-2 s-1','f')
    phyto(SMALL)%id_jexuloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dynamic silicate diagnostics
    !
    vardesc_temp = vardesc("nlg_diatoms","large phytoplankton nitrogen from diatoms",&
                           'h','L','s','mol kg-1','f')
    cobalt%id_nlg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("q_si_2_n_lg_diatoms","Si:N ratio in large diatoms",&
                           'h','L','s','mol Si mol N','f')
    cobalt%id_q_si_2_n_lg_diatoms = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register Phytoplankton Production Diagnostics
    !

    vardesc_temp = vardesc("juptake_n2_Di","Nitrogen fixation layer integral",'h','L','s','mol N m-2 s-1','f')
    phyto(DIAZO)%id_juptake_n2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Di","Diaz. phyto. Fed uptake layer integral",'h','L','s','mol Fe m-2 s-1','f')
    phyto(DIAZO)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Lg","Large phyto. Fed uptake layer integral",'h','L','s','mol Fe m-2 s-1','f')
    phyto(LARGE)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_fe_Sm","Small phyto. Fed uptake layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    phyto(SMALL)%id_juptake_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Di","Diaz. phyto. NH4 uptake layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Lg","Large phyto. NH4 uptake layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    phyto(LARGE)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_nh4_Sm","Small phyto. NH4 uptake layer integral",&
                           'h','L','s','mol NH4 m-2 s-1','f')
    phyto(SMALL)%id_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Di","Diaz. phyto. NO3 uptake layer integral",'h','L','s','mol NO3 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Lg","Large phyto. NO3 uptake layer integral",'h','L','s','mol NO3 m-2 s-1','f')
    phyto(LARGE)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_no3_Sm","Small phyto. NO3 uptake layer integral",&
                           'h','L','s','mol NO3 m-2 s-1','f')
    phyto(SMALL)%id_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Di","Diaz. phyto. PO4 uptake layer integral",'h','L','s','mol PO4 m-2 s-1','f')
    phyto(DIAZO)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Lg","Large phyto. PO4 uptake layer integral",'h','L','s','mol PO4 m-2 s-1','f')
    phyto(LARGE)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_po4_Sm","Small phyto. PO4 uptake layer integral",&
                           'h','L','s','mol PO4 m-2 s-1','f')
    phyto(SMALL)%id_juptake_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_sio4_Lg","Large phyto. SiO4 uptake layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_juptake_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi","Diazotroph Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp","Small phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp","Large phyto. Nitrogen production layer integral",'h','L','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton diagnostics, starting with losses of zooplankton to ingestion by zooplankton
    !

    vardesc_temp = vardesc("jzloss_n_Smz","Small zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Mdz","Medium-sized zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_n_Lgz","Large zooplankton nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for zooplankton loss terms: higher predators
    !

    vardesc_temp = vardesc("jhploss_n_Smz","Small zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_n_Mdz","Medium-sized zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_n_Lgz","Large zooplankton nitrogen loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton ingestion rates
    !

    vardesc_temp = vardesc("jingest_n_Smz","Ingestion of nitrogen by small zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_Mdz","Ingestion of nitrogen by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_Lgz","Ingestion of nitrogen by large zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Smz","Ingestion of phosphorous by small zooplankton, layer integral", &
                           'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Mdz","Ingestion of phosphorous by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_p_Lgz","Ingestion of phosphorous by large zooplankton, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Smz","Ingestion of sio2 by small zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(1)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Mdz","Ingestion of sio2 by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(2)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_sio2_Lgz","Ingestion of sio2 by large zooplankton, layer integral",&
                           'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(3)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Smz","Ingestion of Fe by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Mdz","Ingestion of Fe by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_fe_Lgz","Ingestion of Fe by large zooplankton, layer integral",&
                           'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register detrital production terms for zooplankton
    !

    vardesc_temp = vardesc("jprod_ndet_Smz","Production of nitrogen detritus by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_Mdz","Production of nitrogen detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_Lgz","Production of nitrogen detritus by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Smz","Production of phosphorous detritus by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Mdz","Production of phosphorous detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_pdet_Lgz","Production of phosphorous detritus by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Smz","Production of opal detritus by small zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(1)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Mdz","Production of opal detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(2)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_Lgz","Production of opal detritus by large zooplankton, layer integral",&
                   'h','L','s','mol SiO2 m-2 s-1','f')
    zoo(3)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Smz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(1)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Mdz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(2)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sio4_Lgz","Production of sio4 through grazing/dissolution, layer integral",&
                   'h','L','s','mol SiO4 m-2 s-1','f')
    zoo(3)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Smz","Production of iron detritus by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Mdz","Production of iron detritus by medium zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fedet_Lgz","Production of iron detritus by large zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dissolved organic/inorganic production terms for zooplankton
    !
    ! Labile dissolved organic nitrogen 
    vardesc_temp = vardesc("jprod_ldon_Smz","Production of labile dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldon_Mdz","Production of labile dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldon_Lgz","Production of labile dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_ldop_Smz","Production of labile dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldop_Mdz","Production of labile dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ldop_Lgz","Production of labile dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Refractory dissolved organic nitrogen
    vardesc_temp = vardesc("jprod_srdon_Smz","Production of semi-refractory dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdon_Mdz","Production of semi-refractory dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdon_Lgz","Production of semi-refractory dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_srdop_Smz","Production of semi-refractory dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdop_Mdz","Production of semi-refractory dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_srdop_Lgz","Production of semi-refractory dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! semi-labile dissolved organic nitrogen 
    vardesc_temp = vardesc("jprod_sldon_Smz","Production of semi-labile dissolved organic nitrogen by small zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldon_Mdz","Production of semi-labile dissolved organic nitrogen by medium zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldon_Lgz","Production of semi-labile dissolved organic nitrogen by large zooplankton, layer integral",&
                   'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! semi-labile dissolved organic phosphorous
    vardesc_temp = vardesc("jprod_sldop_Smz","Production of semi-labile dissolved organic phosphorous by small zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(1)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldop_Mdz","Production of semi-labile dissolved organic phosphorous by medium zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(2)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sldop_Lgz","Production of semi-labile dissolved organic phosphorous by large zooplankton, layer integral",&
                   'h','L','s','mol P m-2 s-1','f')
    zoo(3)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! dissolved iron
    vardesc_temp = vardesc("jprod_fed_Smz","Production of dissolved iron by small zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(1)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fed_Mdz","Production of dissolved iron by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(2)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_fed_Lgz","Production of dissolved iron by large zooplankton, layer integral",&
                   'h','L','s','mol Fe m-2 s-1','f')
    zoo(3)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! phosphate
    vardesc_temp = vardesc("jprod_po4_Smz","Production of phosphate by small zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(1)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Mdz","Production of phosphate by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(2)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4_Lgz","Production of phosphate by large zooplankton, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    zoo(3)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! ammonia
    vardesc_temp = vardesc("jprod_nh4_Smz","Production of ammonia by small zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(1)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Mdz","Production of ammonia by medium-sized zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(2)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nh4_Lgz","Production of ammonia by large zooplankton, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    zoo(3)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register zooplankton production terms 
    !

    vardesc_temp = vardesc("jprod_nsmz","Production of new biomass (nitrogen) by small zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nmdz","Production of new biomass (nitrogen) by medium-sized zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgz","Production of new biomass (nitrogen) by large zooplankton, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Smz","Oxygen limitation of small zooplankton",'h','L','s','dimensionless','f')
    zoo(1)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Mdz","Oxygen limitation of medium-sized zooplankton",'h','L','s','dimensionless','f')
    zoo(2)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Lgz","Oxygen limitation of large zooplankton",'h','L','s','dimensionless','f')
    zoo(3)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Smz","Temperature limitation of small zooplankton",'h','L','s','dimensionless','f')
    zoo(1)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Mdz","Temperature limitation of medium-sized zooplankton",'h','L','s','dimensionless','f')
    zoo(2)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Lgz","Temperature limitation of large zooplankton",'h','L','s','dimensionless','f')
    zoo(3)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial diagnostics, starting with losses of bacteria to ingestion by zooplankton
    ! CAS: limit loss terms to N

    vardesc_temp = vardesc("jzloss_n_Bact","Bacterial nitrogen loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register diagnostics for bacteria loss terms: viruses 
    !

    vardesc_temp = vardesc("jvirloss_n_Bact","Bacterial nitrogen loss to viruses layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jvirloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial uptake terms
    !

    vardesc_temp = vardesc("juptake_ldon","Bacterial uptake of labile dissolved organic nitrogen",'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_juptake_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
         
    vardesc_temp = vardesc("juptake_ldop","Bacterial uptake of labile dissolved organic phosphorous",'h','L','s','mol P m-2 s-1','f')
    bact(1)%id_juptake_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register dissolved inorganic production terms for bacteria 
    !
    ! phosphate
    vardesc_temp = vardesc("jprod_po4_Bact","Production of phosphate by bacteria, layer integral",&
                   'h','L','s','mol PO4 m-2 s-1','f')
    bact(1)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! ammonia
    vardesc_temp = vardesc("jprod_nh4_Bact","Production of ammonia by bacteria, layer integral",&
                   'h','L','s','mol NH4 m-2 s-1','f')
    bact(1)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register bacterial production terms
    !

    vardesc_temp = vardesc("jprod_nbact","Production of new biomass (nitrogen) by bacteria, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    bact(1)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("o2lim_Bact","Oxygen limitation of bacteria",'h','L','s','dimensionless','f')
    bact(1)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ldonlim_Bact","ldon limitation of bacteria",'h','L','s','dimensionless','f')
    bact(1)%id_ldonlim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("temp_lim_Bact","Temperature limitation of bacteria",'h','L','s','dimensionless','f')
    bact(1)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Register general COBALT diagnostics
    !

    vardesc_temp = vardesc("co3_sol_arag","Carbonate Ion Solubility for Aragonite",'h','L','s','mol kg-1','f')
    cobalt%id_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("co3_sol_calc","Carbonate Ion Solubility for Calcite",'h','L','s','mol kg-1','f')
    cobalt%id_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_arag","Carbonate Ion Saturation State for Aragonite",'h','L','s','mol kg-1','f')
    cobalt%id_omega_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("omega_calc","Carbonate Ion Saturation State for Calcite",'h','L','s','mol kg-1','f')
    cobalt%id_omega_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! A few overall production diagnostics
    !

    vardesc_temp = vardesc("jprod_cadet_arag","Aragonite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_calc","Calcite CaCO3 production layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_lithdet","Lithogenic detritus production layer integral",'h','L','s','g m-2 s-1','f')
    cobalt%id_jprod_lithdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_sidet","opal detritus production layer integral",'h','L','s','mol SiO2 m-2 s-1','f')
    cobalt%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_sio4","sio4 production layer integral",'h','L','s','mol SiO2 m-2 s-1','f')
    cobalt%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_fedet","Detrital Fedet production layer integral",'h','L','s','mol Fe m-2 s-1','f')
    cobalt%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_ndet","Detrital PON production layer integral",'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_pdet","Detrital phosphorus production layer integral",'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_ldon","labile dissolved organic nitrogen production layer integral",&
                            'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     
    vardesc_temp = vardesc("jprod_ldop","labile dissolved organic phosphorous production layer integral",&
                            'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_srdon","refractory dissolved organic nitrogen production layer integral",&
                            'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_srdop","refractory dissolved organic phosphorous production layer integral",&
                            'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_sldon","semi-labile dissolved organic nitrogen production layer integral",&
                            'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_sldop","semi-labile dissolved organic phosphorous production layer integral",&
                            'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_po4","phosphate production layer integral",&
                            'h','L','s','mol PO4 m-2 s-1','f')
    cobalt%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("jprod_nh4","NH4 production layer integral",'h','L','s','mol NH4 m-2 s-1','f')
    cobalt%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
! CAS added a "plus_btm" version of jprod_nh4 to use for the remoc CMIP variable
    vardesc_temp = vardesc("jprod_nh4_plus_btm","NH4 production layer integral plus bottom fluxes",'h','L','s','mol NH4 m-2 s-1','f')
    cobalt%id_jprod_nh4_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
   
    !
    ! loss diagnostics: detrital loss terms
    !
    
    vardesc_temp = vardesc("det_jzloss_n","nitrogen detritus loss to zooplankton layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    cobalt%id_det_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    vardesc_temp = vardesc("det_jhploss_n","nitrogen detritus loss to higher predators layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    cobalt%id_det_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
       init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    
    !
    ! Loss diagnostics: dissolution and remineralization
    !

    vardesc_temp = vardesc("jdiss_sidet","SiO2 detritus dissolution, layer integral",&
                           'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdiss_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdiss_cadet_arag","CaCO3 detritus dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! CAS added diagnostic including bottom fluxes for cmip5
    vardesc_temp = vardesc("jdiss_cadet_arag_plus_btm","CaCO3 detritus dissolution plus bottom dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_arag_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdiss_cadet_calc","CaCO3 detritus dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! CAS added diagnostic including bottom fluxes for cmip5
    vardesc_temp = vardesc("jdiss_cadet_calc_plus_btm","CaCO3 detritus dissolution plus bottom dissolution, layer integral", &
                           'h','L','s','mol CaCO3 m-2 s-1','f')
    cobalt%id_jdiss_cadet_calc_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_ndet","Nitrogen detritus remineralization, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    cobalt%id_jremin_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_pdet","Phosphorous detritus remineralization, layer integral",&
                           'h','L','s','mol P m-2 s-1','f')
    cobalt%id_jremin_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_fedet","Iron detritus remineralization, layer integral",&
                           'h','L','s','mol m-2 s-1','f')
    cobalt%id_jremin_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! iron cycling diagnostics 
    !

    vardesc_temp = vardesc("jprod_fed","dissolved iron production layer integral",&
                            'h','L','s','mol Fe m-2 s-1','f')
    cobalt%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
        
    vardesc_temp = vardesc("jfed","Dissolved Iron Change layer integral",'h','L','s','mol Fe m-2 s-1','f')
    cobalt%id_jfed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfe_ads","Iron adsorption layer integral",'h','L','s','mol Fe m-2 s-1','f')
    cobalt%id_jfe_ads = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jfe_coast","Coastal iron efflux layer integral",'h','L','s','mol Fe m-2 s-1','f')
    cobalt%id_jfe_coast = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("kfe_eq_lig","Effective ligand binding strength",'h','L','s','kg mol-1','f')
    cobalt%id_kfe_eq_lig = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("feprime","Free iron concentration",'h','L','s','mol kg-1','f')
    cobalt%id_feprime = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ligand","ligand concentration",'h','L','s','mol kg-1','f')
    cobalt%id_ligand = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fe_sol","iron solubility",'h','L','s','mol kg-1','f')
    cobalt%id_fe_sol = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Temperature limitation diagnostics
    !

    vardesc_temp = vardesc("expkT","Eppley temperature limitation factor",'h','L','s','dimensionless','f')
    cobalt%id_expkT = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("expkreminT","Detritus remineralization temperature limitation factor",'h','L','s','dimensionless','f')
    cobalt%id_expkreminT = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("hp_o2lim","Oxygen limitation of higher predators",'h','L','s','dimensionless','f')
    cobalt%id_hp_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("hp_temp_lim","Temperature limitation of higher predators",'h','L','s','dimensionless','f')
    cobalt%id_hp_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Some additional light field diagnostics
    !

    vardesc_temp = vardesc("irr_inst","Instantaneous Light",'h','L','s','W m-2','f')
    cobalt%id_irr_inst = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("irr_mix","Light averaged over mixing layer",'h','L','s','W m-2','f')
    cobalt%id_irr_mix = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Nitrification/Denitrification diagnostics
    !

    vardesc_temp = vardesc("jnitrif","Nitrification layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jnitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jno3denit_wc","Water column Denitrification layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jno3denit_wc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Track total aerobic respiration in the water column
    !

    vardesc_temp = vardesc("jo2resp_wc","Water column aerobic respiration layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jo2resp_wc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Some useful total layer integrals
    !

    vardesc_temp = vardesc("nphyto_tot","Total NO3: Di+Lg+Sm",'h','L','s','mol m-2 s-1','f')
    cobalt%id_nphyto_tot = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_c","Total Carbon (DIC+OC+IC) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_fe","Total Iron (Fed_OFe) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_n","Total Nitrogen (NO3+NH4+ON) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_p","Total Phosphorus (PO4+OP) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_si","Total Silicon (SiO4+SiO2) boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_si = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_o2","Total oxygen boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("tot_layer_int_alk","Total alkalinity boxwise layer integral",'h','L','s','mol m-2','f')
    cobalt%id_tot_layer_int_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("total_filter_feeding","Total filter feeding by large organisms",'h','L','s','mol N m-2 s-1','f')
    cobalt%id_total_filter_feeding = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    !  Save river, depositon and bulk elemental fluxes
    !

    vardesc_temp = vardesc("dep_dry_fed","Dry Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_lith","Dry Deposition of Lithogenic Material",'h','1','s','g m-2 s-1','f')
    cobalt%id_dep_dry_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_nh4","Dry Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_no3","Dry Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_dry_po4","Dry Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_dry_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_fed","Wet Deposition of Iron to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_lith","Wet Deposition of Lithogenic Material",'h','1','s','g m-2 s-1','f')
    cobalt%id_dep_wet_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_nh4","Wet Deposition of Ammonia to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_no3","Wet Deposition of Nitrate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("dep_wet_po4","Wet Deposition of Phosphate to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_dep_wet_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("pka_nh3","pKa of NH3",'h','1','s','','f')
    cobalt%id_pka_nh3     = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_alk","Alkalinity runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_dic","Dissolved Inorganic Carbon runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_di14c","Dissolved Inorganic Carbon 14 runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_di14c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_fed","Iron runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_lith","Lithogenic runoff flux to the ocean",'h','1','s','g m-2 s-1','f')
    cobalt%id_runoff_flux_lith = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_no3","Nitrate runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ldon","LDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_sldon","SLDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_srdon","SRDON runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ndet","NDET runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_po4","PO4 runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_ldop","LDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_sldop","SLDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("runoff_flux_srdop","SRDOP runoff flux to the ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_runoff_flux_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 3D sinking information 
    ! CAS: Note that passing 1D axesTi acts as a flag to save the tracers on the grid interface
    !      as per the note from Niki in the Bling code where saving on interfaces was developed:
    ! Niki: The register_diag_field interface needs to be extended to take the MOM6 axes_grp as argument
    !      instead of this integer array axes_grp%handle
    !      Currently the actual MOM6 diag axes is chosen to be T or Tl based on the size of the axes argument, 2 or 3.
    !      The actual values of these axes argument are not used, only their size is checked to determine the diag axes!
    !      This is not correct since axesTi and axesTl are both of size 3, likewise there are many axes of size 2.
    !      To accomodate axesTi with the least amount of code modification we can set and check for an input array of size 1.    

    vardesc_temp = vardesc("fcadet_arag","CaCO3 sinking flux",'h','1','s','mol m-2 s-1','f')
    axesTi(:)=0
    cobalt%id_fcadet_arag = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc","CaCO3 sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet","fedet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet","lithdet sinking flux",'h','1','s','g m-2 s-1','f')
    cobalt%id_flithdet = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet","ndet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet","pdet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet","sidet sinking flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet = register_diag_field(package_name, vardesc_temp%name, axesTi(1:1),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 2D sinking, bottom source/sink and burial diagnostics
    !

    vardesc_temp = vardesc("fcadet_arag_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_arag_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_btm","CaCO3 sinking flux at bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_burial","CaCO3 permanent burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_redis","CaCO3 redissolution from sediments",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_redis = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcased_redis_surfresp","CaCO3 redissolution rom sediments, surfresp",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcased_redis_surfresp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("cased_redis_coef","CaCO3 redissolution from sediments, deepresp coefficient, ",'h','1','s','s-1','f')
    cobalt%id_cased_redis_coef = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("cased_redis_delz","CaCO3 redissolution from sediments, effective depth",'h','1','s','none (0-1)','f')
    cobalt%id_cased_redis_delz = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_btm","fedet sinking flux burial",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffe_sed","Sediment iron efflux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffe_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffe_geotherm","Geothermal iron efflux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffe_geotherm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffe_iceberg","iceberg iron efflux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffe_iceberg = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_btm","Lithogenic detrital sinking flux burial",'h','1','s','g m-2 s-1','f')
    cobalt%id_flithdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_btm","ndet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnfeso4red_sed","Sediment Ndet Fe and SO4 reduction flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fnfeso4red_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fno3denit_sed","Sediment denitrification flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fno3denit_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fnoxic_sed","Sediment oxic Ndet remineralization flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fnoxic_sed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_btm","pdet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_btm","sidet sinking flux to bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("frac_burial","fraction of organic matter buried",'h','1','s','dimensionless','f')
    cobalt%id_frac_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fndet_burial","ndet burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_burial","pdet burial flux",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_burial = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Surface Diagnostics
    !

    vardesc_temp = vardesc("pco2surf","Oceanic pCO2",'h','1','s','uatm','f')
    cobalt%id_pco2surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("pnh3surf","Oceanic pNH3",'h','1','s','uatm','f')
    cobalt%id_pnh3surf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_alk","Surface Alkalinity",'h','1','s','eq kg-1','f')
    cobalt%id_sfc_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_arag","Surface Detrital Aragonite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_cadet_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_cadet_calc","Surface Detrital Calcite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_cadet_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_dic","Surface Dissolved Inorganic Carbon",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_fed","Surface Dissolved Iron",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ldon","Surface Labile Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_sldon","Surface semi-labile Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_srdon","Surface semi-refractory Dissolved Organic Nitrogen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3","Surface NO3",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4","Surface NH4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4","Surface PO4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_sio4","Surface SiO4",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_htotal","Surface Htotal",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_o2","Surface Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl","Surface Chl",'h','1','s','ug kg-1','f')
    cobalt%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr","Surface Irradiance",'h','1','s','W m-2','f')
    cobalt%id_sfc_irr = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irr_mem","Surface Irradiance memory",'h','1','s','W m-2','f')
    cobalt%id_sfc_irr_mem = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_temp","Surface Temperature",'h','1','s','deg C','f')
    cobalt%id_sfc_temp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_temp","Bottom Temperature",'h','1','s','deg C','f')
    cobalt%id_btm_temp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_o2","Bottom Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_btm_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_htotal","Bottom Htotal",'h','1','s','mol kg-1','f')
    cobalt%id_btm_htotal = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_co3_ion","Bottom Carbonate Ion",'h','1','s','mol kg-1','f')
    cobalt%id_btm_co3_ion = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_co3_sol_arag","Bottom Aragonite Solubility",'h','1','s','mol kg-1','f')
    cobalt%id_btm_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("btm_co3_sol_calc","Bottom Calcite Solubility",'h','1','s','mol kg-1','f')
    cobalt%id_btm_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("cased_2d","calcium carbonite in sediment",'h','1','s','mol m-3','f')
    cobalt%id_cased_2d = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_ion","Surface Carbonate Ion",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_ion = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_arag","Surface Carbonate Ion Solubility for Aragonite",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_sol_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_co3_sol_calc","Surface Carbonate Ion Solubility for Calcite ",'h','1','s','mol kg-1','f')
    cobalt%id_sfc_co3_sol_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nsmp","Surface small phyto. nitrogen",'h','1','s','mol kg-1','f')
    phyto(SMALL)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nlgp","Surface large phyto. nitrogen",'h','1','s','mol kg-1','f')
    phyto(LARGE)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_ndi","Surface diazotroph nitrogen",'h','1','s','mol kg-1','f')
    phyto(DIAZO)%id_sfc_f_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_smp","Surface small phyto. chlorophyll",'h','1','s','ug kg-1','f')
    phyto(SMALL)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_lgp","Surface large phyto. chlorophyll",'h','1','s','ug kg-1','f')
    phyto(LARGE)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_chl_di","Surface diazotroph chlorophyll",'h','1','s','mol kg-1','f')
    phyto(DIAZO)%id_sfc_chl = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_smp","Surface small phyto. iron deficiency",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_lgp","Surface large phyto. iron deficiency",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_def_fe_di","Surface diazotroph iron deficiency",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_def_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_smp","Surface small phyto. iron uptake limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_lgp","Surface large phyto. iron uptake limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_felim_di","Surface diazotroph iron uptake limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_felim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_di","Surface diazotroph iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(DIAZO)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_smp","Surface small phyto. iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(SMALL)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_q_fe_2_n_lgp","Surface large phyto. iron:nitrogen",'h','1','s','moles Fe (moles N)-1','f')
    phyto(LARGE)%id_sfc_q_fe_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_smp","Surface small phyto. light limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_lgp","Surface large phyto. light limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_irrlim_di","Surface diazotroph light limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_irrlim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_smp","Surface small phyto. Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(SMALL)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_lgp","Surface large phyto. Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(LARGE)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_theta_di","Surface diazotroph Chl:C",'h','1','s','g Chl (g C)-1','f')
    phyto(DIAZO)%id_sfc_theta = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

   vardesc_temp = vardesc("sfc_mu_smp","Surface small phyto. Chl:C",'h','1','s','sec-1','f')
    phyto(SMALL)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_mu_lgp","Surface large phyto. Chl:C",'h','1','s','sec-1','f')
    phyto(LARGE)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_mu_di","Surface diazotroph growth rate",'h','1','s','sec-1','f')
    phyto(DIAZO)%id_sfc_mu = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3lim_smp","Surface small phyto. nitrate limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3lim_lgp","Surface large phyto. nitrate limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_no3lim_di","Surface diazotroph nitrate limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_no3lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4lim_smp","Surface small phyto. ammonia limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4lim_lgp","Surface large phyto. ammonia limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_nh4lim_di","Surface diazotroph ammonia limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_nh4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_smp","Surface small phyto. phosphate limitation",'h','1','s','dimensionsless','f')
    phyto(SMALL)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_lgp","Surface large phyto. phosphate limitation",'h','1','s','dimensionless','f')
    phyto(LARGE)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("sfc_po4lim_di","Surface diazotroph phosphate limitation",'h','1','s','dimensionless','f')
    phyto(DIAZO)%id_sfc_po4lim = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 100m integrated fluxes
    !

    vardesc_temp = vardesc("jprod_allphytos_100","Total Nitrogen prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_allphytos_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

! CAS: Added diagnostic for diatom NPP in top 100m for CMIP
    vardesc_temp = vardesc("jprod_diat_100","Diatom prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_diat_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_100","Diazotroph nitrogen prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp_100","Small phyto. nitrogen  prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp_100","Large phyto. nitrogen  prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_new_100","Diazotroph new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmp_new_100","Small phyto. new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgp_new_100","Large phyto. new (NO3-based) prim. prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jprod_n_new_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndi_n2_100","Diazotroph nitrogen fixation in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jprod_n_n2_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_ndi_100","Diazotroph nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nsmp_100","Small phyto. nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nlgp_100","Large phyto. nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_nsmp_100","Small phyto. nitrogen aggregation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jaggloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jaggloss_nlgp_100","Large phyto. nitrogen aggregation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jaggloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_nsmp_100","Small phyto. nitrogen virus loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jvirloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_ndi_100","Diazotroph nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(DIAZO)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_nsmp_100","Small phyto. nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(SMALL)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jexuloss_nlgp_100","Large phyto. nitrogen exudation loss integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    phyto(LARGE)%id_jexuloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nsmz_100","Small zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nmdz_100","Medium zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nlgz_100","Large zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nsmz_100","Small zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nmdz_100","Medium zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_nlgz_100","Large zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nsmz_100","Small zooplankton nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nmdz_100","Medium zooplankton nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_nmdz_100","Medium zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jhploss_nlgz_100","Large zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_nmdz_100","Medium zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_nlgz_100","Large zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_don_nsmz_100","Small zooplankton dissolved org. nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jprod_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_don_nmdz_100","Medium zooplankton dissolved org. nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jprod_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nsmz_100","Small zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(1)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nmdz_100","Medium zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(2)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nlgz_100","Large zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(3)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_hp_100","Higher predator nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jingest_n_hp_100","Higher predator ingestion of nitrogen integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_ndet_hp_100","Higher predator nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_hp_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_nbact_100","Bacteria nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jzloss_nbact_100","Bacteria nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jvirloss_nbact_100","Bacteria nitrogen loss to viruses integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jvirloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_n_nbact_100","Bacteria nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("juptake_ldon_nbact_100","Bacterial uptake of labile dissolved org. nitrogen in upper 100m",'h','1','s','mol m-2 s-1','f')
    bact(1)%id_juptake_ldon_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_lithdet_100","Lithogenic detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_lithdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_sidet_100","Silica detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_sidet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_calc_100","Calcite detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_cadet_arag_100","Aragonite detritus production integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_cadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jremin_ndet_100","Remineralization of nitrogen detritus integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jremin_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jprod_mesozoo_200","Mesozooplankton Production, 200m integration",'h','1','s','mol m-2 s-1','f')
    cobalt%id_jprod_mesozoo_200 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Water column integrated fluxes
    !

    vardesc_temp = vardesc("wc_vert_int_jdiss_sidet","Water column silica dissolution vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jdiss_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jdiss_cadet","Water column calcium carbonate dissolution vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jdiss_cadet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jo2resp","Water column oxygen respired vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jo2resp = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jprod_cadet","Water column calcium carbonate production vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jprod_cadet = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jno3denit","Water column denitrification vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jno3denit = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jnitrif","Water column nitrification vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jnitrif = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_juptake_nh4"," Water column ammonia based NPP vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_juptake_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_jprod_nh4"," Water column ammonia production vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_juptake_no3","Water column nitrate based NPP, vertical integral",'h','1','s','mol m-2 s-1','f')
    cobalt%id_wc_vert_int_juptake_no3 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_nfix","Nitrogen fixation vertical integral",'h','1','s','mol N m-2 s-1','f')
    cobalt%id_wc_vert_int_nfix = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! 100m integrated biomass
    !

    vardesc_temp = vardesc("nsmp_100","Small phytoplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(SMALL)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nlgp_100","Large phytoplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(LARGE)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ndi_100","Diazotroph nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    phyto(DIAZO)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nsmz_100","Small zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(1)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nmdz_100","Medium zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(2)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nlgz_100","Large zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(3)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("nbact_100","Bacterial nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    bact(1)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("silgp_100","Large phytoplankton silicon biomass in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_silg_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ndet_100","Nitrogen detritus biomass in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("don_100","Dissolved organic nitrogen (sr+sl+l) in upper 100m",'h','1','s','mol m-2','f')
    cobalt%id_f_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("mesozoo_200","Mesozooplankton biomass, 200m integral",'h','1','s','mol m-2','f')
    cobalt%id_f_mesozoo_200 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Water column integrated tracers
    !

    vardesc_temp = vardesc("wc_vert_int_c","Total carbon (DIC+OC+IC) vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_dic","Dissolved inorganic carbon vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_doc","Total dissolved organic carbon  vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_doc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_poc","Total particulate organic carbon vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_poc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_n","Total nitrogen vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_n = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_p","Total phosphorus vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_p = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_fe","Total iron vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_si","Total silicon vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_si = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_o2","Total oxygen vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_o2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("wc_vert_int_alk","Total alkalinity vertical integral",'h','1','s','mol m-2','f')
    cobalt%id_wc_vert_int_alk = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     
    !
    ! sinking flux = 100m
    !

    vardesc_temp = vardesc("fndet_100","Nitrogen detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fpdet_100","Phosphorous detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fpdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("ffedet_100","Iron detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ffedet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fsidet_100","Silicon detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsidet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_calc_100","Calcite detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_calc_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("fcadet_arag_100","Aragonite detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fcadet_arag_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("flithdet_100","Lithogenic detritus sinking flux @ 100m",'h','1','s','mol m-2 s-1','f')
    cobalt%id_flithdet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    ! Oxygen minima (value and location
    !

    vardesc_temp = vardesc("o2min","Minimum Oxygen",'h','1','s','mol kg-1','f')
    cobalt%id_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("z_o2min","Depth of Oxygen minimum",'h','1','s','m','f')
    cobalt%id_z_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    !
    ! Calcite and aragonite saturation depths
    !

    vardesc_temp = vardesc("z_sat_arag","Depth of Aragonite Saturation",'h','1','s','m','f')
    cobalt%id_z_sat_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)

    vardesc_temp = vardesc("z_sat_calc","Depth of Calcite Saturation",'h','1','s','m','f')
    cobalt%id_z_sat_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         mask_variant=.TRUE.)

      if (do_14c) then                                        !<<RADIOCARBON
    vardesc_temp = vardesc&
    ("b_di14c","Bottom flux of DI14C into sediment",'h','1','s','mol m-2 s-1','f')
    cobalt%id_b_di14c = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14_2_n","Ratio of DI14C to N-nutrients",'h','L','s','mol kg-1','f')
    cobalt%id_c14_2_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14o2_alpha","Saturation surface 14CO2* per uatm",'h','1','s','mol kg-1 atm-1','f')
    cobalt%id_c14o2_alpha = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("c14o2_csurf","14CO2* concentration at surface",'h','1','s','mol kg-1','f')
    cobalt%id_c14o2_csurf = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("fpo14c","PO14C sinking flux at layer bottom",'h','L','s','mol m-2 s-1','f')
    cobalt%id_fpo14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_decay_dic","DI14C radioactive decay layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_j14c_decay_dic = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_decay_doc","DO14C radioactive decay layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_j14c_decay_doc = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("j14c_reminp","Sinking PO14C remineralization layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_j14c_reminp = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jdi14c","DI14C source layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdi14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc&
    ("jdo14c","DO14C source layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdo14c = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
      endif                                                   !RADIOCARBON>>


    !
    ! Additional diagnostics added for debugging jgj 2015/10/26
    !

    vardesc_temp = vardesc("jalk","Alkalinity source layer integral",'h','L','s','eq m-2 s-1','f')
    cobalt%id_jalk = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jalk_plus_btm","Alkalinity source plus btm layer integral",'h','L','s','eq m-2 s-1','f')
    cobalt%id_jalk_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdic","Dissolved Inorganic Carbon source layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdic = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jdic_plus_btm","Dissolved Inorganic Carbon source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jdic_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jnh4","NH4 source layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jnh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jndet","NDET source layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jnh4_plus_btm","NH4 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jnh4_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

    vardesc_temp = vardesc("jo2_plus_btm","O2 source plus btm layer integral",'h','L','s','mol m-2 s-1','f')
    cobalt%id_jo2_plus_btm = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)

!==============================================================================================================
! 2016/07/05 jgj register and send temperature as a test

    vardesc_temp = vardesc("temp","Potential Temperature",'h','L','s','Celsius','f')
    cobalt%id_thetao = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="thetao", cmor_units="C",                          &
         cmor_standard_name="sea_water_potential_temperature",              &
         cmor_long_name ="Sea Water Potential Temperature")

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Oyr/Omon/day: Marine Biogeochemical Fields

    vardesc_temp = vardesc("dissic_raw","Dissolved Inorganic Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dissic = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissic", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_in_sea_water",  &
         cmor_long_name="Dissolved Inorganic Carbon Concentration")

    vardesc_temp = vardesc("dissicnat_raw","Natural Dissolved Inorganic Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dissicnat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissicnat", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_natural_analogue_in_sea_water",  &
         cmor_long_name="Natural Dissolved Inorganic Carbon Concentration")

    vardesc_temp = vardesc("dissicabio_raw","Abiotic Dissolved Inorganic Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dissicabio = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissicabio", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_abiotic_analogue_in_sea_water",  &
         cmor_long_name="Abiotic Dissolved Inorganic Carbon Concentration")

! 2017/12/28 jgj Updated standard name in data request: added _abiotic_analogue
    vardesc_temp = vardesc("dissi14cabio_raw","Abiotic Dissolved Inorganic 14Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dissi14cabio = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissi14cabio", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon14_abiotic_analogue_in_sea_water", &
         cmor_long_name="Abiotic Dissolved Inorganic 14Carbon Concentration")

    vardesc_temp = vardesc("dissoc_raw","Dissolved Organic Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dissoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissoc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_organic_carbon_in_sea_water",  &
         cmor_long_name="Dissolved Organic Carbon Concentration")

    vardesc_temp = vardesc("phyc_raw","Phytoplankton Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_phyc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phyc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Phytoplankton Carbon Concentration")

    vardesc_temp = vardesc("zooc_raw","Zooplankton Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_zooc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zooc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Zooplankton Carbon Concentration")

    vardesc_temp = vardesc("bacc_raw","Bacterial Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_bacc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bacc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_bacteria_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Bacterial Carbon Concentration")

    vardesc_temp = vardesc("detoc_raw","Detrital Organic Carbon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_detoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="detoc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_organic_detritus_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Detrital Organic Carbon Concentration")

    vardesc_temp = vardesc("calc_raw","Calcite Concentration",'h','L','s','mol m-3','f')
    cobalt%id_calc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="calc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_calcite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Calcite Concentration")

    vardesc_temp = vardesc("arag_raw","Aragonite Concentration",'h','L','s','mol m-3','f')
    cobalt%id_arag = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="arag", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_aragonite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Aragonite Concentration")

    vardesc_temp = vardesc("phydiat_raw","Mole Concentration of Diatoms expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phydiat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phydiat", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_diatoms_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Diatoms expressed as Carbon in sea water")

    vardesc_temp = vardesc("phydiaz_raw","Mole Concentration of Diazotrophs expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phydiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phydiaz", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_diazotrophs_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Diazotrophs expressed as Carbon in sea water")

    vardesc_temp = vardesc("phypico_raw","Mole Concentration of Picophytoplankton expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phypico = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phypico", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_picophytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Picophytoplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("phymisc_raw","Mole Concentration of Miscellaneous Phytoplankton expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phymisc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phymisc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_miscellaneous_phytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Miscellaneous Phytoplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("zmicro_raw","Mole Concentration of Microzooplankton expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_zmicro = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zmicro", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_microzooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Microzooplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("zmeso_raw","Mole Concentration of Mesozooplankton expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_zmeso = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zmeso", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Mesozooplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("talk_raw","Total Alkalinity",'h','L','s','mol m-3','f')
    cobalt%id_talk = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="talk", cmor_units="mol m-3",                          &
         cmor_standard_name="sea_water_alkalinity_expressed_as_mole_equivalent", &
         cmor_long_name="Total Alkalinity")

    vardesc_temp = vardesc("talknat_raw","Natural Total Alkalinity",'h','L','s','mol m-3','f')
    cobalt%id_talknat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="talknat", cmor_units="mol m-3",                          &
         cmor_standard_name="sea_water_alkalinity_natural_analogue_expressed_as_mole_equivalent", &
         cmor_long_name="Natural Total Alkalinity")

    vardesc_temp = vardesc("ph_raw","pH",'h','L','s','1','f')
    cobalt%id_ph = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ph", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_reported_on_total_scale", &
         cmor_long_name="pH")

    vardesc_temp = vardesc("phnat_raw","Natural pH",'h','L','s','1','f')
    cobalt%id_phnat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phnat", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_natural_analogue_reported_on_total_scale", &
         cmor_long_name="Natural pH")

    vardesc_temp = vardesc("phabio_raw","Abiotic pH",'h','L','s','1','f')
    cobalt%id_phabio = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phabio", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_abiotic_analogue_reported_on_total_scale", &
         cmor_long_name="Abiotic pH")

!! same name in model and CMOR, but different units - use _cmip for now
    vardesc_temp = vardesc("o2_raw","Dissolved Oxygen Concentration",'h','L','s','mol m-3','f')
    cobalt%id_o2_cmip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="o2_cmip", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_molecular_oxygen_in_sea_water", &
         cmor_long_name="Dissolved Oxygen Concentration")

    vardesc_temp = vardesc("o2sat_raw","Dissolved Oxygen Concentration at Saturation",'h','L','s','mol m-3','f')
    cobalt%id_o2sat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="o2sat", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_at_saturation", &
         cmor_long_name="Dissolved Oxygen Concentration at Saturation")

!! same name in model and CMOR, but different units - use _cmip for now
    vardesc_temp = vardesc("no3_raw","Dissolved Nitrate Concentration",'h','L','s','mol m-3','f')
    cobalt%id_no3_cmip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="no3_cmip", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_nitrate_in_sea_water", &
         cmor_long_name="Dissolved Nitrate Concentration")

!! same name in model and CMOR, but different units - use for now
    vardesc_temp = vardesc("nh4_raw","Dissolved Ammonium Concentration",'h','L','s','mol m-3','f')
    cobalt%id_nh4_cmip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="nh4_cmip", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_ammonium_in_sea_water", &
         cmor_long_name="Dissolved Ammonium Concentration")

!! same name in model and CMOR, but different units - use _cmip for now
    vardesc_temp = vardesc("po4_raw","Total Dissolved Inorganic Phosphorus Concentration",'h','L','s','mol m-3','f')
    cobalt%id_po4_cmip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="po4_cmip", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_phosphorus_in_sea_water", &
         cmor_long_name="Total Dissolved Inorganic Phosphorus Concentration")

    vardesc_temp = vardesc("dfe_raw","Dissolved Iron Concentration",'h','L','s','mol m-3','f')
    cobalt%id_dfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dfe", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_iron_in_sea_water", &
         cmor_long_name="Dissolved Iron Concentration")

    vardesc_temp = vardesc("si_raw","Total Dissolved Inorganic Silicon Concentration",'h','L','s','mol m-3','f')
    cobalt%id_si = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="si", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_silicon_in_sea_water", &
         cmor_long_name="Total Dissolved Inorganic Silicon Concentration")

!! same name in model and CMOR, but different units - use _cmip for now
    vardesc_temp = vardesc("chl_raw","Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water",'h','L','s','kg m-3','f')
    cobalt%id_chl_cmip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chl_cmip", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chldiat_raw","Mass Concentration of Diatoms expressed as Chlorophyll in sea water",'h','L','s','kg m-3','f')
    cobalt%id_chldiat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chldiat", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_diatoms_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Mass Concentration of Diatoms expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chldiaz_raw","Mass Concentration of Diazotrophs expressed as Chlorophyll in sea water",'h','L','s','kg m-3','f')
    cobalt%id_chldiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chldiaz", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_diazotrophs_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Mass Concentration of Diazotrophs expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chlpico_raw","Mass Concentration of Picophytoplankton expressed as Chlorophyll in sea water",'h','L','s','kg m-3','f')
    cobalt%id_chlpico = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chlpico", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_picophytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Mass Concentration of Picophytoplankton expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chlmisc_raw","Mass Concentration of Other Phytoplankton expressed as Chlorophyll in sea water",'h','L','s','kg m-3','f')
    cobalt%id_chlmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chlmisc", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_miscellaneous_phytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Mass Concentration of Other Phytoplankton expressed as Chlorophyll in sea water")

! 2017/11/27 not in data request
    vardesc_temp = vardesc("poc_raw","Mole Concentration of Particulate Organic Matter expressed as Carbon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_poc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="poc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Mole Concentration of Particulate Organic Matter expressed as Carbon in sea water")

    vardesc_temp = vardesc("pon_raw","Mole Concentration of Particulate Organic Matter expressed as Nitrogen in sea water",'h','L','s','mol N m-3','f')
    cobalt%id_pon = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pon", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_nitrogen_in_sea_water", &
         cmor_long_name="Mole Concentration of Particulate Organic Matter expressed as Nitrogen in sea water")

    vardesc_temp = vardesc("pop_raw","Mole Concentration of Particulate Organic Matter expressed as Phosphorus in sea water",'h','L','s','mol m-3','f')
    cobalt%id_pop = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pop", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_phosphorus_in_sea_water", &
         cmor_long_name="Mole Concentration of Particulate Organic Matter expressed as Phosphorus in sea water")

    vardesc_temp = vardesc("bfe_raw","Mole Concentration of Particulate Organic Matter expressed as Iron in sea water",'h','L','s','mol m-3','f')
    cobalt%id_bfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bfe", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_iron_in_sea_water", &
         cmor_long_name="Mole Concentration of Particulate Organic Matter expressed as Iron in sea water")

! CHECK3:
! 2017/12/28 jgj Updated standard name in data request to include _organic
    vardesc_temp = vardesc("bsi_raw","Mole Concentration of Particulate Organic Matter expressed as Silicon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_bsi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bsi", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_silicon_in_sea_water", &
         cmor_long_name="Mole Concentration of Particulate Organic Matter expressed as Silicon in sea water")

    vardesc_temp = vardesc("phyn_raw","Mole Concentration of Total Phytoplankton expressed as Nitrogen in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phyn = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phyn", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_nitrogen_in_sea_water", &
         cmor_long_name="Mole Concentration of Total Phytoplankton expressed as Nitrogen in sea water")

    vardesc_temp = vardesc("phyp_raw","Mole Concentration of Total Phytoplankton expressed as Phosphorus in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phyp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phyp", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_phosphorus_in_sea_water", &
         cmor_long_name="Mole Concentration of Total Phytoplankton expressed as Phosphorus in sea water")

    vardesc_temp = vardesc("phyfe_raw","Mole Concentration of Total Phytoplankton expressed as Iron in sea water",'h','L','s','mol m-3','f')
    cobalt%id_phyfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phyfe", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_iron_in_sea_water", &
         cmor_long_name="Mole Concentration of Total Phytoplankton expressed as Iron in sea water")

    vardesc_temp = vardesc("physi_raw","Mole Concentration of Total Phytoplankton expressed as Silicon in sea water",'h','L','s','mol m-3','f')
    cobalt%id_physi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="physi", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_silicon_in_sea_water", &
         cmor_long_name="Mole Concentration of Total Phytoplankton expressed as Silicon in sea water")

    vardesc_temp = vardesc("co3_raw","Carbonate Ion Concentration",'h','L','s','mol m-3','f')
    cobalt%id_co3 = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Carbonate Ion Concentration")

    vardesc_temp = vardesc("co3nat_raw","Natural Carbonate Ion Concentration",'h','L','s','mol m-3','f')
    cobalt%id_co3nat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3nat", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_natural_analogue_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Natural Carbonate Ion Concentration")

    vardesc_temp = vardesc("co3abio_raw","Abiotic Carbonate Ion Concentration",'h','L','s','mol m-3','f')
    cobalt%id_co3abio = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3abio", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_abiotic_analogue_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Abiotic Carbonate Ion Concentration")

! 2018/01/17 jgj Updated standard name in data request to match equivalent surface variable
    vardesc_temp = vardesc("co3satcalc_raw","Mole Concentration of Carbonate Ion in Equilibrium with Pure Calcite in sea water",'h','L','s','mol m-3','f')
    cobalt%id_co3satcalc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3satcalc", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_at_equilibrium_with_pure_calcite_in_sea_water", &
         cmor_long_name="Mole Concentration of Carbonate Ion in Equilibrium with Pure Calcite in sea water")

! 2018/01/17 jgj Updated standard name in data request to match equivalent surface variable
    vardesc_temp = vardesc("co3satarag_raw","Mole Concentration of Carbonate Ion in Equilibrium with Pure Aragonite in sea water",'h','L','s','mol m-3','f')
    cobalt%id_co3satarag = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3satarag", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_at_equilibrium_with_pure_aragonite_in_sea_water", &
         cmor_long_name="Mole Concentration of Carbonate Ion in Equilibrium with Pure Aragonite in sea water")

!------------------------------------------------------------------------------------------------------------------
! 3-D rates
! CHECK3: all GFDL and CMOR units

    vardesc_temp = vardesc("pp_raw","Primary Carbon Production by Total Phytoplankton",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pp", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_net_primary_production", &
         cmor_long_name="Primary Carbon Production by Total Phytoplankton")

    vardesc_temp = vardesc("pnitrate_raw","Primary Carbon Production by Phytoplankton due to Nitrate Uptake Alone",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pnitrate = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pnitrate", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_nitrate_utilization", &
         cmor_long_name="Primary Carbon Production by Phytoplankton due to Nitrate Uptake Alone")

! Not requested
!    vardesc_temp = vardesc("pphosphate_raw","Primary Carbon Production by Phytoplankton due to Phosphorus",'h','L','s','mol m-3 s-1','f')
!    cobalt%id_pphosphate = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
!         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
!         cmor_field_name="pphosphate", cmor_units="mol m-3 s-1",                          &
!         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_phosphorus", &
!         cmor_long_name="Primary Carbon Production by Phytoplankton due to Phosphorus")

    vardesc_temp = vardesc("pbfe_raw","Biogenic Iron Production",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pbfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pbfe", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_iron_in_sea_water_due_to_biological_production", &
         cmor_long_name="Biogenic Iron Production")

    vardesc_temp = vardesc("pbsi_raw","Biogenic Silicon Production",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pbsi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pbsi", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_silicon_in_sea_water_due_to_biological_production", &
         cmor_long_name="Biogenic Silicon Production")

! Oyr only
    vardesc_temp = vardesc("pcalc_raw","Calcite Production",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pcalc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pcalc", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_calcite_expressed_as_carbon_in_sea_water_due_to_biological_production", &
         cmor_long_name="Calcite Production")

! Oyr only
    vardesc_temp = vardesc("parag_raw","Aragonite Production",'h','L','s','mol m-3 s-1','f')
    cobalt%id_parag = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="parag", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_aragonite_expressed_as_carbon_in_sea_water_due_to_biological_production", &
         cmor_long_name="Aragonite Production")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello, volume:volcello
    vardesc_temp = vardesc("expc_raw","Sinking Particulate Organic Carbon Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expc", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Sinking Particulate Organic Carbon Flux")

    vardesc_temp = vardesc("expn_raw","Sinking Particulate Organic Nitrogen Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expn = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expn", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_nitrogen_in_sea_water", &
         cmor_long_name="Sinking Particulate Organic Nitrogen Flux")

    vardesc_temp = vardesc("expp_raw","Sinking Particulate Organic Phosphorus Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expp = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expp", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_phosphorus_in_sea_water", &
         cmor_long_name="Sinking Particulate Organic Phosphorus Flux")

    vardesc_temp = vardesc("expfe_raw","Sinking Particulate Iron Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expfe = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expfe", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_iron_in_sea_water", &
         cmor_long_name="Sinking Particulate Iron Flux")

    vardesc_temp = vardesc("expsi_raw","Sinking Particulate Silicon Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expsi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expsi", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_silicon_in_sea_water", &
         cmor_long_name="Sinking Particulate Silicon Flux")

    vardesc_temp = vardesc("expcalc_raw","Sinking Calcite Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_expcalc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="expcalc", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_calcite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Sinking Calcite Flux")

    vardesc_temp = vardesc("exparag_raw","Sinking Aragonite Flux",'h','L','s','mol m-2 s-1','f')
    cobalt%id_exparag = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="exparag", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_aragonite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Sinking Aragonite Flux")

    vardesc_temp = vardesc("remoc_raw","Remineralization of Organic Carbon",'h','L','s','mol m-3 s-1','f')
    cobalt%id_remoc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="remoc", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic matter_expressed_as_carbon_in_sea_water_due_to_remineralization", &
         cmor_long_name="Remineralization of Organic Carbon")

    vardesc_temp = vardesc("dcalc_raw","Calcite Dissolution",'h','L','s','mol m-3 s-1','f')
    cobalt%id_dcalc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dcalc", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_calcite_expressed_as_carbon_in_sea_water_due_to_dissolution", &
         cmor_long_name="Calcite Dissolution")

    vardesc_temp = vardesc("darag_raw","Aragonite Dissolution",'h','L','s','mol m-3 s-1','f')
    cobalt%id_darag = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="darag", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_aragonite_expressed_as_carbon_in_sea_water_due_to_dissolution", &
         cmor_long_name="Aragonite Dissolution")

! CHECK3:
    vardesc_temp = vardesc("ppdiat_raw","Net Primary Organic Carbon Production by Diatoms",'h','L','s','mol m-3 s-1','f')
    cobalt%id_ppdiat = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ppdiat", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_net_primary_production_by_diatoms", &
         cmor_long_name="Net Primary Organic Carbon Production by Diatoms")

! CHECK3:
    vardesc_temp = vardesc("ppdiaz_raw","Net Primary Mole Productivity of Carbon by Diazotrophs",'h','L','s','mol m-3 s-1','f')
    cobalt%id_ppdiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ppdiaz", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_net_primary_production_by_diazotrophs", &
         cmor_long_name="Net Primary Mole Productivity of Carbon by Diazotrophs")

! CHECK3:
    vardesc_temp = vardesc("pppico_raw","Net Primary Mole Productivity of Carbon by Picophytoplankton",'h','L','s','mol m-3 s-1','f')
    cobalt%id_pppico = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="pppico", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_net_primary_production_by_picophytoplankton", &
         cmor_long_name="Net Primary Mole Productivity of Carbon by Picophytoplankton")

! CHECK3:
    vardesc_temp = vardesc("ppmisc_raw","Net Primary Organic Carbon Production by Other Phytoplankton",'h','L','s','mol m-3 s-1','f')
    cobalt%id_ppmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ppmisc", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_net_primary_production_by_miscellaneous_phytoplankton", &
         cmor_long_name="Net Primary Organic Carbon Production by Other Phytoplankton")

    vardesc_temp = vardesc("bddtdic_raw","Rate of Change of Dissolved Inorganic Carbon due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtdic = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtdic", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_inorganic_carbon_in_sea_water_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Carbon due to Biological Activity")

    vardesc_temp = vardesc("bddtdin_raw","Rate of Change of Nitrogen Nutrient due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtdin = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtdin", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_inorganic_nitrogen_in_sea_water_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Nitrogen Nutrient due to Biological Activity")

    vardesc_temp = vardesc("bddtdip_raw","Rate of Change of Dissolved Phosphorus due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtdip = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtdip", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_inorganic_phosphorus_in_sea_water_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Phosphorus due to Biological Activity")

    vardesc_temp = vardesc("bddtdife_raw","Rate of Change of Dissolved Inorganic Iron due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtdife = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtdife", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_inorganic_iron_in_sea_water_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Iron due to Biological Activity")

    vardesc_temp = vardesc("bddtdisi_raw","Rate of Change of Dissolved Inorganic Silicon due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtdisi = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtdisi", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_inorganic_silicon_in_sea_water_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Silicon due to Biological Activity")

    vardesc_temp = vardesc("bddtalk_raw","Rate of Change of Alkalinity due to Biological Activity",'h','L','s','mol m-3 s-1','f')
    cobalt%id_bddtalk = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bddtalk", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_sea_water_alkalinity_expressed_as_mole_equivalent_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Alkalinity due to Biological Activity")

    vardesc_temp = vardesc("fescav_raw","Nonbiogenic Iron Scavenging",'h','L','s','mol m-3 s-1','f')
    cobalt%id_fescav = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fescav", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_iron_in_sea_water_due_to_scavenging_by_inorganic_particles", &
         cmor_long_name="Nonbiogenic Iron Scavenging")

    vardesc_temp = vardesc("fediss_raw","Particle Source of Dissolved Iron",'h','L','s','mol m-3 s-1','f')
    cobalt%id_fediss = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fediss", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_dissolved_iron_in_sea_water_due_to_dissolution_from_inorganic_particles", &
         cmor_long_name="Particle Source of Dissolved Iron")

! CHECK3:
! 2017/08/04 jgj: CMOR requires area:areacello, volume:volcello
    vardesc_temp = vardesc("graz_raw","Total Grazing of Phytoplankton by Zooplankton",'h','L','s','mol m-3 s-1','f')
    cobalt%id_graz = register_diag_field(package_name, vardesc_temp%name, axes(1:3), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="graz", cmor_units="mol m-3 s-1",                          &
         cmor_standard_name="tendency_of_mole_concentration_of_particulate_organic_matter_expressed_as_carbon_in_sea_water_due_to_grazing_of_phytoplankton", &
         cmor_long_name="Total Grazing of Phytoplankton by Zooplankton")

!------------------------------------------------------------------------------------------------------------------
! 3-D Limitation terms
! 2018/07/18 limitation terms are now 2-D

    vardesc_temp = vardesc("limndiat_raw","Nitrogen Limitation of Diatoms",'h','1','s','1','f')
    cobalt%id_limndiat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limndiat", cmor_units="1",                          &
         cmor_standard_name="nitrogen_growth_limitation_of_diatoms", &
         cmor_long_name="Nitrogen Limitation of Diatoms")

    vardesc_temp = vardesc("limndiaz_raw","Nitrogen Limitation of Diazotrophs",'h','1','s','1','f')
    cobalt%id_limndiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limndiaz", cmor_units="1",                          &
         cmor_standard_name="nitrogen_growth_limitation_of_diazotrophs", &
         cmor_long_name="Nitrogen Limitation of Diazotrophs")

    vardesc_temp = vardesc("limnpico_raw","Nitrogen Limitation of Picophytoplankton",'h','1','s','1','f')
    cobalt%id_limnpico = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limnpico", cmor_units="1",                          &
         cmor_standard_name="nitrogen_growth_limitation_of_picophytoplankton", &
         cmor_long_name="Nitrogen Limitation of Picophytoplankton")

    vardesc_temp = vardesc("limnmisc_raw","Nitrogen Limitation of Other Phytoplankton",'h','1','s','1','f')
    cobalt%id_limnmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limnmisc", cmor_units="1",                          &
         cmor_standard_name="nitrogen_growth_limitation_of_miscellaneous_phytoplankton", &
         cmor_long_name="Nitrogen Limitation of Other Phytoplankton")

    vardesc_temp = vardesc("limirrdiat_raw","Irradiance Limitation of Diatoms",'h','1','s','1','f')
    cobalt%id_limirrdiat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limirrdiat", cmor_units="1",                          &
         cmor_standard_name="growth_limitation_of_diatoms_due_to_solar_irradiance", &
         cmor_long_name="Irradiance Limitation of Diatoms")

    vardesc_temp = vardesc("limirrdiaz_raw","Irradiance Limitation of Diazotrophs",'h','1','s','1','f')
    cobalt%id_limirrdiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limirrdiaz", cmor_units="1",                          &
         cmor_standard_name="growth_limitation_of_diazotrophs_due_to_solar_irradiance", &
         cmor_long_name="Irradiance Limitation of Diazotrophs")

    vardesc_temp = vardesc("limirrpico_raw","Irradiance Limitation of Picophytoplankton",'h','1','s','1','f')
    cobalt%id_limirrpico = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limirrpico", cmor_units="1",                          &
         cmor_standard_name="growth_limitation_of_picophytoplankton_due_to_solar_irradiance", &
         cmor_long_name="Irradiance Limitation of Picophytoplankton")

    vardesc_temp = vardesc("limirrmisc_raw","Irradiance Limitation of Other Phytoplankton",'h','1','s','1','f')
    cobalt%id_limirrmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limirrmisc", cmor_units="1",                          &
         cmor_standard_name="growth_limitation_of_miscellaneous_phytoplankton_due_to_solar_irradiance", &
         cmor_long_name="Irradiance Limitation of Other Phytoplankton")

    vardesc_temp = vardesc("limfediat_raw","Iron Limitation of Diatoms",'h','1','s','1','f')
    cobalt%id_limfediat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limfediat", cmor_units="1",                          &
         cmor_standard_name="iron_growth_limitation_of_diatoms", &
         cmor_long_name="Iron Limitation of Diatoms")

    vardesc_temp = vardesc("limfediaz_raw","Iron Limitation of Diazotrophs",'h','1','s','1','f')
    cobalt%id_limfediaz = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limfediaz", cmor_units="1",                          &
         cmor_standard_name="iron_growth_limitation_of_diazotrophs", &
         cmor_long_name="Iron Limitation of Diazotrophs")

    vardesc_temp = vardesc("limfepico_raw","Iron Limitation of Picophytoplankton",'h','1','s','1','f')
    cobalt%id_limfepico = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limfepico", cmor_units="1",                          &
         cmor_standard_name="iron_growth_limitation_of_picophytoplankton", &
         cmor_long_name="Iron Limitation of Picophytoplankton")

    vardesc_temp = vardesc("limfemisc_raw","Iron Limitation of Other Phytoplankton",'h','1','s','1','f')
    cobalt%id_limfemisc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limfemisc", cmor_units="1",                          &
         cmor_standard_name="iron_growth_limitation_of_miscellaneous_phytoplankton", &
         cmor_long_name="Iron Limitation of Other Phytoplankton")

!-- added by JGJ - not requested by CMIP6/OMIP
    vardesc_temp = vardesc("limpdiat_raw","Phosphorus Limitation of Diatoms",'h','1','s','1','f')
    cobalt%id_limpdiat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limpdiat", cmor_units="1",                          &
         cmor_standard_name="phosphorus_growth_limitation_of_diatoms", &
         cmor_long_name="Phosphorus Limitation of Diatoms")

    vardesc_temp = vardesc("limpdiaz_raw","Phosphorus Limitation of Diazotrophs",'h','1','s','1','f')
    cobalt%id_limpdiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limpdiaz", cmor_units="1",                          &
         cmor_standard_name="phosphorus_growth_limitation_of_diazotrophs", &
         cmor_long_name="Phosphorus Limitation of Diazotrophs")

    vardesc_temp = vardesc("limppico_raw","Phosphorus Limitation of Picophytoplankton",'h','1','s','1','f')
    cobalt%id_limppico = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limppico", cmor_units="1",                          &
         cmor_standard_name="phosphorus_growth_limitation_of_picophytoplankton", &
         cmor_long_name="Phosphorus Limitation of Picophytoplankton")

    vardesc_temp = vardesc("limpmisc_raw","Phosphorus Limitation of Other Phytoplankton",'h','1','s','1','f')
    cobalt%id_limpmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="limpmisc", cmor_units="1",                          &
         cmor_standard_name="phosphorus_growth_limitation_of_miscellaneous_phytoplankton", &
         cmor_long_name="Phosphorus Limitation of Other Phytoplankton")

!------------------------------------------------------------------------------------------------------------------
! 2-D fields
! sfc tracers

    vardesc_temp = vardesc("dissicos_raw","Surface Dissolved Inorganic Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dissicos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissicos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_in_sea_water",  &
         cmor_long_name="Surface Dissolved Inorganic Carbon Concentration")

    vardesc_temp = vardesc("dissicnatos_raw","Surface Natural Dissolved Inorganic Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dissicnatos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissicnatos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_natural_analogue_in_sea_water",  &
         cmor_long_name="Surface Natural Dissolved Inorganic Carbon Concentration")

    vardesc_temp = vardesc("dissicabioos_raw","Surface Abiotic Dissolved Inorganic Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dissicabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissicabioos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon_abiotic_analogue_in_sea_water",  &
         cmor_long_name="Surface Abiotic Dissolved Inorganic Carbon Concentration")

    vardesc_temp = vardesc("dissi14cabioos_raw","Surface Abiotic Dissolved Inorganic 14Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dissi14cabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissi14cabioos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_carbon14_abiotic_analogue_in_sea_water", &
         cmor_long_name="Surface Abiotic Dissolved Inorganic 14Carbon Concentration")

    vardesc_temp = vardesc("dissocos_raw","Surface Dissolved Organic Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dissocos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dissocos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_organic_carbon_in_sea_water",  &
         cmor_long_name="Surface Dissolved Organic Carbon Concentration")

! also in Oday (updated to match Omon long_name, standard_name)
    vardesc_temp = vardesc("phycos_raw","Surface Phytoplankton Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_phycos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phycos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Phytoplankton Carbon Concentration")

    vardesc_temp = vardesc("zoocos_raw","Surface Zooplankton Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_zoocos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zoocos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_zooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Zooplankton Carbon Concentration")

    vardesc_temp = vardesc("baccos_raw","Surface Bacterial Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_baccos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="baccos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_bacteria_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Bacterial Carbon Concentration")

    vardesc_temp = vardesc("detocos_raw","Surface Detrital Organic Carbon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_detocos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="detocos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_organic_detritus_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Detrital Organic Carbon Concentration")

    vardesc_temp = vardesc("calcos_raw","Surface Calcite Concentration",'h','1','s','mol m-3','f')
    cobalt%id_calcos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="calcos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_calcite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Calcite Concentration")

    vardesc_temp = vardesc("aragos_raw","Surface Aragonite Concentration",'h','1','s','mol m-3','f')
    cobalt%id_aragos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="aragos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_aragonite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Aragonite Concentration")

    vardesc_temp = vardesc("phydiatos_raw","Surface Mole Concentration of Diatoms expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phydiatos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phydiatos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_diatoms_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Diatoms expressed as Carbon in sea water")

    vardesc_temp = vardesc("phydiazos_raw","Surface Mole Concentration of Diazotrophs expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phydiazos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phydiazos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_diazotrophs_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Diazotrophs expressed as Carbon in sea water")

    vardesc_temp = vardesc("phypicoos_raw","Surface Mole Concentration of Picophytoplankton expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phypicoos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phypicoos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_picophytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Picophytoplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("phymiscos_raw","Surface Mole Concentration of Miscellaneous Phytoplankton expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phymiscos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phymiscos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_miscellaneous_phytoplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Miscellaneous Phytoplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("zmicroos_raw","Surface Mole Concentration of Microzooplankton expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_zmicroos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zmicroos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_microzooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Microzooplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("zmesoos_raw","Surface Mole Concentration of Mesozooplankton expressed as Carbon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_zmesoos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zmesoos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_mesozooplankton_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Mesozooplankton expressed as Carbon in sea water")

    vardesc_temp = vardesc("talkos_raw","Surface Total Alkalinity",'h','1','s','mol m-3','f')
    cobalt%id_talkos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="talkos", cmor_units="mol m-3",                          &
         cmor_standard_name="sea_water_alkalinity_expressed_as_mole_equivalent", &
         cmor_long_name="Surface Total Alkalinity")

    vardesc_temp = vardesc("talknatos_raw","Surface Natural Total Alkalinity",'h','1','s','mol m-3','f')
    cobalt%id_talknatos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="talknatos", cmor_units="mol m-3",                          &
         cmor_standard_name="sea_water_alkalinity_natural_analogue_expressed_as_mole_equivalent", &
         cmor_long_name="Surface Natural Total Alkalinity")

    vardesc_temp = vardesc("phos_raw","Surface pH",'h','1','s','1','f')
    cobalt%id_phos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phos", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_reported_on_total_scale", &
         cmor_long_name="Surface pH")

    vardesc_temp = vardesc("phnatos_raw","Surface Natural pH",'h','1','s','1','f')
    cobalt%id_phnatos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phnatos", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_natural_analogue_reported_on_total_scale", &
         cmor_long_name="Surface Natural pH")

    vardesc_temp = vardesc("phabioos_raw","Surface Abiotic pH",'h','1','s','1','f')
    cobalt%id_phabioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phabioos", cmor_units="1",                          &
         cmor_standard_name="sea_water_ph_abiotic_analogue_reported_on_total_scale", &
         cmor_long_name="Surface Abiotic pH")

!! jgj 2017/08/04 removed _cmip in cmor_field_name - update diag table
    vardesc_temp = vardesc("o2os_raw","Surface Dissolved Oxygen Concentration",'h','1','s','mol m-3','f')
    cobalt%id_o2os = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="o2os", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_molecular_oxygen_in_sea_water", &
         cmor_long_name="Surface Dissolved Oxygen Concentration")

! CHECK3: need 3-D field
    vardesc_temp = vardesc("o2satos_raw","Surface Dissolved Oxygen Concentration at Saturation",'h','1','s','mol m-3','f')
    cobalt%id_o2satos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="o2satos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_at_saturation", &
         cmor_long_name="Surface Dissolved Oxygen Concentration at Saturation")

!! jgj 2017/08/04 removed _cmip in cmor_field_name - update diag table
    vardesc_temp = vardesc("no3os_raw","Surface Dissolved Nitrate Concentration",'h','1','s','mol m-3','f')
    cobalt%id_no3os = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="no3os", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_nitrate_in_sea_water", &
         cmor_long_name="Surface Dissolved Nitrate Concentration")

    vardesc_temp = vardesc("nh4os_raw","Surface Dissolved Ammonium Concentration",'h','1','s','mol m-3','f')
    cobalt%id_nh4os = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="nh4os", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_ammonium_in_sea_water", &
         cmor_long_name="Surface Dissolved Ammonium Concentration")

! CHECK3:
! 2018/01/12 standard name is  "surface_mole_concentration_of_dissolved_inorganic_phosphorous_in_sea_water" in data request
! 2018/01/17 removed surface_ from CF standard name in data request
    vardesc_temp = vardesc("po4os_raw","Surface Total Dissolved Inorganic Phosphorus Concentration",'h','1','s','mol m-3','f')
    cobalt%id_po4os = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="po4os", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_phosphorus_in_sea_water", &
         cmor_long_name="Surface Total Dissolved Inorganic Phosphorus Concentration")

    vardesc_temp = vardesc("dfeos_raw","Surface Dissolved Iron Concentration",'h','1','s','mol m-3','f')
    cobalt%id_dfeos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dfeos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_iron_in_sea_water", &
         cmor_long_name="Surface Dissolved Iron Concentration")

    vardesc_temp = vardesc("sios_raw","Surface Total Dissolved Inorganic Silicon Concentration",'h','1','s','mol m-3','f')
    cobalt%id_sios = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="sios", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_inorganic_silicon_in_sea_water", &
         cmor_long_name="Surface Total Dissolved Inorganic Silicon Concentration")

!! jgj 2017/08/04 removed _cmip in cmor_field_name - update diag table
! also in Oday (updated to match Omon long_name, standard_name)
    vardesc_temp = vardesc("chlos_raw","Surface Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water",'h','1','s','kg m-3','f')
    cobalt%id_chlos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chlos", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_phytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Surface Mass Concentration of Total Phytoplankton expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chldiatos_raw","Surface Mass Concentration of Diatoms expressed as Chlorophyll in sea water",'h','1','s','kg m-3','f')
    cobalt%id_chldiatos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chldiatos", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_diatoms_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Surface Mass Concentration of Diatoms expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chldiazos_raw","Surface Mass Concentration of Diazotrophs expressed as Chlorophyll in sea water",'h','1','s','kg m-3','f')
    cobalt%id_chldiazos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chldiazos", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_diazotrophs_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Surface Mass Concentration of Diazotrophs expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chlpicoos_raw","Surface Mass Concentration of Picophytoplankton expressed as Chlorophyll in sea water",'h','1','s','kg m-3','f')
    cobalt%id_chlpicoos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chlpicoos", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_picophytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Surface Mass Concentration of Picophytoplankton expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("chlmiscos_raw","Surface Mass Concentration of Other Phytoplankton expressed as Chlorophyll in sea water",'h','1','s','kg m-3','f')
    cobalt%id_chlmiscos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="chlmiscos", cmor_units="kg m-3",                          &
         cmor_standard_name="mass_concentration_of_miscellaneous_phytoplankton_expressed_as_chlorophyll_in_sea_water", &
         cmor_long_name="Surface Mass Concentration of Other Phytoplankton expressed as Chlorophyll in sea water")

    vardesc_temp = vardesc("ponos_raw","Surface Mole Concentration of Particulate Organic Matter expressed as Nitrogen in sea water",'h','1','s','mol N m-3','f')
    cobalt%id_ponos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ponos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_nitrogen_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Particulate Organic Matter expressed as Nitrogen in sea water")

    vardesc_temp = vardesc("popos_raw","Surface Mole Concentration of Particulate Organic Matter expressed as Phosphorus in sea water",'h','1','s','mol m-3','f')
    cobalt%id_popos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="popos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_phosphorus_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Particulate Organic Matter expressed as Phosphorus in sea water")

    vardesc_temp = vardesc("bfeos_raw","Surface Mole Concentration of Particulate Organic Matter expressed as Iron in sea water",'h','1','s','mol m-3','f')
    cobalt%id_bfeos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bfeos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_iron_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Particulate Organic Matter expressed as Iron in sea water")

    vardesc_temp = vardesc("bsios_raw","Surface Mole Concentration of Particulate Organic Matter expressed as Silicon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_bsios = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="bsios", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_particulate_organic_matter_expressed_as_silicon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Particulate Organic Matter expressed as Silicon in sea water")

    vardesc_temp = vardesc("phynos_raw","Surface Mole Concentration of Phytoplankton Nitrogen in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phynos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phynos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_nitrogen_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Phytoplankton Nitrogen in sea water")

    vardesc_temp = vardesc("phypos_raw","Surface Mole Concentration of Total Phytoplankton expressed as Phosphorus in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phypos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phypos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_phosphorus_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Total Phytoplankton expressed as Phosphorus in sea water")

! 2017/12/28 jgj Updated long name in data request: Surface Mole Concentration of Total Phytoplankton expressed as Iron in sea water
    vardesc_temp = vardesc("phyfeos_raw","Surface Mole Concentration of Total Phytoplankton expressed as Iron in sea water",'h','1','s','mol m-3','f')
    cobalt%id_phyfeos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="phyfeos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_iron_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Total Phytoplankton expressed as Iron in sea water")

    vardesc_temp = vardesc("physios_raw","Surface Mole Concentration of Total Phytoplankton expressed as Silicon in sea water",'h','1','s','mol m-3','f')
    cobalt%id_physios = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="physios", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_phytoplankton_expressed_as_silicon_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Total Phytoplankton expressed as Silicon in sea water")

    vardesc_temp = vardesc("co3os_raw","Surface Carbonate Ion Concentration",'h','1','s','mol m-3','f')
    cobalt%id_co3os = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3os", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Carbonate Ion Concentration")

    vardesc_temp = vardesc("co3natos_raw","Surface Natural Carbonate Ion Concentration",'h','1','s','mol m-3','f')
    cobalt%id_co3natos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3natos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_ion_natural_analogue_in_sea_water", &
         cmor_long_name="Surface Natural Carbonate Ion Concentration")

    vardesc_temp = vardesc("co3abioos_raw","Surface Abiotic Carbonate Ion Concentration",'h','1','s','mol m-3','f')
    cobalt%id_co3abioos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3abioos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_abiotic_analogue_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Surface Abiotic Carbonate Ion Concentration")

    vardesc_temp = vardesc("co3satcalcos_raw","Surface Mole Concentration of Carbonate Ion in Equilibrium with Pure Calcite in sea water",'h','1','s','mol m-3','f')
    cobalt%id_co3satcalcos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3satcalcos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_at_equilibrium_with_pure_calcite_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Carbonate Ion in Equilibrium with Pure Calcite in sea water")

    vardesc_temp = vardesc("co3sataragos_raw","Surface Mole Concentration of Carbonate Ion in Equilibrium with Pure Aragonite in sea water",'h','1','s','mol m-3','f')
    cobalt%id_co3sataragos = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="co3sataragos", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_carbonate_expressed_as_carbon_at_equilibrium_with_pure_aragonite_in_sea_water", &
         cmor_long_name="Surface Mole Concentration of Carbonate Ion in Equilibrium with Pure Aragonite in sea water")

!------------------------------------------------------------------------------------------------------------------
! 2-D fields (from Omon)

    vardesc_temp = vardesc("intpp_raw","Primary Organic Carbon Production by All Types of Phytoplankton",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpp = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpp", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_phytoplankton", &
         cmor_long_name="Primary Organic Carbon Production by All Types of Phytoplankton")

    vardesc_temp = vardesc("intppnitrate_raw","Primary Organic Carbon Production by Phytoplankton Based on Nitrate Uptake Alone",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intppnitrate = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intppnitrate", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_due_to_nitrate_utilization", &
         cmor_long_name="Primary Organic Carbon Production by Phytoplankton Based on Nitrate Uptake Alone")

    vardesc_temp = vardesc("intppdiat_raw","Net Primary Organic Carbon Production by Diatoms",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intppdiat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intppdiat", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_diatoms", &
         cmor_long_name="Net Primary Organic Carbon Production by Diatoms")

    vardesc_temp = vardesc("intppdiaz_raw","Net Primary Mole Productivity of Carbon by Diazotrophs",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intppdiaz = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intppdiaz", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_diazotrophs", &
         cmor_long_name="Net Primary Mole Productivity of Carbon by Diazotrophs")

    vardesc_temp = vardesc("intpppico_raw","Net Primary Mole Productivity of Carbon by Picophytoplankton",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpppico = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpppico", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_picophytoplankton", &
         cmor_long_name="Net Primary Mole Productivity of Carbon by Picophytoplankton")

    vardesc_temp = vardesc("intppmisc_raw","Net Primary Organic Carbon Production by Other Phytoplankton",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intppmisc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intppmisc", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="net_primary_mole_productivity_of_biomass_expressed_as_carbon_by_miscellaneous_phytoplankton", &
         cmor_long_name="Net Primary Organic Carbon Production by Other Phytoplankton")

    vardesc_temp = vardesc("intpbn_raw","Nitrogen Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpbn = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpbn", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_nitrogen_due_to_biological_production", &
         cmor_long_name="Nitrogen Production")

    vardesc_temp = vardesc("intpbp_raw","Phosphorus Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpbp = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpbp", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_phosphorus_due_to_biological_production", &
         cmor_long_name="Phosphorus Production")

    vardesc_temp = vardesc("intpbfe_raw","Iron Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpbfe = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpbfe", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_iron_due_to_biological_production", &
         cmor_long_name="Iron Production")

    vardesc_temp = vardesc("intpbsi_raw","Silicon Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpbsi = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpbsi", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_silicon_due_to_biological_production", &
         cmor_long_name="Silicon Production")

    vardesc_temp = vardesc("intpcalcite_raw","Calcite Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpcalcite = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpcalcite", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_calcite_expressed_as_carbon_due_to_biological_production", &
         cmor_long_name="Calcite Production")

    vardesc_temp = vardesc("intparag_raw","Aragonite Production",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intparag = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intparag", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_aragonite_expressed_as_carbon_due_to_biological_production", &
         cmor_long_name="Aragonite Production")

! CHECK3: these should be AT 100m  - check how to output to conform to CMOR/data request
    vardesc_temp = vardesc("epc100_raw","Downward Flux of Particulate Organic Carbon",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epc100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epc100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_matter_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Downward Flux of Particulate Organic Carbon")

    vardesc_temp = vardesc("epn100_raw","Downward Flux of Particulate Nitrogen",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epn100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epn100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_nitrogen_in_sea_water", &
         cmor_long_name="Downward Flux of Particulate Nitrogen")

    vardesc_temp = vardesc("epp100_raw","Downward Flux of Particulate Phosphorus",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epp100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epp100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_organic_phosphorus_in_sea_water", &
         cmor_long_name="Downward Flux of Particulate Phosphorus")

    vardesc_temp = vardesc("epfe100_raw","Downward Flux of Particulate Iron",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epfe100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epfe100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_iron_in_sea_water", &
         cmor_long_name="Downward Flux of Particulate Iron")

    vardesc_temp = vardesc("epsi100_raw","Downward Flux of Particulate Silicon",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epsi100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epsi100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_particulate_silicon_in_sea_water", &
         cmor_long_name="Downward Flux of Particulate Silicon")

    vardesc_temp = vardesc("epcalc100_raw","Downward Flux of Calcite",'h','1','s','mol m-2 s-1','f')
    cobalt%id_epcalc100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="epcalc100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_calcite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Downward Flux of Calcite")

    vardesc_temp = vardesc("eparag100_raw","Downward Flux of Aragonite",'h','1','s','mol m-2 s-1','f')
    cobalt%id_eparag100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="eparag100", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="sinking_mole_flux_of_aragonite_expressed_as_carbon_in_sea_water", &
         cmor_long_name="Downward Flux of Aragonite")

! vertically integrated
    vardesc_temp = vardesc("intdic_raw","Dissolved Inorganic Carbon Content",'h','1','s','kg m-2','f')
    cobalt%id_intdic = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intdic", cmor_units="kg m-2",                          &
         cmor_standard_name="ocean_mass_content_of_dissolved_inorganic_carbon", &
         cmor_long_name="Dissolved Inorganic Carbon Content")

    vardesc_temp = vardesc("intdoc_raw","Dissolved Organic Carbon Content",'h','1','s','kg m-2','f')
    cobalt%id_intdoc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intdoc", cmor_units="kg m-2",                          &
         cmor_standard_name="ocean_mass_content_of_dissolved_organic_carbon", &
         cmor_long_name="Dissolved Organic Carbon Content")

    vardesc_temp = vardesc("intpoc_raw","Particulate Organic Carbon Content",'h','1','s','kg m-2','f')
    cobalt%id_intpoc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpoc", cmor_units="kg m-2",                          &
         cmor_standard_name="ocean_mass_content_of_particulate_organic_matter_expressed_as_carbon", &
         cmor_long_name="Particulate Organic Carbon Content")

    vardesc_temp = vardesc("spco2_raw","Surface Aqueous Partial Pressure of CO2",'h','1','s','Pa','f')
    cobalt%id_spco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="spco2", cmor_units="Pa",                          &
         cmor_standard_name="surface_partial_pressure_of_carbon_dioxide_in_sea_water", &
         cmor_long_name="Surface Aqueous Partial Pressure of CO2")

    vardesc_temp = vardesc("spco2nat_raw","Natural Surface Aqueous Partial Pressure of CO2",'h','1','s','Pa','f')
    cobalt%id_spco2nat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="spco2nat", cmor_units="Pa",                          &
         cmor_standard_name="surface_partial_pressure_of_carbon_dioxide_natural_analogue_in_sea_water", &
         cmor_long_name="Natural Surface Aqueous Partial Pressure of CO2")

    vardesc_temp = vardesc("spco2abio_raw","Abiotic Surface Aqueous Partial Pressure of CO2",'h','1','s','Pa','f')
    cobalt%id_spco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="spco2abio", cmor_units="Pa",                          &
         cmor_standard_name="surface_partial_pressure_of_carbon_dioxide_abiotic_analogue_in_sea_water", &
         cmor_long_name="Abiotic Surface Aqueous Partial Pressure of CO2")

    vardesc_temp = vardesc("dpco2_raw","Delta PCO2",'h','1','s','Pa','f')
    cobalt%id_dpco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dpco2", cmor_units="Pa",                          &
         cmor_standard_name="surface_carbon_dioxide_partial_pressure_difference_between_sea_water_and_air", &
         cmor_long_name="Delta PCO2")

    vardesc_temp = vardesc("dpco2nat_raw","Natural Delta PCO2",'h','1','s','Pa','f')
    cobalt%id_dpco2nat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dpco2nat", cmor_units="Pa",                          &
         cmor_standard_name="surface_carbon_dioxide_natural_analogue_partial_pressure_difference_between_sea_water_and_air", &
         cmor_long_name="Natural Delta PCO2")

    vardesc_temp = vardesc("dpco2abio_raw","Abiotic Delta PCO2",'h','1','s','Pa','f')
    cobalt%id_dpco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dpco2abio", cmor_units="Pa",                          &
         cmor_standard_name="surface_carbon_dioxide_abiotic_analogue_partial_pressure_difference_between_sea_water_and_air", &
         cmor_long_name="Abiotic Delta PCO2")

    vardesc_temp = vardesc("dpo2_raw","Delta PO2",'h','1','s','Pa','f')
    cobalt%id_dpo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="dpo2", cmor_units="Pa",                          &
         cmor_standard_name="surface_molecular_oxygen_partial_pressure_difference_between_sea_water_and_air", &
         cmor_long_name="Delta PO2")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello
    vardesc_temp = vardesc("fgco2_raw","Surface Downward Flux of Total CO2",'h','1','s','kg m-2 s-1','f')
    cobalt%id_fgco2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fgco2", cmor_units="kg m-2 s-1",                          &
         cmor_standard_name="surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon", &
         cmor_long_name="Surface Downward Flux of Total CO2")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello
    vardesc_temp = vardesc("fgco2nat_raw","Surface Downward Flux of Natural CO2",'h','1','s','kg m-2 s-1','f')
    cobalt%id_fgco2nat = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fgco2nat", cmor_units="kg m-2 s-1",                          &
         cmor_standard_name="surface_downward_mass_flux_of_carbon_dioxide_natural_analogue_expressed_as_carbon", &
         cmor_long_name="Surface Downward Flux of Natural CO2")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello
    vardesc_temp = vardesc("fgco2abio_raw","Surface Downward Flux of Abiotic CO2",'h','1','s','kg m-2 s-1','f')
    cobalt%id_fgco2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fgco2abio", cmor_units="kg m-2 s-1",                          &
         cmor_standard_name="surface_downward_mass_flux_of_carbon_dioxide_abiotic_analogue_expressed_as_carbon", &
         cmor_long_name="Surface Downward Flux of Abiotic CO2")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello
    vardesc_temp = vardesc("fg14co2abio_raw","Surface Downward Flux of Abiotic 14CO2",'h','1','s','kg m-2 s-1','f')
    cobalt%id_fg14co2abio = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fg14co2abio", cmor_units="kg m-2 s-1",                          &
         cmor_standard_name="surface_downward_mass_flux_of_carbon14_dioxide_abiotic_analogue_expressed_as_carbon", &
         cmor_long_name="Surface Downward Flux of Abiotic 14CO2")

! CHECK3:
! 2017/08/04 jgj: CMOR requires positive down, area:areacello
    vardesc_temp = vardesc("fgo2_raw","Surface Downward Flux of O2",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fgo2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fgo2", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="surface_downward_mole_flux_of_molecular_oxygen", &
         cmor_long_name="Surface Downward Flux of O2")

    vardesc_temp = vardesc("icfriver_raw","Flux of Inorganic Carbon Into Ocean Surface by Runoff",'h','1','s','mol m-2 s-1','f')
    cobalt%id_icfriver = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="icfriver", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_inorganic_carbon_due_to_runoff_and_sediment_dissolution", &
         cmor_long_name="Flux of Inorganic Carbon Into Ocean Surface by Runoff")

    vardesc_temp = vardesc("fric_raw","Downward Inorganic Carbon Flux at Ocean Bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fric = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fric", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_inorganic_carbon_due_to_sedimentation", &
         cmor_long_name="Downward Inorganic Carbon Flux at Ocean Bottom")

    vardesc_temp = vardesc("ocfriver_raw","Flux of Organic Carbon Into Ocean Surface by Runoff",'h','1','s','mol m-2 s-1','f')
    cobalt%id_ocfriver = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="ocfriver", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_organic_carbon_due_to_runoff_and_sediment_dissolution", &
         cmor_long_name="Flux of Organic Carbon Into Ocean Surface by Runoff")

    vardesc_temp = vardesc("froc_raw","Downward Organic Carbon Flux at Ocean Bottom",'h','1','s','mol m-2 s-1','f')
    cobalt%id_froc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="froc", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_organic_carbon_due_to_sedimentation", &
         cmor_long_name="Downward Organic Carbon Flux at Ocean Bottom")

    vardesc_temp = vardesc("intpn2_raw","Nitrogen Fixation Rate in Ocean",'h','1','s','mol m-2 s-1','f')
    cobalt%id_intpn2 = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="intpn2", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_elemental_nitrogen_due_to_fixation", &
         cmor_long_name="Nitrogen Fixation Rate in Ocean")

    vardesc_temp = vardesc("fsn_raw","Surface Downward Net Flux of Nitrogen",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsn = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fsn", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_elemental_nitrogen_due_to_deposition_and_fixation_and_runoff", &
         cmor_long_name="Surface Downward Net Flux of Nitrogen")

    vardesc_temp = vardesc("frn_raw","Nitrogen Loss to Sediments and through Denitrification",'h','1','s','mol m-2 s-1','f')
    cobalt%id_frn = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="frn", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_elemental_nitrogen_due_to_denitrification_and_sedimentation", &
         cmor_long_name="Nitrogen Loss to Sediments and through Denitrification")

    vardesc_temp = vardesc("fsfe_raw","Surface Downward Net Flux of Iron",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fsfe = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fsfe", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_iron_due_to_deposition_and_runoff_and_sediment_dissolution", &
         cmor_long_name="Surface Downward Net Flux of Iron")

    vardesc_temp = vardesc("frfe_raw","Iron Loss to Sediments",'h','1','s','mol m-2 s-1','f')
    cobalt%id_frfe = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="frfe", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_iron_due_to_sedimentation", &
         cmor_long_name="Iron Loss to Sediments")

    vardesc_temp = vardesc("o2min_raw","Oxygen Minimum Concentration",'h','1','s','mol m-3','f')
    cobalt%id_o2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="o2min", cmor_units="mol m-3",                          &
         cmor_standard_name="mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_at_shallowest_local_minimum_in_vertical_profile", &
         cmor_long_name="Oxygen Minimum Concentration")

    vardesc_temp = vardesc("zo2min_raw","Depth of Oxygen Minimum Concentration",'h','1','s','m','f')
    cobalt%id_zo2min = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zo2min", cmor_units="m",                          &
         cmor_standard_name="depth_at_shallowest_local_minimum_in_vertical_profile_of_mole_concentration_of_dissolved_molecular_oxygen_in_sea_water", &
         cmor_long_name="Depth of Oxygen Minimum Concentration")

    vardesc_temp = vardesc("zsatcalc_raw","Calcite Saturation Depth",'h','1','s','m','f')
    cobalt%id_zsatcalc = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zsatcalc", cmor_units="m",                          &
         cmor_standard_name="minimum_depth_of_calcite_undersaturation_in_sea_water", &
         cmor_long_name="Calcite Saturation Depth")

    vardesc_temp = vardesc("zsatarag_raw","Aragonite Saturation Depth",'h','1','s','m','f')
    cobalt%id_zsatarag = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="zsatarag", cmor_units="m",                          &
         cmor_standard_name="minimum_depth_of_aragonite_undersaturation_in_sea_water", &
         cmor_long_name="Aragonite Saturation Depth")

    vardesc_temp = vardesc("fddtdic_raw","Rate of Change of Net Dissolved Inorganic Carbon",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtdic = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtdic", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_carbon", &
         cmor_long_name="Rate of Change of Net Dissolved Inorganic Carbon")

    vardesc_temp = vardesc("fddtdin_raw","Rate of Change of Net Dissolved Inorganic Nitrogen",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtdin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtdin", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_nitrogen", &
         cmor_long_name="Rate of Change of Net Dissolved Inorganic Nitrogen")

    vardesc_temp = vardesc("fddtdip_raw","Rate of Change of Net Dissolved Inorganic Phosphorus",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtdip = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtdip", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_phosphorus", &
         cmor_long_name="Rate of Change of Net Dissolved Inorganic Phosphorus")

    vardesc_temp = vardesc("fddtdife_raw","Rate of Change of Net Dissolved Inorganic Iron",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtdife = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtdife", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_iron", &
         cmor_long_name="Rate of Change of Net Dissolved Inorganic Iron")

    vardesc_temp = vardesc("fddtdisi_raw","Rate of Change of Net Dissolved Inorganic Silicon",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtdisi = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtdisi", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_silicon", &
         cmor_long_name="Rate of Change of Net Dissolved Inorganic Silicon")

    vardesc_temp = vardesc("fddtalk_raw","Rate of Change of Total Alkalinity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fddtalk = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fddtalk", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="integral_wrt_depth_of_tendency_of_sea_water_alkalinity_expressed_as_mole_equivalent", &
         cmor_long_name="Rate of Change of Total Alkalinity")

    vardesc_temp = vardesc("fbddtdic_raw","Rate of Change of Dissolved Inorganic Carbon due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtdic = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtdic", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_carbon_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Carbon due to Biological Activity")

    vardesc_temp = vardesc("fbddtdin_raw","Rate of Change of Dissolved Inorganic Nitrogen due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtdin = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtdin", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_nitrogen_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Nitrogen due to Biological Activity")

    vardesc_temp = vardesc("fbddtdip_raw","Rate of Change of Dissolved Inorganic Phosphorus due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtdip = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtdip", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_phosphorus_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Phosphorus due to Biological Activity")

    vardesc_temp = vardesc("fbddtdife_raw","Rate of Change of Dissolved Inorganic Iron due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtdife = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtdife", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_iron_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Iron due to Biological Activity")

    vardesc_temp = vardesc("fbddtdisi_raw","Rate of Change of Dissolved Inorganic Silicon due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtdisi = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtdisi", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="tendency_of_ocean_mole_content_of_dissolved_inorganic_silicon_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Dissolved Inorganic Silicon due to Biological Activity")

    vardesc_temp = vardesc("fbddtalk_raw","Rate of Change of Biological Alkalinity due to Biological Activity",'h','1','s','mol m-2 s-1','f')
    cobalt%id_fbddtalk = register_diag_field(package_name, vardesc_temp%name, axes(1:2), &
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1, &
         cmor_field_name="fbddtalk", cmor_units="mol m-2 s-1",                          &
         cmor_standard_name="integral_wrt_depth_of_tendency_of_sea_water_alkalinity_expressed_as_mole_equivalent_due_to_biological_processes", &
         cmor_long_name="Rate of Change of Biological Alkalinity due to Biological Activity")

!------------------------------------------------------------------------------------------------------------------
! 2-D fields (from Oday)  

! previously defined above


!==============================================================================================================

  end subroutine generic_COBALT_register_diag

  !
  !   This is an internal sub, not a public interface.
  !   Add all the parameters to be used in this module. 
  !
  subroutine user_add_params

    !Specify all parameters used in this modules.
    !===============================@===============================
    !User adds one call for each parameter below!
    !User also adds the definition of each parameter in generic_COBALT_params type
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
    call g_tracer_add_param('init', cobalt%init, .false. )

    call g_tracer_add_param('htotal_scale_lo', cobalt%htotal_scale_lo, 0.01)
    call g_tracer_add_param('htotal_scale_hi', cobalt%htotal_scale_hi, 100.0)

    !  Rho_0 is used in the Boussinesq
    !  approximation to calculations of pressure and
    !  pressure gradients, in units of kg m-3.
    call g_tracer_add_param('RHO_0', cobalt%Rho_0, 1035.0)
    call g_tracer_add_param('NKML' , cobalt%nkml, 1)
    !-----------------------------------------------------------------------
    !       coefficients for O2 saturation
    !-----------------------------------------------------------------------
    call g_tracer_add_param('a_0', cobalt%a_0, 2.00907)
    call g_tracer_add_param('a_1', cobalt%a_1, 3.22014)
    call g_tracer_add_param('a_2', cobalt%a_2, 4.05010)
    call g_tracer_add_param('a_3', cobalt%a_3, 4.94457)
    call g_tracer_add_param('a_4', cobalt%a_4, -2.56847e-01)
    call g_tracer_add_param('a_5', cobalt%a_5, 3.88767)
    call g_tracer_add_param('b_0', cobalt%b_0, -6.24523e-03)
    call g_tracer_add_param('b_1', cobalt%b_1, -7.37614e-03)
    call g_tracer_add_param('b_2', cobalt%b_2, -1.03410e-02 )
    call g_tracer_add_param('b_3', cobalt%b_3, -8.17083e-03)
    call g_tracer_add_param('c_0', cobalt%c_0, -4.88682e-07)
    !-----------------------------------------------------------------------
    !     Schmidt number coefficients
    !-----------------------------------------------------------------------
    if (trim(as_param_cobalt) == 'W92') then
        !  Compute the Schmidt number of CO2 in seawater using the 
        !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
        !  7373-7382).
        call g_tracer_add_param('a1_co2', cobalt%a1_co2,  2068.9)
        call g_tracer_add_param('a2_co2', cobalt%a2_co2, -118.63)
        call g_tracer_add_param('a3_co2', cobalt%a3_co2,  2.9311)
        call g_tracer_add_param('a4_co2', cobalt%a4_co2, -0.027)
        call g_tracer_add_param('a5_co2', cobalt%a5_co2,  0.0)     ! Not used for W92
        !  Compute the Schmidt number of O2 in seawater using the 
        !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
        !  Cycles, 12, 141-163).
        call g_tracer_add_param('a1_o2', cobalt%a1_o2, 1929.7)
        call g_tracer_add_param('a2_o2', cobalt%a2_o2, -117.46)
        call g_tracer_add_param('a3_o2', cobalt%a3_o2, 3.116)
        call g_tracer_add_param('a4_o2', cobalt%a4_o2, -0.0306)
        call g_tracer_add_param('a5_o2', cobalt%a5_o2, 0.0)       ! Not used for W92
        if (is_root_pe()) call mpp_error(NOTE,'generic_cobalt: Using Schmidt number coefficients for W92')
    else if ((trim(as_param_cobalt) == 'W14') .or. (trim(as_param_cobalt) == 'gfdl_cmip6')) then
        !  Compute the Schmidt number of CO2 in seawater using the 
        !  formulation presented by Wanninkhof 
        !  (2014, Limnol. Oceanogr., 12, 351-362)
        call g_tracer_add_param('a1_co2', cobalt%a1_co2,    2116.8)
        call g_tracer_add_param('a2_co2', cobalt%a2_co2,   -136.25)
        call g_tracer_add_param('a3_co2', cobalt%a3_co2,    4.7353)
        call g_tracer_add_param('a4_co2', cobalt%a4_co2, -0.092307)
        call g_tracer_add_param('a5_co2', cobalt%a5_co2, 0.0007555)
        !  Compute the Schmidt number of O2 in seawater using the 
        !  formulation presented by Wanninkhof 
        !  (2014, Limnol. Oceanogr., 12, 351-362)
        call g_tracer_add_param('a1_o2', cobalt%a1_o2, 1920.4)
        call g_tracer_add_param('a2_o2', cobalt%a2_o2, -135.6)
        call g_tracer_add_param('a3_o2', cobalt%a3_o2, 5.2122)
        call g_tracer_add_param('a4_o2', cobalt%a4_o2, -0.10939)
        call g_tracer_add_param('a5_o2', cobalt%a5_o2, 0.00093777)
        if (is_root_pe()) call mpp_error(NOTE,'generic_cobalt: Using Schmidt number coefficients for W14')
    else
        call mpp_error(FATAL,'generic_cobalt: unable to set Schmidt number coefficients for as_param '//trim(as_param_cobalt))
    endif
    !-----------------------------------------------------------------------
    ! Stoichiometry
    !-----------------------------------------------------------------------
    !
    !
    ! Values taken from OCMIP-II Biotic protocols after Anderson (1995) for an
    ! organic material of C106H172O38N16(H3PO4) which gives an average oxidation state
    ! for carbon of (3*16+2*38-172)/106 = -0.4528.  These calculations ignore H3PO4.
    !
    ! Nitrate Production:
    !   16*H+ + 16*NO3- + 106*CO2 + 78*H2O <-> C106H172O38N16 + 150*O2
    !   Effect is to increase alkalinity by 16 NO3 equivalents.
    !
    ! Ammonia Production (and reverse for remineralization):
    !   16*NH4+ + 106*CO2 + 62*H2O <-> C106H172O38N16 + 118*O2 + 16*H+
    !   Effect is to decrease (increase for remineralization) alkalinity by 16 NH4 equivalents.
    !
    ! N2 Production:
    !   8*N2 + 106*CO2 + 86*H2O <-> C106H172O38N16 + 130*O2
    !   No effect on alkalinity.
    !
    ! Nitrification:
    !   NH4+ + 2*O2 <-> NO3- + H2O + 2*H+
    !   Effect is to decrease alkalinity by 2 NH4 equivalents.
    !
    ! Denitrification:
    !   C106H172O38N16 + 472/5*NO3- + 552/5*H+ <-> 106*CO2 + 16*NH4+ + 236/5*N2 + 546/5*H2O
    !   Effect is to increase alkalinity by 552/472 = 1.169 NO3 equivalents.
    !
    call g_tracer_add_param('n_2_n_denit', cobalt%n_2_n_denit, 472.0/(5.0*16.0))             ! mol N NO3 mol N org-1
    call g_tracer_add_param('o2_2_nfix', cobalt%o2_2_nfix, 130.0/16.0)                       ! mol O2 mol N-1
!    call g_tracer_add_param('o2_2_nfix', cobalt%o2_2_nfix, (118.0+3.0/(5.0+3.0)*(150.0-118.0))/16.0) ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nh4', cobalt%o2_2_nh4, 118.0 / 16.0)                       ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_nitrif', cobalt%o2_2_nitrif, 2.0)                          ! mol O2 mol N-1
    call g_tracer_add_param('o2_2_no3', cobalt%o2_2_no3, 150.0 / 16.0)                       ! mol O2 mol N-1
    !
    !-----------------------------------------------------------------------
    ! Nutrient Limitation Parameters (phytoplankton) 
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('k_fed_Di', phyto(DIAZO)%k_fed, 5.0e-10)                   ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Lg', phyto(LARGE)%k_fed, 5.0e-10)                   ! mol Fed kg-1
    call g_tracer_add_param('k_fed_Sm', phyto(SMALL)%k_fed,  1.0e-10)                 ! mol Fed kg-1
    call g_tracer_add_param('k_nh4_Lg', phyto(LARGE)%k_nh4,  k_nh4_large)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_nh4_Sm', phyto(SMALL)%k_nh4,  k_nh4_small)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_nh4_Di', phyto(DIAZO)%k_nh4,  k_nh4_diazo)                  ! mol NH4 kg-1
    call g_tracer_add_param('k_no3_Lg', phyto(LARGE)%k_no3,  k_no3_large)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_no3_Sm', phyto(SMALL)%k_no3,  k_no3_small)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_no3_Di', phyto(DIAZO)%k_no3,  k_no3_diazo)                  ! mol NO3 kg-1
    call g_tracer_add_param('k_po4_Di', phyto(DIAZO)%k_po4,  5.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Lg', phyto(LARGE)%k_po4,  5.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_po4_Sm', phyto(SMALL)%k_po4,  1.0e-8)                  ! mol PO4 kg-1
    call g_tracer_add_param('k_sio4_Lg',phyto(LARGE)%k_sio4, 2.0e-6)                        ! mol SiO4 kg-1
    call g_tracer_add_param('k_fe_2_n_Di', phyto(DIAZO)%k_fe_2_n, 12.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Lg', phyto(LARGE)%k_fe_2_n, 6.0e-6 * 106.0 / 16.0)   ! mol Fe mol N-1
    call g_tracer_add_param('k_fe_2_n_Sm',phyto(SMALL)%k_fe_2_n, 3.0e-6*106.0/16.0)        ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Sm',phyto(SMALL)%fe_2_n_max, 50.e-6*106.0/16.0)     ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Lg', phyto(LARGE)%fe_2_n_max, 500.0e-6*106.0/16.0)  ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_max_Di', phyto(DIAZO)%fe_2_n_max, 500.0e-6*106.0/16.0)  ! mol Fe mol N-1
    call g_tracer_add_param('fe_2_n_upt_fac', cobalt%fe_2_n_upt_fac, 15.0e-6)               ! mol Fe mol N-1
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton light limitation/growth rate
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('alpha_Di', phyto(DIAZO)%alpha,  0.8e-5 * 2.77e18 / 6.022e17)  ! g C g Chl-1 m2 J-1
    call g_tracer_add_param('alpha_Lg', phyto(LARGE)%alpha,  0.8e-5 * 2.77e18 / 6.022e17)  ! g C g Chl-1 m2 J-1 
    call g_tracer_add_param('alpha_Sm', phyto(SMALL)%alpha,  2.4e-5*2.77e18/6.022e17)      ! g C g Chl-1 m-2 J-1
    call g_tracer_add_param('kappa_eppley', cobalt%kappa_eppley, 0.063)                    ! deg C-1
    call g_tracer_add_param('P_C_max_Di', phyto(DIAZO)%P_C_max, 0.50/sperd)                ! s-1
    ! Uncomment for "no mass change" check
    ! call g_tracer_add_param('P_C_max_Di', phyto(DIAZO)%P_C_max, 0.01/sperd)              ! s-1
    call g_tracer_add_param('P_C_max_Lg', phyto(LARGE)%P_C_max, 1.25/sperd)                ! s-1
    call g_tracer_add_param('P_C_max_Sm', phyto(SMALL)%P_C_max, 1.25/sperd)               ! s-1
    call g_tracer_add_param('thetamax_Di', phyto(DIAZO)%thetamax, 0.03)                    ! g Chl g C-1
    call g_tracer_add_param('thetamax_Lg', phyto(LARGE)%thetamax, 0.05)                    ! g Chl g C-1
    call g_tracer_add_param('thetamax_Sm', phyto(SMALL)%thetamax, 0.03)                    ! g Chl g C-1
    call g_tracer_add_param('bresp_Di', phyto(DIAZO)%bresp,0.05/sperd)                     ! sec-1 
    call g_tracer_add_param('bresp_Lg', phyto(LARGE)%bresp,0.05/sperd)                     ! sec-1 
    call g_tracer_add_param('bresp_Sm', phyto(SMALL)%bresp,0.03/sperd)                     ! sec-1 
    call g_tracer_add_param('thetamin', cobalt%thetamin, 0.002)                            ! g Chl g C-1
    call g_tracer_add_param('thetamin_nolim', cobalt%thetamin_nolim, 0.0)                  ! g Chl g C-1
    call g_tracer_add_param('zeta', cobalt%zeta, 0.05)                                     ! dimensionless
    call g_tracer_add_param('gamma_irr_mem', cobalt%gamma_irr_mem, 1.0 / sperd)            ! s-1
    call g_tracer_add_param('gamma_mu_mem', cobalt%gamma_mu_mem, 1.0 / sperd)              ! s-1
    call g_tracer_add_param('refuge_conc', cobalt%refuge_conc, 1.0e-9)                     ! moles N kg-1
    !
    !-----------------------------------------------------------------------
    ! Nitrogen fixation inhibition parameters
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('o2_inhib_Di_pow', cobalt%o2_inhib_Di_pow, 4.0)                 ! mol O2-1 m3
    call g_tracer_add_param('o2_inhib_Di_sat', cobalt%o2_inhib_Di_sat, 3.0e-4)              ! mol O2 kg-1
    !
    !-----------------------------------------------------------------------
    ! Other stoichiometry
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('p_2_n_static', cobalt%p_2_n_static, .true. )
    call g_tracer_add_param('c_2_n', cobalt%c_2_n, 106.0 / 16.0)
    call g_tracer_add_param('alk_2_n_denit', cobalt%alk_2_n_denit, 552.0/472.0)             ! eq. alk mol NO3-1
    call g_tracer_add_param('p_2_n_static_Di', phyto(DIAZO)%p_2_n_static,1.0/40.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_static_Lg', phyto(LARGE)%p_2_n_static,1.0/12.0 )         ! mol P mol N-1
    call g_tracer_add_param('p_2_n_static_Sm', phyto(SMALL)%p_2_n_static,1.0/20.0 )         ! mol P mol N-1
    call g_tracer_add_param('si_2_n_static_Lg', phyto(LARGE)%si_2_n_static, 2.0)            ! mol Si mol N-1
    call g_tracer_add_param('si_2_n_max_Lg', phyto(LARGE)%si_2_n_max, 3.0)                  ! mol Si mol N-1
    call g_tracer_add_param('ca_2_n_arag', cobalt%ca_2_n_arag, 0.030 * 106.0 / 16.0)        ! mol Ca mol N-1
    call g_tracer_add_param('ca_2_n_calc', cobalt%ca_2_n_calc, 0.013 * 106.0 / 16.0)        ! mol Ca mol N-1
    call g_tracer_add_param('caco3_sat_max', cobalt%caco3_sat_max,10.0)                     ! dimensionless
    !
    !-----------------------------------------------------------------------
    ! Zooplankton Stoichiometry - presently static
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('q_p_2_n_smz',zoo(1)%q_p_2_n, 1.0/18.0)          ! mol P mol N-1 
    call g_tracer_add_param('q_p_2_n_mdz',zoo(2)%q_p_2_n, 1.0/16.0)          ! mol P mol N-1 
    call g_tracer_add_param('q_p_2_n_lgz',zoo(3)%q_p_2_n, 1.0/16.0)          ! mol P mol N-1 
    !
    !-----------------------------------------------------------------------
    ! Bacteria Stoichiometry - presently static
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('q_p_2_n_bact',bact(1)%q_p_2_n, 1.0/16.0)        ! mol P mol N-1
    !
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton aggregation
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('agg_Sm',phyto(SMALL)%agg,0.1*1e6 / sperd)           ! s-1 (mole N kg)-1
    call g_tracer_add_param('agg_Di',phyto(DIAZO)%agg,  0.0    / sperd)          ! s-1 (mole N kg)-1
    call g_tracer_add_param('agg_Lg',phyto(LARGE)%agg,0.3*1e6/ sperd)            ! s-1 (mole N kg)-1
    call g_tracer_add_param('frac_mu_agg_Sm',phyto(SMALL)%frac_mu_agg,0.25)      ! none
    call g_tracer_add_param('frac_mu_agg_Di',phyto(DIAZO)%frac_mu_agg,0.25)      ! none
    call g_tracer_add_param('frac_mu_agg_Lg',phyto(LARGE)%frac_mu_agg,0.25)      ! none 
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton and bacterial losses to viruses
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('vir_Sm',phyto(SMALL)%vir, 0.20*1e6/sperd )  ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Di',phyto(DIAZO)%vir, 0.0 )        ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Lg',phyto(LARGE)%vir, 0.0 )        ! s-1 (mole N kg)-1
    call g_tracer_add_param('vir_Bact',bact(1)%vir,   0.20*1e6/sperd)   ! s-1 (mole N kg)-1
    call g_tracer_add_param('ktemp_vir',cobalt%vir_ktemp, 0.063)       ! C-1
    !
    !-----------------------------------------------------------------------
    ! Phytoplankton losses to exudation
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('exu_Sm',phyto(SMALL)%exu, 0.13) ! dimensionless (fraction of NPP)
    call g_tracer_add_param('exu_Di',phyto(DIAZO)%exu, 0.13)       ! dimensionless (fraction of NPP) 
    call g_tracer_add_param('exu_Lg',phyto(LARGE)%exu, 0.13)       ! dimensionless (fraction of NPP) 
    !
    !-----------------------------------------------------------------------
    ! Zooplankton ingestion parameterization and temperature dependence
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('imax_smz',zoo(1)%imax, 0.9*1.42 / sperd)              ! s-1
    call g_tracer_add_param('imax_mdz',zoo(2)%imax, 0.57 / sperd)              ! s-1
    call g_tracer_add_param('imax_lgz',zoo(3)%imax, 0.23 / sperd)              ! s-1
    call g_tracer_add_param('ki_smz',zoo(1)%ki, 1.25e-6)                       ! moles N kg-1
    call g_tracer_add_param('ki_mdz',zoo(2)%ki, 1.25e-6)                       ! moles N kg-1
    call g_tracer_add_param('ki_lgz',zoo(3)%ki, 1.25e-6)                       ! moles N kg-1
    call g_tracer_add_param('ktemp_smz',zoo(1)%ktemp, 0.063)                   ! C-1
    call g_tracer_add_param('ktemp_mdz',zoo(2)%ktemp, 0.063)                   ! C-1
    call g_tracer_add_param('ktemp_lgz',zoo(3)%ktemp, 0.063)                   ! C-1
    !
    !-----------------------------------------------------------------------
    ! Bacterial growth and uptake parameters
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('mu_max_bact',bact(1)%mu_max, 1.0/sperd )          ! s-1 
    call g_tracer_add_param('k_ldon_bact', bact(1)%k_ldon,  5.0e-7)            ! mol ldon kg-1
    call g_tracer_add_param('ktemp_bact', bact(1)%ktemp, 0.063)                ! C-1
    call g_tracer_add_param('gge_max_bact',bact(1)%gge_max,0.4)                ! dimensionless
    call g_tracer_add_param('bresp_bact',bact(1)%bresp, 0.0075/sperd)          ! s-1
    !
    !-----------------------------------------------------------------------
    ! Zooplankton switching and prey preference parameters
    !-----------------------------------------------------------------------
    !
    ! parameters controlling the extent of biomass-based switching between
    ! multiple prey options 
    call g_tracer_add_param('nswitch_smz',zoo(1)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('nswitch_mdz',zoo(2)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('nswitch_lgz',zoo(3)%nswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_smz',zoo(1)%mswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_mdz',zoo(2)%mswitch, 2.0)          ! dimensionless
    call g_tracer_add_param('mswitch_lgz',zoo(3)%mswitch, 2.0)          ! dimensionless
    ! innate prey availability for small zooplankton 
    call g_tracer_add_param('smz_ipa_smp',zoo(1)%ipa_smp, 1.0)    ! dimensionless
    call g_tracer_add_param('smz_ipa_lgp',zoo(1)%ipa_lgp, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_diaz',zoo(1)%ipa_diaz,0.0)         ! dimensionless
    call g_tracer_add_param('smz_ipa_smz',zoo(1)%ipa_smz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_mdz',zoo(1)%ipa_mdz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_lgz',zoo(1)%ipa_lgz, 0.0)          ! dimensionless
    call g_tracer_add_param('smz_ipa_bact',zoo(1)%ipa_bact,0.25)         ! dimensionless
    call g_tracer_add_param('smz_ipa_det',zoo(1)%ipa_det, 0.0)          ! dimensionless
    ! innate prey availability for large zooplankton 
    call g_tracer_add_param('mdz_ipa_smp',zoo(2)%ipa_smp, 0.0)    ! dimensionless
    call g_tracer_add_param('mdz_ipa_lgp',zoo(2)%ipa_lgp, 1.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_diaz',zoo(2)%ipa_diaz,1.0)         ! dimensionless
    call g_tracer_add_param('mdz_ipa_smz',zoo(2)%ipa_smz, 1.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_mdz',zoo(2)%ipa_mdz, 0.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_lgz',zoo(2)%ipa_lgz, 0.0)          ! dimensionless
    call g_tracer_add_param('mdz_ipa_bact',zoo(2)%ipa_bact, 0.0)        ! dimensionless
    call g_tracer_add_param('mdz_ipa_det',zoo(2)%ipa_det, 0.0)          ! dimensionless
    ! innate prey availability large predatory zooplankton/krill
    call g_tracer_add_param('lgz_ipa_smp',zoo(3)%ipa_smp, 0.0)   ! dimensionless
    call g_tracer_add_param('lgz_ipa_lgp',zoo(3)%ipa_lgp, 1.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_diaz',zoo(3)%ipa_diaz, 1.0)       ! dimensionless
    call g_tracer_add_param('lgz_ipa_smz',zoo(3)%ipa_smz, 0.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_mdz',zoo(3)%ipa_mdz, 1.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_lgz',zoo(3)%ipa_lgz, 0.0)         ! dimensionless
    call g_tracer_add_param('lgz_ipa_bact',zoo(3)%ipa_bact, 0.0)       ! dimensionless
    call g_tracer_add_param('lgz_ipa_det',zoo(3)%ipa_det, 0.0)         ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Zooplankton bioenergetics
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('gge_max_smz',zoo(1)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('gge_max_mdz',zoo(2)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('gge_max_lgz',zoo(3)%gge_max, 0.4)                   ! dimensionless
    call g_tracer_add_param('bresp_smz',zoo(1)%bresp, 0.9*0.020 / sperd)        ! s-1
    call g_tracer_add_param('bresp_mdz',zoo(2)%bresp, 0.008 / sperd)        ! s-1
    call g_tracer_add_param('bresp_lgz',zoo(3)%bresp, 0.0032 / sperd)       ! s-1
    !
    !----------------------------------------------------------------------
    ! Bacterial bioenergetics
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('gge_max_bact',bact(1)%gge_max,0.4)              ! dimensionless
    call g_tracer_add_param('bresp_bact',bact(1)%bresp, 0.0075/sperd)         ! s-1
    !
    !----------------------------------------------------------------------
    ! Partitioning of zooplankton ingestion to other compartments
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('phi_det_smz',zoo(1)%phi_det, 0.10)            ! dimensionless
    call g_tracer_add_param('phi_det_mdz',zoo(2)%phi_det, 0.20)            ! dimensionless
    call g_tracer_add_param('phi_det_lgz',zoo(3)%phi_det, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_ldon_smz',zoo(1)%phi_ldon, 0.7*0.20)    ! dimensionless
    call g_tracer_add_param('phi_ldon_mdz',zoo(2)%phi_ldon, 0.7*0.10)    ! dimensionless
    call g_tracer_add_param('phi_ldon_lgz',zoo(3)%phi_ldon, 0.7*0.0)     ! dimensionless
    call g_tracer_add_param('phi_ldop_smz',zoo(1)%phi_ldop, 0.65*0.20)     ! dimensionless
    call g_tracer_add_param('phi_ldop_mdz',zoo(2)%phi_ldop, 0.65*0.10)     ! dimensionless
    call g_tracer_add_param('phi_ldop_lgz',zoo(3)%phi_ldop, 0.65*0.0)      ! dimensionless
    call g_tracer_add_param('phi_srdon_smz',zoo(1)%phi_srdon, 0.1*0.20)    ! dimensionless
    call g_tracer_add_param('phi_srdon_mdz',zoo(2)%phi_srdon, 0.1*0.10)    ! dimensionless
    call g_tracer_add_param('phi_srdon_lgz',zoo(3)%phi_srdon, 0.1*0.0)     ! dimensionless
    call g_tracer_add_param('phi_srdop_smz',zoo(1)%phi_srdop, 0.15*0.20)   ! dimensionless
    call g_tracer_add_param('phi_srdop_mdz',zoo(2)%phi_srdop, 0.15*0.10)   ! dimensionless
    call g_tracer_add_param('phi_srdop_lgz',zoo(3)%phi_srdop, 0.15*0.0)    ! dimensionless
    call g_tracer_add_param('phi_sldon_smz',zoo(1)%phi_sldon, 0.2*0.20)    ! dimensionless
    call g_tracer_add_param('phi_sldon_mdz',zoo(2)%phi_sldon, 0.2*0.10)    ! dimensionless
    call g_tracer_add_param('phi_sldon_lgz',zoo(3)%phi_sldon, 0.2*0.0)     ! dimensionless
    call g_tracer_add_param('phi_sldop_smz',zoo(1)%phi_sldop, 0.2*0.20)    ! dimensionless
    call g_tracer_add_param('phi_sldop_mdz',zoo(2)%phi_sldop, 0.2*0.10)    ! dimensionless
    call g_tracer_add_param('phi_sldop_lgz',zoo(3)%phi_sldop, 0.2*0.0)     ! dimensionless
    call g_tracer_add_param('phi_nh4_smz',zoo(1)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_nh4_mdz',zoo(2)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_nh4_lgz',zoo(3)%phi_nh4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_smz',zoo(1)%phi_po4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_mdz',zoo(2)%phi_po4, 0.30)            ! dimensionless
    call g_tracer_add_param('phi_po4_lgz',zoo(3)%phi_po4, 0.30)            ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Partitioning of viral losses to various dissolved pools
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('phi_ldon_vir',cobalt%lysis_phi_ldon, 0.7)    ! dimensionless
    call g_tracer_add_param('phi_srdon_vir',cobalt%lysis_phi_srdon, 0.1)  ! dimensionless
    call g_tracer_add_param('phi_sldon_vir',cobalt%lysis_phi_sldon, 0.2)  ! dimensionless
    call g_tracer_add_param('phi_ldop_vir',cobalt%lysis_phi_ldop, 0.65)    ! dimensionless
    call g_tracer_add_param('phi_srdop_vir',cobalt%lysis_phi_srdop, 0.15)  ! dimensionless
    call g_tracer_add_param('phi_sldop_vir',cobalt%lysis_phi_sldop, 0.20)  ! dimensionless
    ! 
    !----------------------------------------------------------------------
    ! Parameters for unresolved higher predators
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('imax_hp',     cobalt%imax_hp, 0.09/sperd)     ! s-1 
    call g_tracer_add_param('ki_hp',       cobalt%ki_hp, 1.25e-6)          ! mol N kg-1
    call g_tracer_add_param('coef_hp',     cobalt%coef_hp, 2.0)            ! dimensionless
    call g_tracer_add_param('ktemp_hp',    cobalt%ktemp_hp, 0.063)         ! C-1 
    call g_tracer_add_param('nswitch_hp',  cobalt%nswitch_hp, 2.0)         ! dimensionless
    call g_tracer_add_param('mswitch_hp',  cobalt%mswitch_hp, 2.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_smp',  cobalt%hp_ipa_smp, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_lgp',  cobalt%hp_ipa_lgp, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_diaz', cobalt%hp_ipa_diaz, 0.0)        ! dimensionless
    call g_tracer_add_param('hp_ipa_smz',  cobalt%hp_ipa_smz, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_mdz',  cobalt%hp_ipa_mdz, 1.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_lgz',  cobalt%hp_ipa_lgz, 1.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_bact', cobalt%hp_ipa_bact,0.0)         ! dimensionless
    call g_tracer_add_param('hp_ipa_det',  cobalt%hp_ipa_det, 0.0)         ! dimensionless
    call g_tracer_add_param('hp_phi_det',  cobalt%hp_phi_det, 0.35)        ! dimensionless
    !
    !----------------------------------------------------------------------
    ! Iron chemistry
    !----------------------------------------------------------------------
    !
    call g_tracer_add_param('felig_bkg', cobalt%felig_bkg, 0.5e-9)                          ! mol Fe kg-1
    call g_tracer_add_param('felig_2_don', cobalt%felig_2_don, 0.5e-3)                          ! mol lig mol N-1
    call g_tracer_add_param('fe_2_n_sed', cobalt%fe_2_n_sed, 100.0e-5 * 106 / 16)            ! mol Fe mol N-1
    call g_tracer_add_param('ffe_sed_max', cobalt%ffe_sed_max, 170.0/1.0e6/sperd)            ! mol Fe m-2 s-1 
    call g_tracer_add_param('ffe_geotherm_ratio', cobalt%ffe_geotherm_ratio,2.0e-12)         ! mol Fe m-2 s-1 (watt m-2)-1
    call g_tracer_add_param('ffe_iceberg_ratio', cobalt%ffe_iceberg_ratio,1.0e-7)            ! mol Fe kg-1 ice melt
    call g_tracer_add_param('fe_coast', cobalt%fe_coast,0.0 )                                ! mol Fe m kg-1 s-1
    call g_tracer_add_param('alpha_fescav',cobalt%alpha_fescav, 0.0/spery)                   ! sec-1
    call g_tracer_add_param('beta_fescav',cobalt%beta_fescav, 2.5e9/spery )                  ! sec-1 (mole ndet kg-1)-1
    call g_tracer_add_param('remin_eff_fedet',cobalt%remin_eff_fedet, 0.25)                  ! unitless 
    call g_tracer_add_param('io_fescav',cobalt%io_fescav, 10.0 )                             ! watts m-2
    call g_tracer_add_param('kfe_eq_lig_ll',cobalt%kfe_eq_lig_ll, 1.0e12)                    ! mol lig-1 kg
    call g_tracer_add_param('kfe_eq_lig_hl',cobalt%kfe_eq_lig_hl, 1.0e9)                     ! mol lig-1 kg

    ! Radiocarbon
    call g_tracer_add_param('half_life_14c', cobalt%half_life_14c, 5730.0 )                  ! s
    call g_tracer_add_param('lambda_14c', cobalt%lambda_14c, log(2.0) / (cobalt%half_life_14c * spery)) ! s-1
    !-------------------------------------------------------------------------
    ! Remineralization
    !-------------------------------------------------------------------------
    !
    call g_tracer_add_param('k_o2', cobalt%k_o2, 8.0e-6)                                     ! mol O2 kg-1
    call g_tracer_add_param('o2_min', cobalt%o2_min, 0.8e-6 )                                ! mol O2 kg-1
    call g_tracer_add_param('k_o2_nit', cobalt%k_o2_nit, k_o2_nit)                           ! mol O2 kg-1
    call g_tracer_add_param('o2_min_nit', cobalt%o2_min_nit, o2_min_nit )                    ! mol O2 kg-1
    call g_tracer_add_param('kappa_remin', cobalt%kappa_remin, 0.063 )                       ! deg C-1
    call g_tracer_add_param('remin_ramp_scale', cobalt%remin_ramp_scale, 50.0 )              ! m
    call g_tracer_add_param('rpcaco3', cobalt%rpcaco3, 0.070/12.0*16.0/106.0*100.0)          ! mol N mol Ca-1
    call g_tracer_add_param('rplith',  cobalt%rplith,  0.065/12.0*16.0/106.0)                ! mol N g lith-1
    call g_tracer_add_param('rpsio2',  cobalt%rpsio2,  0.026/12.0*16.0/106.0*60.0)           ! mol N mol Si-1
    call g_tracer_add_param('gamma_ndet',  cobalt%gamma_ndet, cobalt%wsink / 350.0 )         ! s-1
    call g_tracer_add_param('gamma_cadet_arag',cobalt%gamma_cadet_arag,cobalt%wsink/760.0)   ! s-1
    call g_tracer_add_param('gamma_cadet_calc',cobalt%gamma_cadet_calc,cobalt%wsink/1343.0)  ! s-1
    call g_tracer_add_param('kappa_sidet',  cobalt%kappa_sidet, 0.063 )                      ! deg C -1 
    call g_tracer_add_param('gamma_sidet',  cobalt%gamma_sidet, cobalt%wsink / 1.0e4 )       ! s-1
    call g_tracer_add_param('phi_lith' ,  cobalt%phi_lith, 0.002)                            ! dimensionless 
    call g_tracer_add_param('k_lith',  cobalt%k_lith, 0.5/spery )                            ! s-1
    call g_tracer_add_param('z_sed',  cobalt%z_sed, 0.1 )                                    ! m
    call g_tracer_add_param('k_no3_denit',cobalt%k_no3_denit,1.0e-6)                         ! mol NO3 kg-1
    call g_tracer_add_param('z_burial',cobalt%z_burial,50.0)                                 ! m
    !
    !-----------------------------------------------------------------------
    ! Calcium carbonate in sediments (see Dunne et al., 2012, GBC, 26)
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('cased_steady', cobalt%cased_steady, .false. )
    call g_tracer_add_param('phi_surfresp_cased',cobalt%phi_surfresp_cased,0.14307)      ! none 
    call g_tracer_add_param('phi_deepresp_cased',cobalt%phi_deepresp_cased,4.1228)       ! none
    call g_tracer_add_param('alpha_cased',cobalt%alpha_cased,2.7488)                     ! none
    call g_tracer_add_param('beta_cased',cobalt%beta_cased,-2.2185)                      ! none
    call g_tracer_add_param('gamma_cased',cobalt%gamma_cased,0.03607/spery)              ! sec-1
    call g_tracer_add_param('Co_cased',cobalt%Co_cased,8.1e3)                            ! moles m-3 
    !
    !-----------------------------------------------------------------------
    ! Dissolved Organic Material
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('gamma_srdon',  cobalt%gamma_srdon, 1.0 / (10.0 * spery))          ! s-1
    call g_tracer_add_param('gamma_srdop',  cobalt%gamma_srdop, 1.0 / (4.0 * spery))           ! s-1
    call g_tracer_add_param('gamma_sldon',  cobalt%gamma_sldon, 1.0 / (90 * sperd))           ! s-1
    call g_tracer_add_param('gamma_sldop',  cobalt%gamma_sldop, 1.0 / (90 * sperd))           ! s-1
    ! 2016/08/24 jgj add parameter for background dissolved organic material
    ! For the oceanic carbon budget, a constant 42 uM of dissolved organic
    ! carbon is added to represent the refractory component.
    ! For the oceanic nitrogen budget, a constant 2 uM of dissolved organic
    ! nitrogen is added to represent the refractory component.
    ! 2016/09/22 jgj changed background DOC to 4.0e-5 per agreement with CAS, JPD
    !
    call g_tracer_add_param('doc_background',  cobalt%doc_background, 4.0e-5)    ! uM

    !---------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    ! Nitrification
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('gamma_nitrif',  cobalt%gamma_nitrif, gamma_nitrif / (30.0 * sperd))     ! s-1
    call g_tracer_add_param('gamma_nitrif',  cobalt%k_nh3_nitrif, k_nh3_nitrif )                     ! moles kg-1
    call g_tracer_add_param('irr_inhibit',  cobalt%irr_inhibit, irr_inhibit)                         ! W m-2
    !
    !-----------------------------------------------------------------------
    ! Miscellaneous
    !-----------------------------------------------------------------------
    !
    call g_tracer_add_param('tracer_debug',  cobalt%tracer_debug, .false.)

    call g_tracer_end_param_list(package_name)
    !===========
    !Block Ends: g_tracer_add_param
    !===========

  end subroutine user_add_params

  subroutine user_add_tracers(tracer_list)
    type(g_tracer_type), pointer :: tracer_list
    character(len=fm_string_len), parameter :: sub_name = 'user_add_tracers'
    real :: as_coeff_cobalt

    if ((trim(as_param_cobalt) == 'W92') .or. (trim(as_param_cobalt) == 'gfdl_cmip6')) then
      ! Air-sea gas exchange coefficient presented in OCMIP2 protocol.
      ! Value is 0.337 cm/hr in units of m/s.
      as_coeff_cobalt=9.36e-7
    else
      ! Value is 0.251 cm/hr in units of m/s
      as_coeff_cobalt=6.972e-7
    endif

    !
    !Add here only the parameters that are required at the time of registeration 
    !(to make flux exchanging Ocean tracers known for all PE's) 
    !
    call g_tracer_start_param_list(package_name)
    !
    call g_tracer_add_param('htotal_in', cobalt%htotal_in, 1.0e-08)
    !
    ! Sinking velocity of detritus: a value of 20 m d-1 is consistent with a characteristic sinking
    ! velocity of 100 m d-1 of marine aggregates and a disaggregation rate constant
    ! of 5 d-1 in the surface ocean (Clegg and Whitfield, 1992; Dunne, 1999).  Alternatively, 100 m d-1 
    ! is more in line with the deep water synthesis of Berelson (2002; Particel settling rates increase
    ! with depth in the ocean, DSR-II, 49, 237-252).
    !
    call g_tracer_add_param('wsink',  cobalt%wsink, 100.0 / sperd)                             ! m s-1

    call g_tracer_add_param('ice_restart_file'   , cobalt%ice_restart_file   , 'ice_cobalt.res.nc')
    call g_tracer_add_param('ocean_restart_file' , cobalt%ocean_restart_file , 'ocean_cobalt.res.nc' )
    call g_tracer_add_param('IC_file'       , cobalt%IC_file       , '')
    !
    call g_tracer_end_param_list(package_name)

    ! Set Restart files
    call g_tracer_set_files(ice_restart_file    = cobalt%ice_restart_file,&
         ocean_restart_file  = cobalt%ocean_restart_file )

    !All tracer fields shall be registered for diag output.

    !=====================================================
    !Specify all prognostic tracers of this modules.
    !=====================================================
    !User adds one call for each prognostic tracer below!
    !User should specify if fluxes must be extracted from boundary 
    !by passing one or more of the following methods as .true.  
    !and provide the corresponding parameters array
    !methods: flux_gas,flux_runoff,flux_wetdep,flux_drydep  
    !
    !Pass an init_value arg if the tracers should be initialized to a nonzero value everywhere
    !otherwise they will be initialized to zero.
    !
    !===========================================================
    !Prognostic Tracers
    !===========================================================
    !
    !
    !       ALK (Total carbonate alkalinity)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'alk',         &
         longname   = 'Alkalinity',  &
         units      = 'mol/kg',      &
         prog       = .true.,        &
         flux_runoff= .true.,        &
         flux_param = (/ 1.0e-03 /), &
         flux_bottom= .true.         )
    !
    !       Aragonite (Sinking detrital/particulate CaCO3)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_arag',     &
         longname   = 'Detrital CaCO3', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.         )
    !
    !       Calcite 
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cadet_calc',     &
         longname   = 'Detrital CaCO3', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.         )
    !
    !       DIC (Dissolved inorganic carbon)
    !
    call g_tracer_add(tracer_list,package_name,                       &
         name       = 'dic',                                           &
         longname   = 'Dissolved Inorganic Carbon',                    &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'co2_flux',                                  &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMCO2,                                      &
         flux_gas_param = (/ as_coeff_cobalt, 9.7561e-06 /),           &
         ! Uncomment for "no mass change" check
         ! flux_gas_param = (/ 0.0, 0.0 /),                            &
         flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc',   &
         flux_runoff= .true.,                                          &
         flux_param = (/12.011e-03  /),                                &
         flux_bottom= .true.,                                          &
         init_value = 0.001                                            )
    !
    !       Dissolved Fe (assumed to be all available to phytoplankton)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fed',            &
         longname   = 'Dissolved Iron', & 
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_runoff= .true.,           &
         flux_wetdep= .true.,           &
         flux_drydep= .true.,           &
         flux_param = (/ 55.847e-03 /), &
         flux_bottom= .true.            )
    !
    !    Fedet (Sinking detrital/particulate iron)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fedet',         &
         longname   = 'Detrital Iron', &
         units      = 'mol/kg',        &
         prog       = .true.,          &
         sink_rate  = cobalt%wsink,     &
         btm_reservoir = .true.        )
    !
    !       Diazotroph Fe (Iron in N2-fixing phytoplankton for variable Fe:N ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fedi',            &
         longname   = 'Diazotroph Iron', &
         units      = 'mol/kg',          &
         prog       = .true.             )
    !
    !       Large Fe (Iron in large phytoplankton to allow for variable Fe:N ratios)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'felg',       &
         longname   = 'Large Phytoplankton Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       Small Fe
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'fesm',       &
         longname   = 'Small Phytoplankton Iron', &
         units      = 'mol/kg',     &
         prog       = .true.        )
    !
    !       LDON (Labile dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ldon',           &
         flux_runoff= .true.,         &
         longname   = 'labile DON',     &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)       ) 
    !
    !       LDOP (Labile dissolved organic phosphorous)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ldop',           &
         flux_runoff= .true.,         &
         longname   = 'labile DOP',     &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)       ) 
    !
    !       LITH (Lithogenic aluminosilicate particles)
    !
    call g_tracer_add(tracer_list,package_name,     &
         name       = 'lith',                       &
         longname   = 'Lithogenic Aluminosilicate', &
         units      = 'g/kg',                       &
         prog       = .true.,                       &
!         const_init_value = 0.0 ,                   &
         flux_runoff= .true.,                       &
         flux_wetdep= .true.,                       &
         flux_drydep= .true.,                       &
         flux_param = (/ 1.0e-03 /)                 )
    !
    !     LITHdet (Detrital Lithogenic aluminosilicate particles)  
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'lithdet',   &
         longname   = 'lithdet',   &
         units      = 'g/kg',      &
         prog       = .true.,      &
         sink_rate  = cobalt%wsink, &
         btm_reservoir = .true.    )
    !
    !       NBact: Bacteria nitrogen
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nbact',          &
         longname   = 'bacterial',      &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !    Ndet (Sinking detrital/particulate Nitrogen)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndet',      &
         longname   = 'ndet',      &
         flux_runoff= .true.,      &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         sink_rate  = cobalt%wsink,&
         btm_reservoir = .true.,   &
         flux_param = (/ 1.0e-3 /) ) 
    !
    !    NDi (assumed to be facultative N2-fixers, with a variable N:P ratio
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndi',                 &
         longname   = 'Diazotroph Nitrogen', &
         units      = 'mol/kg',              &
         prog       = .true.                 )

    !
    !    NLg (assumed to be a dynamic combination of diatoms and other 
    !         eukaryotes all effectively greater than 5 um in diameter,
    !         and having a fixed C:N ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nlg',            &
         longname   = 'Large Phytoplankton Nitrogen', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       NSm (Nitrogen in picoplankton and nanoplankton
    !            ~less than 5 um in diameter and having a fixed C:N:P ratio)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nsm',            &
         longname   = 'Small Phytoplankton Nitrogen', &
         units      = 'mol/kg',         &
         prog       = .true.            )
    !
    !       NH4
    !
    if (do_nh3_atm_ocean_exchange) then
       if (mpp_root_pe().eq.mpp_pe()) then
          write(*,*) 'setting up nh3 tracer (ocean)'
       end if
       call g_tracer_add(tracer_list,package_name,&
            name       = 'nh4',             &
            longname   = 'Ammonia',         &
            units      = 'mol/kg',          &
            prog       = .true.,            &
            flux_wetdep= .true.,            &
            flux_drydep= .true.,            &
            flux_gas   = .true.,            &
!            implementation='duce',   &
            implementation='johnson',   &
            flux_gas_name  = 'nh3_flux',    &
            flux_gas_type  = 'air_sea_gas_flux_generic', &
            flux_gas_molwt = WTMN,                       &
            flux_gas_param = (/ 17.,vb_nh3 /), &
            flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc',    &
            flux_param = (/ WTMN*1.e-3/),    &         
            flux_bottom= .true.,            &
            init_value = 0.001              )
    else
       if (mpp_root_pe().eq.mpp_pe()) then
          write(*,*) 'setting up nh3 tracer (ocean - no exchange)'
       end if

       call g_tracer_add(tracer_list,package_name,&                                 
         name       = 'nh4',             &         
         longname   = 'Ammonia',         &         
         units      = 'mol/kg',          &         
         prog       = .true.,            &         
         flux_wetdep= .true.,            &         
         flux_drydep= .true.,            &         
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )  
    end if
    !
    !       NO3
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'no3',             &
         longname   = 'Nitrate',         &
         units      = 'mol/kg',          &
         prog       = .true.,            &
         flux_runoff= .true.,            &
         flux_wetdep= .true.,            &
         flux_drydep= .true.,            &
         flux_param = (/ 14.0067e-03 /), &
         flux_bottom= .true.             )
    !
    !       O2
    !
    call g_tracer_add(tracer_list,package_name,                        &
         name       = 'o2',                                            &
         longname   = 'Oxygen',                                        &
         units      = 'mol/kg',                                        &
         prog       = .true.,                                          &
         flux_gas   = .true.,                                          &
         flux_gas_name  = 'o2_flux',                                   &
         flux_gas_type  = 'air_sea_gas_flux_generic',                  &
         flux_gas_molwt = WTMO2,                                       &
         flux_gas_param = (/ as_coeff_cobalt, 9.7561e-06 /),           &
         flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc',   &
         flux_bottom= .true.             )
    !
    !    Pdet (Sinking detrital/particulate Phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,         &
         name       = 'pdet',                           &
         longname   = 'Detrital Phosphorus',            &
         units      = 'mol/kg',                         &
         prog       = .true.,                           &
         sink_rate  = cobalt%wsink,                      &
         btm_reservoir = .true.    )
    !
    !       PO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'po4',       &
         longname   = 'Phosphate', &
         flux_runoff= .true.,      &
         units      = 'mol/kg',    &
         prog       = .true.,      &
         flux_wetdep= .true.,      &
         flux_drydep= .true.,      &
         flux_bottom= .true.,      &
         flux_param = (/ 1.0e-3 /) )
    !
    !       SRDON (Semi-Refractory dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'srdon',           &
         longname   = 'Semi-Refractory DON', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)  )
    !
    !       SRDOP (Semi-Refractory dissolved organic phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'srdop',           &
         longname   = 'Semi-Refractory DOP', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       SLDON (Semilabile dissolved organic nitrogen)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sldon',           &
         longname   = 'Semilabile DON', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       SLDOP (Semilabile dissolved organic phosphorus)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sldop',           &
         longname   = 'Semilabile DOP', &
         flux_runoff= .true.,           &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         flux_param = (/ 1.0e-3 /)      )
    !
    !       Sidet (Sinking detrital/particulate Silicon)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sidet',            &
         longname   = 'Detrital Silicon', &
         units      = 'mol/kg',           &
         prog       = .true.,             &
         sink_rate  = cobalt%wsink,        &
         btm_reservoir = .true.    )
    !
    !    SiLg (Silicon in large phytoplankton for variable Si:N ratios
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'silg',          &
         longname   = 'Large Phytoplankton Silicon', &
         units      = 'mol/kg',        &
         prog       = .true.           )
    !
    !       SiO4
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'sio4',     &
         longname   = 'Silicate', &
         units      = 'mol/kg',   &
         prog       = .true.,     &
         flux_bottom= .true.      )

    !
    !     Small zooplankton N  
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nsmz',     &
         longname   = 'Small Zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

    !
    !     Medium-sized zooplankton N
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nmdz',     &
         longname   = 'Medium-sized zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

    !
    !     Large zooplankton N (Pred zoo + krill)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'nlgz',     &
         longname   = 'large Zooplankton Nitrogen', &
         units      = 'mol/kg',   &
         prog       = .true.     ) 

      if (do_14c) then                                        !<<RADIOCARBON
      !       D14IC (Dissolved inorganic radiocarbon)
      !
      call g_tracer_add(tracer_list,package_name,       &
         name       = 'di14c',                          &
         longname   = 'Dissolved Inorganic Radiocarbon',&
         units      = 'mol/kg',                         &
         prog       = .true.,                           &
         flux_gas       = .true.,                       &
         flux_gas_name  = 'c14o2_flux',                 &
         flux_gas_type  = 'air_sea_gas_flux',           &
         flux_gas_molwt = WTMCO2,                       &
         flux_gas_param = (/ as_coeff_cobalt, 9.7561e-06 /),   &
         flux_gas_restart_file  = 'ocean_cobalt_airsea_flux.res.nc', &
         flux_runoff= .true.,                           &
         flux_param     = (/14.e-03  /),                &
         flux_bottom    = .true.,                       &
         init_value     = 0.001)
      !
      !       DO14C (Dissolved organic radiocarbon)
      !
      call g_tracer_add(tracer_list,package_name, &
         name       = 'do14c',                    &
         longname   = 'DO14C',                    &
         units      = 'mol/kg',                   &
         prog       = .true.)
      endif                                                   !RADIOCARBON>>

    !===========================================================
    !Diagnostic Tracers
    !===========================================================
    !
    !    Cased (CaCO3 in sediments)   
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'cased',          &
         longname   = 'Sediment CaCO3', &
         units      = 'mol m-3',       &
         prog       = .false.           )
    !
    !       Chl (Chlorophyll)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'chl',         &
         longname   = 'Chlorophyll', &
         units      = 'ug kg-1',     &
         prog       = .false.,       &
         init_value = 0.08           )
    !
    !       CO3_ion (Carbonate ion)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'co3_ion',       &
         longname   = 'Carbonate ion', &
         units      = 'mol/kg',        &
         prog       = .false.          )
    !
    !      cadet_arag_btf (Aragonite flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'cadet_arag_btf',          &
         longname   = 'aragonite flux to Sediments', &
         units      = 'mol m-2 s-1',             &
         prog       = .false.                    )
    !
    !      cadet_calc_btf (Calcite flux to sediments) 
    !
    call g_tracer_add(tracer_list,package_name,  &
         name       = 'cadet_calc_btf',          &
         longname   = 'calcite flux to Sediments', &
         units      = 'mol m-2 s-1',             &
         prog       = .false.                    )
    !
    !      lithdet_btf (Lithogenic flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'lithdet_btf',            &
         longname   = 'Lith flux to Sediments', &
         units      = 'g m-2 s-1',              &
         prog       = .false.                   )
    !
    !      ndet_btf (N flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'ndet_btf',            &
         longname   = 'N flux to Sediments', &
         units      = 'mol m-2 s-1',         &
         prog       = .false.                )
    !
    !      pdet_btf (P flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'pdet_btf',            &
         longname   = 'P flux to Sediments', &
         units      = 'mol m-2 s-1',         &
         prog       = .false.                )
    !
    !      sidet_btf (SiO2 flux to sediments)
    !
    call g_tracer_add(tracer_list,package_name, &
         name       = 'sidet_btf',              &
         longname   = 'SiO2 flux to Sediments', &
         units      = 'mol m-2 s-1',            &
         prog       = .false.                   )
    !
    !      fedet_btf (Fe flux to sediments)
    !      (only used in "no mass change" check)
    call g_tracer_add(tracer_list,package_name, &
         name       = 'fedet_btf',              &
         longname   = 'Fe flux to Sediments',   &
         units      = 'mol m-2 s-1',            &
         prog       = .false.                   )
    !
    !       htotal (H+ ion concentration)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'htotal',               &
         longname   = 'H+ ion concentration', &
         units      = 'mol/kg',               &
         prog       = .false.,                &
         init_value = cobalt%htotal_in         )
    !
    !       Irr_mem (Irradiance Memory)
    !
    call g_tracer_add(tracer_list,package_name,&
         name       = 'irr_mem',           &
         longname   = 'Irradiance memory', &
         units      = 'Watts/m^2',         &
         prog       = .false.              )

    call g_tracer_add(tracer_list,package_name,&
         name       = 'mu_mem_ndi',           &
         longname   = 'Growth memory', &
         units      = 'sec-1',         &
         prog       = .false.              )

    call g_tracer_add(tracer_list,package_name,&
         name       = 'mu_mem_nlg',           &
         longname   = 'Growth memory', &
         units      = 'sec-1',         &
         prog       = .false.              )

    call g_tracer_add(tracer_list,package_name,&
         name       = 'mu_mem_nsm',           &
         longname   = 'Growth memory', &
         units      = 'sec-1',         &
         prog       = .false.              )

    if (do_nh3_atm_ocean_exchange .or. scheme_nitrif.eq.2 .or. scheme_nitrif.eq.3) then
       call g_tracer_add(tracer_list,package_name,&
            name       = 'nh3',         &
            longname   = 'NH3', &
            units      = 'mol/kg',     &
            prog       = .false.,       &
            init_value = 1.e-10           )
    end if

  end subroutine user_add_tracers


  ! <SUBROUTINE NAME="generic_COBALT_update_from_coupler">
  !  <OVERVIEW>
  !   Modify the values obtained from the coupler if necessary.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Some tracer fields need to be modified after values are obtained from the coupler.
  !   This subroutine is the place for specific tracer manipulations.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_coupler(tracer_list) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_coupler(tracer_list)
    type(g_tracer_type), pointer :: tracer_list

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_update_from_copler'

    real, dimension(:,:)  ,pointer    :: stf_alk,dry_no3,wet_no3

    !
    ! NO3 has deposition, river flux, and negative deposition contribution to alkalinity
    !
    call g_tracer_get_pointer(tracer_list,'no3','drydep',dry_no3)
    call g_tracer_get_pointer(tracer_list,'no3','wetdep',wet_no3)

    call g_tracer_get_pointer(tracer_list,'alk','stf',stf_alk)

    stf_alk = stf_alk - dry_no3 - wet_no3 ! update 'tracer%stf' thru pointer

    return
  end subroutine generic_COBALT_update_from_coupler

  ! <SUBROUTINE NAME="generic_COBALT_update_from_bottom">
  ! 
  !  <OVERVIEW>
  !   Set values of bottom fluxes and reservoirs
  !  </OVERVIEW>
  !
  !  <DESCRIPTION>
  !
  !   This routine calculates bottom fluxes for tracers with bottom reservoirs.
  !   It is called near the end of the time step, meaning that the fluxes 
  !   calculated pertain to the next time step.  
  ! 
  !  </DESCRIPTION>
  !
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_bottom(tracer_list,dt, tau, model_time) 
  !  </TEMPLATE>
  !
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment 
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index to be used for %field
  !  </IN>
  !
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_bottom(tracer_list, dt, tau, model_time)
    type(g_tracer_type), pointer :: tracer_list
    real,               intent(in) :: dt
    integer,            intent(in) :: tau
    type(time_type),    intent(in) :: model_time
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau
    logical :: used
    real, dimension(:,:,:),pointer :: grid_tmask
    real, dimension(:,:,:),pointer :: temp_field

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask) 

    !
    ! The bottom reservoirs of aragonite and calcite are immediately redistributed to the
    ! water column as a bottom flux (btf) where they impact the alkalinity and DIC
    !
    call g_tracer_get_values(tracer_list,'cadet_arag','btm_reservoir',cobalt%fcadet_arag_btm,isd,jsd)
    cobalt%fcadet_arag_btm = cobalt%fcadet_arag_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_arag_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_arag_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_arag','btm_reservoir',0.0)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_arag_btm,cobalt%fcadet_arag_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'cadet_calc','btm_reservoir',cobalt%fcadet_calc_btm,isd,jsd)
    cobalt%fcadet_calc_btm = cobalt%fcadet_calc_btm/dt
    call g_tracer_get_pointer(tracer_list,'cadet_calc_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fcadet_calc_btm(:,:)
    call g_tracer_set_values(tracer_list,'cadet_calc','btm_reservoir',0.0)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_calc_btm, cobalt%fcadet_calc_btm, &
         model_time, rmask = grid_tmask(:,:,1), &
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Iron is buried, but can re-enter the water column in association with
    ! organic matter degradation (see ffe_sed in update_from_source)
    !
    call g_tracer_get_values(tracer_list,'fedet','btm_reservoir',cobalt%ffedet_btm,isd,jsd)
    cobalt%ffedet_btm = cobalt%ffedet_btm/dt
    ! uncomment for "no mass change check"
    !call g_tracer_get_pointer(tracer_list,'fedet_btf','field',temp_field)
    !temp_field(:,:,1) = cobalt%ffedet_btm(:,:)
    call g_tracer_set_values(tracer_list,'fedet','btm_reservoir',0.0)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_ffedet_btm, cobalt%ffedet_btm, &
         model_time, rmask = grid_tmask(:,:,1), & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! Lithogenic material is buried
    !
    call g_tracer_get_values(tracer_list,'lithdet','btm_reservoir',cobalt%flithdet_btm,isd,jsd)
    cobalt%flithdet_btm = cobalt%flithdet_btm /dt
    call g_tracer_get_pointer(tracer_list,'lithdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%flithdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'lithdet','btm_reservoir',0.0)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_flithdet_btm, cobalt%flithdet_btm, &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    !
    ! N, P, and Si detritus that hits the bottom is re-entered as a bottom source of 
    ! nh4, po4, and SiO4 respectively
    !
    call g_tracer_get_values(tracer_list,'ndet','btm_reservoir',cobalt%fndet_btm,isd,jsd)
    cobalt%fndet_btm = cobalt%fndet_btm/dt
    call g_tracer_get_pointer(tracer_list,'ndet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fndet_btm(:,:)
    call g_tracer_set_values(tracer_list,'ndet','btm_reservoir',0.0)
    if (cobalt%id_fndet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fndet_btm,cobalt%fndet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'pdet','btm_reservoir',cobalt%fpdet_btm,isd,jsd)
    cobalt%fpdet_btm = cobalt%fpdet_btm/dt
    call g_tracer_get_pointer(tracer_list,'pdet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fpdet_btm(:,:)
    call g_tracer_set_values(tracer_list,'pdet','btm_reservoir',0.0)
    if (cobalt%id_fpdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fpdet_btm,cobalt%fpdet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    call g_tracer_get_values(tracer_list,'sidet','btm_reservoir',cobalt%fsidet_btm,isd,jsd)
    cobalt%fsidet_btm = cobalt%fsidet_btm/dt
    call g_tracer_get_pointer(tracer_list,'sidet_btf','field',temp_field)
    temp_field(:,:,1) = cobalt%fsidet_btm(:,:)
    call g_tracer_set_values(tracer_list,'sidet','btm_reservoir',0.0)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fsidet_btm,    cobalt%fsidet_btm,          &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

  end subroutine generic_COBALT_update_from_bottom

  ! <SUBROUTINE NAME="generic_COBALT_update_from_source">
  !  <OVERVIEW>
  !   Update tracer concentration fields due to the source/sink contributions.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This is the subroutine to contain most of the biogeochemistry for calculating the 
  !   interaction of tracers with each other and with outside forcings.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_update_from_source(tracer_list,Temp,Salt,dzt,hblt_depth,&
  !                                         ilb,jlb,tau,dt, grid_dat,sw_pen,opacity) 
  !  </TEMPLATE>
  !  <IN NAME="tracer_list" TYPE="type(g_tracer_type), pointer">
  !   Pointer to the head of generic tracer list.
  !  </IN>
  !  <IN NAME="ilb,jlb" TYPE="integer">
  !   Lower bounds of x and y extents of input arrays on data domain
  !  </IN>
  !  <IN NAME="Temp" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean temperature   
  !  </IN>
  !  <IN NAME="Salt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean salinity
  !  </IN>
  !  <IN NAME="dzt" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean layer thickness (meters)
  !  </IN>
  !  <IN NAME="opacity" TYPE="real, dimension(ilb:,jlb:,:)">
  !   Ocean opacity
  !  </IN>
  !  <IN NAME="sw_pen" TYPE="real, dimension(ilb:,jlb:)">
  !   Shortwave peneteration
  !  </IN>
  !  <IN NAME="hblt_depth" TYPE="real, dimension(ilb:,jlb:)">
  !   
  !  </IN>
  !  <IN NAME="grid_dat" TYPE="real, dimension(ilb:,jlb:)">
  !   Grid area
  !  </IN>
  !  <IN NAME="tau" TYPE="integer">
  !   Time step index of %field
  !  </IN>
  !  <IN NAME="dt" TYPE="real">
  !   Time step increment
  !  </IN>
  ! </SUBROUTINE>
  subroutine generic_COBALT_update_from_source(tracer_list,Temp,Salt,rho_dzt,dzt,hblt_depth,&
       ilb,jlb,tau,dt,grid_dat,model_time,nbands,max_wavelength_band,sw_pen_band,opacity_band,internal_heat,frunoff)

    type(g_tracer_type),            pointer    :: tracer_list
    real, dimension(ilb:,jlb:,:),   intent(in) :: Temp,Salt,rho_dzt,dzt
    real, dimension(ilb:,jlb:),     intent(in) :: hblt_depth
    integer,                        intent(in) :: ilb,jlb,tau
    real,                           intent(in) :: dt
    real, dimension(ilb:,jlb:),     intent(in) :: grid_dat
    type(time_type),                intent(in) :: model_time

    integer,                        intent(in) :: nbands
    real, dimension(:),             intent(in) :: max_wavelength_band
    real, dimension(:,ilb:,jlb:),   intent(in) :: sw_pen_band
    real, dimension(:,ilb:,jlb:,:), intent(in) :: opacity_band
    real, dimension(ilb:,jlb:),     intent(in) :: internal_heat
    real, dimension(ilb:,jlb:),     intent(in) :: frunoff 

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_update_from_source'
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau, i, j, k , kblt, m, n, k_100, k_200, kbot
    real, dimension(:,:,:) ,pointer :: grid_tmask
    integer, dimension(:,:),pointer :: mask_coast,grid_kmt
    !
    !------------------------------------------------------------------------
    ! Local Variables
    !------------------------------------------------------------------------
    !
    logical :: used, first
    integer :: nb
    real :: r_dt
    real :: feprime_temp
    real :: juptake_di_tot2nterm
    real :: log_btm_flx
    real :: P_C_m
    real :: p_lim_nhet
    real :: TK, PRESS, PKSPA, PKSPC
    real :: tmp_hblt, tmp_irrad, tmp_irrad_ML,tmp_opacity,tmp_mu_ML
    real :: drho_dzt
    real, dimension(:), Allocatable   :: tmp_irr_band
    real, dimension(:,:), Allocatable :: rho_dzt_100, rho_dzt_200
    real, dimension(:,:,:), Allocatable :: z_remin_ramp
    real,dimension(1:NUM_ZOO,1:NUM_PREY) :: ipa_matrix,pa_matrix,ingest_matrix
    real,dimension(1:NUM_PREY) :: hp_ipa_vec,hp_pa_vec,hp_ingest_vec
    real,dimension(1:NUM_PREY) :: prey_vec,prey_p2n_vec,prey_fe2n_vec,prey_si2n_vec
    real,dimension(1:NUM_ZOO)  :: tot_prey
    real :: tot_prey_hp, sw_fac_denom, assim_eff
    real :: bact_uptake_ratio, vmax_bact, growth_ratio
    real :: fpoc_btm, log_fpoc_btm
    real :: fe_salt
    real :: sal,tt,tkb,ts,ts2,ts3,ts4,ts5

    logical ::  phos_nh3_override

    real, dimension(:,:,:), Allocatable :: pre_totn, net_srcn, post_totn
    real, dimension(:,:,:), Allocatable :: pre_totp, post_totp
    real, dimension(:,:,:), Allocatable :: pre_totsi, post_totsi
    real, dimension(:,:,:), Allocatable :: pre_totfe, net_srcfe, post_totfe
    real, dimension(:,:,:), Allocatable :: pre_totc, net_srcc, post_totc
    real, dimension(:,:),   Allocatable :: pka_nh3,phos_nh3_exchange

    real :: tr,ltr
    real :: imbal
    integer :: stdoutunit, imbal_flag, outunit


    r_dt = 1.0 / dt

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
         grid_tmask=grid_tmask,grid_mask_coast=mask_coast,grid_kmt=grid_kmt)

    call mpp_clock_begin(id_clock_carbon_calculations)
    !Get necessary fields
    call g_tracer_get_values(tracer_list,'htotal','field', cobalt%f_htotal,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'po4'   ,'field', cobalt%f_po4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'sio4'  ,'field', cobalt%f_sio4,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'alk'   ,'field', cobalt%f_alk,isd,jsd,ntau=tau)
    call g_tracer_get_values(tracer_list,'dic'   ,'field', cobalt%f_dic  ,isd,jsd,ntau=tau)
    if (do_nh3_diag) then
       allocate(pka_nh3(isd:ied,jsd:jed))
       pka_nh3 = 0.
       call g_tracer_get_values(tracer_list,'nh4'   ,'field', cobalt%f_nh4  ,isd,jsd,ntau=tau)
    end if
    allocate(phos_nh3_exchange(isd:ied,jsd:jed))

    !---------------------------------------------------------------------
    !Calculate co3_ion
    !Also calculate co2 fluxes csurf and alpha for the next round of exchange
    !---------------------------------------------------------------------
   
    cobalt%zt = 0.0
    cobalt%zm = 0.0
    do j = jsc, jec ; do i = isc, iec   !{
       cobalt%zt(i,j,1) = dzt(i,j,1)
       cobalt%zm(i,j,1) = 0.5*dzt(i,j,1)
    enddo; enddo !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%zt(i,j,k) = cobalt%zt(i,j,k-1) + dzt(i,j,k)
       cobalt%zm(i,j,k) = cobalt%zm(i,j,k-1) + dzt(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    k=1
    do j = jsc, jec ; do i = isc, iec  !{
       cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * cobalt%f_htotal(i,j,k)
       cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * cobalt%f_htotal(i,j,k)
    enddo; enddo ; !} i, j


    call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
         Temp(:,:,k), Salt(:,:,k),                    &
         cobalt%f_dic(:,:,k),                          &
         cobalt%f_po4(:,:,k),                          &  
         cobalt%f_sio4(:,:,k),                         &
         cobalt%f_alk(:,:,k),                          &
         cobalt%htotallo, cobalt%htotalhi,&
                                !InOut
         cobalt%f_htotal(:,:,k),                       & 
                                !Optional In
         co2_calc=trim(co2_calc),                      & 
         zt=cobalt%zt(:,:,k),                          & 
                                !OUT
         co2star=cobalt%co2_csurf(:,:), alpha=cobalt%co2_alpha(:,:), &
         pCO2surf=cobalt%pco2_csurf(:,:), &
         co3_ion=cobalt%f_co3_ion(:,:,k), &
         omega_arag=cobalt%omegaa(:,:,k), &
         omega_calc=cobalt%omegac(:,:,k))

    do k = 2, nk
       do j = jsc, jec ; do i = isc, iec  !{
          cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * cobalt%f_htotal(i,j,k)
          cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * cobalt%f_htotal(i,j,k)
       enddo; enddo ; !} i, j
  
       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,k),&
            Temp(:,:,k), Salt(:,:,k),                    &
            cobalt%f_dic(:,:,k),                          &
            cobalt%f_po4(:,:,k),                          &  
            cobalt%f_sio4(:,:,k),                         &
            cobalt%f_alk(:,:,k),                          &
            cobalt%htotallo, cobalt%htotalhi,&
                                !InOut
            cobalt%f_htotal(:,:,k),                       & 
                                !Optional In
            co2_calc=trim(co2_calc),                      & 
            zt=cobalt%zt(:,:,k),                          & 
                                !OUT
            co3_ion=cobalt%f_co3_ion(:,:,k), &
            omega_arag=cobalt%omegaa(:,:,k), &
            omega_calc=cobalt%omegac(:,:,k))
    enddo

    call g_tracer_set_values(tracer_list,'htotal','field',cobalt%f_htotal  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'dic','alpha',cobalt%co2_alpha    ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','csurf',cobalt%co2_csurf    ,isd,jsd)

    call mpp_clock_end(id_clock_carbon_calculations)

    if (do_nh3_atm_ocean_exchange) then
       !to override pH used for ocean nh3 exchange
       call data_override('OCN', 'phos_nh3_exchange', phos_nh3_exchange(isc:iec,jsc:jec), model_time,override=phos_nh3_override)       

       do j = jsc, jec ; do i = isc, iec 
          pka_nh3(i,j)   = calc_pka_nh3(temp(i,j,1),salt(i,j,1))*grid_tmask(i,j,1)
          tr             = 298.15/(temp(i,j,1)+273.15)-1.
          ltr            = -tr+log(298.15/(temp(i,j,1)+273.15))       
          !f1p
          !henry's constant is from Mark Jacobson's book "Fundamental of Atmospheric Modeling".
          !I used this expression to be consistent with the cloud chemistry module
          !the units are in mol/kg/atm. However it's probably in mol/kg(pure water)/atm
          !the density of pure water is ~997 kg/m3 (25C). alpha will then be scaled by the density of seawater, which for some reason I don't quite understand is always set to 1035. 
          !to be consistent, I am scaling the Jacobson's number by 997/1035. This decreases the solubility of NH3 by less than 4%. For references, Sander's estimate is alpha(298.15)=0.59*101.63=59.96 mol/kg(pure water)/atm, about 4% greater than Jacobson's.
          cobalt%nh3_alpha(i,j)   =  5.76e1*exp(13.79*tr-5.39*ltr)*997. 
          !apply salinity correction
          cobalt%nh3_alpha(i,j)   =  cobalt%nh3_alpha(i,j)/saltout_correction(101325./(1.e-3*rdgas*wtmair*(temp(i,j,1)+273.15)*cobalt%nh3_alpha(i,j)),vb_nh3,salt(i,j,1))* 1./cobalt%Rho_0 !mol/kg/atm  
          if (phos_nh3_override) then
             cobalt%nh3_csurf(i,j)   =  cobalt%f_nh4(i,j,1)/(1.+10**(pka_nh3(i,j)-max(min(phos_nh3_exchange(i,j),11.),3.))) !in mol/kg
          else
             cobalt%nh3_csurf(i,j)   =  cobalt%f_nh4(i,j,1)/(1.+10**(pka_nh3(i,j)+log10(min(max(cobalt%f_htotal(i,j,1),1e-11),1e-3)))) !in mol/kg
          end if
          cobalt%pnh3_csurf(i,j)  =  cobalt%nh3_csurf(i,j)/cobalt%nh3_alpha(i,j)*1.e6 !in uatm
       enddo; enddo ; !

       call g_tracer_set_values(tracer_list,'nh4','alpha',cobalt%nh3_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'nh4','csurf',cobalt%nh3_csurf    ,isd,jsd)         
    end if

      if (do_14c) then                                        !<<RADIOCARBON
      
      ! Normally, the alpha would be multiplied by the atmospheric 14C/12C ratio. However,
      ! here that is set to 1, so that alpha_14C = alpha_12C. This needs to be changed!

   call g_tracer_get_values(tracer_list,'di14c' ,'field', cobalt%f_di14c,isd,jsd,ntau=tau,positive=.true.)

    ! This is not used until later, but get it now
    call g_tracer_get_values(tracer_list,'do14c' ,'field', cobalt%f_do14c,isd,jsd,ntau=tau,positive=.true.)
    
       do j = jsc, jec ; do i = isc, iec  !{
       cobalt%c14o2_csurf(i,j) =  cobalt%co2_csurf(i,j) *                &
         cobalt%f_di14c(i,j,1) / (cobalt%f_dic(i,j,1) + epsln)
       cobalt%c14o2_alpha(i,j) =  cobalt%co2_alpha(i,j)
       enddo; enddo ; !} i, j

    call g_tracer_set_values(tracer_list,'di14c','alpha',cobalt%c14o2_alpha      ,isd,jsd)
    call g_tracer_set_values(tracer_list,'di14c','csurf',cobalt%c14o2_csurf      ,isd,jsd)

      endif                                                   !RADIOCARBON>>

    !---------------------------------------------------------------------
    ! Get positive tracer concentrations
    !---------------------------------------------------------------------

    call mpp_clock_begin(id_clock_phyto_growth)

    call g_tracer_get_values(tracer_list,'cadet_arag','field',cobalt%f_cadet_arag ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'cadet_calc','field',cobalt%f_cadet_calc ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fed'    ,'field',cobalt%f_fed      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fedet'  ,'field',cobalt%f_fedet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ldon'   ,'field',cobalt%f_ldon     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ldop'   ,'field',cobalt%f_ldop     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lith'   ,'field',cobalt%f_lith     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'lithdet','field',cobalt%f_lithdet  ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndet'   ,'field',cobalt%f_ndet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nh4'    ,'field',cobalt%f_nh4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'no3'    ,'field',cobalt%f_no3      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'o2'     ,'field',cobalt%f_o2       ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'pdet'   ,'field',cobalt%f_pdet     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'po4'    ,'field',cobalt%f_po4      ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'srdon'   ,'field',cobalt%f_srdon   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'srdop'   ,'field',cobalt%f_srdop   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sldon'   ,'field',cobalt%f_sldon   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sldop'   ,'field',cobalt%f_sldop   ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sidet'  ,'field',cobalt%f_sidet    ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'sio4'   ,'field',cobalt%f_sio4     ,isd,jsd,ntau=tau,positive=.true.)
!
    ! phytoplankton fields
    !
    call g_tracer_get_values(tracer_list,'fedi'   ,'field',phyto(DIAZO)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'felg'   ,'field',phyto(LARGE)%f_fe(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'fesm','field',phyto(SMALL)%f_fe(:,:,:),isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'ndi'    ,'field',phyto(DIAZO)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nlg'    ,'field',phyto(LARGE)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nsm' ,'field',phyto(SMALL)%f_n(:,:,:),isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'silg'   ,'field',cobalt%f_silg     ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'mu_mem_ndi' ,'field',phyto(DIAZO)%f_mu_mem,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'mu_mem_nlg' ,'field',phyto(LARGE)%f_mu_mem,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'mu_mem_nsm' ,'field',phyto(SMALL)%f_mu_mem,isd,jsd,ntau=1)
    !
    ! zooplankton fields
    !
    call g_tracer_get_values(tracer_list,'nsmz'    ,'field',zoo(1)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nmdz'    ,'field',zoo(2)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    call g_tracer_get_values(tracer_list,'nlgz'    ,'field',zoo(3)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! bacteria
    !
    call g_tracer_get_values(tracer_list,'nbact'   ,'field',bact(1)%f_n(:,:,:) ,isd,jsd,ntau=tau,positive=.true.)
    !
    ! diagnostic tracers that are passed between time steps (except chlorophyll)
    !
    call g_tracer_get_values(tracer_list,'cased'  ,'field',cobalt%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1,positive=.true.)
    call g_tracer_get_values(tracer_list,'cadet_arag_btf','field',cobalt%f_cadet_arag_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'cadet_calc_btf','field',cobalt%f_cadet_calc_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'lithdet_btf','field',cobalt%f_lithdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'ndet_btf','field',cobalt%f_ndet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'pdet_btf','field',cobalt%f_pdet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'sidet_btf','field',cobalt%f_sidet_btf,isd,jsd,ntau=1)
    ! uncomment for "no mass change" test
    !call g_tracer_get_values(tracer_list,'fedet_btf','field',cobalt%f_fedet_btf,isd,jsd,ntau=1)
    call g_tracer_get_values(tracer_list,'irr_mem','field',cobalt%f_irr_mem ,isd,jsd,ntau=1)

    ! zero out cumulative COBALT-wide production diagnostics
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec
       cobalt%jprod_fed(i,j,k) = 0.0
       cobalt%jprod_fedet(i,j,k) = 0.0
       cobalt%jprod_ndet(i,j,k) = 0.0
       cobalt%jprod_pdet(i,j,k) = 0.0
       cobalt%jprod_sldon(i,j,k) = 0.0
       cobalt%jprod_ldon(i,j,k) = 0.0
       cobalt%jprod_srdon(i,j,k) = 0.0
       cobalt%jprod_sldop(i,j,k) = 0.0
       cobalt%jprod_ldop(i,j,k) = 0.0
       cobalt%jprod_srdop(i,j,k) = 0.0
       cobalt%jprod_sidet(i,j,k) = 0.0
       cobalt%jprod_sio4(i,j,k) = 0.0
       cobalt%jprod_po4(i,j,k) = 0.0
       cobalt%jprod_nh4(i,j,k) = 0.0
       cobalt%jno3denit_wc(i,j,k) = 0.0
! added jgj - do we need it every timestep
       cobalt%jremin_ndet(i,j,k) = 0.0
       cobalt%jo2resp_wc(i,j,k) = 0.0
    enddo;  enddo ;  enddo !} i,j,k
!
!-----------------------------------------------------------------------------------
! 1: Phytoplankton growth and nutrient uptake calculations
!-----------------------------------------------------------------------------------
!
    !
    !-----------------------------------------------------------------------------------
    ! 1.1: Nutrient Limitation Calculations
    !-----------------------------------------------------------------------------------
    !
    ! Calculate iron cell quota
    !
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec
       do n = 1,NUM_PHYTO    !{
          phyto(n)%q_fe_2_n(i,j,k) = max(0.0, phyto(n)%f_fe(i,j,k)/ &
                 max(epsln,phyto(n)%f_n(i,j,k)))
          phyto(n)%q_p_2_n(i,j,k) = phyto(n)%p_2_n_static
       enddo  !} n
       !
       ! N limitation with NH4 inhibition after Frost and Franzen (1992)
       ! (Note nitrate does not limit diazotroph growth but uptake limitation is used 
       ! to determine nitrogen fixation versus facultative no3/nh4 uptake, see Sec. 1.3)
       do n= 1, NUM_PHYTO   !{
          if (scheme_no3_nh4_lim .eq. 1) then
             phyto(n)%no3lim(i,j,k) = cobalt%f_no3(i,j,k) / &
                  ( (phyto(n)%k_no3+cobalt%f_no3(i,j,k)) * (1.0 + cobalt%f_nh4(i,j,k)/phyto(n)%k_nh4) )
             phyto(n)%nh4lim(i,j,k) = cobalt%f_nh4(i,j,k) / (phyto(n)%k_nh4 + cobalt%f_nh4(i,j,k))             
          elseif (scheme_no3_nh4_lim .eq. 2) then
             phyto(n)%no3lim(i,j,k) = cobalt%f_no3(i,j,k) / &
                  ( phyto(n)%k_no3+cobalt%f_no3(i,j,k)+phyto(n)%k_no3/phyto(n)%k_nh4*cobalt%f_nh4(i,j,k) )
             phyto(n)%nh4lim(i,j,k) = cobalt%f_nh4(i,j,k) / &
                  ( phyto(n)%k_nh4+cobalt%f_nh4(i,j,k)+phyto(n)%k_nh4/phyto(n)%k_no3*cobalt%f_no3(i,j,k) )
          end if
       enddo !} n
       !
       ! O2 inhibition term for diazotrophs
       !
       n = DIAZO
       phyto(n)%o2lim(i,j,k) = (1.0 - cobalt%f_o2(i,j,k)**cobalt%o2_inhib_Di_pow / &
         (cobalt%f_o2(i,j,k)**cobalt%o2_inhib_Di_pow+cobalt%o2_inhib_Di_sat**cobalt%o2_inhib_Di_pow))
       !
       ! SiO4, PO4 and Fe uptake limitation with Michaelis-Mentin 
       !
       phyto(LARGE)%silim(i,j,k) = cobalt%f_sio4(i,j,k) / (phyto(LARGE)%k_sio4 + cobalt%f_sio4(i,j,k))
       do n= 1, NUM_PHYTO   !{
          phyto(n)%po4lim(i,j,k) = cobalt%f_po4(i,j,k) / (phyto(n)%k_po4 + cobalt%f_po4(i,j,k))
          phyto(n)%felim(i,j,k)  = cobalt%f_fed(i,j,k) / (phyto(n)%k_fed + cobalt%f_fed(i,j,k))
          phyto(n)%def_fe(i,j,k) = phyto(n)%q_fe_2_n(i,j,k)**2.0 / (phyto(n)%k_fe_2_n**2.0 +  &
               phyto(n)%q_fe_2_n(i,j,k)**2.0)
       enddo !} n
    enddo;  enddo ;  enddo !} i,j,k
    !
    ! Calculate nutrient limitation based on the most limiting nutrient (liebig_lim)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       n=DIAZO
       phyto(n)%liebig_lim(i,j,k) = phyto(n)%o2lim(i,j,k)* &
          min(phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
       do n= 2, NUM_PHYTO   !{
          phyto(n)%liebig_lim(i,j,k) = min(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k),&
             phyto(n)%po4lim(i,j,k), phyto(n)%def_fe(i,j,k))
       enddo !} n
    enddo;  enddo ;  enddo !} i,j,k
    !
    !-----------------------------------------------------------------------
    ! 1.2: Light Limitation/Growth Calculations
    !-----------------------------------------------------------------------
    !
    ! Create relevant light fields based on incident radiation and opacity
    ! information passed from the ocean code
    !
    allocate(tmp_irr_band(nbands))
    do j = jsc, jec ; do i = isc, iec   !{

       do nb=1,nbands !{
          if (max_wavelength_band(nb) .lt. 710.0) then !{
             tmp_irr_band(nb) = sw_pen_band(nb,i,j)
          else
             tmp_irr_band(nb) = 0.0
          endif !}
       enddo !}

       kblt = 0 ; tmp_irrad_ML = 0.0 ; tmp_hblt = 0.0
       do k = 1, nk !{
          tmp_irrad = 0.0
          do nb=1,nbands !{
             tmp_opacity = opacity_band(nb,i,j,k)
             tmp_irrad = tmp_irrad + max(0.0,tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k) * 0.5))
             ! Change tmp_irr_band from being the value at the middle of layer k to the value
             ! at the bottom of layer k.
             tmp_irr_band(nb) = tmp_irr_band(nb) * exp(-tmp_opacity * dzt(i,j,k))
          enddo !}
          cobalt%irr_inst(i,j,k) = tmp_irrad * grid_tmask(i,j,k)
          cobalt%irr_mix(i,j,k) = tmp_irrad * grid_tmask(i,j,k)
          if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
             kblt = kblt+1
             tmp_irrad_ML = tmp_irrad_ML + cobalt%irr_mix(i,j,k) * dzt(i,j,k)
             tmp_hblt = tmp_hblt + dzt(i,j,k)
          endif !}
       enddo !} k-loop
       cobalt%irr_mix(i,j,1:kblt) = tmp_irrad_ML / max(1.0e-6,tmp_hblt)
    enddo;  enddo !} i,j

    deallocate(tmp_irr_band)
    !
    ! Calculate the temperature limitation (expkT) and the time integrated
    ! irradiance (f_irr_mem) to which the Chl:C ratio responds (~24 hours)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        
       cobalt%expkT(i,j,k) = exp(cobalt%kappa_eppley * Temp(i,j,k))
       cobalt%f_irr_mem(i,j,k) = (cobalt%f_irr_mem(i,j,k) + (cobalt%irr_mix(i,j,k) - &
          cobalt%f_irr_mem(i,j,k)) * min(1.0,cobalt%gamma_irr_mem * dt)) * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    !nh3
    if (do_nh3_diag) then
    cobalt%f_nh3(:,:,:) = 0.
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%f_nh3(i,j,k) = cobalt%f_nh4(i,j,k)/(1.+10**(calc_pka_nh3(temp(i,j,k),salt(i,j,k))+log10(min(max(cobalt%f_htotal(i,j,1),1e-10),1e-5)))) * grid_tmask(i,j,k)
    enddo;  enddo ; enddo !} i,j,k
    end if


    !
    ! Phytoplankton growth rate calculation based on Geider et al. (1997)
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%f_chl(i,j,k) = 0.0

       do n = 1, NUM_PHYTO   !{
          P_C_m = max(phyto(n)%liebig_lim(i,j,k)*phyto(n)%P_C_max*cobalt%expkT(i,j,k),epsln)
          phyto(n)%theta(i,j,k) = (phyto(n)%thetamax-cobalt%thetamin) / (1.0 +                   &
             phyto(n)%thetamax*phyto(n)%alpha*cobalt%f_irr_mem(i,j,k)*0.5 /  &
             P_C_m) + cobalt%thetamin
          cobalt%f_chl(i,j,k) = cobalt%f_chl(i,j,k)+cobalt%c_2_n*12.0e6*phyto(n)%theta(i,j,k)*   &
             phyto(n)%f_n(i,j,k)
          phyto(n)%irrlim(i,j,k) = (1.0-exp(-phyto(n)%alpha*cobalt%irr_inst(i,j,k)*              &
             phyto(n)%theta(i,j,k)/P_C_m))

          ! calculate the growth rate
          phyto(n)%mu(i,j,k) = P_C_m / (1.0 + cobalt%zeta) * phyto(n)%irrlim(i,j,k) - &
             cobalt%expkT(i,j,k)*phyto(n)%bresp*                                      &
             phyto(n)%f_n(i,j,k)/(cobalt%refuge_conc + phyto(n)%f_n(i,j,k))

          ! calculate net production by phytoplankton group
          phyto(n)%jprod_n(i,j,k) = phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)

          phyto(n)%mu_mix(i,j,k) = phyto(n)%mu(i,j,k)

       enddo !} n

    enddo;  enddo ; enddo !} i,j,k

    do j = jsc, jec ; do i = isc, iec ; do n = 1,NUM_PHYTO !{
       kblt = 0 ; tmp_mu_ML = 0.0 ; tmp_hblt = 0.0
       do k = 1, nk !{
          if ((k == 1) .or. (tmp_hblt .lt. hblt_depth(i,j))) then !{
             kblt = kblt+1
             tmp_mu_ML = tmp_mu_ML + phyto(n)%mu_mix(i,j,k) * dzt(i,j,k)
             tmp_hblt = tmp_hblt + dzt(i,j,k)
          endif !}
       enddo !} k-loop
       phyto(n)%mu_mix(i,j,1:kblt) = tmp_mu_ML / max(epsln,tmp_hblt)
    enddo;  enddo; enddo !} i,j,n

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; do n = 1,NUM_PHYTO !{        
       phyto(n)%f_mu_mem(i,j,k) = phyto(n)%f_mu_mem(i,j,k) + (phyto(n)%mu_mix(i,j,k) - &
             phyto(n)%f_mu_mem(i,j,k))*min(1.0,cobalt%gamma_mu_mem*dt)*grid_tmask(i,j,k)
    enddo; enddo ; enddo; enddo !} i,j,k,n

    !-----------------------------------------------------------------------
    ! 1.3: Nutrient uptake calculations 
    !-----------------------------------------------------------------------
    !
    ! Uptake of nitrate and ammonia
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       n = DIAZO
       phyto(n)%juptake_n2(i,j,k) =  max(0.0,(1.0 - phyto(n)%no3lim(i,j,k) - phyto(n)%nh4lim(i,j,k))* &
          phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       phyto(n)%juptake_nh4(i,j,k) = max(0.0,phyto(n)%nh4lim(i,j,k)* phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       phyto(n)%juptake_no3(i,j,k) = max(0.0,phyto(n)%no3lim(i,j,k)* phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)) 
       ! uncomment for "no mass change" test (next 2 lines)
       ! phyto(n)%juptake_nh4(i,j,k) = phyto(n)%juptake_nh4(i,j,k) + phyto(n)%juptake_n2(i,j,k)
       ! phyto(n)%juptake_n2(i,j,k) = 0.0

       ! If growth is negative, results in net respiration and production of nh4, aerobic loss in all cases
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) - min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
       cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) - min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))*cobalt%o2_2_nh4
       do n = 2, NUM_PHYTO !{
          phyto(n)%juptake_no3(i,j,k) = max( 0.0, phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*   & 
             phyto(n)%no3lim(i,j,k)/(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k)+epsln) )
          phyto(n)%juptake_nh4(i,j,k) = max( 0.0, phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k)*   & 
             phyto(n)%nh4lim(i,j,k)/(phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k)+epsln) )
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) - min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))
          cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) - min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))*cobalt%o2_2_nh4
       enddo !} n
    enddo;  enddo ; enddo !} i,j,k
    !
    ! Phosphorous uptake
    ! 
    do k = 1, nk  ;    do j = jsc, jec ;      do i = isc, iec   !{
       n=DIAZO
       phyto(n)%juptake_po4(i,j,k) = (phyto(n)%juptake_n2(i,j,k)+phyto(n)%juptake_nh4(i,j,k) + &
          phyto(n)%juptake_no3(i,j,k))*phyto(n)%p_2_n_static
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) - &
          min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))*phyto(n)%p_2_n_static
       do n = 2, NUM_PHYTO
          phyto(n)%juptake_po4(i,j,k) = (phyto(n)%juptake_no3(i,j,k)+   &
                  phyto(n)%juptake_nh4(i,j,k)) * phyto(n)%p_2_n_static
          cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) - &
                  min(0.0,phyto(n)%mu(i,j,k)*phyto(n)%f_n(i,j,k))*phyto(n)%p_2_n_static
       enddo !} n
    enddo; enddo ; enddo !} i,j,k
    !
    ! Iron uptake
    ! 
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       do n = 1, NUM_PHYTO  !{
          if (phyto(n)%q_fe_2_n(i,j,k).lt.phyto(n)%fe_2_n_max) then
             phyto(n)%juptake_fe(i,j,k) = phyto(n)%P_C_max*cobalt%expkT(i,j,k)*phyto(n)%f_n(i,j,k)* &
                phyto(n)%felim(i,j,k)*cobalt%fe_2_n_upt_fac
          else 
             phyto(n)%juptake_fe(i,j,k) = 0.0
          endif
       enddo   !} n
    enddo; enddo ; enddo !} i,j,k
    !
    ! Silicate uptake
    !
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%nlg_diatoms(i,j,k)=phyto(LARGE)%f_n(i,j,k)*phyto(LARGE)%silim(i,j,k)
       cobalt%q_si_2_n_lg_diatoms(i,j,k)= cobalt%f_silg(i,j,k)/ &
             (cobalt%nlg_diatoms(i,j,k) + epsln)
       phyto(LARGE)%juptake_sio4(i,j,k) = &
             max(phyto(LARGE)%juptake_no3(i,j,k)+phyto(LARGE)%juptake_nh4(i,j,k),0.0)*phyto(LARGE)%silim(i,j,k)* &
             phyto(LARGE)%silim(i,j,k)*phyto(LARGE)%si_2_n_max 

       ! Note that this is si_2_n in large phytoplankton pool, not in diatoms themselves (q_si_2_n_lg_diatoms) 
       phyto(LARGE)%q_si_2_n(i,j,k) = cobalt%f_silg(i,j,k)/(phyto(LARGE)%f_n(i,j,k)+epsln)

    enddo; enddo ; enddo !} i,j,k
!
    call mpp_clock_end(id_clock_phyto_growth)
!
!-----------------------------------------------------------------------
! 2: Bacterial Growth and Uptake Calculations 
!-----------------------------------------------------------------------
!
    !
    ! calculate an effective maximum ldon uptake rate (at 0 deg. C) for bacteria
    ! from specified values of bact(1)%gge_max, bact(1)%mu_max and bact(1)%bresp
    !

    call mpp_clock_begin(id_clock_bacteria_growth)
    vmax_bact = (1.0/bact(1)%gge_max)*(bact(1)%mu_max + bact(1)%bresp)
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec   !{
       bact(1)%temp_lim(i,j,k) = exp(bact(1)%ktemp*Temp(i,j,k))
       bact(1)%ldonlim(i,j,k) = cobalt%f_ldon(i,j,k)/(bact(1)%k_ldon + cobalt%f_ldon(i,j,k))
       bact(1)%o2lim(i,j,k) = max(cobalt%f_o2(i,j,k),cobalt%o2_min)/  &
                              (cobalt%k_o2 + max(cobalt%f_o2(i,j,k),cobalt%o2_min))
       bact(1)%juptake_ldon(i,j,k) = vmax_bact*bact(1)%temp_lim(i,j,k)*bact(1)%ldonlim(i,j,k)* &
                                     bact(1)%o2lim(i,j,k)*bact(1)%f_n(i,j,k)
       bact_uptake_ratio = ( cobalt%f_ldop(i,j,k)/max(cobalt%f_ldon(i,j,k),epsln) )
       bact(1)%juptake_ldop(i,j,k) = bact(1)%juptake_ldon(i,j,k)*bact_uptake_ratio
       ! calculate bacteria production if N-limited, adjust down if P-limited
       bact(1)%jprod_n(i,j,k) = bact(1)%gge_max*bact(1)%juptake_ldon(i,j,k) - &
          bact(1)%f_n(i,j,k)/(cobalt%refuge_conc + bact(1)%f_n(i,j,k)) *      &
          bact(1)%temp_lim(i,j,k)*bact(1)%bresp*bact(1)%f_n(i,j,k)
       bact(1)%jprod_n(i,j,k) = min(bact(1)%jprod_n(i,j,k), &
                                    bact(1)%juptake_ldop(i,j,k)/bact(1)%q_p_2_n)
       ! remineralization of oragnic N to nh4 = difference between uptake and production
       bact(1)%jprod_nh4(i,j,k) = bact(1)%juptake_ldon(i,j,k) - max(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + bact(1)%jprod_nh4(i,j,k)

       if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
          ! aerobic remineralization, nh4 production, o2 respired
          cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) + bact(1)%jprod_nh4(i,j,k)*cobalt%o2_2_nh4
       else
          ! low o2 leads to water column denitrification. nh4 is created, but no o2 is used
          cobalt%jno3denit_wc(i,j,k) = cobalt%jno3denit_wc(i,j,k) + & 
                                       bact(1)%jprod_nh4(i,j,k)*cobalt%n_2_n_denit
          ! uncomment for "no mass change" test
          ! cobalt%jno3denit_wc(i,j,k) = 0.0
       endif  !}

       ! produce phosphate at the same rate regardless of whether aerobic/anaerobic
       bact(1)%jprod_po4(i,j,k) = bact(1)%juptake_ldop(i,j,k) - max(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + bact(1)%jprod_po4(i,j,k)

    enddo; enddo ; enddo !} i,j,k
!
    call mpp_clock_end(id_clock_bacteria_growth)
!
!-----------------------------------------------------------------------
! 3: Plankton foodweb dynamics
!-----------------------------------------------------------------------
!
    !
    ! 3.1 Plankton foodweb dynamics: consumption by zooplankton and higher predators
    !

    call mpp_clock_begin(id_clock_zooplankton_calculations)

    !
    ! Set-up local matrices for calculating zooplankton ingestion of
    ! multiple prey types.  The rows are consumers (i.e., NUM_ZOO zooplankton
    ! groups), the columns are food sources (i.e., NUM_PREY potential food sources)
    !
    ! ipa_matrix = the innate prey availability matrix
    ! pa_matrix = prey availability matrix after accounting for switching 
    ! ingest_matrix = ingestion matrix
    ! tot_prey = total prey available to predator m 
    !
    ! The definition of predator-prey matrices is intended to allow for
    ! efficient experimentation with predator-prey interconnections.
    ! However, we are still working to reduce the runtime required to
    ! include this feature.  The matrix structures are thus included,
    ! but the standard COBALT interactions have been hard-coded.  This
    ! makes the code faster, but adding consumer-resource linkages 
    ! requires new code rather than just changing parameters
    !
    ! Note that the primary ingestion calculations allow for variable 
    ! stoichiometry.
    !

    do m = 1,NUM_ZOO !{
       ipa_matrix(m,1) = zoo(m)%ipa_diaz
       ipa_matrix(m,2) = zoo(m)%ipa_lgp
       ipa_matrix(m,3) = zoo(m)%ipa_smp
       ipa_matrix(m,4) = zoo(m)%ipa_bact
       ipa_matrix(m,5) = zoo(m)%ipa_smz
       ipa_matrix(m,6) = zoo(m)%ipa_mdz
       ipa_matrix(m,7) = zoo(m)%ipa_lgz
       ipa_matrix(m,8) = zoo(m)%ipa_det
       tot_prey(m) = 0.0
       do n = 1,NUM_PREY !{
           pa_matrix(m,n) = 0.0
           ingest_matrix(m,n) = 0.0
       enddo !} n
    enddo !} m

    !
    ! Set-up local matrices for calculating higher predator ingestion
    ! of multiple prey types
    !

    hp_ipa_vec(1) = cobalt%hp_ipa_diaz
    hp_ipa_vec(2) = cobalt%hp_ipa_lgp
    hp_ipa_vec(3) = cobalt%hp_ipa_smp
    hp_ipa_vec(4) = cobalt%hp_ipa_bact
    hp_ipa_vec(5) = cobalt%hp_ipa_smz
    hp_ipa_vec(6) = cobalt%hp_ipa_mdz
    hp_ipa_vec(7) = cobalt%hp_ipa_lgz
    hp_ipa_vec(8) = cobalt%hp_ipa_det
    tot_prey_hp = 0.0
    do n = 1,NUM_PREY  !{  
       hp_pa_vec(n) = 0.0                  
       hp_ingest_vec(n) = 0.0              
    enddo !} n

    ! 
    ! Set all static stoichiometric ratios outside k,j,i loop
    !

    prey_p2n_vec(1) = phyto(DIAZO)%p_2_n_static
    prey_p2n_vec(2) = phyto(LARGE)%p_2_n_static
    prey_p2n_vec(3) = phyto(SMALL)%p_2_n_static
    prey_p2n_vec(4) = bact(1)%q_p_2_n
    prey_p2n_vec(5) = zoo(1)%q_p_2_n
    prey_p2n_vec(6) = zoo(2)%q_p_2_n
    prey_p2n_vec(7) = zoo(3)%q_p_2_n

    prey_fe2n_vec(4) = 0.0
    prey_fe2n_vec(5) = 0.0
    prey_fe2n_vec(6) = 0.0
    prey_fe2n_vec(7) = 0.0

    prey_si2n_vec(1) = 0.0
    prey_si2n_vec(3) = 0.0
    prey_si2n_vec(4) = 0.0
    prey_si2n_vec(5) = 0.0
    prey_si2n_vec(6) = 0.0
    prey_si2n_vec(7) = 0.0

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{

       !
       ! 3.1.1: Calculate zooplankton ingestion fluxes
       !

       ! Calculate the temperature and oxygen limitations, no ingestion
       ! in low o2 environments
       do m = 1,3  !{
          zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*Temp(i,j,k))
          zoo(m)%o2lim(i,j,k) = max((cobalt%f_o2(i,j,k) - cobalt%o2_min),0.0)/ & 
                                (cobalt%k_o2 + max(cobalt%f_o2(i,j,k)-cobalt%o2_min,0.0))
       enddo  !}  m
       cobalt%hp_temp_lim(i,j,k) = exp(cobalt%ktemp_hp*Temp(i,j,k))
       cobalt%hp_o2lim(i,j,k) = max((cobalt%f_o2(i,j,k) - cobalt%o2_min),0.0)/ &
                                (cobalt%k_o2 + max(cobalt%f_o2(i,j,k)-cobalt%o2_min,0.0))

       ! Prey vectors for ingestion and loss calculations 
       ! (note: ordering of phytoplankton must be consistent with
       !  DIAZO, LARGE, SMALL ordering inherited from TOPAZ)
       !
       prey_vec(1) = max(phyto(DIAZO)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(2) = max(phyto(LARGE)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(3) = max(phyto(SMALL)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(4) = max(bact(1)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(5) = max(zoo(1)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(6) = max(zoo(2)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(7) = max(zoo(3)%f_n(i,j,k) - cobalt%refuge_conc,0.0)
       prey_vec(8) = max(cobalt%f_ndet(i,j,k) - cobalt%refuge_conc,0.0)
       ! 
       ! Set dynamic stoichiometric rations inside k,j,i loop
       prey_p2n_vec(8) = cobalt%f_pdet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_fe2n_vec(1) = phyto(DIAZO)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(2) = phyto(LARGE)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(3) = phyto(SMALL)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(8) = cobalt%f_fedet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_si2n_vec(2) = phyto(LARGE)%q_si_2_n(i,j,k)
       prey_si2n_vec(8) = cobalt%f_sidet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)

       !
       ! Calculate zooplankton ingestion
       !
       ! Small zooplankton (m = 1) consuming small phytoplankton (3) and
       ! bacteria (4).  sw_fac_denom is the denominator of the abundance-
       ! based switching factor, tot_prey is the total available prey 
       ! after accounting for switching.
       !
       ! CAS: speed up code by using integer "switch" terms and sqrt?

       m = 1 
       sw_fac_denom = (ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch + &
                      (ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch
       pa_matrix(m,3) = ipa_matrix(m,3)* &
                        ( (ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,4) = ipa_matrix(m,4)* &
                        ( (ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,3)*prey_vec(3) + pa_matrix(m,4)*prey_vec(4)
       ingest_matrix(m,3) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                 pa_matrix(m,3)*prey_vec(3)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,4) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                 pa_matrix(m,4)*prey_vec(4)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,3) + ingest_matrix(m,4)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,3)*prey_p2n_vec(3) + &
                                 ingest_matrix(m,4)*prey_p2n_vec(4)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,3)*prey_fe2n_vec(3)

       !
       ! Medium zooplankton (m = 2) consuming diazotrophs (1), large
       ! phytoplankton (2), and small zooplankton (5) 
       !

       m = 2 
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch + &
                      (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch + &
                      (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)* &
                        ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)* & 
                        ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,5) = ipa_matrix(m,5)* & 
                        ( (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,1)*prey_vec(1) + pa_matrix(m,2)*prey_vec(2) + &
                     pa_matrix(m,5)*prey_vec(5)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,1)*prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,2)*prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,5) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,5)*prey_vec(5)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1) + ingest_matrix(m,2) + &
                                 ingest_matrix(m,5)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_p2n_vec(2) + &
                                 ingest_matrix(m,5)*prey_p2n_vec(5)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       !
       ! Large zooplankton (m = 3) consuming diazotrophs (2), large phytoplankton (2)
       ! and medium zooplankton (6)
       !

       m = 3
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch + &
                      (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch + &
                      (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)* &
                        ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)* &
                        ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,6) = ipa_matrix(m,6)* &
                        ( (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch / &
                          (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) = pa_matrix(m,1)*prey_vec(1) + pa_matrix(m,2)*prey_vec(2) + &
                     pa_matrix(m,6)*prey_vec(6)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,1)*prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,2)*prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,6) = zoo(m)%temp_lim(i,j,k)*zoo(m)%o2lim(i,j,k)*zoo(m)%imax* &
                     pa_matrix(m,6)*prey_vec(6)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1) + ingest_matrix(m,2) + &
                                 ingest_matrix(m,6)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_p2n_vec(2) + &
                                 ingest_matrix(m,6)*prey_p2n_vec(6)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) + &
                                 ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       cobalt%total_filter_feeding(i,j,k) = ingest_matrix(2,1) + ingest_matrix(2,2) + &
          ingest_matrix(2,3) + ingest_matrix(3,1) + ingest_matrix(3,2) + & 
          ingest_matrix(3,3)

       !
       ! Calculate losses to zooplankton
       !

       do n = 1,NUM_PHYTO
          phyto(n)%jzloss_n(i,j,k) = 0.0
       enddo

       do m = 1,NUM_ZOO !{
          phyto(DIAZO)%jzloss_n(i,j,k) = phyto(DIAZO)%jzloss_n(i,j,k) + ingest_matrix(m,DIAZO)
          phyto(LARGE)%jzloss_n(i,j,k) = phyto(LARGE)%jzloss_n(i,j,k) + ingest_matrix(m,LARGE)
          phyto(SMALL)%jzloss_n(i,j,k) = phyto(SMALL)%jzloss_n(i,j,k) + ingest_matrix(m,SMALL)
       enddo !} m

       do n = 1,NUM_PHYTO !{
          phyto(n)%jzloss_p(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_p2n_vec(n)
          phyto(n)%jzloss_fe(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_fe2n_vec(n)
          phyto(n)%jzloss_sio2(i,j,k) = phyto(n)%jzloss_n(i,j,k)*prey_si2n_vec(n)  
       enddo !} n

       !
       ! losses of bacteria to zooplankton 
       !

       bact(1)%jzloss_n(i,j,k) = 0.0
       do m = 1,NUM_ZOO !{
          bact(1)%jzloss_n(i,j,k) = bact(1)%jzloss_n(i,j,k) + ingest_matrix(m,4)
       enddo !} m
       bact(1)%jzloss_p(i,j,k) = bact(1)%jzloss_n(i,j,k)*prey_p2n_vec(4)

       !
       ! losses of zooplankton to zooplankton
       !

       do n = 1,NUM_ZOO !{
          zoo(n)%jzloss_n(i,j,k) = 0.0
          do m = 1,NUM_ZOO !{
             zoo(n)%jzloss_n(i,j,k) = zoo(n)%jzloss_n(i,j,k) + ingest_matrix(m,NUM_PHYTO+1+n)
          enddo !} m
          zoo(n)%jzloss_p(i,j,k) = zoo(n)%jzloss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+n)
       enddo !} n

       !
       ! 3.1.2 Calculate ingestion by higher predators
       !

       ! The higher-predator ingestion calculations mirror those used for zooplankton
       !
       sw_fac_denom = (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp + &
                      (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp
       hp_pa_vec(6) = hp_ipa_vec(6)* &
                      ( (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp / &
                        (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       hp_pa_vec(7) = hp_ipa_vec(7)* &
                      ( (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp / &
                        (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       tot_prey_hp = hp_pa_vec(6)*prey_vec(6) + hp_pa_vec(7)*prey_vec(7)
       hp_ingest_vec(6) = cobalt%hp_temp_lim(i,j,k)*cobalt%hp_o2lim(i,j,k)*cobalt%imax_hp* &
                          hp_pa_vec(6)*prey_vec(6)*tot_prey_hp**(cobalt%coef_hp-1.0)/ &
                            (cobalt%ki_hp+tot_prey_hp)
       hp_ingest_vec(7) = cobalt%hp_temp_lim(i,j,k)*cobalt%hp_o2lim(i,j,k)*cobalt%imax_hp* &
                          hp_pa_vec(7)*prey_vec(7)*tot_prey_hp**(cobalt%coef_hp-1.0)/ &
                            (cobalt%ki_hp+tot_prey_hp)
       cobalt%hp_jingest_n(i,j,k) = hp_ingest_vec(6) + hp_ingest_vec(7)
       cobalt%hp_jingest_p(i,j,k) = hp_ingest_vec(6)*prey_p2n_vec(6) + &
                                    hp_ingest_vec(7)*prey_p2n_vec(7)
       !
       ! Calculate losses to higher predators
       !

       do n = 1,NUM_ZOO !{
         zoo(n)%jhploss_n(i,j,k) = hp_ingest_vec(NUM_PHYTO+1+n)
         zoo(n)%jhploss_p(i,j,k) = zoo(n)%jhploss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+n)
       enddo !} n

    enddo; enddo; enddo  !} i,j,k
    call mpp_clock_end(id_clock_zooplankton_calculations)

    !
    ! 3.2: Plankton foodweb dynamics: Other mortality and loss terms
    !

    call mpp_clock_begin(id_clock_other_losses)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{

       !  
       ! 3.2.1 Calculate losses of phytoplankton to aggregation 
       !

       do n = 1,NUM_PHYTO !{
            growth_ratio = min(phyto(n)%f_mu_mem(i,j,k)/ &
                           (phyto(n)%frac_mu_agg*phyto(n)%P_C_max*cobalt%expkT(i,j,k)),1.0)
            phyto(n)%agg_lim(i,j,k) = (1.0-growth_ratio)**2
            phyto(n)%jaggloss_n(i,j,k) = phyto(n)%agg_lim(i,j,k)*phyto(n)%agg*phyto(n)%f_n(i,j,k)**2.0 
            phyto(n)%jaggloss_p(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k)
            phyto(n)%jaggloss_fe(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k)
            phyto(n)%jaggloss_sio2(i,j,k) = phyto(n)%jaggloss_n(i,j,k)*phyto(n)%q_si_2_n(i,j,k)
       enddo !} n

       !
       ! 3.2.2 Calculate phytoplankton and bacterial losses to viruses
       !

       do n = 1,NUM_PHYTO !{
          phyto(n)%jvirloss_n(i,j,k) = bact(1)%temp_lim(i,j,k)*phyto(n)%vir*phyto(n)%f_n(i,j,k)**2.0 
          phyto(n)%jvirloss_p(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_p_2_n(i,j,k)
          phyto(n)%jvirloss_fe(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_fe_2_n(i,j,k)
          phyto(n)%jvirloss_sio2(i,j,k) = phyto(n)%jvirloss_n(i,j,k)*phyto(n)%q_si_2_n(i,j,k)
       enddo !} n

       bact(1)%jvirloss_n(i,j,k) = bact(1)%temp_lim(i,j,k)*bact(1)%vir*bact(1)%f_n(i,j,k)**2.0
       bact(1)%jvirloss_p(i,j,k) = bact(1)%jvirloss_n(i,j,k)*bact(1)%q_p_2_n

       !
       ! 3.2.3 Calculate losses to exudation
       !

       n = DIAZO
       phyto(n)%jexuloss_n(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_no3(i,j,k)+ &
                                    phyto(n)%juptake_nh4(i,j,k)+phyto(n)%juptake_n2(i,j,k),0.0)
       phyto(n)%jexuloss_p(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_po4(i,j,k),0.0)
       phyto(n)%jexuloss_fe(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_fe(i,j,k),0.0)
       do n = 2,NUM_PHYTO !{
          phyto(n)%jexuloss_n(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_no3(i,j,k)+phyto(n)%juptake_nh4(i,j,k),0.0)
          phyto(n)%jexuloss_p(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_po4(i,j,k),0.0)
          phyto(n)%jexuloss_fe(i,j,k) = phyto(n)%exu*max(phyto(n)%juptake_fe(i,j,k),0.0)
       enddo

    enddo; enddo; enddo  !} i,j,k
    call mpp_clock_end(id_clock_other_losses)

    !
    ! 3.3: Plankton foodweb dynamics: Production calculations
    !

    call mpp_clock_begin(id_clock_production_loop)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

       !
       ! 3.3.1: Calculate the production of detritus and dissolved organic material
       !
       !
       ! Production of detritus and dissolved organic material from zooplankton egestion 
       !   

       do m = 1,NUM_ZOO
           zoo(m)%jprod_ndet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_pdet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_sldon(i,j,k) = zoo(m)%phi_sldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_ldon(i,j,k) = zoo(m)%phi_ldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_srdon(i,j,k) = zoo(m)%phi_srdon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_sldop(i,j,k) = zoo(m)%phi_sldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_ldop(i,j,k) = zoo(m)%phi_ldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_srdop(i,j,k) = zoo(m)%phi_srdop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_fedet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_fe(i,j,k)
           zoo(m)%jprod_sidet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_sio2(i,j,k)


           ! augment cumulative production with zooplankton terms
           cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + zoo(m)%jprod_ndet(i,j,k)
           cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + zoo(m)%jprod_pdet(i,j,k)
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + zoo(m)%jprod_sldon(i,j,k)
           cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + zoo(m)%jprod_ldon(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + zoo(m)%jprod_srdon(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + zoo(m)%jprod_sldop(i,j,k)
           cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + zoo(m)%jprod_ldop(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + zoo(m)%jprod_srdop(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + zoo(m)%jprod_fedet(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + zoo(m)%jprod_sidet(i,j,k)
       enddo !} m

       !
       ! Production of detritus and dissolved organic material from higher predator egestion 
       ! (did not track individual terms, just add to cumulative total)
       !

       cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_sio2(i,j,k)
       
       !
       ! Sources from phytoplankton aggregation
       !

       do m = 1,NUM_PHYTO
           cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) + phyto(m)%jaggloss_n(i,j,k)
           cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) + phyto(m)%jaggloss_p(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + phyto(m)%jaggloss_fe(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + phyto(m)%jaggloss_sio2(i,j,k)
       enddo !} m

       !
       ! Sources from viral lysis of phytoplankton (0 in default formulation) and exudation
       !

       do m = 1,NUM_PHYTO
           cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + cobalt%lysis_phi_ldon*phyto(m)%jvirloss_n(i,j,k) + &
                                      phyto(m)%jexuloss_n(i,j,k) 
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%lysis_phi_sldon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%lysis_phi_srdon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + cobalt%lysis_phi_ldop*phyto(m)%jvirloss_p(i,j,k) + &
                                      phyto(m)%jexuloss_p(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%lysis_phi_sldop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%lysis_phi_srdop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_fed(i,j,k)   = cobalt%jprod_fed(i,j,k)   + phyto(m)%jvirloss_fe(i,j,k) + &
                                       phyto(m)%jexuloss_fe(i,j,k) 
           cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + phyto(m)%jvirloss_sio2(i,j,k)
       enddo !} m


       !
       ! Sources of dissolved organic material from viral lysis due to bacteria 
       !

       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + cobalt%lysis_phi_ldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%lysis_phi_sldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%lysis_phi_srdon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + cobalt%lysis_phi_ldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%lysis_phi_sldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%lysis_phi_srdop*bact(1)%jvirloss_p(i,j,k)

       !
       ! Sources of dissolved organic material from bacterial mortality (metabolic costs higher than food uptake).
       ! These conditions are assumed to lead to a lysis-like redistribution of bacteria organic matter.
       !

       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) - cobalt%lysis_phi_ldon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) - cobalt%lysis_phi_sldon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) - cobalt%lysis_phi_srdon* &
                                  min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) - cobalt%lysis_phi_ldop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) - cobalt%lysis_phi_sldop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - cobalt%lysis_phi_srdop* &
                                  min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       !
       ! 3.3.2: Zooplankton production and excretion calculations
       !

       do m = 1,NUM_ZOO
          assim_eff = 1.0-zoo(m)%phi_det-zoo(m)%phi_ldon-zoo(m)%phi_sldon-zoo(m)%phi_srdon
          zoo(m)%jprod_n(i,j,k) = zoo(m)%gge_max*zoo(m)%jingest_n(i,j,k) - &
                                     zoo(m)%f_n(i,j,k)/(cobalt%refuge_conc + zoo(m)%f_n(i,j,k))* &
                                     zoo(m)%temp_lim(i,j,k)*zoo(m)%bresp*zoo(m)%f_n(i,j,k)
          zoo(m)%jprod_n(i,j,k) = min(zoo(m)%jprod_n(i,j,k), &
                                      assim_eff*zoo(m)%jingest_p(i,j,k)/zoo(m)%q_p_2_n)

          !
          ! Ingested material that does not go to zooplankton production, detrital production
          ! or production of dissolved organic material is excreted as nh4 or po4.  If production
          ! is negative, zooplankton are lost to large detritus 
          !
          if (zoo(m)%jprod_n(i,j,k) .gt. 0.0) then 
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k) - zoo(m)%jprod_ndet(i,j,k) -  &
                                        zoo(m)%jprod_n(i,j,k) - zoo(m)%jprod_ldon(i,j,k) - &
                                        zoo(m)%jprod_sldon(i,j,k) - zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k) - zoo(m)%jprod_pdet(i,j,k) - & 
                                        zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n - zoo(m)%jprod_ldop(i,j,k) -  &
                                        zoo(m)%jprod_sldop(i,j,k) - zoo(m)%jprod_srdop(i,j,k)
          else
             ! None of the ingestion material goes to zooplankton production
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k) - zoo(m)%jprod_ndet(i,j,k) - & 
                                        zoo(m)%jprod_ldon(i,j,k) - zoo(m)%jprod_sldon(i,j,k) - & 
                                        zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k) - zoo(m)%jprod_pdet(i,j,k) - &
                                        zoo(m)%jprod_ldop(i,j,k) - zoo(m)%jprod_sldop(i,j,k) - & 
                                        zoo(m)%jprod_srdop(i,j,k)

             ! The negative production (i.e., mortality) is lost to large detritus. Update values
             ! for zooplankton and for total.

             zoo(m)%jprod_ndet(i,j,k) = zoo(m)%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             zoo(m)%jprod_pdet(i,j,k) = zoo(m)%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n
             cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n 
          endif

          ! cumulative production of inorganic nutrients 
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + zoo(m)%jprod_nh4(i,j,k)
          cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + zoo(m)%jprod_po4(i,j,k)
          cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) + zoo(m)%jprod_nh4(i,j,k)*cobalt%o2_2_nh4

          !
          ! Any ingested iron that is not allocated to detritus is routed back to the
          ! dissolved pool.       
          !
          zoo(m)%jprod_fed(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_fe(i,j,k)
          cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + zoo(m)%jprod_fed(i,j,k)
          !
          ! Any ingested opal that is not allocated to detritus is assumed to undergo
          ! rapid dissolution to dissolved silica
          !
          zoo(m)%jprod_sio4(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_sio2(i,j,k)
          cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + zoo(m)%jprod_sio4(i,j,k)
 
       enddo !} m

       !
       ! Excretion by higher predators
       !
       cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_sio2(i,j,k)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_p(i,j,k)
       cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) + (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_n(i,j,k)*cobalt%o2_2_nh4

    enddo; enddo ; enddo !} i,j,k
    call mpp_clock_end(id_clock_production_loop)

    call mpp_clock_begin(id_clock_ballast_loops)
!
!
!------------------------------------------------------------------------------------
! 4: Production of calcium carbonate (Calcite and Aragonite) and lithogenic material
!------------------------------------------------------------------------------------
!
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

    !
    ! 4.1: Calculate aragonite and calcite saturation states
    !
       if (trim(co2_calc) == "ocmip2") then
         TK = Temp(i,j,k) + 273.15
         PRESS = 0.1016 * cobalt%zt(i,j,k) + 1.013
         PKSPA = 171.945 + 0.077993 * TK - 2903.293 / TK - 71.595 * log10(TK) - (-0.068393 + 1.7276e-3 * &
            TK + 88.135 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.10018 * max(epsln, Salt(i,j,k)) -      &
            5.9415e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 2.8 - 0.5304 * Temp(i,j,k)) *   &
            (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
            (PRESS - 1.013) / (382.92 * TK)
         cobalt%co3_sol_arag(i,j,k) = 10**(-PKSPA) / (2.937d-4 * max(5.0, Salt(i,j,k)))
         cobalt%omega_arag(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_arag(i,j,k)
         PKSPC = 171.9065 + 0.077993 * TK - 2839.319 / TK - 71.595 * log10(TK) - (-0.77712 + 2.8426e-3 * &
            TK + 178.34 / TK) * sqrt(max(epsln, Salt(i,j,k))) + 0.07711 * max(epsln, Salt(i,j,k)) -      &
            4.1249e-3 * max(epsln, Salt(i,j,k))**(1.5) - 0.02 - (48.76 - 0.5304 * Temp(i,j,k)) *         &
            (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 - 0.3692 * Temp(i,j,k))) * (PRESS - 1.013) *&
            (PRESS - 1.013) / (382.92 * TK)
         cobalt%co3_sol_calc(i,j,k) = 10**(-PKSPC) / (2.937d-4 * max(5.0, Salt(i,j,k)))
         cobalt%omega_calc(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_calc(i,j,k)
      else if (trim(co2_calc) == "mocsy") then
         cobalt%omega_arag(i,j,k) = cobalt%omegaa(i,j,k)  ! from Mocsy
         cobalt%omega_calc(i,j,k) = cobalt%omegac(i,j,k)  ! from Mocsy
         cobalt%co3_sol_arag(i,j,k) = cobalt%f_co3_ion(i,j,k) / max(cobalt%omega_arag(i,j,k),epsln)
         cobalt%co3_sol_calc(i,j,k) = cobalt%f_co3_ion(i,j,k) / max(cobalt%omega_calc(i,j,k),epsln)
      else
        call mpp_error(FATAL,"Unable to compute aragonite and calcite saturation states")
      endif

    enddo; enddo ; enddo !} i,j,k

    !
    ! 4.2: Calculate the production rate of aragonite and calcite detritus 
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
        cobalt%jprod_cadet_arag(i,j,k) = (zoo(2)%jzloss_n(i,j,k)*zoo(3)%phi_det + &
                       (zoo(2)%jhploss_n(i,j,k) + zoo(3)%jhploss_n(i,j,k))*cobalt%hp_phi_det)* &
                       cobalt%ca_2_n_arag*min(cobalt%caco3_sat_max, max(0.0,cobalt%omega_arag(i,j,k) - 1.0)) + epsln 
        cobalt%jprod_cadet_calc(i,j,k) = (zoo(1)%jzloss_n(i,j,k)*zoo(2)%phi_det + & 
                       phyto(SMALL)%jzloss_n(i,j,k)*zoo(1)%phi_det + phyto(LARGE)%jzloss_n(i,j,k)*zoo(3)%phi_det + &  
                       phyto(SMALL)%jaggloss_n(i,j,k) + phyto(LARGE)%jaggloss_n(i,j,k))*cobalt%ca_2_n_calc* &
                       min(cobalt%caco3_sat_max, max(0.0, cobalt%omega_calc(i,j,k) - 1.0)) + epsln
    enddo; enddo ; enddo !} i,j,k

    !
    ! 4.3: Lithogenic detritus production (repackaged from f_lith during filter feeding)
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%jprod_lithdet(i,j,k)=( cobalt%total_filter_feeding(i,j,k)/ &
                                   ( phyto(LARGE)%f_n(i,j,k) + phyto(DIAZO)%f_n(i,j,k) + epsln ) * &  
                                    cobalt%phi_lith + cobalt%k_lith ) * cobalt%f_lith(i,j,k)
    enddo; enddo ; enddo !} i,j,k

!
!---------------------------------------------------------------------------------------------------------
! 5: Detrital dissolution and remineralization calculation
!---------------------------------------------------------------------------------------------------------
!

    !
    ! 5.1: Dissolution of aragonite, calcite and opal detrital particles
    !

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%jdiss_cadet_arag(i,j,k) = cobalt%gamma_cadet_arag * & 
         max(0.0, 1.0 - cobalt%omega_arag(i,j,k)) * cobalt%f_cadet_arag(i,j,k)
       cobalt%jdiss_cadet_calc(i,j,k) = cobalt%gamma_cadet_calc * &
         max(0.0, 1.0 - cobalt%omega_calc(i,j,k)) * cobalt%f_cadet_calc(i,j,k)
       !cobalt%jdiss_sidet(i,j,k) = cobalt%gamma_sidet * cobalt%f_sidet(i,j,k)
       cobalt%jdiss_sidet(i,j,k) = cobalt%gamma_sidet * exp(cobalt%kappa_sidet * &
          Temp(i,j,k)) * cobalt%f_sidet(i,j,k)
       cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + cobalt%jdiss_sidet(i,j,k)
    enddo; enddo ; enddo !} i,j,k

    !
    ! 5.2: Remineralization of nitrogen, phosphorous and iron detritus accounting for oxygen 
    !      and mineral protection 
    !

    ! Calculate the depth for scaling of remineralization near the surface
    allocate(z_remin_ramp(isd:ied,jsd:jed,1:nk)); z_remin_ramp = 0.0
    z_remin_ramp(:,:,1) = dzt(:,:,1)
    do k = 2,nk !{
      z_remin_ramp(:,:,k) = z_remin_ramp(:,:,k-1) + dzt(:,:,k)
    enddo !}k 
!
!---------------------------------------------------------------------------------------------------------

    do k=1,nk ; do j=jsc,jec ; do i=isc,iec  !{
       cobalt%expkreminT(i,j,k) = exp(cobalt%kappa_remin * Temp(i,j,k))
       !cobalt%jno3denit_wc(i,j,k) = 0.0
       !
       !   Under oxic conditions
       !
       if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet * cobalt%expkreminT(i,j,k) * &
               z_remin_ramp(i,j,k)/(z_remin_ramp(i,j,k) + cobalt%remin_ramp_scale) * cobalt%f_o2(i,j,k) / & 
               ( cobalt%k_o2 + cobalt%f_o2(i,j,k) )*max( 0.0, cobalt%f_ndet(i,j,k) - &
               cobalt%rpcaco3*(cobalt%f_cadet_arag(i,j,k) + cobalt%f_cadet_calc(i,j,k)) - & 
               cobalt%rplith*cobalt%f_lithdet(i,j,k) - cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + cobalt%jremin_ndet(i,j,k)
          cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) + cobalt%jremin_ndet(i,j,k)*cobalt%o2_2_nh4
       !
       ! Under sub-oxic conditions
       !
       else !}{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet * cobalt%o2_min / &
               (cobalt%k_o2 + cobalt%o2_min)* &
               cobalt%f_no3(i,j,k) / (phyto(SMALL)%k_no3 + cobalt%f_no3(i,j,k))* &
               max(0.0, cobalt%f_ndet(i,j,k) - &
               cobalt%rpcaco3*(cobalt%f_cadet_arag(i,j,k) + cobalt%f_cadet_calc(i,j,k)) - &
               cobalt%rplith*cobalt%f_lithdet(i,j,k) - cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
          cobalt%jno3denit_wc(i,j,k) = cobalt%jno3denit_wc(i,j,k) + cobalt%jremin_ndet(i,j,k) * cobalt%n_2_n_denit
          ! uncomment for "no mass change" test
          ! cobalt%jno3denit_wc(i,j,k) = 0.0
          ! using TOPAZ stoichiometry, denitrification produces ammonia.
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + cobalt%jremin_ndet(i,j,k)
       endif !}
       !
       ! P and Fe assumed to be protected similarly to N
       !
       cobalt%jremin_pdet(i,j,k) = cobalt%jremin_ndet(i,j,k)/(cobalt%f_ndet(i,j,k) + epsln)* &
         cobalt%f_pdet(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + cobalt%jremin_pdet(i,j,k)
       cobalt%jremin_fedet(i,j,k) = cobalt%jremin_ndet(i,j,k) / (cobalt%f_ndet(i,j,k) + epsln) * &
         cobalt%remin_eff_fedet*cobalt%f_fedet(i,j,k)
       cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + cobalt%jremin_fedet(i,j,k) 
    enddo; enddo; enddo  !} i,j,k

    deallocate(z_remin_ramp)
!
!
!--------------------------------------------------------------------------------------------
! 6: Nitrification
!--------------------------------------------------------------------------------------------
!

       !
       !  Nitrification
       !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%jnitrif(i,j,k) = 0.0
       if (scheme_nitrif .eq. 2 .or. scheme_nitrif .eq. 3) then             
          if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min_nit) then  !{
             cobalt%jnitrif(i,j,k) = cobalt%gamma_nitrif * &
                  cobalt%f_nh3(i,j,k)/(cobalt%f_nh3(i,j,k)+cobalt%k_nh3_nitrif) *  &
                  (1.-cobalt%f_irr_mem(i,j,k)/(cobalt%irr_inhibit+cobalt%f_irr_mem(i,j,k))) * &
                  cobalt%f_o2(i,j,k)/(cobalt%k_o2_nit+cobalt%f_o2(i,j,k)) * cobalt%f_nh4(i,j,k)**2
             
             if (scheme_nitrif .eq. 3) then
                cobalt%jnitrif(i,j,k) = cobalt%jnitrif(i,j,k)*cobalt%expkT(i,j,k)
             end if
          end if
       elseif (scheme_nitrif .eq. 1) then
          if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
             cobalt%jnitrif(i,j,k) = cobalt%gamma_nitrif * cobalt%expkT(i,j,k) * cobalt%f_nh4(i,j,k) * &
                  phyto(SMALL)%nh4lim(i,j,k) * (1.0 - cobalt%f_irr_mem(i,j,k) / &
                  (cobalt%irr_inhibit + cobalt%f_irr_mem(i,j,k))) * cobalt%f_o2(i,j,k) / &
                  ( cobalt%k_o2 + cobalt%f_o2(i,j,k) )
          end if
       end if
       cobalt%jo2resp_wc(i,j,k) = cobalt%jo2resp_wc(i,j,k) + cobalt%jnitrif(i,j,k)*cobalt%o2_2_nitrif
    enddo; enddo; enddo  !} i,j,k

!
!
!--------------------------------------------------------------------------------------------
! 7: Iron
!--------------------------------------------------------------------------------------------
!

    ! 
    ! Iron scavenging and coastal sources
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%kfe_eq_lig(i,j,k) = min(cobalt%kfe_eq_lig_ll, 10.0**( log10(cobalt%kfe_eq_lig_hl) + &
          max(0.0,log10(cobalt%io_fescav/max(epsln,cobalt%irr_inst(i,j,k)))) ) ) 

       cobalt%ligand(i,j,k) = cobalt%felig_bkg + cobalt%felig_2_don*(cobalt%f_sldon(i,j,k) + &
            cobalt%f_srdon(i,j,k) + cobalt%f_ldon(i,j,k))
       feprime_temp = 1.0 + cobalt%kfe_eq_lig(i,j,k) * (cobalt%ligand(i,j,k) - cobalt%f_fed(i,j,k))
       cobalt%feprime(i,j,k) = (-feprime_temp + (feprime_temp * feprime_temp + 4.0 * cobalt%kfe_eq_lig(i,j,k) * &
            cobalt%f_fed(i,j,k))**(0.5)) / (2.0 * max(epsln,cobalt%kfe_eq_lig(i,j,k)))

       fe_salt = 19.922*Salt(i,j,k)/(1000.0 - 1.005*Salt(i,j,k))
       cobalt%fe_sol(i,j,k) = 10**(-10.53 + 322.5/(Temp(i,j,k)+273.15) - 2.524*sqrt(fe_salt) + &
                              2.921*fe_salt)

       !
       ! Iron adsorption to detrital particles
       !
       if (cobalt%feprime(i,j,k).lt.cobalt%fe_sol(i,j,k)) then
         cobalt%jfe_ads(i,j,k) = cobalt%alpha_fescav*cobalt%feprime(i,j,k) + &
                                 cobalt%beta_fescav*cobalt%feprime(i,j,k)*cobalt%f_ndet(i,j,k)
       else
         cobalt%jfe_ads(i,j,k) = 10.0*(cobalt%alpha_fescav*cobalt%feprime(i,j,k) + &
                                 cobalt%beta_fescav*cobalt%feprime(i,j,k)*cobalt%f_ndet(i,j,k))
       endif
       ! uncomment if running "no mass change" test
       !cobalt%jfe_ads(i,j,k) = 0.0

    enddo; enddo; enddo  !} i,j,k

!
!-------------------------------------------------------------------------------------------------
! 8: Sedimentary/coastal fluxes/transformations
!-------------------------------------------------------------------------------------------------
!

    !
    ! Coastal iron input (default is 0)
    !
    !do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
    !   cobalt%jfe_coast(i,j,k) = cobalt%fe_coast * mask_coast(i,j) * grid_tmask(i,j,k) / &
    !        sqrt(grid_dat(i,j))
    !     ! uncomment if running "no mass change" test
    !     !cobalt%jfe_coast(i,j,k) = 0.0
    !enddo; enddo; enddo  !} i,j,k

    do j = jsc, jec; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{
          !
          ! Nitrogen flux from the sediments
          ! 
          if (cobalt%f_ndet_btf(i,j,1) .gt. 0.0) then !{
             ! fpoc_bottom in mmoles C m-2 day-1 for burial relationship
             fpoc_btm = cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n*sperd*1000.0
             !cobalt%frac_burial(i,j) = (0.013 + 0.53*fpoc_btm**2.0)/((7.0+fpoc_btm)**2.0)
             cobalt%frac_burial(i,j) = 0.013 + 0.53*fpoc_btm**2.0/((7.0+fpoc_btm)**2.0) * &
                  cobalt%zt(i,j,k) / (cobalt%z_burial + cobalt%zt(i,j,k))
             ! uncomment for "no mass change" test
             !cobalt%frac_burial(i,j) = 0.0
             cobalt%fndet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_ndet_btf(i,j,1)
             cobalt%fpdet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_pdet_btf(i,j,1)
             ! fpoc_bottom in micromoles C cm-2 day-1 for denitrification relationship, cap at 43
             ! to prevent anomalous extrapolation of the relationship
             log_fpoc_btm = log(min(43.0,0.1*fpoc_btm))
             cobalt%fno3denit_sed(i,j) = min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*r_dt,  &      
                  min((cobalt%f_ndet_btf(i,j,1)-cobalt%fndet_burial(i,j))*cobalt%n_2_n_denit, & 
                  10.0**(-0.9543+0.7662*log_fpoc_btm - 0.235*log_fpoc_btm**2.0)/(cobalt%c_2_n*sperd*100.0)* &
                  cobalt%n_2_n_denit*cobalt%f_no3(i,j,k)/(cobalt%k_no3_denit + cobalt%f_no3(i,j,k)))) * &
                  cobalt%zt(i,j,k) / (cobalt%z_burial + cobalt%zt(i,j,k))
             ! uncomment "no mass change" test 
             !cobalt%fno3denit_sed(i,j) = 0.0             
             if (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) then  !{
                cobalt%fnoxic_sed(i,j) = max(0.0, min(cobalt%f_o2(i,j,k)*cobalt%Rho_0*r_dt*(1.0/cobalt%o2_2_nh4), &
                                         cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j) - &
                                         cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit))
             else
                cobalt%fnoxic_sed(i,j) = 0.0
             endif !}
             cobalt%fno3denit_sed(i,j) = cobalt%fno3denit_sed(i,j) + &
                                         min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*r_dt-cobalt%fno3denit_sed(i,j), &
                                         (cobalt%f_ndet_btf(i,j,1)-cobalt%fnoxic_sed(i,j)-cobalt%fndet_burial(i,j) - &
                                         cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)*cobalt%n_2_n_denit)
             cobalt%fnfeso4red_sed(i,j) = max(0.0, cobalt%f_ndet_btf(i,j,1)-cobalt%fnoxic_sed(i,j)- &
                                          cobalt%fndet_burial(i,j)-cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)
          else
             cobalt%fnfeso4red_sed(i,j) = 0.0
             cobalt%fno3denit_sed(i,j) = 0.0
             cobalt%fnoxic_sed(i,j) = 0.0
          endif !}

          ! iron from sediment (Elrod) 
          !cobalt%ffe_sed(i,j) = cobalt%fe_2_n_sed * cobalt%f_ndet_btf(i,j,1)
          ! iron from sediment (Dale)
          cobalt%ffe_sed(i,j) = cobalt%ffe_sed_max * tanh( (cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n*sperd*1.0e3)/ &
                                max(cobalt%f_o2(i,j,k)*1.0e6,epsln) )

          cobalt%ffe_geotherm(i,j) = cobalt%ffe_geotherm_ratio*internal_heat(i,j)*4184.0/dt
          ! default for icebergs: 40 nanomoles fe dissolved per kg of icemelt
          ! sediments: Raiswell et al., 2008: 0.5 kg sed per m-3 of iceberg; 0.1% mean Fe, 5-10% soluble
          ! ~500 nanomoles Fe per kg-1 icemelt 
          cobalt%ffe_iceberg(i,j) = cobalt%ffe_iceberg_ratio*max(frunoff(i,j),0.0)
          cobalt%jprod_fed(i,j,1) = cobalt%jprod_fed(i,j,1) + cobalt%ffe_iceberg(i,j)/rho_dzt(i,j,1) 

          !
          ! Calcium carbonate flux and burial
          ! 2015/11/18 JGJ: fix from JPD to cap the absolute cased dissolution rate to 10 mmol m-2 d-1
          ! Calcite cycling in the sediments is based on the model of   Dunne et al., 2012.
          !
          ! phi_surfresp_cased = 0.14307   ! const for enhanced diss., surf sed respiration (dimensionless)
          ! phi_deepresp_cased = 4.1228    ! const for enhanced diss., deep sed respiration (dimensionless)
          ! alpha_cased = 2.7488 ! exponent controlling non-linearity of deep dissolution           
          ! beta_cased = -2.2185 ! exponent controlling non-linearity of effective thickness
          ! gamma_cased = 0.03607/spery   ! dissolution rate constant
          ! Co_cased = 8.1e3        ! moles CaCo3 m-3 for pure calcite sediment with porosity = 0.7
          !
          ! if cased_steady is true, burial is calculated from Dunne's eq. (2) assuming dcased/dt = 0.
          ! This ensures that all the calcite bottom flux is partitioned between burial and redissolution.     
          ! The steady state cased value of cased is calculated to reflect the changing bottom conditions.
          ! This influences the the partitioning of burial and redissolution over time, but there are
          ! no alkalinity changes/drifts associated with the long-term evolution of cased
          ! 
          ! If cased_steady is false, calcite is partitioned between dissolution, burial and evolving
          ! cased as described in Dunne et al. (2012).  The multi-century scale evolution of cased
          ! impacts alkalinity, but care must to ensure that cased starts in equilibrium with the 
          ! mean ocean state to avoid unrealistic drifts. 

          ! Enhanced dissolution by fast respiration near the sediment surface, proportional 
          ! to organic flux, moles Ca m-2 s-1, limited to a max 1/2 the instantaneous calcite flux
          cobalt%fcased_redis_surfresp(i,j)=min(0.5*cobalt%f_cadet_calc_btf(i,j,1), &
            cobalt%phi_surfresp_cased*cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n)
          ! Ca-specific dissolution coeficient, depends on calcite saturation state and is enhanced by
          ! respiration deep in the sediment (s-1), non-linearity controlled by alpha_cased
          cobalt%cased_redis_coef(i,j) = cobalt%gamma_cased*max(0.0,1.0-cobalt%omega_calc(i,j,k)+ &
            cobalt%phi_deepresp_cased*cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n*spery)**cobalt%alpha_cased
          ! Effective thickness term that enhances burial of calcite when total sediment accumulation is high
          ! dimensionless value between 0 and 1
          cobalt%cased_redis_delz(i,j) = max(1.0, &
            cobalt%f_lithdet_btf(i,j,1)*spery+cobalt%f_cadet_calc_btf(i,j,1)*100.0*spery)**cobalt%beta_cased  
          ! calculate the sediment redissolution rate (moles Ca m-2 sec-1). This calculation is subject to
          ! three limiters: a) a maximum of 1/2 of the total cased over one time step; b) a maximum of 0.01
          ! moles Ca per day; and c) a minimum of 0.0
          cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd, min(0.5*cobalt%f_cased(i,j,1)*r_dt,  &
            cobalt%fcased_redis_surfresp(i,j)+cobalt%cased_redis_coef(i,j)*cobalt%cased_redis_delz(i,j)*cobalt%f_cased(i,j,1))) ) 
          !
          ! Old expression
          !
          !cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd,min(0.5 * cobalt%f_cased(i,j,1) * r_dt, min(0.5 *       &                          
          !   cobalt%f_cadet_calc_btf(i,j,1), 0.14307 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n) +        &
          !   0.03607 / spery * max(0.0, 1.0 - cobalt%omega_calc(i,j,k) +   &
          !   4.1228 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n * spery)**(2.7488) *                        &
          !   max(1.0, cobalt%f_lithdet_btf(i,j,1) * spery + cobalt%f_cadet_calc_btf(i,j,1) * 100.0 *  &
          !   spery)**(-2.2185) * cobalt%f_cased(i,j,1))))*grid_tmask(i,j,k)

          if (cobalt%cased_steady) then
            cobalt%fcased_burial(i,j) = cobalt%f_cadet_calc_btf(i,j,1) - cobalt%fcased_redis(i,j)
            cobalt%f_cased(i,j,1) = cobalt%fcased_burial(i,j)*cobalt%Co_cased/cobalt%f_cadet_calc_btf(i,j,1)
          else
            cobalt%fcased_burial(i,j) = max(0.0, cobalt%f_cadet_calc_btf(i,j,1) * cobalt%f_cased(i,j,1) / &
              cobalt%Co_cased)
            cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1) + (cobalt%f_cadet_calc_btf(i,j,1) -            &
              cobalt%fcased_redis(i,j) - cobalt%fcased_burial(i,j)) / cobalt%z_sed * dt *                &
              grid_tmask(i,j,k)
          endif

          ! uncomment for "no mass change" test (next 3 lines)
          !cobalt%fcased_redis(i,j) = cobalt%f_cadet_calc_btf(i,j,1)
          !cobalt%fcased_burial(i,j) = 0.0
          !cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1)
          !
          ! Bottom flux boundaries passed to the vertical mixing routine 
          !

          cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
             cobalt%fnoxic_sed(i,j) - cobalt%fno3denit_sed(i,j)*cobalt%alk_2_n_denit
          cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) -       &
             (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n
          ! uncomment for "no mass change" test (next 2 lines)
          !cobalt%b_dic(i,j) =  - cobalt%f_cadet_calc_btf(i,j,1)  - cobalt%f_cadet_arag_btf(i,j,1) -            &
          !   (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n 
          cobalt%b_fed(i,j) = - cobalt%ffe_sed(i,j) - cobalt%ffe_geotherm(i,j)
          ! uncomment for "no mass change" test (next line)
          !cobalt%b_fed(i,j) = - cobalt%f_fedet_btf(i,j,1)
          cobalt%b_nh4(i,j) = - cobalt%f_ndet_btf(i,j,1) + cobalt%fndet_burial(i,j)
          cobalt%b_no3(i,j) = cobalt%fno3denit_sed(i,j)
          cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnfeso4red_sed(i,j))
          cobalt%b_po4(i,j) = - cobalt%f_pdet_btf(i,j,1) + cobalt%fpdet_burial(i,j)
          cobalt%b_sio4(i,j)= - cobalt%f_sidet_btf(i,j,1)

       endif !}
    enddo; enddo  !} i, j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%f_cased(i,j,k) = 0.0
    enddo; enddo ; enddo  !} i,j,k

    call mpp_clock_end(id_clock_ballast_loops)

    call g_tracer_set_values(tracer_list,'alk',  'btf', cobalt%b_alk ,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic',  'btf', cobalt%b_dic ,isd,jsd)
    call g_tracer_set_values(tracer_list,'fed',  'btf', cobalt%b_fed ,isd,jsd)
    call g_tracer_set_values(tracer_list,'nh4',  'btf', cobalt%b_nh4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'no3',  'btf', cobalt%b_no3 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'o2',   'btf', cobalt%b_o2  ,isd,jsd)
    call g_tracer_set_values(tracer_list,'po4',  'btf', cobalt%b_po4 ,isd,jsd)
    call g_tracer_set_values(tracer_list,'sio4', 'btf', cobalt%b_sio4,isd,jsd)
!
    call mpp_clock_begin(id_clock_source_sink_loop1)
!
!-----------------------------------------------------------------------
! 8: Source/sink calculations 
!-----------------------------------------------------------------------
!
    !  
    !-------------------------------------------------------------------
    ! 4.1: Update the prognostics tracer fields via their pointers.
    !-------------------------------------------------------------------
    !
    call g_tracer_get_pointer(tracer_list,'alk'    ,'field',cobalt%p_alk    )
    call g_tracer_get_pointer(tracer_list,'cadet_arag','field',cobalt%p_cadet_arag)
    call g_tracer_get_pointer(tracer_list,'cadet_calc','field',cobalt%p_cadet_calc)
    call g_tracer_get_pointer(tracer_list,'dic'    ,'field',cobalt%p_dic    )
    call g_tracer_get_pointer(tracer_list,'fed'    ,'field',cobalt%p_fed    )
    call g_tracer_get_pointer(tracer_list,'fedi'   ,'field',cobalt%p_fedi   )
    call g_tracer_get_pointer(tracer_list,'felg'   ,'field',cobalt%p_felg   )
    call g_tracer_get_pointer(tracer_list,'fesm'   ,'field',cobalt%p_fesm )
    call g_tracer_get_pointer(tracer_list,'fedet'  ,'field',cobalt%p_fedet  )
    call g_tracer_get_pointer(tracer_list,'ldon'   ,'field',cobalt%p_ldon   )
    call g_tracer_get_pointer(tracer_list,'ldop'   ,'field',cobalt%p_ldop   )
    call g_tracer_get_pointer(tracer_list,'lith'   ,'field',cobalt%p_lith   )
    call g_tracer_get_pointer(tracer_list,'lithdet','field',cobalt%p_lithdet)
    call g_tracer_get_pointer(tracer_list,'nbact'  ,'field',cobalt%p_nbact  )
    call g_tracer_get_pointer(tracer_list,'ndet'   ,'field',cobalt%p_ndet   )
    call g_tracer_get_pointer(tracer_list,'ndi'    ,'field',cobalt%p_ndi    )
    call g_tracer_get_pointer(tracer_list,'nlg'    ,'field',cobalt%p_nlg    )
    call g_tracer_get_pointer(tracer_list,'nsm' ,'field',cobalt%p_nsm )
    call g_tracer_get_pointer(tracer_list,'nh4'    ,'field',cobalt%p_nh4    )
    call g_tracer_get_pointer(tracer_list,'no3'    ,'field',cobalt%p_no3    )
    call g_tracer_get_pointer(tracer_list,'o2'     ,'field',cobalt%p_o2     )
    call g_tracer_get_pointer(tracer_list,'pdet'   ,'field',cobalt%p_pdet   )
    call g_tracer_get_pointer(tracer_list,'po4'    ,'field',cobalt%p_po4    )
    call g_tracer_get_pointer(tracer_list,'srdon'   ,'field',cobalt%p_srdon   )
    call g_tracer_get_pointer(tracer_list,'srdop'   ,'field',cobalt%p_srdop   )
    call g_tracer_get_pointer(tracer_list,'sldon'   ,'field',cobalt%p_sldon   )
    call g_tracer_get_pointer(tracer_list,'sldop'   ,'field',cobalt%p_sldop   )
    call g_tracer_get_pointer(tracer_list,'sidet'  ,'field',cobalt%p_sidet  )
    call g_tracer_get_pointer(tracer_list,'silg'   ,'field',cobalt%p_silg   )
    call g_tracer_get_pointer(tracer_list,'sio4'   ,'field',cobalt%p_sio4   )
    call g_tracer_get_pointer(tracer_list,'nsmz'   ,'field',cobalt%p_nsmz   )
    call g_tracer_get_pointer(tracer_list,'nmdz'   ,'field',cobalt%p_nmdz   )
    call g_tracer_get_pointer(tracer_list,'nlgz'   ,'field',cobalt%p_nlgz   )

    if (do_14c) then
       call g_tracer_get_pointer(tracer_list,'di14c','field',cobalt%p_di14c)
       call g_tracer_get_pointer(tracer_list,'do14c','field',cobalt%p_do14c)
    endif 

    ! CAS calculate total N and P before source/sink
    ! calculate internal sources (those not applied as air-sea or benthos
    ! exchanges) to close the balance
    allocate(pre_totn(isc:iec,jsc:jec,1:nk))
    allocate(pre_totc(isc:iec,jsc:jec,1:nk))
    allocate(net_srcn(isc:iec,jsc:jec,1:nk))
    allocate(pre_totp(isc:iec,jsc:jec,1:nk))
    allocate(pre_totfe(isc:iec,jsc:jec,1:nk))
    allocate(net_srcfe(isc:iec,jsc:jec,1:nk))
    allocate(pre_totsi(isc:iec,jsc:jec,1:nk))
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
         pre_totn(i,j,k) = (cobalt%p_no3(i,j,k,tau) + cobalt%p_nh4(i,j,k,tau) + & 
                    cobalt%p_ndi(i,j,k,tau) + cobalt%p_nlg(i,j,k,tau) + &
                    cobalt%p_nsm(i,j,k,tau) + cobalt%p_nbact(i,j,k,tau) + &
                    cobalt%p_ldon(i,j,k,tau) + cobalt%p_sldon(i,j,k,tau) + &
                    cobalt%p_srdon(i,j,k,tau) +  cobalt%p_ndet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau) + cobalt%p_nmdz(i,j,k,tau) + &
                    cobalt%p_nlgz(i,j,k,tau))*grid_tmask(i,j,k)
         net_srcn(i,j,k) = (phyto(DIAZO)%juptake_n2(i,j,k) - cobalt%jno3denit_wc(i,j,k))* &
                    dt*grid_tmask(i,j,k)
         pre_totc(i,j,k) = (cobalt%p_dic(i,j,k,tau) + &
                    cobalt%p_cadet_arag(i,j,k,tau) + cobalt%p_cadet_calc(i,j,k,tau) + &
                    cobalt%c_2_n*(cobalt%p_ndi(i,j,k,tau) + cobalt%p_nlg(i,j,k,tau) + &
                    cobalt%p_nsm(i,j,k,tau) + cobalt%p_nbact(i,j,k,tau) + &
                    cobalt%p_ldon(i,j,k,tau) + cobalt%p_sldon(i,j,k,tau) + &
                    cobalt%p_srdon(i,j,k,tau) +  cobalt%p_ndet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau) + cobalt%p_nmdz(i,j,k,tau) + &
                    cobalt%p_nlgz(i,j,k,tau)))*grid_tmask(i,j,k)
         pre_totp(i,j,k) = (cobalt%p_po4(i,j,k,tau) + cobalt%p_ndi(i,j,k,tau)*phyto(1)%p_2_n_static + &
                    cobalt%p_nlg(i,j,k,tau)*phyto(2)%p_2_n_static + &
                    cobalt%p_nsm(i,j,k,tau)*phyto(3)%p_2_n_static + &
                    cobalt%p_ldop(i,j,k,tau) + cobalt%p_sldop(i,j,k,tau) + &
                    cobalt%p_srdop(i,j,k,tau) +  cobalt%p_pdet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau)*zoo(1)%q_p_2_n + &
                    cobalt%p_nmdz(i,j,k,tau)*zoo(2)%q_p_2_n + &
                    cobalt%p_nlgz(i,j,k,tau)*zoo(3)%q_p_2_n + &
                    bact(1)%q_p_2_n*cobalt%p_nbact(i,j,k,tau))*grid_tmask(i,j,k)
         pre_totfe(i,j,k) = (cobalt%p_fed(i,j,k,tau) + cobalt%p_fedi(i,j,k,tau) + &
                    cobalt%p_felg(i,j,k,tau) + cobalt%p_fesm(i,j,k,tau) + & 
                    cobalt%p_fedet(i,j,k,tau))*grid_tmask(i,j,k)
         net_srcfe(i,j,k) = cobalt%jfe_coast(i,j,k)*dt*grid_tmask(i,j,k)
         pre_totsi(i,j,k) = (cobalt%p_sio4(i,j,k,tau) + cobalt%p_silg(i,j,k,tau) + &
                    cobalt%p_sidet(i,j,k,tau))*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    do j = jsc, jec ; do i = isc, iec  !{
      net_srcfe(i,j,1) = net_srcfe(i,j,1)+cobalt%ffe_iceberg(i,j)*dt*grid_tmask(i,j,1)/rho_dzt(i,j,1)
    enddo; enddo

    if (cobalt%id_no3_in_source .gt. 0)                &
         used = g_send_data(cobalt%id_no3_in_source,         cobalt%f_no3,          &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    call mpp_clock_end(id_clock_source_sink_loop1)
    !
    !-----------------------------------------------------------------------
    ! 4.2: Source sink calculations
    !-----------------------------------------------------------------------
    !
    !     Phytoplankton Nitrogen and Phosphorus
    !
    call mpp_clock_begin(id_clock_source_sink_loop2)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Diazotrophic Phytoplankton Nitrogen
       !
       cobalt%jndi(i,j,k) = phyto(DIAZO)%mu(i,j,k)*phyto(DIAZO)%f_n(i,j,k) - &
                            phyto(DIAZO)%jzloss_n(i,j,k) -       &
                            phyto(DIAZO)%jhploss_n(i,j,k) - phyto(DIAZO)%jaggloss_n(i,j,k) -       &
                            phyto(DIAZO)%jvirloss_n(i,j,k) - phyto(DIAZO)%jexuloss_n(i,j,k)
       cobalt%p_ndi(i,j,k,tau) = cobalt%p_ndi(i,j,k,tau) + cobalt%jndi(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Nitrogen
       !
       cobalt%jnlg(i,j,k) = phyto(LARGE)%mu(i,j,k)*phyto(LARGE)%f_n(i,j,k) -    &
                            phyto(LARGE)%jzloss_n(i,j,k) - phyto(LARGE)%jhploss_n(i,j,k) -         &
                            phyto(LARGE)%jaggloss_n(i,j,k) - phyto(LARGE)%jvirloss_n(i,j,k) -      &
                            phyto(LARGE)%jexuloss_n(i,j,k)
       cobalt%p_nlg(i,j,k,tau) = cobalt%p_nlg(i,j,k,tau) + cobalt%jnlg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Nitrogen
       !
       cobalt%jnsm(i,j,k) = phyto(SMALL)%mu(i,j,k)*phyto(SMALL)%f_n(i,j,k) -    &
                            phyto(SMALL)%jzloss_n(i,j,k) - phyto(SMALL)%jhploss_n(i,j,k) -         &
                            phyto(SMALL)%jaggloss_n(i,j,k) - phyto(SMALL)%jvirloss_n(i,j,k) -      &
                            phyto(SMALL)%jexuloss_n(i,j,k)                                         
       cobalt%p_nsm(i,j,k,tau) = cobalt%p_nsm(i,j,k,tau) + cobalt%jnsm(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
!
    call mpp_clock_end(id_clock_source_sink_loop2)
    !
    !     Phytoplankton Silicon and Iron
    !
    call mpp_clock_begin(id_clock_source_sink_loop3)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Large Phytoplankton Silicon
       !
       cobalt%jsilg(i,j,k) = phyto(LARGE)%juptake_sio4(i,j,k) - & 
                             phyto(LARGE)%jzloss_sio2(i,j,k) - phyto(LARGE)%jhploss_sio2(i,j,k) - &
                             phyto(LARGE)%jaggloss_sio2(i,j,k) - phyto(LARGE)%jvirloss_sio2(i,j,k)
       cobalt%p_silg(i,j,k,tau) = cobalt%p_silg(i,j,k,tau) + cobalt%jsilg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Diazotrophic Phytoplankton Iron
       !
       cobalt%jfedi(i,j,k) = phyto(DIAZO)%juptake_fe(i,j,k) - &
                             phyto(DIAZO)%jzloss_fe(i,j,k) - &
                             phyto(DIAZO)%jhploss_fe(i,j,k) - phyto(DIAZO)%jaggloss_fe(i,j,k) - &
                             phyto(DIAZO)%jvirloss_fe(i,j,k) - phyto(DIAZO)%jexuloss_fe(i,j,k)
       cobalt%p_fedi(i,j,k,tau) = cobalt%p_fedi(i,j,k,tau) + cobalt%jfedi(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large Phytoplankton Iron
       !
       cobalt%jfelg(i,j,k) = phyto(LARGE)%juptake_fe(i,j,k) - & 
                             phyto(LARGE)%jzloss_fe(i,j,k) - &
                             phyto(LARGE)%jhploss_fe(i,j,k) - phyto(LARGE)%jaggloss_fe(i,j,k) - &
                             phyto(LARGE)%jvirloss_fe(i,j,k) - phyto(LARGE)%jexuloss_fe(i,j,k)
       cobalt%p_felg(i,j,k,tau) = cobalt%p_felg(i,j,k,tau) + cobalt%jfelg(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Small Phytoplankton Iron
       !
       cobalt%jfesm(i,j,k) = phyto(SMALL)%juptake_fe(i,j,k) - &
                                phyto(SMALL)%jzloss_fe(i,j,k) - &
                                phyto(SMALL)%jhploss_fe(i,j,k) - phyto(SMALL)%jaggloss_fe(i,j,k) - &
                                phyto(SMALL)%jvirloss_fe(i,j,k) - phyto(SMALL)%jexuloss_fe(i,j,k)
       cobalt%p_fesm(i,j,k,tau) = cobalt%p_fesm(i,j,k,tau) + cobalt%jfesm(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Bacteria
       !
       cobalt%jnbact(i,j,k) = bact(1)%jprod_n(i,j,k) - bact(1)%jzloss_n(i,j,k) - &
                              bact(1)%jvirloss_n(i,j,k) - bact(1)%jhploss_n(i,j,k)  
       cobalt%p_nbact(i,j,k,tau) = cobalt%p_nbact(i,j,k,tau) + cobalt%jnbact(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    call mpp_clock_end(id_clock_source_sink_loop3)
    !
    !    Zooplankton 
    !
    call mpp_clock_begin(id_clock_source_sink_loop4)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Small zooplankton
       !
       cobalt%jnsmz(i,j,k) = zoo(1)%jprod_n(i,j,k) - zoo(1)%jzloss_n(i,j,k) - &
                             zoo(1)%jhploss_n(i,j,k)
       cobalt%p_nsmz(i,j,k,tau) = cobalt%p_nsmz(i,j,k,tau) + cobalt%jnsmz(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Medium zooplankton
       !
       cobalt%jnmdz(i,j,k) = zoo(2)%jprod_n(i,j,k) - zoo(2)%jzloss_n(i,j,k) - &
                             zoo(2)%jhploss_n(i,j,k)
       cobalt%p_nmdz(i,j,k,tau) = cobalt%p_nmdz(i,j,k,tau) + cobalt%jnmdz(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Large zooplankton
       !
       cobalt%jnlgz(i,j,k) = zoo(3)%jprod_n(i,j,k) - zoo(3)%jzloss_n(i,j,k) - &
                             zoo(3)%jhploss_n(i,j,k)
       cobalt%p_nlgz(i,j,k,tau) = cobalt%p_nlgz(i,j,k,tau) + cobalt%jnlgz(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
!
    call mpp_clock_end(id_clock_source_sink_loop4)
    !
    !     NO3
    !
    call mpp_clock_begin(id_clock_source_sink_loop5)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%jno3(i,j,k) =  cobalt%jnitrif(i,j,k) - phyto(DIAZO)%juptake_no3(i,j,k) -  &
                             phyto(LARGE)%juptake_no3(i,j,k) - phyto(SMALL)%juptake_no3(i,j,k) - &
                             cobalt%jno3denit_wc(i,j,k)
       cobalt%p_no3(i,j,k,tau) = cobalt%p_no3(i,j,k,tau) + cobalt%jno3(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     Other nutrients
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! NH4
       !
       cobalt%jnh4(i,j,k) = cobalt%jprod_nh4(i,j,k) - phyto(DIAZO)%juptake_nh4(i,j,k) - &
                            phyto(LARGE)%juptake_nh4(i,j,k) - phyto(SMALL)%juptake_nh4(i,j,k) - &
                            cobalt%jnitrif(i,j,k)
       cobalt%p_nh4(i,j,k,tau) = cobalt%p_nh4(i,j,k,tau) + cobalt%jnh4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! PO4
       !
       cobalt%jpo4(i,j,k) = cobalt%jprod_po4(i,j,k) - phyto(DIAZO)%juptake_po4(i,j,k) - &
                            phyto(LARGE)%juptake_po4(i,j,k) - phyto(SMALL)%juptake_po4(i,j,k)
       cobalt%p_po4(i,j,k,tau) = cobalt%p_po4(i,j,k,tau) + cobalt%jpo4(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! SiO4
       !
       cobalt%jsio4(i,j,k) = cobalt%jprod_sio4(i,j,k) - phyto(LARGE)%juptake_sio4(i,j,k)
       cobalt%p_sio4(i,j,k,tau) = cobalt%p_sio4(i,j,k,tau) + cobalt%jsio4(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    ! 2016/06/13 JGJ: keep original Fed calculation
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
          !
          ! Fed
          ! use original code to compute jprod_fed, jfed and p_fed
          !
       cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + cobalt%jfe_coast(i,j,k) 
       cobalt%jfed(i,j,k) = cobalt%jprod_fed(i,j,k) - phyto(DIAZO)%juptake_fe(i,j,k) - &
                            phyto(LARGE)%juptake_fe(i,j,k) -  phyto(SMALL)%juptake_fe(i,j,k) - &
                            cobalt%jfe_ads(i,j,k)
       cobalt%p_fed(i,j,k,tau) = cobalt%p_fed(i,j,k,tau) + cobalt%jfed(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo; enddo  !} i,j,k

    call mpp_clock_end(id_clock_source_sink_loop5)
    !
    !-----------------------------------------------------------------------
    !     Detrital Components
    !-----------------------------------------------------------------------
    !
    call mpp_clock_begin(id_clock_source_sink_loop6)
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Cadet_arag
       !
       cobalt%jcadet_arag(i,j,k) = cobalt%jprod_cadet_arag(i,j,k) - cobalt%jdiss_cadet_arag(i,j,k) 
       cobalt%p_cadet_arag(i,j,k,tau) = cobalt%p_cadet_arag(i,j,k,tau) + cobalt%jcadet_arag(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Cadet_calc
       !
       cobalt%jcadet_calc(i,j,k) = cobalt%jprod_cadet_calc(i,j,k) - cobalt%jdiss_cadet_calc(i,j,k)
       cobalt%p_cadet_calc(i,j,k,tau) = cobalt%p_cadet_calc(i,j,k,tau) + cobalt%jcadet_calc(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Lithdet
       !
       cobalt%jlithdet(i,j,k) = cobalt%jprod_lithdet(i,j,k) 
       cobalt%p_lithdet(i,j,k,tau) = cobalt%p_lithdet(i,j,k,tau) + cobalt%jlithdet(i,j,k) * dt *  &
                                     grid_tmask(i,j,k)
       !
       ! Ndet
       !
       cobalt%jndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - cobalt%jremin_ndet(i,j,k) - &
                             cobalt%det_jzloss_n(i,j,k) - cobalt%det_jhploss_n(i,j,k)
       cobalt%p_ndet(i,j,k,tau) = cobalt%p_ndet(i,j,k,tau) + cobalt%jndet(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Pdet
       !
       cobalt%jpdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - cobalt%jremin_pdet(i,j,k) - &
                             cobalt%det_jzloss_p(i,j,k) - cobalt%det_jhploss_p(i,j,k)         
       cobalt%p_pdet(i,j,k,tau) = cobalt%p_pdet(i,j,k,tau) + cobalt%jpdet(i,j,k)*dt*grid_tmask(i,j,k)
       !
       ! Sidet
       !
       cobalt%jsidet(i,j,k) = cobalt%jprod_sidet(i,j,k) - & 
                              cobalt%jdiss_sidet(i,j,k) - cobalt%det_jzloss_si(i,j,k) - &
                              cobalt%det_jhploss_si(i,j,k)
       cobalt%p_sidet(i,j,k,tau) = cobalt%p_sidet(i,j,k,tau) + cobalt%jsidet(i,j,k)*dt*grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k

    ! 2016/06/13 JGJ: keep original jfedet calculation
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
          !
          ! Fedet
          ! use original code to compute fedet
          !
       cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%jfe_ads(i,j,k)
       cobalt%jfedet(i,j,k) = cobalt%jprod_fedet(i,j,k) - &
                              cobalt%jremin_fedet(i,j,k) - cobalt%det_jzloss_fe(i,j,k) - & 
                              cobalt%det_jhploss_fe(i,j,k)
       cobalt%p_fedet(i,j,k,tau) = cobalt%p_fedet(i,j,k,tau) + cobalt%jfedet(i,j,k)*dt*grid_tmask(i,j,k) 
    enddo; enddo; enddo  !} i,j,k
    !
    !     Dissolved Organic Matter
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Labile Dissolved Organic Nitrogen
       !
       cobalt%jldon(i,j,k) = cobalt%jprod_ldon(i,j,k) + &
                             cobalt%gamma_sldon*cobalt%f_sldon(i,j,k) + &
                             cobalt%gamma_srdon*cobalt%f_srdon(i,j,k) - bact(1)%juptake_ldon(i,j,k)
       cobalt%p_ldon(i,j,k,tau) = cobalt%p_ldon(i,j,k,tau) +  cobalt%jldon(i,j,k)*dt*               &
            grid_tmask(i,j,k)
       !
       ! Labile Dissolved Organic Phosphorous
       !
       cobalt%jldop(i,j,k) = cobalt%jprod_ldop(i,j,k) + &
                             cobalt%gamma_sldop*cobalt%f_sldop(i,j,k) + &
                             cobalt%gamma_srdop*cobalt%f_srdop(i,j,k) - bact(1)%juptake_ldop(i,j,k)
       cobalt%p_ldop(i,j,k,tau) = cobalt%p_ldop(i,j,k,tau) +  cobalt%jldop(i,j,k)*dt*               &
                             grid_tmask(i,j,k)
       !
       ! Semilabile Dissolved Organic Nitrogen
       !
       cobalt%jsldon(i,j,k) = cobalt%jprod_sldon(i,j,k) - &
                              cobalt%gamma_sldon*cobalt%f_sldon(i,j,k)
       cobalt%p_sldon(i,j,k,tau) = cobalt%p_sldon(i,j,k,tau) +  cobalt%jsldon(i,j,k) * dt *               &
            grid_tmask(i,j,k)
       !
       ! Semilabile dissolved organic phosphorous  
       !
       cobalt%jsldop(i,j,k) = cobalt%jprod_sldop(i,j,k) - &
                              cobalt%gamma_sldop*cobalt%f_sldop(i,j,k)
       cobalt%p_sldop(i,j,k,tau) = cobalt%p_sldop(i,j,k,tau) + cobalt%jsldop(i,j,k) * dt *                &
                                  grid_tmask(i,j,k)
       !
       ! Refractory Dissolved Organic Nitrogen
       ! 
       cobalt%jsrdon(i,j,k) = cobalt%jprod_srdon(i,j,k) -  cobalt%gamma_srdon * cobalt%f_srdon(i,j,k)
       cobalt%p_srdon(i,j,k,tau) = cobalt%p_srdon(i,j,k,tau) +  cobalt%jsrdon(i,j,k) * dt *               &
            grid_tmask(i,j,k)
       !
       ! Refractory dissolved organic phosphorous
       !
       cobalt%jsrdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - cobalt%gamma_srdop * cobalt%f_srdop(i,j,k)
       cobalt%p_srdop(i,j,k,tau) = cobalt%p_srdop(i,j,k,tau) + cobalt%jsrdop(i,j,k) * dt *                &
                                  grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     O2
    !
    do k = 1, nk ; do j =jsc, jec ; do i = isc, iec  !{
       cobalt%jo2(i,j,k) = (cobalt%o2_2_no3 * (phyto(DIAZO)%juptake_no3(i,j,k) +   &
            phyto(LARGE)%juptake_no3(i,j,k) + phyto(SMALL)%juptake_no3(i,j,k)) + & 
             cobalt%o2_2_nh4 *       &
            (phyto(DIAZO)%juptake_nh4(i,j,k) + phyto(LARGE)%juptake_nh4(i,j,k) +      &
            phyto(SMALL)%juptake_nh4(i,j,k) + &  
            phyto(DIAZO)%juptake_n2(i,j,k))) * grid_tmask(i,j,k)
       cobalt%jo2(i,j,k) = cobalt%jo2(i,j,k) - cobalt%jo2resp_wc(i,j,k)
       cobalt%p_o2(i,j,k,tau) = cobalt%p_o2(i,j,k,tau) + cobalt%jo2(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    !
    !     The Carbon system
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       !
       ! Alkalinity
       ! CAS: remove o2 removal via nitrification from the total o2 respired 
       !      to isolate the change in alkalinity due to organic matter remineralization
       !
       cobalt%jalk(i,j,k) = 2.0 * (cobalt%jdiss_cadet_arag(i,j,k) +        &
          cobalt%jdiss_cadet_calc(i,j,k) - cobalt%jprod_cadet_arag(i,j,k) - &
          cobalt%jprod_cadet_calc(i,j,k)) + phyto(DIAZO)%juptake_no3(i,j,k) + &
          phyto(LARGE)%juptake_no3(i,j,k) + phyto(SMALL)%juptake_no3(i,j,k) + &
          (cobalt%jo2resp_wc(i,j,k)-cobalt%jnitrif(i,j,k)*cobalt%o2_2_nitrif)/cobalt%o2_2_nh4 + &
          cobalt%alk_2_n_denit*cobalt%jno3denit_wc(i,j,k) - & 
          phyto(DIAZO)%juptake_nh4(i,j,k) - phyto(LARGE)%juptake_nh4(i,j,k) - &  
          phyto(SMALL)%juptake_nh4(i,j,k) - 2.0 * cobalt%jnitrif(i,j,k)

       cobalt%p_alk(i,j,k,tau) = cobalt%p_alk(i,j,k,tau) + cobalt%jalk(i,j,k) * dt * grid_tmask(i,j,k)
       !
       ! Dissolved Inorganic Carbon
       !
       cobalt%jdic(i,j,k) =(cobalt%c_2_n * (cobalt%jno3(i,j,k) + &
          cobalt%jnh4(i,j,k) + cobalt%jno3denit_wc(i,j,k) - phyto(DIAZO)%juptake_n2(i,j,k)) + &
          cobalt%jdiss_cadet_arag(i,j,k) + cobalt%jdiss_cadet_calc(i,j,k) - &
          cobalt%jprod_cadet_arag(i,j,k) - cobalt%jprod_cadet_calc(i,j,k))
       cobalt%p_dic(i,j,k,tau) = cobalt%p_dic(i,j,k,tau) + cobalt%jdic(i,j,k) * dt * grid_tmask(i,j,k)
    enddo; enddo ; enddo !} i,j,k
!
       
    if (do_14c) then                                        !<<RADIOCARBON

         do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{

          cobalt%c14_2_n(i,j,k) = cobalt%c_2_n *                             &
            cobalt%f_di14c(i,j,k) / (epsln + cobalt%f_dic(i,j,k))

         enddo; enddo ; enddo !} i,j,k

      ! Sinking particulate 14C is generated in the local ratio of 14C/12C
      ! to sinking 12C, which itself is strictly tied to P through a fixed
      ! C:P. Therefore, jpop can be used to calculate fpo14c.

      do j = jsc, jec ;      do i = isc, iec   !{
        cobalt%fpo14c(i,j,1) =  (cobalt%jprod_ndet(i,j,1) - (cobalt%jremin_ndet(i,j,1) +          &
                             cobalt%det_jzloss_n(i,j,1) + cobalt%det_jhploss_n(i,j,1))) *         &
                             cobalt%c14_2_n(i,j,1) * rho_dzt(i,j,1) 
        cobalt%j14c_reminp(i,j,1) = (-1) * cobalt%fpo14c(i,j,1) / rho_dzt(i,j,1)
      enddo; enddo !} i,j

      do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
        cobalt%fpo14c(i,j,k) = max(0.,cobalt%fpo14c(i,j,k-1) +          &
                               (cobalt%jprod_ndet(i,j,k) * cobalt%c14_2_n(i,j,k) - (cobalt%jremin_ndet(i,j,k) + &
                               cobalt%det_jzloss_n(i,j,k) + cobalt%det_jhploss_n(i,j,k)) *                      &
                               cobalt%fpo14c(i,j,k-1) / max(epsln,cobalt%f_ndet(i,j,k-1) * cobalt%Rho_0 *       &
                               cobalt%wsink)) * rho_dzt(i,j,k))
 
         cobalt%j14c_reminp(i,j,k) = (cobalt%fpo14c(i,j,k-1) - cobalt%fpo14c(i,j,k)) / rho_dzt(i,j,k)
      enddo; enddo ; enddo !} i,j,k

     ! Decay the radiocarbon in both DIC and DOC

      do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{        

        cobalt%j14c_decay_dic(i,j,k) = cobalt%f_di14c(i,j,k) *               &
          cobalt%lambda_14c 

        cobalt%j14c_decay_doc(i,j,k) = cobalt%f_do14c(i,j,k) *               &
          cobalt%lambda_14c 

      enddo; enddo ; enddo !} i,j,k

      do j = jsc, jec ; do i = isc, iec  !{
         k = grid_kmt(i,j)
         if (k .gt. 0) then !{
           cobalt%b_di14c(i,j) = - cobalt%fpo14c(i,j,k)- cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) 
         endif  
      enddo; enddo  !} i, j

     call g_tracer_set_values(tracer_list,'di14c','btf',cobalt%b_di14c,isd,jsd)
!
! Include only 14C in the semirefractory component of DOC
!
     do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       cobalt%jdo14c(i,j,k) = cobalt%jprod_srdon(i,j,k) * cobalt%c14_2_n(i,j,k) - &
           cobalt%gamma_srdon * cobalt%f_do14c(i,j,k)
         
       cobalt%p_do14c(i,j,k,tau) = cobalt%p_do14c(i,j,k,tau) +               &
         (cobalt%jdo14c(i,j,k) - cobalt%j14c_decay_doc(i,j,k)) * dt          &
         * grid_tmask(i,j,k)
!
! Use the DIC budget except remove the srdon component and sinking detritus components which are treated separately
!
       cobalt%jdi14c(i,j,k) =(cobalt%c14_2_n(i,j,k) * (cobalt%jno3(i,j,k) + &
          cobalt%jnh4(i,j,k) + cobalt%jno3denit_wc(i,j,k) - phyto(DIAZO)%juptake_n2(i,j,k)) + &
          cobalt%jsrdon(i,j,k)) + cobalt%jdiss_cadet_arag(i,j,k) + cobalt%jdiss_cadet_calc(i,j,k) - &
          cobalt%jprod_cadet_arag(i,j,k) - cobalt%jprod_cadet_calc(i,j,k) -&
          cobalt%jdo14c(i,j,k) + cobalt%j14c_reminp(i,j,k) 

       cobalt%p_di14c(i,j,k,tau) = cobalt%p_di14c(i,j,k,tau) +               &
         (cobalt%jdi14c(i,j,k) - cobalt%j14c_decay_dic(i,j,k)) * dt          &
         * grid_tmask(i,j,k)
     enddo; enddo ; enddo !} i,j,k
    endif                                                   !RADIOCARBON>>
    !
    !-----------------------------------------------------------------------
    !     Lithogenic aluminosilicate particulates
    !-----------------------------------------------------------------------
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
       cobalt%p_lith(i,j,k,tau) = cobalt%p_lith(i,j,k,tau) - cobalt%jlithdet(i,j,k) * dt *        &
            grid_tmask(i,j,k)
    enddo; enddo ; enddo  !} i,j,k
    call mpp_clock_end(id_clock_source_sink_loop6)
    call mpp_clock_begin(id_clock_cobalt_calc_diagnostics)
    !
    !Set the diagnostics tracer fields.
    !
    call g_tracer_set_values(tracer_list,'cased',  'field',cobalt%f_cased    ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'chl',    'field',cobalt%f_chl      ,isd,jsd,ntau=1)
    if (do_nh3_diag) call g_tracer_set_values(tracer_list,'nh3',    'field',cobalt%f_nh3      ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'co3_ion','field',cobalt%f_co3_ion  ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'irr_mem' ,'field',cobalt%f_irr_mem ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'mu_mem_ndi' ,'field',phyto(DIAZO)%f_mu_mem ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'mu_mem_nlg' ,'field',phyto(LARGE)%f_mu_mem ,isd,jsd,ntau=1)
    call g_tracer_set_values(tracer_list,'mu_mem_nsm' ,'field',phyto(SMALL)%f_mu_mem ,isd,jsd,ntau=1)

    ! CAS calculate totals after source/sinks have been applied
    imbal_flag = 0;
    stdoutunit = stdout();
    allocate(post_totn(isc:iec,jsc:jec,1:nk))
    allocate(post_totc(isc:iec,jsc:jec,1:nk))
    allocate(post_totp(isc:iec,jsc:jec,1:nk))
    allocate(post_totsi(isc:iec,jsc:jec,1:nk))
    allocate(post_totfe(isc:iec,jsc:jec,1:nk))
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
         post_totn(i,j,k) = (cobalt%p_no3(i,j,k,tau) + cobalt%p_nh4(i,j,k,tau) + &
                    cobalt%p_ndi(i,j,k,tau) + cobalt%p_nlg(i,j,k,tau) + &
                    cobalt%p_nsm(i,j,k,tau) + cobalt%p_nbact(i,j,k,tau) + &
                    cobalt%p_ldon(i,j,k,tau) + cobalt%p_sldon(i,j,k,tau) + &
                    cobalt%p_srdon(i,j,k,tau) +  cobalt%p_ndet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau) + cobalt%p_nmdz(i,j,k,tau) + &
                    cobalt%p_nlgz(i,j,k,tau))*grid_tmask(i,j,k)
         imbal = (post_totn(i,j,k) - pre_totn(i,j,k) - net_srcn(i,j,k))*86400.0/dt*1.03e6
         if (abs(imbal).gt.imbalance_tolerance) then
           call mpp_error(FATAL,&
           '==>biological source/sink imbalance (generic_COBALT_update_from_source): Nitrogen')
         endif

         post_totc(i,j,k) = (cobalt%p_dic(i,j,k,tau) + &
                    cobalt%p_cadet_arag(i,j,k,tau) + cobalt%p_cadet_calc(i,j,k,tau) + &
                    cobalt%c_2_n*(cobalt%p_ndi(i,j,k,tau) + cobalt%p_nlg(i,j,k,tau) + &
                    cobalt%p_nsm(i,j,k,tau) + cobalt%p_nbact(i,j,k,tau) + &
                    cobalt%p_ldon(i,j,k,tau) + cobalt%p_sldon(i,j,k,tau) + &
                    cobalt%p_srdon(i,j,k,tau) +  cobalt%p_ndet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau) + cobalt%p_nmdz(i,j,k,tau) + &
                    cobalt%p_nlgz(i,j,k,tau)))*grid_tmask(i,j,k)
        imbal = (post_totc(i,j,k) - pre_totc(i,j,k))*86400.0/dt*1.03e6
         if (abs(imbal).gt.imbalance_tolerance) then
           call mpp_error(FATAL,&
           '==>biological source/sink imbalance (generic_COBALT_update_from_source): Carbon')
         endif

         post_totp(i,j,k) = (cobalt%p_po4(i,j,k,tau) + cobalt%p_ndi(i,j,k,tau)*phyto(1)%p_2_n_static + &
                    cobalt%p_nlg(i,j,k,tau)*phyto(2)%p_2_n_static + &
                    cobalt%p_nsm(i,j,k,tau)*phyto(3)%p_2_n_static + &
                    cobalt%p_ldop(i,j,k,tau) + cobalt%p_sldop(i,j,k,tau) + &
                    cobalt%p_srdop(i,j,k,tau) +  cobalt%p_pdet(i,j,k,tau) + &
                    cobalt%p_nsmz(i,j,k,tau)*zoo(1)%q_p_2_n + &
                    cobalt%p_nmdz(i,j,k,tau)*zoo(2)%q_p_2_n + &
                    cobalt%p_nlgz(i,j,k,tau)*zoo(3)%q_p_2_n + &
                    bact(1)%q_p_2_n*cobalt%p_nbact(i,j,k,tau))*grid_tmask(i,j,k)
         imbal = (post_totp(i,j,k) - pre_totp(i,j,k))*86400.0/dt*1.03e6
         if (abs(imbal).gt.imbalance_tolerance) then
           call mpp_error(FATAL,&
           '==>biological source/sink imbalance (generic_COBALT_update_from_source): Phosphorus')
         endif

         post_totfe(i,j,k) = (cobalt%p_fed(i,j,k,tau) + cobalt%p_fedi(i,j,k,tau) + &
                    cobalt%p_felg(i,j,k,tau) + cobalt%p_fesm(i,j,k,tau) + &
                    cobalt%p_fedet(i,j,k,tau))*grid_tmask(i,j,k)
         imbal = (post_totfe(i,j,k) - pre_totfe(i,j,k) - net_srcfe(i,j,k))*86400.0/dt*1.03e6
         if (abs(imbal).gt.imbalance_tolerance) then
           call mpp_error(FATAL,&
           '==>biological source/sink imbalance (generic_COBALT_update_from_source): Iron')
         endif

         post_totsi(i,j,k) = (cobalt%p_sio4(i,j,k,tau) + cobalt%p_silg(i,j,k,tau) + &
                    cobalt%p_sidet(i,j,k,tau))*grid_tmask(i,j,k)
         imbal = (post_totsi(i,j,k) - pre_totsi(i,j,k))*86400.0/dt*1.03e6
         if (abs(imbal).gt.imbalance_tolerance) then
           call mpp_error(FATAL,&
           '==>biological source/sink imbalance (generic_COBALT_update_from_source): Silica')
         endif
    enddo; enddo ; enddo  !} i,j,k
 

    !
    !
    !-----------------------------------------------------------------------
    !       Save variables for diagnostics
    !-----------------------------------------------------------------------
    !

    do j = jsc, jec ; do i = isc, iec  !{
      if (grid_kmt(i,j) .gt. 0) then !{
        cobalt%o2min(i,j)=cobalt%p_o2(i,j,1,tau)
        cobalt%z_o2min(i,j)=cobalt%zt(i,j,1)
        cobalt%z_sat_arag(i,j)=missing_value1
        cobalt%z_sat_calc(i,j)=missing_value1
        cobalt%mask_z_sat_arag(i,j) = .FALSE.
        cobalt%mask_z_sat_calc(i,j) = .FALSE.
        if (cobalt%omega_arag(i,j,1) .le. 1.0) cobalt%z_sat_arag(i,j)=0.0
        if (cobalt%omega_calc(i,j,1) .le. 1.0) cobalt%z_sat_calc(i,j)=0.0
      endif !}
    enddo ; enddo  !} i,j,k
    do j = jsc, jec ; do i = isc, iec  !{
    first = .true.
      do k = 2, nk
         if (k .le. grid_kmt(i,j) .and. first) then !{
           if (cobalt%p_o2(i,j,k,tau) .lt. cobalt%p_o2(i,j,k-1,tau)) then
             cobalt%o2min(i,j)=cobalt%p_o2(i,j,k,tau)
             cobalt%z_o2min(i,j)=cobalt%zt(i,j,k)
           else
             first = .false.
           endif !}
         endif !}
      enddo;
    enddo ; enddo  !} i,j

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec  !{
      if (k .le. grid_kmt(i,j)) then !{
        if (cobalt%omega_arag(i,j,k) .le. 1.0 .and. cobalt%z_sat_arag(i,j) .lt. 0.0) then
          cobalt%z_sat_arag(i,j)=cobalt%zt(i,j,k)
          cobalt%mask_z_sat_arag(i,j) = .TRUE.
        endif
        if (cobalt%omega_calc(i,j,k) .le. 1.0 .and. cobalt%z_sat_calc(i,j) .lt. 0.0) then
          cobalt%z_sat_calc(i,j)=cobalt%zt(i,j,k)
          cobalt%mask_z_sat_calc(i,j) = .TRUE.
        endif
      endif !}
    enddo; enddo ; enddo  !} i,j,k

    ! O2 saturation
    do k = 1, nk  ; do j = jsc, jec ; do i = isc, iec
       sal = min(42.0,max(0.0,Salt(i,j,k)))
       tt = 298.15 - min(40.0,max(0.0,Temp(i,j,k)))
       tkb = 273.15 + min(40.0,max(0.0,Temp(i,j,k)))
       ts = log(tt / tkb)
       ts2 = ts  * ts
       ts3 = ts2 * ts
       ts4 = ts3 * ts
       ts5 = ts4 * ts

       !The atmospheric code needs solubilities in units of mol/m3/atm
       cobalt%o2sat(i,j,k) = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
            exp( cobalt%a_0 + cobalt%a_1*ts + cobalt%a_2*ts2 + cobalt%a_3*ts3 + cobalt%a_4*ts4 + cobalt%a_5*ts5 + &
            (cobalt%b_0 + cobalt%b_1*ts + cobalt%b_2*ts2 + cobalt%b_3*ts3 + cobalt%c_0*sal)*sal)

    enddo; enddo ; enddo  !} i,j,k

    !
    !---------------------------------------------------------------------
    ! Calculate total carbon  = Dissolved Inorganic Carbon + Phytoplankton Carbon
    !   + Dissolved Organic Carbon (including refractory) + Heterotrophic Biomass
    !   + Detrital Orgainc and Inorganic Carbon
    ! For the oceanic carbon budget, a constant 42 uM of dissolved organic
    ! carbon is added to represent the refractory component.
    ! For the oceanic nitrogen budget, a constant 2 uM of dissolved organic
    ! nitrogen is added to represent the refractory component.
    !---------------------------------------------------------------------
    !
    cobalt%tot_layer_int_c(:,:,:) = (cobalt%p_dic(:,:,:,tau) + cobalt%doc_background + cobalt%p_cadet_arag(:,:,:,tau) +&
         cobalt%p_cadet_calc(:,:,:,tau) + cobalt%c_2_n * (cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) +      &
         cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) + &
         cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) + cobalt%p_srdon(:,:,:,tau) +  &
         cobalt%p_ndet(:,:,:,tau) + cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + &
         cobalt%p_nlgz(:,:,:,tau))) * rho_dzt(:,:,:)

    cobalt%tot_layer_int_fe(:,:,:) = (cobalt%p_fed(:,:,:,tau) + cobalt%p_fedi(:,:,:,tau) + &
         cobalt%p_felg(:,:,:,tau) + cobalt%p_fesm(:,:,:,tau) + & 
         cobalt%p_fedet(:,:,:,tau)) * rho_dzt(:,:,:) 

    cobalt%tot_layer_int_n(:,:,:) = (cobalt%p_no3(:,:,:,tau) + &
         cobalt%p_nh4(:,:,:,tau) + cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) + &
         cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) + &
         cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) + cobalt%p_srdon(:,:,:,tau) +  cobalt%p_ndet(:,:,:,tau) + &
         cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + cobalt%p_nlgz(:,:,:,tau)) * & 
         rho_dzt(:,:,:)

    cobalt%tot_layer_int_p(:,:,:) = (cobalt%p_po4(:,:,:,tau) + &
         cobalt%p_ndi(:,:,:,tau)*phyto(1)%p_2_n_static + &
         cobalt%p_nlg(:,:,:,tau)*phyto(2)%p_2_n_static + &
         cobalt%p_nsm(:,:,:,tau)*phyto(3)%p_2_n_static + &
         cobalt%p_ldop(:,:,:,tau) + cobalt%p_sldop(:,:,:,tau) + &
         cobalt%p_srdop(:,:,:,tau) + cobalt%p_pdet(:,:,:,tau) + &
         bact(1)%q_p_2_n*cobalt%p_nbact(:,:,:,tau) + zoo(1)%q_p_2_n*cobalt%p_nsmz(:,:,:,tau) +  &
         zoo(2)%q_p_2_n*cobalt%p_nmdz(:,:,:,tau) + zoo(3)%q_p_2_n*cobalt%p_nlgz(:,:,:,tau))  &
         * rho_dzt(:,:,:)

    cobalt%tot_layer_int_si(:,:,:) = (cobalt%p_sio4(:,:,:,tau) + cobalt%p_silg(:,:,:,tau) +   &
         cobalt%p_sidet(:,:,:,tau)) * rho_dzt(:,:,:)

    cobalt%tot_layer_int_o2(:,:,:) = cobalt%p_o2(:,:,:,tau)*rho_dzt(:,:,:)

    cobalt%tot_layer_int_alk(:,:,:) = cobalt%p_alk(:,:,:,tau)*rho_dzt(:,:,:)

! CHECK3
    !add background of 42 uM (as in other parts of cobalt)- may need to change to 3.8e-5 per JPD
    cobalt%tot_layer_int_dic(:,:,:) = cobalt%p_dic(:,:,:,tau)*rho_dzt(:,:,:)

! CHECK3
! CAS: spreadsheet indicates no background, should we remove it?
    cobalt%tot_layer_int_doc(:,:,:) = cobalt%doc_background +  cobalt%c_2_n * (cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) +  &
         cobalt%p_srdon(:,:,:,tau)) * rho_dzt(:,:,:)

   cobalt%tot_layer_int_poc(:,:,:) = (cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) + cobalt%p_nsm(:,:,:,tau) + &
         cobalt%p_nbact(:,:,:,tau) + cobalt%p_ndet(:,:,:,tau) + cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + &
         cobalt%p_nlgz(:,:,:,tau))*cobalt%c_2_n*rho_dzt(:,:,:)


    !
    !---------------------------------------------------------------------
    ! calculate water column vertical integrals for diagnostics
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec !{
       cobalt%wc_vert_int_c(i,j) = 0.0
       cobalt%wc_vert_int_dic(i,j) = 0.0
       cobalt%wc_vert_int_doc(i,j) = 0.0
       cobalt%wc_vert_int_poc(i,j) = 0.0
       cobalt%wc_vert_int_n(i,j) = 0.0
       cobalt%wc_vert_int_p(i,j) = 0.0
       cobalt%wc_vert_int_fe(i,j) = 0.0
       cobalt%wc_vert_int_si(i,j) = 0.0
       cobalt%wc_vert_int_o2(i,j) = 0.0
       cobalt%wc_vert_int_alk(i,j) = 0.0
       cobalt%wc_vert_int_jdiss_sidet(i,j) = 0.0
       cobalt%wc_vert_int_jdiss_cadet(i,j) = 0.0
       cobalt%wc_vert_int_jo2resp(i,j) = 0.0
       cobalt%wc_vert_int_jprod_cadet(i,j) = 0.0
       cobalt%wc_vert_int_jno3denit(i,j) = 0.0
       cobalt%wc_vert_int_jnitrif(i,j) = 0.0
       cobalt%wc_vert_int_juptake_nh4(i,j) = 0.0
       cobalt%wc_vert_int_jprod_nh4(i,j) = 0.0
       cobalt%wc_vert_int_juptake_no3(i,j) = 0.0
       cobalt%wc_vert_int_nfix(i,j) = 0.0
    enddo; enddo !} i,j
    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
          cobalt%wc_vert_int_c(i,j) = cobalt%wc_vert_int_c(i,j) + cobalt%tot_layer_int_c(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_dic(i,j) = cobalt%wc_vert_int_dic(i,j) + cobalt%tot_layer_int_dic(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_doc(i,j) = cobalt%wc_vert_int_doc(i,j) + cobalt%tot_layer_int_doc(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_poc(i,j) = cobalt%wc_vert_int_poc(i,j) + cobalt%tot_layer_int_poc(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_n(i,j) = cobalt%wc_vert_int_n(i,j) + cobalt%tot_layer_int_n(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_p(i,j) = cobalt%wc_vert_int_p(i,j) + cobalt%tot_layer_int_p(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_fe(i,j) = cobalt%wc_vert_int_fe(i,j) + cobalt%tot_layer_int_fe(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_si(i,j) = cobalt%wc_vert_int_si(i,j) + cobalt%tot_layer_int_si(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_o2(i,j) = cobalt%wc_vert_int_o2(i,j) + cobalt%tot_layer_int_o2(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_alk(i,j) = cobalt%wc_vert_int_alk(i,j) + cobalt%tot_layer_int_alk(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_jdiss_sidet(i,j) = cobalt%wc_vert_int_jdiss_sidet(i,j) +                 &
             cobalt%jdiss_sidet(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          cobalt%wc_vert_int_jdiss_cadet(i,j) = cobalt%wc_vert_int_jdiss_cadet(i,j) +                 &
             (cobalt%jdiss_cadet_calc(i,j,k)+cobalt%jdiss_cadet_arag(i,j,k))*rho_dzt(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_jo2resp(i,j) = cobalt%wc_vert_int_jo2resp(i,j) +                         &
             cobalt%jo2resp_wc(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          cobalt%wc_vert_int_jprod_cadet(i,j) = cobalt%wc_vert_int_jprod_cadet(i,j) +                 &
             (cobalt%jprod_cadet_calc(i,j,k)+cobalt%jprod_cadet_arag(i,j,k))*rho_dzt(i,j,k)*grid_tmask(i,j,k)
          cobalt%wc_vert_int_jno3denit(i,j) = cobalt%wc_vert_int_jno3denit(i,j) +                     &
             cobalt%jno3denit_wc(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          cobalt%wc_vert_int_jnitrif(i,j) = cobalt%wc_vert_int_jnitrif(i,j) +                         &
             cobalt%jnitrif(i,j,k) * rho_dzt(i,j,k) * grid_tmask(i,j,k)
          cobalt%wc_vert_int_juptake_nh4(i,j) = cobalt%wc_vert_int_juptake_nh4(i,j) +                     &
             (phyto(1)%juptake_nh4(i,j,k)+phyto(2)%juptake_nh4(i,j,k)+phyto(3)%juptake_nh4(i,j,k))* & 
             rho_dzt(i,j,k) * grid_tmask(i,j,k)

          cobalt%wc_vert_int_jprod_nh4(i,j) = cobalt%wc_vert_int_jprod_nh4(i,j) +                     &
             cobalt%jprod_nh4(i,j,k)*rho_dzt(i,j,k) * grid_tmask(i,j,k)

          cobalt%wc_vert_int_juptake_no3(i,j) = cobalt%wc_vert_int_juptake_no3(i,j) +                     &
             (phyto(1)%juptake_no3(i,j,k)+phyto(2)%juptake_no3(i,j,k)+phyto(3)%juptake_no3(i,j,k))* &
             rho_dzt(i,j,k) * grid_tmask(i,j,k)
          cobalt%wc_vert_int_nfix(i,j) = cobalt%wc_vert_int_nfix(i,j) + phyto(DIAZO)%juptake_n2(i,j,k) *&
             rho_dzt(i,j,k) * grid_tmask(i,j,k)
    enddo; enddo; enddo  !} i,j,k
    !
    !---------------------------------------------------------------------
    ! Add external bottom fluxes to specific rates  
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec ; do k = 1, nk  !{
       ! CAS added calcite and aragonite redisolution terms 
       cobalt%jdiss_cadet_calc_plus_btm(i,j,k)  = cobalt%jdiss_cadet_calc(i,j,k)
       cobalt%jdiss_cadet_arag_plus_btm(i,j,k)  = cobalt%jdiss_cadet_arag(i,j,k)
       ! CAS added a jprod_nh4_plus_btm for remoc CMIP variable
       cobalt%jprod_nh4_plus_btm(i,j,k) = cobalt%jprod_nh4(i,j,k)
       cobalt%jalk_plus_btm(i,j,k)  = cobalt%jalk(i,j,k) 
       cobalt%jdic_plus_btm(i,j,k)  = cobalt%jdic(i,j,k) 
       cobalt%jfed_plus_btm(i,j,k)  = cobalt%jfed(i,j,k) 
       cobalt%jnh4_plus_btm(i,j,k)  = cobalt%jnh4(i,j,k) 
       cobalt%jno3_plus_btm(i,j,k)  = cobalt%jno3(i,j,k) 
       cobalt%jo2_plus_btm(i,j,k)   = cobalt%jo2(i,j,k) 
       cobalt%jpo4_plus_btm(i,j,k)  = cobalt%jpo4(i,j,k) 
       cobalt%jsio4_plus_btm(i,j,k) = cobalt%jsio4(i,j,k) 
       cobalt%jdin_plus_btm(i,j,k)  = cobalt%jno3(i,j,k) + cobalt%jnh4(i,j,k)
    enddo; enddo; enddo  !} i,j,k

    do j = jsc, jec ; do i = isc, iec  !{
       k = grid_kmt(i,j)
       if (k .gt. 0) then !{

          ! CAS added calcite and aragonite redissolution terms
          cobalt%jdiss_cadet_calc_plus_btm(i,j,k)  = cobalt%jdiss_cadet_calc(i,j,k) +  &
             cobalt%fcased_redis(i,j) / rho_dzt(i,j,k)
          cobalt%jdiss_cadet_arag_plus_btm(i,j,k)  = cobalt%jdiss_cadet_arag(i,j,k) +  &
             cobalt%f_cadet_arag_btf(i,j,1) / rho_dzt(i,j,k)

          ! CAS added for remoc calculation
          cobalt%jprod_nh4_plus_btm(i,j,k)  = cobalt%jprod_nh4(i,j,k) + (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) / rho_dzt(i,j,k)


          cobalt%jalk_plus_btm(i,j,k)  = cobalt%jalk(i,j,k) +                       &
            (2.0 * (cobalt%fcased_redis(i,j) + cobalt%f_cadet_arag_btf(i,j,1)) +    &
             cobalt%f_ndet_btf(i,j,1) + cobalt%alk_2_n_denit * cobalt%fno3denit_sed(i,j)) / rho_dzt(i,j,k)
! updated
          cobalt%jdic_plus_btm(i,j,k)  = cobalt%jdic(i,j,k) +                       &
            (cobalt%fcased_redis(i,j) + cobalt%f_cadet_arag_btf(i,j,1) +            &
            ((cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n)) / rho_dzt(i,j,k)

! CAS is ffe_sed a biogenic source or is it similar to coast/atmosphere source?
          cobalt%jfed_plus_btm(i,j,k)  = cobalt%jfed(i,j,k) + (cobalt%ffe_sed(i,j)+cobalt%ffe_geotherm(i,j)) / rho_dzt(i,j,k)
! updated
! CAS: fixed parentheses (commented out old for comparison, think rho_dzt should only divide bottom fluxes)
          cobalt%jnh4_plus_btm(i,j,k)  = cobalt%jnh4(i,j,k) + (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) / rho_dzt(i,j,k)     
!          cobalt%jnh4_plus_btm(i,j,k)  = (cobalt%jnh4(i,j,k) + cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) / rho_dzt(i,j,k)

! NOTE: should fno3denit_sed be SUBTRACTED ?
! CAS: yes, I think so, I've made the change
          cobalt%jno3_plus_btm(i,j,k)  = cobalt%jno3(i,j,k) - cobalt%fno3denit_sed(i,j) / rho_dzt(i,j,k)

          cobalt%jo2_plus_btm(i,j,k)   = cobalt%jo2(i,j,k) +                        &
            (cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnfeso4red_sed(i,j))) / rho_dzt(i,j,k)

! updated
! CAS: fixed parentheses to not include jpo4 as in jnh4 example above
          cobalt%jpo4_plus_btm(i,j,k)  = cobalt%jpo4(i,j,k) + (cobalt%f_pdet_btf(i,j,1) - cobalt%fpdet_burial(i,j)) / rho_dzt(i,j,k)     

          cobalt%jsio4_plus_btm(i,j,k) = cobalt%jsio4(i,j,k) + cobalt%f_sidet_btf(i,j,1) / rho_dzt(i,j,k)   

          cobalt%jdin_plus_btm(i,j,k)  = cobalt%jno3_plus_btm(i,j,k) + cobalt%jnh4_plus_btm(i,j,k) 
          
       endif !}
    enddo; enddo  !} i, j

!
! CHECK3 Remineralization of Organic Carbon, remoc=(jprod_nh4*c_2_n/dht) for k=1,kbot-1 + (jprod_nh4+f_ndet_btf-fndet_burial)*c_2_n/dht for k=kbot
! CAS: will code this up when I address denitrification issue
    do j = jsc, jec ; do i = isc, iec  !{
       kbot = grid_kmt(i,j)
       if (kbot .gt. 0) then !{
       do k = 1, kbot-1  !{
          cobalt%remoc(i,j,k) = cobalt%jprod_nh4(i,j,k) * cobalt%c_2_n / dzt(i,j,k)
       enddo  !} k
       cobalt%remoc(i,j,kbot) = (cobalt%jprod_nh4(i,j,kbot) + cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n / dzt(i,j,kbot)
       endif !}
    enddo; enddo  !} i, j


!****************************************************************************************************

    allocate(rho_dzt_100(isc:iec,jsc:jec))
    !
    !---------------------------------------------------------------------
    ! calculate upper 100 m vertical integrals
    !---------------------------------------------------------------------
    !
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_100(i,j) = rho_dzt(i,j,1)
       cobalt%f_alk_int_100(i,j) = cobalt%p_alk(i,j,1,tau) * rho_dzt(i,j,1)
       cobalt%f_dic_int_100(i,j) = cobalt%p_dic(i,j,1,tau) * rho_dzt(i,j,1)
       cobalt%f_din_int_100(i,j) = (cobalt%p_no3(i,j,1,tau) + cobalt%p_nh4(i,j,1,tau)) * rho_dzt(i,j,1)
       cobalt%f_fed_int_100(i,j) = cobalt%p_fed(i,j,1,tau) * rho_dzt(i,j,1)
       cobalt%f_po4_int_100(i,j) = cobalt%p_po4(i,j,1,tau) * rho_dzt(i,j,1)
       cobalt%f_sio4_int_100(i,j) = cobalt%p_sio4(i,j,1,tau) * rho_dzt(i,j,1)
       cobalt%jalk_100(i,j) = cobalt%jalk(i,j,1) * rho_dzt(i,j,1)
       cobalt%jdic_100(i,j) = cobalt%jdic(i,j,1) * rho_dzt(i,j,1)
       cobalt%jdin_100(i,j) = (cobalt%jno3(i,j,1) + cobalt%jnh4(i,j,1)) * rho_dzt(i,j,1)
       cobalt%jfed_100(i,j) = cobalt%jfed(i,j,1) * rho_dzt(i,j,1)
       cobalt%jpo4_100(i,j) = cobalt%jpo4(i,j,1) * rho_dzt(i,j,1)
       cobalt%jsio4_100(i,j) = cobalt%jsio4(i,j,1) * rho_dzt(i,j,1)
!       cobalt%jprod_ptot_100(i,j) = (phyto(DIAZO)%jprod_po4(i,j,1) + phyto(LARGE)%jprod_po4(i,j,1) + &
!          phyto(SMALL)%jprod_po4(i,j,1)) * rho_dzt(i,j,1)
! previously computed in COBALT
       cobalt%jprod_ptot_100(i,j) = cobalt%jprod_po4(i,j,1) * rho_dzt(i,j,1)
       do n = 1, NUM_PHYTO  !{
          phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jprod_n_new_100(i,j) = phyto(n)%juptake_no3(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n(i,j,1) * rho_dzt(i,j,1)
          phyto(n)%f_n_100(i,j) = phyto(n)%f_n(i,j,1) * rho_dzt(i,j,1)
! added juptake_fe_100 
          phyto(n)%juptake_fe_100(i,j) = phyto(n)%juptake_fe(i,j,1) * rho_dzt(i,j,1)
! CAS: added juptake_po4_100
          phyto(n)%juptake_po4_100(i,j) = phyto(n)%juptake_po4(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n
       phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%juptake_n2(i,j,1) * rho_dzt(i,j,1)
       phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n(i,j,1) * rho_dzt(i,j,1)
       phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n(i,j,1) * rho_dzt(i,j,1)
       phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n(i,j,1) * rho_dzt(i,j,1)
! CAS: added diagnotistic for depth integrated diatom production
       cobalt%jprod_diat_100(i,j) = phyto(LARGE)%jprod_n(i,j,1)*phyto(LARGE)%silim(i,j,1)*rho_dzt(i,j,1)
! added juptake_sio4_100 (large only)
       phyto(LARGE)%juptake_sio4_100(i,j) = phyto(LARGE)%juptake_sio4(i,j,1) * rho_dzt(i,j,1)
       do n = 1, NUM_ZOO  !{
          zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jremin_n_100(i,j) = zoo(n)%jprod_nh4(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%f_n_100(i,j) = zoo(n)%f_n(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n

       do n = 1,2  !{
          zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jprod_don_100(i,j) = (zoo(n)%jprod_ldon(i,j,1) + zoo(n)%jprod_sldon(i,j,1) + &
             zoo(n)%jprod_srdon(i,j,1))  * rho_dzt(i,j,1)
       enddo   !} n

       do n = 2,3  !{
          zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n(i,j,1) * rho_dzt(i,j,1)
          zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet(i,j,1) * rho_dzt(i,j,1)
       enddo   !} n

       cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)
       cobalt%hp_jremin_n_100(i,j) =  cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)*(1.0-cobalt%hp_phi_det)
       cobalt%hp_jprod_ndet_100(i,j) =  cobalt%hp_jingest_n(i,j,1)*rho_dzt(i,j,1)*cobalt%hp_phi_det

       bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n(i,j,1) * rho_dzt(i,j,1)
       bact(1)%jremin_n_100(i,j) = bact(1)%jprod_nh4(i,j,1) * rho_dzt(i,j,1)
       bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon(i,j,1) * rho_dzt(i,j,1)
       bact(1)%f_n_100(i,j) = bact(1)%f_n(i,j,1) * rho_dzt(i,j,1)

       cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc(i,j,1) * rho_dzt(i,j,1)
       cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag(i,j,1) * rho_dzt(i,j,1)
       cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet(i,j,1) * rho_dzt(i,j,1)

       cobalt%f_ndet_100(i,j) = cobalt%f_ndet(i,j,1)*rho_dzt(i,j,1)
       cobalt%f_don_100(i,j) = (cobalt%f_ldon(i,j,1)+cobalt%f_sldon(i,j,1)+cobalt%f_srdon(i,j,1))* &
           rho_dzt(i,j,1)
       cobalt%f_silg_100(i,j) = cobalt%f_silg(i,j,1)*rho_dzt(i,j,1)

       cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,1) * cobalt%Rho_0 * cobalt%wsink
       cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,1) * cobalt%Rho_0 * cobalt%wsink
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; !{
       k_100 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_100(i,j) .lt. cobalt%Rho_0 * 100.0) then 
             k_100 = k
             rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
             cobalt%f_alk_int_100(i,j) = cobalt%f_alk_int_100(i,j) + cobalt%p_alk(i,j,k,tau) * rho_dzt(i,j,k)
             cobalt%f_dic_int_100(i,j) = cobalt%f_dic_int_100(i,j) + cobalt%p_dic(i,j,k,tau) * rho_dzt(i,j,k)
             cobalt%f_din_int_100(i,j) = cobalt%f_din_int_100(i,j) + (cobalt%p_no3(i,j,k,tau) +        &
                cobalt%p_nh4(i,j,k,tau)) * rho_dzt(i,j,k)
             cobalt%f_fed_int_100(i,j) = cobalt%f_fed_int_100(i,j) + cobalt%p_fed(i,j,k,tau) * rho_dzt(i,j,k)
             cobalt%f_po4_int_100(i,j) = cobalt%f_po4_int_100(i,j) + cobalt%p_po4(i,j,k,tau) * rho_dzt(i,j,k)
             cobalt%f_sio4_int_100(i,j) = cobalt%f_sio4_int_100(i,j) + cobalt%p_sio4(i,j,k,tau) *  rho_dzt(i,j,k)
             cobalt%jalk_100(i,j) = cobalt%jalk_100(i,j) + cobalt%jalk(i,j,k) * rho_dzt(i,j,k)
             cobalt%jdic_100(i,j) = cobalt%jdic_100(i,j) + cobalt%jdic(i,j,k) * rho_dzt(i,j,k)
             cobalt%jdin_100(i,j) = cobalt%jdin_100(i,j) + (cobalt%jno3(i,j,k) + cobalt%jnh4(i,j,k)) * rho_dzt(i,j,k)
             cobalt%jfed_100(i,j) = cobalt%jfed_100(i,j) + cobalt%jfed(i,j,k) * rho_dzt(i,j,k)
             cobalt%jpo4_100(i,j) = cobalt%jpo4_100(i,j) + cobalt%jpo4(i,j,k) * rho_dzt(i,j,k)
             cobalt%jsio4_100(i,j) = cobalt%jsio4_100(i,j) + cobalt%jsio4(i,j,k) * rho_dzt(i,j,k)
!            cobalt%jprod_ptot_100(i,j) = cobalt%jprod_ptot_100(i,j) + (phyto(DIAZO)%jprod_po4(i,j,k) &
!               + phyto(LARGE)%jprod_po4(i,j,k) + phyto(SMALL)%jprod_po4(i,j,k)) * rho_dzt(i,j,k)
! previously computed in COBALT
             cobalt%jprod_ptot_100(i,j) = cobalt%jprod_ptot_100(i,j) + cobalt%jprod_po4(i,j,k) * rho_dzt(i,j,k)

             do n = 1, NUM_PHYTO !{
                phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k)* & 
                   rho_dzt(i,j,k)
                phyto(n)%jprod_n_new_100(i,j) = phyto(n)%jprod_n_new_100(i,j) + phyto(n)%juptake_no3(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n_100(i,j) + phyto(n)%jzloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n_100(i,j) + phyto(n)%jexuloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                phyto(n)%f_n_100(i,j) = phyto(n)%f_n_100(i,j) + phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k) 
! added juptake_fe_100 
                phyto(n)%juptake_fe_100(i,j) = phyto(n)%juptake_fe_100(i,j) + phyto(n)%juptake_fe(i,j,k)*rho_dzt(i,j,k) 
! CAS: added juptake_po4_100
                phyto(n)%juptake_po4_100(i,j) = phyto(n)%juptake_po4_100(i,j) + phyto(n)%juptake_po4(i,j,k)*rho_dzt(i,j,k)
             enddo !} n
             phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%jprod_n_n2_100(i,j) + &
                 phyto(DIAZO)%juptake_n2(i,j,k)*rho_dzt(i,j,k)
             phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n_100(i,j) + &
                 phyto(SMALL)%jvirloss_n(i,j,k)*rho_dzt(i,j,k)
             phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n_100(i,j) + &
                 phyto(SMALL)%jaggloss_n(i,j,k)*rho_dzt(i,j,k)
             phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n_100(i,j) + &
                 phyto(LARGE)%jaggloss_n(i,j,k)*rho_dzt(i,j,k)
! CAS: added diagnotistic for depth integrated diatom production
             cobalt%jprod_diat_100(i,j) = cobalt%jprod_diat_100(i,j) + & 
               phyto(LARGE)%jprod_n(i,j,k)*phyto(LARGE)%silim(i,j,k)*rho_dzt(i,j,k)
! added juptake_sio4_100 (large only)
             phyto(LARGE)%juptake_sio4_100(i,j) = phyto(LARGE)%juptake_sio4_100(i,j) + &
                 phyto(LARGE)%juptake_sio4(i,j,k)*rho_dzt(i,j,k)

             do n = 1, NUM_ZOO !{
                zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n_100(i,j) + zoo(n)%jprod_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n_100(i,j) + zoo(n)%jingest_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jremin_n_100(i,j) = zoo(n)%jremin_n_100(i,j) + zoo(n)%jprod_nh4(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%f_n_100(i,j) = zoo(n)%f_n_100(i,j) + zoo(n)%f_n(i,j,k)*rho_dzt(i,j,k)
             enddo !} n

             do n = 1,2 !{
                zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n_100(i,j) + zoo(n)%jzloss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jprod_don_100(i,j) = zoo(n)%jprod_don_100(i,j) + (zoo(n)%jprod_ldon(i,j,k) + &
                   zoo(n)%jprod_sldon(i,j,k) + zoo(n)%jprod_srdon(i,j,k))*rho_dzt(i,j,k)
             enddo !} n

             do n = 2,3 !{
                zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n_100(i,j) + zoo(n)%jhploss_n(i,j,k)* &
                   rho_dzt(i,j,k)
                zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet_100(i,j) + zoo(n)%jprod_ndet(i,j,k)* &
                   rho_dzt(i,j,k)
             enddo !} n

             cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 rho_dzt(i,j,k)
             cobalt%hp_jremin_n_100(i,j) = cobalt%hp_jremin_n_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 (1.0-cobalt%hp_phi_det)*rho_dzt(i,j,k)
             cobalt%hp_jprod_ndet_100(i,j) = cobalt%hp_jprod_ndet_100(i,j) + cobalt%hp_jingest_n(i,j,k)* &
                 cobalt%hp_phi_det*rho_dzt(i,j,k)

             bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n_100(i,j) + bact(1)%jprod_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n_100(i,j) + bact(1)%jzloss_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n_100(i,j) + bact(1)%jvirloss_n(i,j,k) * rho_dzt(i,j,k)
             bact(1)%jremin_n_100(i,j) = bact(1)%jremin_n_100(i,j) + bact(1)%jprod_nh4(i,j,k) * rho_dzt(i,j,k)
             bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon_100(i,j) + bact(1)%juptake_ldon(i,j,k) * rho_dzt(i,j,k)
             bact(1)%f_n_100(i,j) = bact(1)%f_n_100(i,j) + bact(1)%f_n(i,j,k)*rho_dzt(i,j,k)

             cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet_100(i,j) + cobalt%jprod_lithdet(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet_100(i,j) + cobalt%jprod_sidet(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc_100(i,j) + cobalt%jprod_cadet_calc(i,j,k) * rho_dzt(i,j,k)
             cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag_100(i,j) + cobalt%jprod_cadet_arag(i,j,k) * rho_dzt(i,j,k)
             cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet_100(i,j) + cobalt%jremin_ndet(i,j,k) * rho_dzt(i,j,k)
             cobalt%f_ndet_100(i,j) = cobalt%f_ndet_100(i,j) + cobalt%f_ndet(i,j,k)*rho_dzt(i,j,k)
             cobalt%f_don_100(i,j) = cobalt%f_don_100(i,j) + (cobalt%f_ldon(i,j,k) + cobalt%f_sldon(i,j,k) + &
                cobalt%f_srdon(i,j,k))*rho_dzt(i,j,k)
             cobalt%f_silg_100(i,j) = cobalt%f_silg_100(i,j) + cobalt%f_silg(i,j,k)*rho_dzt(i,j,k)

             cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,k) * cobalt%Rho_0 * cobalt%wsink
             cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,k) * cobalt%Rho_0 * cobalt%wsink

          endif
       enddo  !} k

       if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
          drho_dzt = cobalt%Rho_0 * 100.0 - rho_dzt_100(i,j)
          cobalt%f_alk_int_100(i,j) = cobalt%f_alk_int_100(i,j) + cobalt%p_alk(i,j,k_100,tau) * drho_dzt
          cobalt%f_dic_int_100(i,j) = cobalt%f_dic_int_100(i,j) + cobalt%p_dic(i,j,k_100,tau) * drho_dzt
          cobalt%f_din_int_100(i,j) = cobalt%f_din_int_100(i,j) + (cobalt%p_no3(i,j,k_100,tau) +       &
             cobalt%p_nh4(i,j,k_100,tau)) * drho_dzt
          cobalt%f_fed_int_100(i,j) = cobalt%f_fed_int_100(i,j) + cobalt%p_fed(i,j,k_100,tau) * drho_dzt
          cobalt%f_po4_int_100(i,j) = cobalt%f_po4_int_100(i,j) + cobalt%p_po4(i,j,k_100,tau) * drho_dzt
          cobalt%f_sio4_int_100(i,j) = cobalt%f_sio4_int_100(i,j) + cobalt%p_sio4(i,j,k_100,tau) * drho_dzt
          cobalt%jalk_100(i,j) = cobalt%jalk_100(i,j) + cobalt%jalk(i,j,k_100) * drho_dzt
          cobalt%jdic_100(i,j) = cobalt%jdic_100(i,j) + cobalt%jdic(i,j,k_100) * drho_dzt
          cobalt%jdin_100(i,j) = cobalt%jdin_100(i,j) + (cobalt%jno3(i,j,k_100) +  cobalt%jnh4(i,j,k_100)) * drho_dzt
          cobalt%jfed_100(i,j) = cobalt%jfed_100(i,j) + cobalt%jfed(i,j,k_100) * drho_dzt
          cobalt%jpo4_100(i,j) = cobalt%jpo4_100(i,j) + cobalt%jpo4(i,j,k_100) * drho_dzt
          cobalt%jsio4_100(i,j) = cobalt%jsio4_100(i,j) + cobalt%jsio4(i,j,k_100) * drho_dzt
!         cobalt%jprod_ptot_100(i,j) = cobalt%jprod_ptot_100(i,j) +                                   &
!            (phyto(DIAZO)%jprod_po4(i,j,k_100) + phyto(LARGE)%jprod_po4(i,j,k_100) +               &
!            phyto(SMALL)%jprod_po4(i,j,k_100)) * drho_dzt
! previously computed in COBALT
          cobalt%jprod_ptot_100(i,j) = cobalt%jprod_ptot_100(i,j) + cobalt%jprod_po4(i,j,k_100) * drho_dzt

          do n = 1, NUM_PHYTO !{
              phyto(n)%jprod_n_100(i,j) = phyto(n)%jprod_n_100(i,j) + phyto(n)%jprod_n(i,j,k_100)* &
                 drho_dzt
              phyto(n)%jprod_n_new_100(i,j) = phyto(n)%jprod_n_new_100(i,j) + phyto(n)%juptake_no3(i,j,k_100)* &
                 drho_dzt
              phyto(n)%jzloss_n_100(i,j) = phyto(n)%jzloss_n_100(i,j) + phyto(n)%jzloss_n(i,j,k_100)* &
                 drho_dzt
             phyto(n)%jexuloss_n_100(i,j) = phyto(n)%jexuloss_n_100(i,j) + phyto(n)%jexuloss_n(i,j,k_100)* &
                 drho_dzt
              phyto(n)%f_n_100(i,j) = phyto(n)%f_n_100(i,j) + phyto(n)%f_n(i,j,k_100)*drho_dzt
! added juptake_fe_100 
              phyto(n)%juptake_fe_100(i,j) = phyto(n)%juptake_fe_100(i,j) + phyto(n)%juptake_fe(i,j,k_100)*drho_dzt
! CAS: added juptake_po4_100
              phyto(n)%juptake_po4_100(i,j) = phyto(n)%juptake_po4_100(i,j) + phyto(n)%juptake_po4(i,j,k_100)*drho_dzt
           enddo !} n
           phyto(DIAZO)%jprod_n_n2_100(i,j) = phyto(DIAZO)%jprod_n_n2_100(i,j) + &
               phyto(DIAZO)%juptake_n2(i,j,k_100)*drho_dzt
           phyto(SMALL)%jvirloss_n_100(i,j) = phyto(SMALL)%jvirloss_n_100(i,j) + &
               phyto(SMALL)%jvirloss_n(i,j,k_100)*drho_dzt
           phyto(SMALL)%jaggloss_n_100(i,j) = phyto(SMALL)%jaggloss_n_100(i,j) + &
               phyto(SMALL)%jaggloss_n(i,j,k_100)*drho_dzt
           phyto(LARGE)%jaggloss_n_100(i,j) = phyto(LARGE)%jaggloss_n_100(i,j) + &
               phyto(LARGE)%jaggloss_n(i,j,k_100)*drho_dzt
! CAS: added diagnotistic for depth integrated diatom production
           cobalt%jprod_diat_100(i,j) = cobalt%jprod_diat_100(i,j) + &
               phyto(LARGE)%jprod_n(i,j,k_100)*phyto(LARGE)%silim(i,j,k_100)*drho_dzt
! added juptake_sio4_100 (large only)
           phyto(LARGE)%juptake_sio4_100(i,j) = phyto(LARGE)%juptake_sio4_100(i,j) + &
               phyto(LARGE)%juptake_sio4(i,j,k_100)*drho_dzt

           do n = 1, NUM_ZOO !{
               zoo(n)%jprod_n_100(i,j) = zoo(n)%jprod_n_100(i,j) + zoo(n)%jprod_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jingest_n_100(i,j) = zoo(n)%jingest_n_100(i,j) + zoo(n)%jingest_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jremin_n_100(i,j) = zoo(n)%jremin_n_100(i,j) + zoo(n)%jprod_nh4(i,j,k_100)* &
                 drho_dzt
               zoo(n)%f_n_100(i,j) = zoo(n)%f_n_100(i,j) + zoo(n)%f_n(i,j,k_100)*drho_dzt
           enddo !} n

           do n = 1,2 !{
               zoo(n)%jzloss_n_100(i,j) = zoo(n)%jzloss_n_100(i,j) + zoo(n)%jzloss_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jprod_don_100(i,j) = zoo(n)%jprod_don_100(i,j) + (zoo(n)%jprod_ldon(i,j,k_100) + &
                 zoo(n)%jprod_sldon(i,j,k_100) + zoo(n)%jprod_srdon(i,j,k_100))*drho_dzt
           enddo !} n

           do n = 2,3 !{
               zoo(n)%jhploss_n_100(i,j) = zoo(n)%jhploss_n_100(i,j) + zoo(n)%jhploss_n(i,j,k_100)* &
                 drho_dzt
               zoo(n)%jprod_ndet_100(i,j) = zoo(n)%jprod_ndet_100(i,j) + zoo(n)%jprod_ndet(i,j,k_100)* &
                 drho_dzt
           enddo !} n

           cobalt%hp_jingest_n_100(i,j) = cobalt%hp_jingest_n_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               drho_dzt
           cobalt%hp_jremin_n_100(i,j) = cobalt%hp_jremin_n_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               (1.0-cobalt%hp_phi_det)*drho_dzt
           cobalt%hp_jprod_ndet_100(i,j) = cobalt%hp_jprod_ndet_100(i,j) + cobalt%hp_jingest_n(i,j,k_100)* &
               cobalt%hp_phi_det*drho_dzt

           bact(1)%jprod_n_100(i,j) = bact(1)%jprod_n_100(i,j) + bact(1)%jprod_n(i,j,k_100)* &
                drho_dzt
           bact(1)%jzloss_n_100(i,j) = bact(1)%jzloss_n_100(i,j) + bact(1)%jzloss_n(i,j,k_100)* & 
                drho_dzt
           bact(1)%jvirloss_n_100(i,j) = bact(1)%jvirloss_n_100(i,j) + bact(1)%jvirloss_n(i,j,k_100)* & 
                drho_dzt
           bact(1)%jremin_n_100(i,j) = bact(1)%jremin_n_100(i,j) + bact(1)%jprod_nh4(i,j,k_100)* & 
                drho_dzt
           bact(1)%juptake_ldon_100(i,j) = bact(1)%juptake_ldon_100(i,j) + bact(1)%juptake_ldon(i,j,k_100)* &
                drho_dzt
           bact(1)%f_n_100(i,j) = bact(1)%f_n_100(i,j) + bact(1)%f_n(i,j,k_100)*drho_dzt

           cobalt%jprod_lithdet_100(i,j) = cobalt%jprod_lithdet_100(i,j) + cobalt%jprod_lithdet(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_sidet_100(i,j) = cobalt%jprod_sidet_100(i,j) + cobalt%jprod_sidet(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_cadet_calc_100(i,j) = cobalt%jprod_cadet_calc_100(i,j) + cobalt%jprod_cadet_calc(i,j,k_100)* &
                drho_dzt
           cobalt%jprod_cadet_arag_100(i,j) = cobalt%jprod_cadet_arag_100(i,j) + cobalt%jprod_cadet_arag(i,j,k_100)* &
                drho_dzt
           cobalt%jremin_ndet_100(i,j) = cobalt%jremin_ndet_100(i,j) + cobalt%jremin_ndet(i,j,k_100)* &
                drho_dzt

           cobalt%f_ndet_100(i,j) = cobalt%f_ndet_100(i,j) + cobalt%f_ndet(i,j,k_100)*drho_dzt
           cobalt%f_don_100(i,j) = cobalt%f_don_100(i,j) + (cobalt%f_ldon(i,j,k_100) + cobalt%f_sldon(i,j,k_100) + &
              cobalt%f_srdon(i,j,k_100))*drho_dzt
           cobalt%f_silg_100(i,j) = cobalt%f_silg_100(i,j) + cobalt%f_silg(i,j,k_100)*drho_dzt

           cobalt%fndet_100(i,j) = cobalt%f_ndet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fpdet_100(i,j) = cobalt%f_pdet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%ffedet_100(i,j) = cobalt%f_fedet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%flithdet_100(i,j) = cobalt%f_lithdet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fsidet_100(i,j) = cobalt%f_sidet(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fcadet_arag_100(i,j) = cobalt%f_cadet_arag(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
           cobalt%fcadet_calc_100(i,j) = cobalt%f_cadet_calc(i,j,k_100) * cobalt%Rho_0 * cobalt%wsink
       endif

       cobalt%jprod_allphytos_100(i,j) = phyto(SMALL)%jprod_n_100(i,j) + phyto(LARGE)%jprod_n_100(i,j) + &
          phyto(DIAZO)%jprod_n_100(i,j) 
    enddo ; enddo  !} i,j

    !
    ! Calculate biomass-weighted nutrient and light limitation terms for CMIP6  
    ! Note: these need to be done after other 100m integrals because the biomass
    ! weighting requires the pre-calculated 100m biomass
    !
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_100(i,j) = rho_dzt(i,j,1)
       do n = 1,NUM_PHYTO
          phyto(n)%nlim_bw_100(i,j) = (phyto(n)%no3lim(i,j,1)+phyto(n)%nh4lim(i,j,1))* &
                phyto(n)%f_n(i,j,1)*rho_dzt(i,j,1)/(phyto(n)%f_n_100(i,j)+epsln)
          phyto(n)%plim_bw_100(i,j) = phyto(n)%po4lim(i,j,1)* &
                phyto(n)%f_n(i,j,1)*rho_dzt(i,j,1)/(phyto(n)%f_n_100(i,j)+epsln)
          phyto(n)%def_fe_bw_100(i,j) = phyto(n)%def_fe(i,j,1)* &
                phyto(n)%f_n(i,j,1)*rho_dzt(i,j,1)/(phyto(n)%f_n_100(i,j)+epsln)
          phyto(n)%irrlim_bw_100(i,j) = phyto(n)%irrlim(i,j,1)* &
                phyto(n)%f_n(i,j,1)*rho_dzt(i,j,1)/(phyto(n)%f_n_100(i,j)+epsln)
       enddo   !} n
    enddo; enddo  !} i, j

    do j = jsc, jec ; do i = isc, iec ; !{
       k_100 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_100(i,j) .lt. cobalt%Rho_0 * 100.0) then 
             k_100 = k
             rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
             do n = 1,NUM_PHYTO
                phyto(n)%nlim_bw_100(i,j) = phyto(n)%nlim_bw_100(i,j) + & 
                   (phyto(n)%no3lim(i,j,k)+phyto(n)%nh4lim(i,j,k))* &
                   phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k)/phyto(n)%f_n_100(i,j)
                phyto(n)%plim_bw_100(i,j) = phyto(n)%plim_bw_100(i,j) + phyto(n)%po4lim(i,j,k)* &
                   phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k)/phyto(n)%f_n_100(i,j)
                phyto(n)%def_fe_bw_100(i,j) = phyto(n)%def_fe_bw_100(i,j) + phyto(n)%def_fe(i,j,k)* &
                   phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k)/phyto(n)%f_n_100(i,j)
                phyto(n)%irrlim_bw_100(i,j) = phyto(n)%irrlim_bw_100(i,j) + phyto(n)%irrlim(i,j,k)* &
                   phyto(n)%f_n(i,j,k)*rho_dzt(i,j,k)/phyto(n)%f_n_100(i,j)
             enddo
          endif
       enddo  !} k

       if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
          drho_dzt = cobalt%Rho_0 * 100.0 - rho_dzt_100(i,j)
          do n = 1,NUM_PHYTO
             phyto(n)%nlim_bw_100(i,j) = phyto(n)%nlim_bw_100(i,j) + &
                (phyto(n)%no3lim(i,j,k_100)+phyto(n)%nh4lim(i,j,k_100))* &
                phyto(n)%f_n(i,j,k_100)*drho_dzt/phyto(n)%f_n_100(i,j)
             phyto(n)%plim_bw_100(i,j) = phyto(n)%plim_bw_100(i,j) + phyto(n)%po4lim(i,j,k_100)* &
                phyto(n)%f_n(i,j,k_100)*drho_dzt/phyto(n)%f_n_100(i,j)
             phyto(n)%def_fe_bw_100(i,j) = phyto(n)%def_fe_bw_100(i,j) + phyto(n)%def_fe(i,j,k_100)* &
                phyto(n)%f_n(i,j,k_100)*drho_dzt/phyto(n)%f_n_100(i,j)
             phyto(n)%irrlim_bw_100(i,j) = phyto(n)%irrlim_bw_100(i,j) + phyto(n)%irrlim(i,j,k_100)* &
                phyto(n)%f_n(i,j,k_100)*drho_dzt/phyto(n)%f_n_100(i,j)
          enddo
        endif
    enddo; enddo  !} i, j
    deallocate(rho_dzt_100)

    do j = jsc, jec ; do i = isc, iec ; !{
      if (grid_kmt(i,j) .gt. 0) then !{
         cobalt%btm_temp(i,j) = TEMP(i,j,grid_kmt(i,j))
         cobalt%btm_o2(i,j) = cobalt%f_o2(i,j,grid_kmt(i,j))
         cobalt%btm_htotal(i,j) = cobalt%f_htotal(i,j,grid_kmt(i,j))
         cobalt%btm_co3_sol_arag(i,j) = cobalt%co3_sol_arag(i,j,grid_kmt(i,j))
         cobalt%btm_co3_sol_calc(i,j) = cobalt%co3_sol_calc(i,j,grid_kmt(i,j))
         cobalt%btm_co3_ion(i,j) = cobalt%f_co3_ion(i,j,grid_kmt(i,j))
         cobalt%cased_2d(i,j) = cobalt%f_cased(i,j,1)      
      endif
    enddo; enddo  !} i, j

    !
    !---------------------------------------------------------------------
    ! calculate upper 200m vertical integrals for mesozooplankton
    ! quantities for comparison with COPEPOD database
    !---------------------------------------------------------------------
    !
    allocate(rho_dzt_200(isc:iec,jsc:jec))
    do j = jsc, jec ; do i = isc, iec !{
       rho_dzt_200(i,j) = rho_dzt(i,j,1)
       cobalt%jprod_mesozoo_200(i,j) = (zoo(2)%jprod_n(i,j,1) + zoo(3)%jprod_n(i,j,1))*rho_dzt(i,j,1)
       cobalt%f_mesozoo_200(i,j) = (zoo(2)%f_n(i,j,1)+zoo(3)%f_n(i,j,1))*rho_dzt(i,j,1)
    enddo; enddo !} i,j

    do j = jsc, jec ; do i = isc, iec ; !{
       k_200 = 1
       do k = 2, grid_kmt(i,j)  !{
          if (rho_dzt_200(i,j) .lt. cobalt%Rho_0 * 200.0) then
             k_200 = k
             rho_dzt_200(i,j) = rho_dzt_200(i,j) + rho_dzt(i,j,k)
             cobalt%jprod_mesozoo_200(i,j) = cobalt%jprod_mesozoo_200(i,j) + &
                (zoo(2)%jprod_n(i,j,k) + zoo(3)%jprod_n(i,j,k))*rho_dzt(i,j,k)
             cobalt%f_mesozoo_200(i,j) = cobalt%f_mesozoo_200(i,j) + &
                (zoo(2)%f_n(i,j,k)+zoo(3)%f_n(i,j,k))*rho_dzt(i,j,k)
          endif
       enddo  !} k

       if (k_200 .gt. 1 .and. k_200 .lt. grid_kmt(i,j)) then
          drho_dzt = cobalt%Rho_0 * 200.0 - rho_dzt_200(i,j)
          cobalt%jprod_mesozoo_200(i,j) = cobalt%jprod_mesozoo_200(i,j) + &
              (zoo(2)%jprod_n(i,j,k_200) + zoo(3)%jprod_n(i,j,k_200))*drho_dzt
          cobalt%f_mesozoo_200(i,j) = cobalt%f_mesozoo_200(i,j) + &
              (zoo(2)%f_n(i,j,k_200)+zoo(3)%f_n(i,j,k_200))*drho_dzt
       endif
    enddo ; enddo  !} i,j

    call g_tracer_get_values(tracer_list,'alk','runoff_tracer_flux',cobalt%runoff_flux_alk,isd,jsd)
    call g_tracer_get_values(tracer_list,'dic','runoff_tracer_flux',cobalt%runoff_flux_dic,isd,jsd)
    if (do_14c) then  !{
      call g_tracer_get_values(tracer_list,'di14c','runoff_tracer_flux',cobalt%runoff_flux_di14c,isd,jsd)
    endif  !}
    call g_tracer_get_values(tracer_list,'fed','runoff_tracer_flux',cobalt%runoff_flux_fed,isd,jsd)
    call g_tracer_get_values(tracer_list,'fed','drydep',cobalt%dry_fed,isd,jsd)
    call g_tracer_get_values(tracer_list,'fed','wetdep',cobalt%wet_fed,isd,jsd)
    call g_tracer_get_values(tracer_list,'lith','runoff_tracer_flux',cobalt%runoff_flux_lith,isd,jsd)
    call g_tracer_get_values(tracer_list,'lith','drydep',cobalt%dry_lith,isd,jsd)
    call g_tracer_get_values(tracer_list,'lith','wetdep',cobalt%wet_lith,isd,jsd)
    call g_tracer_get_values(tracer_list,'no3','runoff_tracer_flux',cobalt%runoff_flux_no3,isd,jsd)
    call g_tracer_get_values(tracer_list,'no3','drydep',cobalt%dry_no3,isd,jsd)
    call g_tracer_get_values(tracer_list,'no3','wetdep',cobalt%wet_no3,isd,jsd)
    call g_tracer_get_values(tracer_list,'nh4','drydep',cobalt%dry_nh4,isd,jsd)
    call g_tracer_get_values(tracer_list,'nh4','wetdep',cobalt%wet_nh4,isd,jsd)
    call g_tracer_get_values(tracer_list,'po4','drydep',cobalt%dry_po4,isd,jsd)
    call g_tracer_get_values(tracer_list,'po4','wetdep',cobalt%wet_po4,isd,jsd)
    call g_tracer_get_values(tracer_list,'ldon','runoff_tracer_flux',cobalt%runoff_flux_ldon,isd,jsd)
    call g_tracer_get_values(tracer_list,'sldon','runoff_tracer_flux',cobalt%runoff_flux_sldon,isd,jsd)
    call g_tracer_get_values(tracer_list,'srdon','runoff_tracer_flux',cobalt%runoff_flux_srdon,isd,jsd)
    call g_tracer_get_values(tracer_list,'ndet','runoff_tracer_flux',cobalt%runoff_flux_ndet,isd,jsd)
    call g_tracer_get_values(tracer_list,'po4','runoff_tracer_flux',cobalt%runoff_flux_po4,isd,jsd)
    call g_tracer_get_values(tracer_list,'ldop','runoff_tracer_flux',cobalt%runoff_flux_ldop,isd,jsd)
    call g_tracer_get_values(tracer_list,'sldop','runoff_tracer_flux',cobalt%runoff_flux_sldop,isd,jsd)
    call g_tracer_get_values(tracer_list,'srdop','runoff_tracer_flux',cobalt%runoff_flux_srdop,isd,jsd)
! JGJ: Added for CMIP6
    call g_tracer_get_values(tracer_list,'dic','stf_gas',cobalt%stf_gas_dic,isd,jsd)
    call g_tracer_get_values(tracer_list,'o2','stf_gas',cobalt%stf_gas_o2,isd,jsd)
    call g_tracer_get_values(tracer_list,'dic','deltap',cobalt%deltap_dic,isd,jsd)
    call g_tracer_get_values(tracer_list,'o2','deltap',cobalt%deltap_o2,isd,jsd)


!---------------------------------------------------------------------
! Add vertical integrals for diagnostics
!---------------------------------------------------------------------
!

    call mpp_clock_end(id_clock_cobalt_calc_diagnostics)
    call mpp_clock_begin(id_clock_cobalt_send_diagnostics)

    if (cobalt%id_pka_nh3 .gt. 0) then
       used = g_send_data(cobalt%id_pka_nh3,  pka_nh3,   &
            model_time, rmask = grid_tmask(:,:,1),&
            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    end if
    if (allocated(pka_nh3)) deallocate(pka_nh3)
    deallocate(phos_nh3_exchange)
!
!---------------------------------------------------------------------
!
! Send phytoplankton diagnostic data

    do n= 1, NUM_PHYTO
       if (phyto(n)%id_def_fe .gt. 0)          &
            used = g_send_data(phyto(n)%id_def_fe,     phyto(n)%def_fe,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_felim .gt. 0)           &
            used = g_send_data(phyto(n)%id_felim,      phyto(n)%felim,            &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_irrlim .gt. 0)          &
            used = g_send_data(phyto(n)%id_irrlim,     phyto(n)%irrlim,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jzloss_n .gt. 0)          &
            used = g_send_data(phyto(n)%id_jzloss_n, phyto(n)%jzloss_n*rho_dzt,      &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jaggloss_n .gt. 0)          &
            used = g_send_data(phyto(n)%id_jaggloss_n, phyto(n)%jaggloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jvirloss_n .gt. 0)          &
            used = g_send_data(phyto(n)%id_jvirloss_n, phyto(n)%jvirloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jexuloss_n .gt. 0)          &
            used = g_send_data(phyto(n)%id_jexuloss_n, phyto(n)%jexuloss_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_fe .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_fe, phyto(n)%juptake_fe*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_nh4 .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_nh4, phyto(n)%juptake_nh4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_no3 .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_no3, phyto(n)%juptake_no3*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_po4 .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_po4, phyto(n)%juptake_po4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_sio4 .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_sio4, phyto(n)%juptake_sio4*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_juptake_n2 .gt. 0)          &
            used = g_send_data(phyto(n)%id_juptake_n2, phyto(n)%juptake_n2*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_jprod_n .gt. 0)          &
            used = g_send_data(phyto(n)%id_jprod_n, phyto(n)%jprod_n*rho_dzt,   &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_liebig_lim .gt. 0)      &
            used = g_send_data(phyto(n)%id_liebig_lim,phyto(n)%liebig_lim,          &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_mu .gt. 0)              &
            used = g_send_data(phyto(n)%id_mu,        phyto(n)%mu,                  &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_nh4lim .gt. 0)          &
            used = g_send_data(phyto(n)%id_nh4lim,     phyto(n)%nh4lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_no3lim .gt. 0)          &
            used = g_send_data(phyto(n)%id_no3lim,     phyto(n)%no3lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_po4lim .gt. 0)          &
            used = g_send_data(phyto(n)%id_po4lim,     phyto(n)%po4lim,             &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_o2lim .gt. 0)          &
            used = g_send_data(phyto(n)%id_o2lim,     phyto(n)%o2lim,             &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_fe_2_n .gt. 0)        &
            used = g_send_data(phyto(n)%id_q_fe_2_n,   phyto(n)%q_fe_2_n,           &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_silim .gt. 0)     &
            used = g_send_data(phyto(n)%id_silim, phyto(n)%silim,       &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_q_si_2_n .gt. 0)     &
            used = g_send_data(phyto(n)%id_q_si_2_n, phyto(n)%q_si_2_n,       &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_theta .gt. 0)           &
            used = g_send_data(phyto(n)%id_theta,      phyto(n)%theta,              &
            model_time, rmask = grid_tmask,& 
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_f_mu_mem .gt. 0)           &
            used = g_send_data(phyto(n)%id_f_mu_mem,      phyto(n)%f_mu_mem,              &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_mu_mix .gt. 0)           &
            used = g_send_data(phyto(n)%id_mu_mix,      phyto(n)%mu_mix,              &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (phyto(n)%id_agg_lim .gt. 0)           &
            used = g_send_data(phyto(n)%id_agg_lim,      phyto(n)%agg_lim,              &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
    !--------------------------------------------------------------------------------------
    ! Send bacterial diagnostic data
    !

    if (bact(1)%id_jzloss_n .gt. 0)          &
       used = g_send_data(bact(1)%id_jzloss_n, bact(1)%jzloss_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jvirloss_n .gt. 0)          &
       used = g_send_data(bact(1)%id_jvirloss_n, bact(1)%jvirloss_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_juptake_ldon .gt. 0)          &
       used = g_send_data(bact(1)%id_juptake_ldon, bact(1)%juptake_ldon*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_juptake_ldop .gt. 0)          &
       used = g_send_data(bact(1)%id_juptake_ldop, bact(1)%juptake_ldop*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_nh4 .gt. 0)          &
       used = g_send_data(bact(1)%id_jprod_nh4, bact(1)%jprod_nh4*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_po4 .gt. 0)          &
       used = g_send_data(bact(1)%id_jprod_po4, bact(1)%jprod_po4*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_jprod_n .gt. 0)          &
       used = g_send_data(bact(1)%id_jprod_n, bact(1)%jprod_n*rho_dzt,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_o2lim .gt. 0)          &
       used = g_send_data(bact(1)%id_o2lim, bact(1)%o2lim,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_ldonlim .gt. 0)          &
       used = g_send_data(bact(1)%id_ldonlim, bact(1)%ldonlim,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (bact(1)%id_temp_lim .gt. 0)          &
       used = g_send_data(bact(1)%id_temp_lim, bact(1)%temp_lim,           &
       model_time, rmask = grid_tmask,&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    !--------------------------------------------------------------------------------------
    ! Send zooplankton diagnostic data
    !

    do n= 1, NUM_ZOO
       if (zoo(n)%id_jzloss_n .gt. 0)          &
            used = g_send_data(zoo(n)%id_jzloss_n, zoo(n)%jzloss_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jhploss_n .gt. 0)          &
            used = g_send_data(zoo(n)%id_jhploss_n, zoo(n)%jhploss_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_n .gt. 0)          &
            used = g_send_data(zoo(n)%id_jingest_n, zoo(n)%jingest_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_p .gt. 0)          &
            used = g_send_data(zoo(n)%id_jingest_p, zoo(n)%jingest_p*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_sio2 .gt. 0)          &
            used = g_send_data(zoo(n)%id_jingest_sio2, zoo(n)%jingest_sio2*rho_dzt,      &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jingest_fe .gt. 0)          &
            used = g_send_data(zoo(n)%id_jingest_fe, zoo(n)%jingest_fe*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ndet .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_ndet, zoo(n)%jprod_ndet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_pdet .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_pdet, zoo(n)%jprod_pdet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ldon .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_ldon, zoo(n)%jprod_ldon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_ldop .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_ldop, zoo(n)%jprod_ldop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sldon .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_sldon, zoo(n)%jprod_sldon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sldop .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_sldop, zoo(n)%jprod_sldop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
      if (zoo(n)%id_jprod_srdon .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_srdon, zoo(n)%jprod_srdon*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_srdop .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_srdop, zoo(n)%jprod_srdop*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_fed .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_fed,  zoo(n)%jprod_fed*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_fedet .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_fedet, zoo(n)%jprod_fedet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sidet .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_sidet, zoo(n)%jprod_sidet*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_sio4 .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_sio4, zoo(n)%jprod_sio4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_po4 .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_po4,  zoo(n)%jprod_po4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_nh4 .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_nh4,  zoo(n)%jprod_nh4*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_jprod_n .gt. 0)          &
            used = g_send_data(zoo(n)%id_jprod_n,   zoo(n)%jprod_n*rho_dzt,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_o2lim .gt. 0)          &
            used = g_send_data(zoo(n)%id_o2lim, zoo(n)%o2lim,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
       if (zoo(n)%id_temp_lim .gt. 0)          &
            used = g_send_data(zoo(n)%id_temp_lim, zoo(n)%temp_lim,           &
            model_time, rmask = grid_tmask,&
            is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    enddo
    !
    ! Production diagnostics
    !
    if (cobalt%id_jprod_cadet_arag .gt. 0)    &
       used = g_send_data(cobalt%id_jprod_cadet_arag, cobalt%jprod_cadet_arag * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_cadet_calc .gt. 0)    &
       used = g_send_data(cobalt%id_jprod_cadet_calc, cobalt%jprod_cadet_calc * rho_dzt, &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_ndet .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_ndet, cobalt%jprod_ndet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_pdet .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_pdet, cobalt%jprod_pdet*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_srdon .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_srdon, cobalt%jprod_srdon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sldon .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_sldon, cobalt%jprod_sldon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_ldon .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_ldon, cobalt%jprod_ldon*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_srdop .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_srdop, cobalt%jprod_srdop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_sldop .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_sldop, cobalt%jprod_sldop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !if (cobalt%id_jprod_ldop .gt. 0)          &
    !    used = g_send_data(cobalt%id_jprod_ldop, cobalt%jprod_ldop*rho_dzt,           &
    !    model_time, rmask = grid_tmask,&
    !    is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_nh4 .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_nh4, cobalt%jprod_nh4*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
   if (cobalt%id_jprod_nh4_plus_btm .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_nh4_plus_btm, cobalt%jprod_nh4_plus_btm*rho_dzt, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_po4 .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_po4, cobalt%jprod_po4*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_fed .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_fed, cobalt%jprod_fed*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_fedet .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_fedet,  cobalt%jprod_fedet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_sidet .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_sidet, cobalt%jprod_sidet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_sio4 .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_sio4, cobalt%jprod_sio4*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jprod_lithdet .gt. 0)          &
        used = g_send_data(cobalt%id_jprod_lithdet, cobalt%jprod_lithdet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_cadet_arag .gt. 0)          &
        used = g_send_data(cobalt%id_jdiss_cadet_arag, cobalt%jdiss_cadet_arag*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_cadet_calc .gt. 0)          &
        used = g_send_data(cobalt%id_jdiss_cadet_calc, cobalt%jdiss_cadet_calc*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdiss_sidet .gt. 0)          &
        used = g_send_data(cobalt%id_jdiss_sidet, cobalt%jdiss_sidet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_ndet .gt. 0)          &
        used = g_send_data(cobalt%id_jremin_ndet, cobalt%jremin_ndet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_pdet .gt. 0)          &
        used = g_send_data(cobalt%id_jremin_pdet, cobalt%jremin_pdet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jremin_fedet .gt. 0)          &
        used = g_send_data(cobalt%id_jremin_fedet, cobalt%jremin_fedet*rho_dzt,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jfed .gt. 0)              &
         used = g_send_data(cobalt%id_jfed,       cobalt%jfed*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jfe_ads .gt. 0)              &
         used = g_send_data(cobalt%id_jfe_ads,       cobalt%jfe_ads*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jfe_coast .gt. 0)            &
         used = g_send_data(cobalt%id_jfe_coast, cobalt%jfe_coast*rho_dzt,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_kfe_eq_lig .gt. 0)            &
         used = g_send_data(cobalt%id_kfe_eq_lig, log10(cobalt%kfe_eq_lig+epsln),         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_feprime .gt. 0)            &
         used = g_send_data(cobalt%id_feprime, cobalt%feprime,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_ligand .gt. 0)            &
         used = g_send_data(cobalt%id_ligand, cobalt%ligand,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fe_sol .gt. 0)            &
         used = g_send_data(cobalt%id_fe_sol, cobalt%fe_sol,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_expkT .gt. 0)              &
         used = g_send_data(cobalt%id_expkT,       cobalt%expkT,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_expkreminT .gt. 0)              &
         used = g_send_data(cobalt%id_expkreminT,       cobalt%expkreminT,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_hp_o2lim .gt. 0)            &
         used = g_send_data(cobalt%id_hp_o2lim, cobalt%hp_o2lim,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_hp_temp_lim .gt. 0)            &
         used = g_send_data(cobalt%id_hp_temp_lim, cobalt%hp_temp_lim,         &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_irr_inst .gt. 0)           &
         used = g_send_data(cobalt%id_irr_inst,      cobalt%irr_inst,              &
         model_time, rmask = grid_tmask(:,:,:),&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)             
    if (cobalt%id_irr_mix .gt. 0)           &
         used = g_send_data(cobalt%id_irr_mix,       cobalt%irr_mix,               &
         model_time, rmask = grid_tmask(:,:,:),&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jno3denit_wc .gt. 0)            &
         used = g_send_data(cobalt%id_jno3denit_wc,  cobalt%jno3denit_wc*rho_dzt,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
   if (cobalt%id_jo2resp_wc .gt. 0)            &
         used = g_send_data(cobalt%id_jo2resp_wc,  cobalt%jo2resp_wc*rho_dzt,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jnitrif .gt. 0)              &
         used = g_send_data(cobalt%id_jnitrif,       cobalt%jnitrif*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_c .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_c, cobalt%tot_layer_int_c,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_fe .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_fe,cobalt%tot_layer_int_fe,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_n .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_n,cobalt%tot_layer_int_n,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_p .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_p,cobalt%tot_layer_int_p,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_si .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_si,cobalt%tot_layer_int_si,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_o2 .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_o2,cobalt%tot_layer_int_o2,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_tot_layer_int_alk .gt. 0)  &
         used = g_send_data(cobalt%id_tot_layer_int_alk,cobalt%tot_layer_int_alk,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_total_filter_feeding .gt. 0)  &
         used = g_send_data(cobalt%id_total_filter_feeding,cobalt%total_filter_feeding,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_nlg_diatoms.gt. 0)  &
         used = g_send_data(cobalt%id_nlg_diatoms,cobalt%nlg_diatoms,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_q_si_2_n_lg_diatoms.gt. 0)  &
         used = g_send_data(cobalt%id_q_si_2_n_lg_diatoms,cobalt%q_si_2_n_lg_diatoms,&
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_co2_csurf .gt. 0)             &
         used = g_send_data(cobalt%id_co2_csurf,      cobalt%co2_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_pco2_csurf .gt. 0)             &
         used = g_send_data(cobalt%id_pco2_csurf,      cobalt%pco2_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_co2_alpha .gt. 0)             &
         used = g_send_data(cobalt%id_co2_alpha,      cobalt%co2_alpha,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_nh3_csurf .gt. 0)             &
         used = g_send_data(cobalt%id_nh3_csurf,      cobalt%nh3_csurf,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_nh3_alpha .gt. 0)             &
         used = g_send_data(cobalt%id_nh3_alpha,      cobalt%nh3_alpha,              &
         model_time, rmask = grid_tmask(:,:,1),& 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_arag_btm .gt. 0)           &  
         used = g_send_data(cobalt%id_fcadet_arag_btm,   cobalt%fcadet_arag_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_calc_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fcadet_calc_btm,   cobalt%fcadet_calc_btm,      &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_ffedet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_ffedet_btm,   cobalt%ffedet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fndet_btm .gt. 0)            &
         used = g_send_data(cobalt%id_fndet_btm,    cobalt%fndet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_btm .gt. 0)            &
         used = g_send_data(cobalt%id_fpdet_btm,    cobalt%fpdet_btm,              &
         model_time, rmask = grid_tmask(:,:,1),&  
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fsidet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_fsidet_btm,   cobalt%fsidet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_flithdet_btm .gt. 0)           &
         used = g_send_data(cobalt%id_flithdet_btm,   cobalt%flithdet_btm,             &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_burial .gt. 0)        &
         used = g_send_data(cobalt%id_fcased_burial, cobalt%fcased_burial,         &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_redis .gt. 0)         &
         used = g_send_data(cobalt%id_fcased_redis,  cobalt%fcased_redis,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fcased_redis_surfresp .gt. 0)         &
         used = g_send_data(cobalt%id_fcased_redis_surfresp,cobalt%fcased_redis_surfresp, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_cased_redis_coef .gt. 0)         &
         used = g_send_data(cobalt%id_cased_redis_coef,  cobalt%cased_redis_coef, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_cased_redis_delz .gt. 0)         &
         used = g_send_data(cobalt%id_cased_redis_delz,  cobalt%cased_redis_delz, &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_ffe_sed .gt. 0)              &
         used = g_send_data(cobalt%id_ffe_sed,       cobalt%ffe_sed,               &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_ffe_geotherm .gt. 0)              &
         used = g_send_data(cobalt%id_ffe_geotherm,  cobalt%ffe_geotherm,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_ffe_iceberg .gt. 0)              &
         used = g_send_data(cobalt%id_ffe_iceberg,  cobalt%ffe_iceberg,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fnfeso4red_sed .gt. 0)           &
         used = g_send_data(cobalt%id_fnfeso4red_sed,cobalt%fnfeso4red_sed,        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fno3denit_sed .gt. 0)           &
         used = g_send_data(cobalt%id_fno3denit_sed, cobalt%fno3denit_sed,         &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fnoxic_sed .gt. 0)           &
         used = g_send_data(cobalt%id_fnoxic_sed,    cobalt%fnoxic_sed,            &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_frac_burial .gt. 0)           &
         used = g_send_data(cobalt%id_frac_burial,    cobalt%frac_burial,          &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fndet_burial .gt. 0)           &
         used = g_send_data(cobalt%id_fndet_burial,    cobalt%fndet_burial,        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_burial .gt. 0)           &
         used = g_send_data(cobalt%id_fpdet_burial,    cobalt%fpdet_burial,            &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_co3_sol_arag .gt. 0)        &
       used = g_send_data(cobalt%id_co3_sol_arag,  cobalt%co3_sol_arag,             &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_co3_sol_calc .gt. 0)        &
       used = g_send_data(cobalt%id_co3_sol_calc,  cobalt%co3_sol_calc,             &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_omega_arag .gt. 0)          &
       used = g_send_data(cobalt%id_omega_arag,  cobalt%omega_arag,                 &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_omega_calc .gt. 0)          &
       used = g_send_data(cobalt%id_omega_calc,  cobalt%omega_calc,                 &
       model_time, rmask = grid_tmask(:,:,:),&
       is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fcadet_arag .gt. 0)               &
         used = g_send_data(cobalt%id_fcadet_arag, cobalt%p_cadet_arag(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:), &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fcadet_calc .gt. 0)               &
         used = g_send_data(cobalt%id_fcadet_calc, cobalt%p_cadet_calc(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink*grid_tmask(:,:,:), &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_ffedet .gt. 0)                &
         used = g_send_data(cobalt%id_ffedet,        cobalt%p_fedet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_flithdet .gt. 0)                &
         used = g_send_data(cobalt%id_flithdet,      cobalt%p_lithdet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fndet .gt. 0)                 &
         used = g_send_data(cobalt%id_fndet,         cobalt%p_ndet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fpdet .gt. 0)                 &
         used = g_send_data(cobalt%id_fpdet,         cobalt%p_pdet(:,:,:,tau) * cobalt%Rho_0 * &
         cobalt%wsink * grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_fsidet .gt. 0)                &
         used = g_send_data(cobalt%id_fsidet,        cobalt%p_sidet(:,:,:,tau)  * cobalt%Rho_0 * &
         cobalt%wsink  *grid_tmask(:,:,:),&
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_nphyto_tot .gt. 0)                &
         used = g_send_data(cobalt%id_nphyto_tot,   (cobalt%p_ndi(:,:,:,tau) +  &
         cobalt%p_nlg(:,:,:,tau) + cobalt%p_nsm(:,:,:,tau)), &
         model_time, rmask = grid_tmask(:,:,:),& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
! Radiocarbon fields
!
    if (cobalt%id_b_di14c .gt. 0)                                                 &
         used = g_send_data(cobalt%id_b_di14c,        cobalt%b_di14c,                &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_c14_2_n .gt. 0)                                                 &
         used = g_send_data(cobalt%id_c14_2_n,        cobalt%c14_2_n,                &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_c14o2_csurf .gt. 0)                                             &
         used = g_send_data(cobalt%id_c14o2_csurf,    cobalt%c14o2_csurf,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_c14o2_alpha .gt. 0)                                             &
         used = g_send_data(cobalt%id_c14o2_alpha,    cobalt%c14o2_alpha,            &
         model_time, rmask = grid_tmask(:,:,1),                                  & 
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_fpo14c .gt. 0)                                                  &
         used = g_send_data(cobalt%id_fpo14c,         cobalt%fpo14c,                 &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_j14c_decay_dic .gt. 0)                                          &
         used = g_send_data(cobalt%id_j14c_decay_dic, cobalt%j14c_decay_dic*rho_dzt, &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_j14c_decay_doc .gt. 0)                                          &
         used = g_send_data(cobalt%id_j14c_decay_doc, cobalt%j14c_decay_doc*rho_dzt, &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_j14c_reminp .gt. 0)                                             &
         used = g_send_data(cobalt%id_j14c_reminp,    cobalt%j14c_reminp*rho_dzt,    &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdi14c .gt. 0)                                                  &
         used = g_send_data(cobalt%id_jdi14c,         cobalt%jdi14c*rho_dzt,         &
         model_time, rmask = grid_tmask,& 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    if (cobalt%id_jdo14c .gt. 0)                                                  &
         used = g_send_data(cobalt%id_jdo14c,         cobalt%jdo14c*rho_dzt,         &
         model_time, rmask = grid_tmask,                                         & 
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
    !
    ! 2D COBALT fields
    !
   if (cobalt%id_pco2surf .gt. 0)              &
         used = g_send_data(cobalt%id_pco2surf,      cobalt%pco2_csurf,           &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_pnh3surf .gt. 0)              &
         used = g_send_data(cobalt%id_pnh3surf,      cobalt%pnh3_csurf,           &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_alk .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_alk,       cobalt%p_alk(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_cadet_arag .gt. 0)      &
       used = g_send_data(cobalt%id_sfc_cadet_arag,cobalt%p_cadet_arag(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_cadet_calc .gt. 0)      &
       used = g_send_data(cobalt%id_sfc_cadet_calc,cobalt%p_cadet_calc(:,:,1,tau),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_dic .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_dic,       cobalt%p_dic(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_fed .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_fed,       cobalt%p_fed(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_ldon .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_ldon,       cobalt%p_ldon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_sldon .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_sldon,       cobalt%p_sldon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_srdon .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_srdon,       cobalt%p_srdon(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_no3 .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_no3,       cobalt%p_no3(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_nh4 .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_nh4,       cobalt%p_nh4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_po4 .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_po4,       cobalt%p_po4(:,:,1,tau),         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_sio4 .gt. 0)            &
       used = g_send_data(cobalt%id_sfc_sio4,      cobalt%p_sio4(:,:,1,tau),        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_htotal .gt. 0)          &
       used = g_send_data(cobalt%id_sfc_htotal,    cobalt%f_htotal(:,:,1),          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_o2 .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_o2,       cobalt%p_o2(:,:,1,tau),           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_chl .gt. 0)             &
       used = g_send_data(cobalt%id_sfc_chl,       cobalt%f_chl(:,:,1),             &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_irr .gt. 0)              &
         used = g_send_data(cobalt%id_sfc_irr,      cobalt%irr_inst(:,:,1),        &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
    if (cobalt%id_sfc_irr_mem .gt. 0)          &
         used = g_send_data(cobalt%id_sfc_irr_mem,  cobalt%f_irr_mem(:,:,1),       &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_temp .gt. 0)            &
       used = g_send_data(cobalt%id_sfc_temp,      Temp(:,:,1),                    &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_temp .gt. 0)            &
       used = g_send_data(cobalt%id_btm_temp,      cobalt%btm_temp,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_o2 .gt. 0)            &
       used = g_send_data(cobalt%id_btm_o2,      cobalt%btm_o2,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_htotal .gt. 0)            &
       used = g_send_data(cobalt%id_btm_htotal,      cobalt%btm_htotal,        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_cased_2d .gt. 0)            &
       used = g_send_data(cobalt%id_cased_2d,      cobalt%cased_2d,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_ion .gt. 0)            &
       used = g_send_data(cobalt%id_sfc_co3_ion, cobalt%f_co3_ion(:,:,1),            &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_sol_arag .gt. 0)            &
       used = g_send_data(cobalt%id_sfc_co3_sol_arag, cobalt%co3_sol_arag(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_sfc_co3_sol_calc .gt. 0)            &
       used = g_send_data(cobalt%id_sfc_co3_sol_calc, cobalt%co3_sol_calc(:,:,1),  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_co3_ion .gt. 0)            &
       used = g_send_data(cobalt%id_btm_co3_ion, cobalt%btm_co3_ion,            &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_co3_sol_arag .gt. 0)            &
       used = g_send_data(cobalt%id_btm_co3_sol_arag, cobalt%btm_co3_sol_arag,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
   if (cobalt%id_btm_co3_sol_calc .gt. 0)            &
       used = g_send_data(cobalt%id_btm_co3_sol_calc, cobalt%btm_co3_sol_calc,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)

    do n= 1, NUM_PHYTO
       if (phyto(n)%id_sfc_f_n .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_f_n, phyto(n)%f_n(:,:,1),            & 
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_chl .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_chl,    cobalt%c_2_n * 12.0e6 *      &
          phyto(n)%theta(:,:,1) * phyto(n)%f_n(:,:,1),                          &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_def_fe .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_def_fe, phyto(n)%def_fe(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (phyto(n)%id_sfc_felim .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_felim, phyto(n)%felim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_sfc_irrlim .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_irrlim, phyto(n)%irrlim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_theta .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_theta, phyto(n)%theta(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_mu .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_mu, phyto(n)%mu(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(n)%id_sfc_po4lim .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_po4lim, phyto(n)%po4lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(n)%id_sfc_q_fe_2_n .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_q_fe_2_n, phyto(n)%q_fe_2_n(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(n)%id_sfc_nh4lim .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_nh4lim, phyto(n)%nh4lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (phyto(n)%id_sfc_no3lim .gt. 0)              &
          used = g_send_data(phyto(n)%id_sfc_no3lim, phyto(n)%no3lim(:,:,1),      &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo

    ! 
    ! Save river, depositon and bulk elemental fluxes
    !
    if (cobalt%id_dep_dry_fed .gt. 0)     &
       used = g_send_data(cobalt%id_dep_dry_fed, cobalt%dry_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_lith .gt. 0)     &
       used = g_send_data(cobalt%id_dep_dry_lith, cobalt%dry_lith,                        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_nh4 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_dry_nh4, cobalt%dry_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_no3 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_dry_no3, cobalt%dry_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_dry_po4 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_dry_po4, cobalt%dry_po4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_fed .gt. 0)     &
       used = g_send_data(cobalt%id_dep_wet_fed, cobalt%wet_fed,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_lith .gt. 0)     &
       used = g_send_data(cobalt%id_dep_wet_lith, cobalt%wet_lith,                        &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_nh4 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_wet_nh4, cobalt%wet_nh4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_no3 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_wet_no3, cobalt%wet_no3,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_dep_wet_po4 .gt. 0)     &
       used = g_send_data(cobalt%id_dep_wet_po4, cobalt%wet_po4,                          &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_alk .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_alk, cobalt%runoff_flux_alk,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_dic .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_dic, cobalt%runoff_flux_dic,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_runoff_flux_di14c .gt. 0)     &
        used = g_send_data(cobalt%id_runoff_flux_di14c, cobalt%runoff_flux_di14c,           &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_fed .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_fed, cobalt%runoff_flux_fed,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_lith .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_lith, cobalt%runoff_flux_lith,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_no3 .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_no3, cobalt%runoff_flux_no3,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ldon .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_ldon, cobalt%runoff_flux_ldon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_sldon .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_sldon, cobalt%runoff_flux_sldon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_srdon .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_srdon, cobalt%runoff_flux_srdon,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ndet .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_ndet, cobalt%runoff_flux_ndet,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_po4 .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_po4, cobalt%runoff_flux_po4,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_ldop .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_ldop, cobalt%runoff_flux_ldop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_sldop .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_sldop, cobalt%runoff_flux_sldop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_runoff_flux_srdop .gt. 0)     &
       used = g_send_data(cobalt%id_runoff_flux_srdop, cobalt%runoff_flux_srdop,           &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    ! Save 100m integral fluxes
    !
    if (cobalt%id_jprod_allphytos_100 .gt. 0)     &
       used = g_send_data(cobalt%id_jprod_allphytos_100, cobalt%jprod_allphytos_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_jprod_diat_100 .gt. 0)     &
       used = g_send_data(cobalt%id_jprod_diat_100, cobalt%jprod_diat_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    do n= 1, NUM_PHYTO  !{
       if (phyto(n)%id_jprod_n_100 .gt. 0)     &
          used = g_send_data(phyto(n)%id_jprod_n_100, phyto(n)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jprod_n_new_100 .gt. 0)     &
          used = g_send_data(phyto(n)%id_jprod_n_new_100, phyto(n)%jprod_n_new_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jzloss_n_100 .gt. 0)     &
          used = g_send_data(phyto(n)%id_jzloss_n_100, phyto(n)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_jexuloss_n_100 .gt. 0)     &
          used = g_send_data(phyto(n)%id_jexuloss_n_100, phyto(n)%jexuloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (phyto(n)%id_f_n_100 .gt. 0)     &
          used = g_send_data(phyto(n)%id_f_n_100, phyto(n)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    enddo !} n
    if (phyto(DIAZO)%id_jprod_n_n2_100 .gt. 0)     &
       used = g_send_data(phyto(DIAZO)%id_jprod_n_n2_100, phyto(DIAZO)%jprod_n_n2_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (phyto(SMALL)%id_jvirloss_n_100 .gt. 0)     &
       used = g_send_data(phyto(SMALL)%id_jvirloss_n_100, phyto(SMALL)%jvirloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   if (phyto(SMALL)%id_jaggloss_n_100 .gt. 0)     &
       used = g_send_data(phyto(SMALL)%id_jaggloss_n_100, phyto(SMALL)%jaggloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   if (phyto(LARGE)%id_jaggloss_n_100 .gt. 0)     &
       used = g_send_data(phyto(LARGE)%id_jaggloss_n_100, phyto(LARGE)%jaggloss_n_100,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     do n= 1, NUM_ZOO  !{
       if (zoo(n)%id_jprod_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jprod_n_100, zoo(n)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jingest_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jingest_n_100, zoo(n)%jingest_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jremin_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jremin_n_100, zoo(n)%jremin_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if (zoo(n)%id_f_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_f_n_100, zoo(n)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     do n= 1,2  !{
       if (zoo(n)%id_jzloss_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jzloss_n_100, zoo(n)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jprod_don_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jprod_don_100, zoo(n)%jprod_don_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     do n= 2,3  !{
       if (zoo(n)%id_jhploss_n_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jhploss_n_100, zoo(n)%jhploss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
       if (zoo(n)%id_jprod_ndet_100 .gt. 0)     &
          used = g_send_data(zoo(n)%id_jprod_ndet_100, zoo(n)%jprod_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     enddo !} n

     if (cobalt%id_hp_jingest_n_100 .gt. 0)     &
        used = g_send_data(cobalt%id_hp_jingest_n_100, cobalt%hp_jingest_n_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_hp_jremin_n_100 .gt. 0)     &
        used = g_send_data(cobalt%id_hp_jremin_n_100, cobalt%hp_jremin_n_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_hp_jprod_ndet_100 .gt. 0)     &
        used = g_send_data(cobalt%id_hp_jprod_ndet_100, cobalt%hp_jprod_ndet_100,         &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (bact(1)%id_jprod_n_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_jprod_n_100, bact(1)%jprod_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jzloss_n_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_jzloss_n_100, bact(1)%jzloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jvirloss_n_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_jvirloss_n_100, bact(1)%jvirloss_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_jremin_n_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_jremin_n_100, bact(1)%jremin_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_juptake_ldon_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_juptake_ldon_100, bact(1)%juptake_ldon_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (bact(1)%id_f_n_100 .gt. 0)     &
          used = g_send_data(bact(1)%id_f_n_100, bact(1)%f_n_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_jprod_lithdet_100 .gt. 0)     &
          used = g_send_data(cobalt%id_jprod_lithdet_100, cobalt%jprod_lithdet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_sidet_100 .gt. 0)     &
          used = g_send_data(cobalt%id_jprod_sidet_100, cobalt%jprod_sidet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_cadet_calc_100 .gt. 0)     &
          used = g_send_data(cobalt%id_jprod_cadet_calc_100, cobalt%jprod_cadet_calc_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_cadet_arag_100 .gt. 0)     &
          used = g_send_data(cobalt%id_jprod_cadet_arag_100, cobalt%jprod_cadet_arag_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jremin_ndet_100 .gt. 0)     &
          used = g_send_data(cobalt%id_jremin_ndet_100, cobalt%jremin_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     if (cobalt%id_jprod_mesozoo_200 .gt. 0)     &
          used = g_send_data(cobalt%id_jprod_mesozoo_200, cobalt%jprod_mesozoo_200,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_f_ndet_100 .gt. 0)     &
          used = g_send_data(cobalt%id_f_ndet_100, cobalt%f_ndet_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_don_100 .gt. 0)     &
          used = g_send_data(cobalt%id_f_don_100, cobalt%f_don_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_silg_100 .gt. 0)     &
          used = g_send_data(cobalt%id_f_silg_100, cobalt%f_silg_100,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_f_mesozoo_200 .gt. 0)     &
          used = g_send_data(cobalt%id_f_mesozoo_200, cobalt%f_mesozoo_200,         &
          model_time, rmask = grid_tmask(:,:,1),&
          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_fndet_100 .gt. 0)           &
       used = g_send_data(cobalt%id_fndet_100,     cobalt%fndet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fpdet_100 .gt. 0)           &
       used = g_send_data(cobalt%id_fpdet_100,     cobalt%fpdet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fsidet_100 .gt. 0)           &
       used = g_send_data(cobalt%id_fsidet_100,     cobalt%fsidet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_flithdet_100 .gt. 0)           &
       used = g_send_data(cobalt%id_flithdet_100,     cobalt%flithdet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_calc_100 .gt. 0)           &
       used = g_send_data(cobalt%id_fcadet_calc_100,     cobalt%fcadet_calc_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_fcadet_arag_100 .gt. 0)           &
       used = g_send_data(cobalt%id_fcadet_arag_100,     cobalt%fcadet_arag_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_ffedet_100 .gt. 0)           &
       used = g_send_data(cobalt%id_ffedet_100,     cobalt%ffedet_100,                &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    !---------------------------------------------------------------------
    ! Save water column vertical integrals
    !---------------------------------------------------------------------
    !
    if (cobalt%id_wc_vert_int_c .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_c,    cobalt%wc_vert_int_c,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_dic .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_dic,    cobalt%wc_vert_int_dic,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_doc .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_doc,    cobalt%wc_vert_int_doc,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_poc .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_poc,    cobalt%wc_vert_int_poc,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_n .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_n,    cobalt%wc_vert_int_n,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_p .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_p,    cobalt%wc_vert_int_p,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_fe .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_fe,    cobalt%wc_vert_int_fe,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_si .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_si,    cobalt%wc_vert_int_si,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_o2 .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_o2,    cobalt%wc_vert_int_o2,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_alk .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_alk,    cobalt%wc_vert_int_alk,         &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jdiss_sidet .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jdiss_sidet,    cobalt%wc_vert_int_jdiss_sidet,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jdiss_cadet .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jdiss_cadet,    cobalt%wc_vert_int_jdiss_cadet,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jo2resp .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jo2resp,    cobalt%wc_vert_int_jo2resp,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jprod_cadet .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jprod_cadet,    cobalt%wc_vert_int_jprod_cadet,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jno3denit .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jno3denit,    cobalt%wc_vert_int_jno3denit,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jnitrif .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jnitrif,    cobalt%wc_vert_int_jnitrif,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_juptake_nh4 .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_juptake_nh4,    cobalt%wc_vert_int_juptake_nh4,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_jprod_nh4 .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_jprod_nh4,    cobalt%wc_vert_int_jprod_nh4,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_juptake_no3 .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_juptake_no3,    cobalt%wc_vert_int_juptake_no3,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_wc_vert_int_nfix .gt. 0)       &
       used = g_send_data(cobalt%id_wc_vert_int_nfix,    cobalt%wc_vert_int_nfix,  &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    !
    !---------------------------------------------------------------------
    ! Save CaCO3 saturation and O2 minimum depths
    !---------------------------------------------------------------------
    !
    if (cobalt%id_o2min .gt. 0)               &
       used = g_send_data(cobalt%id_o2min,         cobalt%o2min,                    &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_o2min .gt. 0)             &
       used = g_send_data(cobalt%id_z_o2min,    cobalt%z_o2min,                     &
       model_time, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_sat_arag .gt. 0)          &
       used = g_send_data(cobalt%id_z_sat_arag,    cobalt%z_sat_arag,               &
       model_time, mask = cobalt%mask_z_sat_arag, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    if (cobalt%id_z_sat_calc .gt. 0)          &
       used = g_send_data(cobalt%id_z_sat_calc,    cobalt%z_sat_calc,               &
       model_time, mask = cobalt%mask_z_sat_calc, rmask = grid_tmask(:,:,1),&
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    !
    !---------------------------------------------------------------------
    ! Send additional diagnostics  jgj 2015/10/26
    !---------------------------------------------------------------------
    !

    if (cobalt%id_jalk .gt. 0)              &
         used = g_send_data(cobalt%id_jalk, cobalt%jalk*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jalk_plus_btm .gt. 0)              &
         used = g_send_data(cobalt%id_jalk_plus_btm, cobalt%jalk_plus_btm*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jdic .gt. 0)              &
         used = g_send_data(cobalt%id_jdic, cobalt%jdic*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jdic_plus_btm .gt. 0)              &
         used = g_send_data(cobalt%id_jdic_plus_btm, cobalt%jdic_plus_btm*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jnh4 .gt. 0)              &
         used = g_send_data(cobalt%id_jnh4, cobalt%jnh4*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

     if (cobalt%id_jndet .gt. 0)              &
          used = g_send_data(cobalt%id_jndet, cobalt%jndet*rho_dzt,       &
          model_time, rmask = grid_tmask,&
          is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jnh4_plus_btm .gt. 0)              &
         used = g_send_data(cobalt%id_jnh4_plus_btm, cobalt%jnh4_plus_btm*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_jo2_plus_btm .gt. 0)              &
         used = g_send_data(cobalt%id_jo2_plus_btm, cobalt%jo2_plus_btm*rho_dzt,       &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!==============================================================================================================
!  2016/07/05 jgj  send temperature as a test

    if (cobalt%id_thetao .gt. 0)            &
        used = g_send_data(cobalt%id_thetao,  Temp,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Oyr and Omon: 3-D Marine Biogeochemical Tracer Fields
!
    if (cobalt%id_dissic .gt. 0)            &
        used = g_send_data(cobalt%id_dissic,  cobalt%p_dic(:,:,:,tau) * cobalt%Rho_0,           &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING:
!    if (cobalt%id_dissicnat .gt. 0)            &
!        used = g_send_data(cobalt%id_dissicnat,                                                    &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!    if (cobalt%id_dissicabio .gt. 0)            
!        used = g_send_data(cobalt%id_dissicabio,                                                    &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!    if (cobalt%id_dissi14cabio .gt. 0)            
!        used = g_send_data(cobalt%id_dissi14cabio,                                                    &
!         model_time, rmask = grid_tmask,&
!         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CHECK3:
    ! CAS comment on spreadsheet implies that this is the explicitly represented
    ! pool and suggests that we shouldn't add the background
    cobalt%dissoc(:,:,:) = cobalt%doc_background +                                                     &
        cobalt%c_2_n * (cobalt%p_ldon(:,:,:,tau) + cobalt%p_sldon(:,:,:,tau) + cobalt%p_srdon(:,:,:,tau) )

    if (cobalt%id_dissoc .gt. 0)            &
        used = g_send_data(cobalt%id_dissoc,  cobalt%dissoc * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phyc .gt. 0)            &
        used = g_send_data(cobalt%id_phyc,  (cobalt%p_nlg(:,:,:,tau) + cobalt%p_nsm(:,:,:,tau) +  &
        cobalt%p_ndi(:,:,:,tau)) * cobalt%c_2_n * cobalt%Rho_0, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_zooc .gt. 0)            &
        used = g_send_data(cobalt%id_zooc,  (cobalt%p_nlgz(:,:,:,tau) + cobalt%p_nsmz(:,:,:,tau) +  &
        cobalt%p_nmdz(:,:,:,tau)) * cobalt%c_2_n * cobalt%Rho_0, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bacc .gt. 0)            &
        used = g_send_data(cobalt%id_bacc,  cobalt%p_nbact(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_detoc .gt. 0)            &
        used = g_send_data(cobalt%id_detoc,  cobalt%p_ndet(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_calc .gt. 0)            &
        used = g_send_data(cobalt%id_calc,  cobalt%p_cadet_calc(:,:,:,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_arag .gt. 0)            &
        used = g_send_data(cobalt%id_arag,  cobalt%p_cadet_arag(:,:,:,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phydiat.gt. 0)  &
         used = g_send_data(cobalt%id_phydiat,  cobalt%nlg_diatoms * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phydiaz.gt. 0)  &
         used = g_send_data(cobalt%id_phydiaz,  cobalt%p_ndi(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phypico.gt. 0)  &
         used = g_send_data(cobalt%id_phypico,  cobalt%p_nsm(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phymisc.gt. 0)  &
         used = g_send_data(cobalt%id_phymisc,  (cobalt%p_nlg(:,:,:,tau)-cobalt%nlg_diatoms) * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_zmicro.gt. 0)  &
         used = g_send_data(cobalt%id_zmicro,  cobalt%p_nsmz(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_zmeso.gt. 0)  &
         used = g_send_data(cobalt%id_zmeso,  (cobalt%p_nlgz(:,:,:,tau)+cobalt%p_nmdz(:,:,:,tau)) * cobalt%c_2_n * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask,&
         is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_talk .gt. 0)            &
        used = g_send_data(cobalt%id_talk,  cobalt%p_alk(:,:,:,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING:
!    if (cobalt%id_talknat .gt. 0)            &
!        used = g_send_data(cobalt%id_talknat,                     
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CHECK3: this is using ntau=1
    if (cobalt%id_ph .gt. 0)            &
        used = g_send_data(cobalt%id_ph,  log10(cobalt%f_htotal+epsln) * (-1.0),       &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING:
!    if (cobalt%id_phnat .gt. 0)            &
!        used = g_send_data(cobalt%id_phnat,  log10(cobalt%f_htotal+epsln) * -1.0,       &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING: 
!    if (cobalt%id_phabio .gt. 0)            &
!        used = g_send_data(cobalt%id_phabio,  log10(cobalt%f_htotal+epsln) * -1.0,       &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_o2_cmip.gt. 0)            &
        used = g_send_data(cobalt%id_o2_cmip,  cobalt%p_o2(:,:,:,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_o2sat .gt. 0)            &
        used = g_send_data(cobalt%id_o2sat,  cobalt%o2sat, & 
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_no3_cmip.gt. 0)            &
        used = g_send_data(cobalt%id_no3_cmip,  cobalt%p_no3(:,:,:,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_nh4_cmip.gt. 0)            &
        used = g_send_data(cobalt%id_nh4_cmip,  cobalt%p_nh4(:,:,:,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_po4_cmip.gt. 0)            &
        used = g_send_data(cobalt%id_po4_cmip,  cobalt%p_po4(:,:,:,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_dfe .gt. 0)            &
        used = g_send_data(cobalt%id_dfe,  cobalt%p_fed(:,:,:,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_si .gt. 0)            &
        used = g_send_data(cobalt%id_si,  cobalt%p_sio4(:,:,:,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_chl_cmip.gt. 0)            &
        used = g_send_data(cobalt%id_chl_cmip,  cobalt%f_chl * cobalt%Rho_0 / 1e9,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_chldiat .gt. 0)            &
        used = g_send_data(cobalt%id_chldiat,  phyto(LARGE)%theta * cobalt%nlg_diatoms * cobalt%c_2_n * cobalt%Rho_0 * 12e-3,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_chldiaz .gt. 0)            &
        used = g_send_data(cobalt%id_chldiaz,  phyto(DIAZO)%theta * cobalt%p_ndi(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0 * 12e-3,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_chlpico .gt. 0)            &
        used = g_send_data(cobalt%id_chlpico,  phyto(SMALL)%theta * cobalt%p_nsm(:,:,:,tau) * cobalt%c_2_n * cobalt%Rho_0 * 12e-3,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_chlmisc .gt. 0)            &
        used = g_send_data(cobalt%id_chlmisc,  phyto(LARGE)%theta * (cobalt%p_nlg(:,:,:,tau)-cobalt%nlg_diatoms) *  &
        cobalt%c_2_n * cobalt%Rho_0 * 12e-3,   &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_poc .gt. 0)            &
        used = g_send_data(cobalt%id_poc,  (cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) + &
        cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) +  cobalt%p_ndet(:,:,:,tau) + &
        cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + cobalt%p_nlgz(:,:,:,tau)) * cobalt%Rho_0 * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_pon .gt. 0)            &
        used = g_send_data(cobalt%id_pon,  (cobalt%p_ndi(:,:,:,tau) + cobalt%p_nlg(:,:,:,tau) + &
        cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) +  cobalt%p_ndet(:,:,:,tau) + &
        cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + cobalt%p_nlgz(:,:,:,tau)) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: added bacteria and more general accomodation of static but different p_2_n ratios
    if (cobalt%id_pop .gt. 0)            &
        used = g_send_data(cobalt%id_pop,  (phyto(DIAZO)%p_2_n_static * cobalt%p_ndi(:,:,:,tau) + &
        phyto(LARGE)%p_2_n_static * cobalt%p_nlg(:,:,:,tau) + phyto(SMALL)%p_2_n_static * cobalt%p_nsm(:,:,:,tau) + &
        cobalt%p_pdet(:,:,:,tau) + zoo(1)%q_p_2_n * cobalt%p_nsmz(:,:,:,tau) + zoo(2)%q_p_2_n * cobalt%p_nmdz(:,:,:,tau) + &
        zoo(3)%q_p_2_n * cobalt%p_nlgz(:,:,:,tau) + bact(1)%q_p_2_n * cobalt%p_nbact(:,:,:,tau)) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: old code, delete once satisfied with new code above
!    if (cobalt%id_pop .gt. 0)            &
!        used = g_send_data(cobalt%id_pop,  ((phyto(DIAZO)%p_2_n_static * cobalt%p_ndi(:,:,:,tau)) + cobalt%p_nlg(:,:,:,tau) + &
!        cobalt%p_nsm(:,:,:,tau) + cobalt%p_nbact(:,:,:,tau) +  cobalt%p_pdet(:,:,:,tau) + &
!        cobalt%p_nsmz(:,:,:,tau) + cobalt%p_nmdz(:,:,:,tau) + cobalt%p_nlgz(:,:,:,tau)) * cobalt%Rho_0,  &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bfe .gt. 0)            &
        used = g_send_data(cobalt%id_bfe,  (cobalt%p_fedi(:,:,:,tau) + cobalt%p_felg(:,:,:,tau) + cobalt%p_fesm(:,:,:,tau) + & 
        cobalt%p_fedet(:,:,:,tau))  * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bsi .gt. 0)            &
        used = g_send_data(cobalt%id_bsi,  (cobalt%p_silg(:,:,:,tau) + cobalt%p_sidet(:,:,:,tau))  * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phyn .gt. 0)            &
        used = g_send_data(cobalt%id_phyn,  (cobalt%p_nlg(:,:,:,tau) + cobalt%p_nsm(:,:,:,tau) +  &
        cobalt%p_ndi(:,:,:,tau)) * cobalt%Rho_0, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: added p_2_n ratios for other phyto groups
    if (cobalt%id_phyp .gt. 0)            &
        used = g_send_data(cobalt%id_phyp,  (phyto(DIAZO)%p_2_n_static * cobalt%p_ndi(:,:,:,tau) + &
        phyto(LARGE)%p_2_n_static * cobalt%p_nlg(:,:,:,tau) + phyto(SMALL)%p_2_n_static * cobalt%p_nsm(:,:,:,tau) )* &
        cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_phyfe .gt. 0)            &
        used = g_send_data(cobalt%id_phyfe,  (cobalt%p_fedi(:,:,:,tau) + cobalt%p_felg(:,:,:,tau) +  &
        cobalt%p_fesm(:,:,:,tau)) * cobalt%Rho_0, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_physi .gt. 0)            &
        used = g_send_data(cobalt%id_physi,  cobalt%p_silg(:,:,:,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_co3 .gt. 0)            &
        used = g_send_data(cobalt%id_co3,  cobalt%f_co3_ion * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING:
!    if (cobalt%id_co3nat .gt. 0)            &
!        used = g_send_data(cobalt%id_co3nat,  cobalt%f_co3_ion * cobalt%Rho_0,  &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! PENDING:
!    if (cobalt%id_co3abio .gt. 0)            &
!        used = g_send_data(cobalt%id_co3abio,  cobalt%f_co3_ion * cobalt%Rho_0,  &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_co3satcalc .gt. 0)            &
        used = g_send_data(cobalt%id_co3satcalc,  cobalt%co3_sol_calc * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_co3satarag .gt. 0)            &
        used = g_send_data(cobalt%id_co3satarag,  cobalt%co3_sol_arag * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Oyr: Marine Biogeochemical 3-D Fields: Rates of Production and Removal
!
! CHECK3: using dzt for layer thickness
! Maybe just use cobalt%Rho_0 instead of rho_dzt / dzt in production terms
!
! also in Omon
    if (cobalt%id_pp .gt. 0)            &
        used = g_send_data(cobalt%id_pp,  (phyto(DIAZO)%jprod_n +  phyto(LARGE)%jprod_n +  &
        phyto(SMALL)%jprod_n) * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! also in Omon
    if (cobalt%id_pnitrate .gt. 0)            &
        used = g_send_data(cobalt%id_pnitrate,  (phyto(DIAZO)%juptake_no3 +  phyto(LARGE)%juptake_no3 +  &
        phyto(SMALL)%juptake_no3) * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! Not requested
!    if (cobalt%id_pphosphate .gt. 0)            &
!        used = g_send_data(cobalt%id_pphosphate,  (phyto(DIAZO)%juptake_po4 +  phyto(LARGE)%juptake_po4 +  &
!        phyto(SMALL)%juptake_po4) * rho_dzt * cobalt%c_2_n / dzt,  &
!        model_time, rmask = grid_tmask,&
!        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! also in Omon
    if (cobalt%id_pbfe .gt. 0)            &
        used = g_send_data(cobalt%id_pbfe,  (phyto(DIAZO)%juptake_fe +  phyto(LARGE)%juptake_fe +  &
        phyto(SMALL)%juptake_fe) * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! also in Omon
    if (cobalt%id_pbsi .gt. 0)            &
        used = g_send_data(cobalt%id_pbsi,  phyto(LARGE)%juptake_sio4 * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_pcalc .gt. 0)            &
        used = g_send_data(cobalt%id_pcalc,  cobalt%jprod_cadet_calc * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_parag .gt. 0)            &
        used = g_send_data(cobalt%id_parag,  cobalt%jprod_cadet_arag * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! also in Omon
    if (cobalt%id_expc .gt. 0)            &
        used = g_send_data(cobalt%id_expc,  cobalt%p_ndet(:,:,:,tau) * cobalt%Rho_0 *cobalt%wsink * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_expn .gt. 0)            &
        used = g_send_data(cobalt%id_expn,  cobalt%p_ndet(:,:,:,tau) * cobalt%Rho_0 *cobalt%wsink,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_expp .gt. 0)            &
        used = g_send_data(cobalt%id_expp,  cobalt%p_pdet(:,:,:,tau) * cobalt%Rho_0 *cobalt%wsink, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_expfe .gt. 0)            &
        used = g_send_data(cobalt%id_expfe,  cobalt%p_fedet(:,:,:,tau) * cobalt%Rho_0 *cobalt%wsink, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_expsi .gt. 0)            &
        used = g_send_data(cobalt%id_expsi,  cobalt%p_sidet(:,:,:,tau) * cobalt%Rho_0 *cobalt%wsink, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_expcalc .gt. 0)            &
        used = g_send_data(cobalt%id_expcalc,  cobalt%p_cadet_calc(:,:,:,tau) * cobalt%Rho_0 * cobalt%wsink,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_exparag .gt. 0)            &
        used = g_send_data(cobalt%id_exparag,  cobalt%p_cadet_arag(:,:,:,tau) * cobalt%Rho_0 * cobalt%wsink,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: added jprod_nh4_plus_btm to calculate
    if (cobalt%id_remoc .gt. 0)            &
        used = g_send_data(cobalt%id_remoc,  cobalt%jprod_nh4_plus_btm*cobalt%c_2_n*rho_dzt/dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: added redissolution from sediment
    if (cobalt%id_dcalc .gt. 0)            &
        used = g_send_data(cobalt%id_dcalc,  cobalt%jdiss_cadet_calc_plus_btm*rho_dzt/dzt, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: added redissolution from sediment
    if (cobalt%id_darag .gt. 0)            &
        used = g_send_data(cobalt%id_darag,  cobalt%jdiss_cadet_arag_plus_btm*rho_dzt/dzt, &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: fixed unit conversion on production from all groups by adding *rho_dzt,
    if (cobalt%id_ppdiat .gt. 0)            &
        used = g_send_data(cobalt%id_ppdiat,  phyto(LARGE)%jprod_n * phyto(LARGE)%silim * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_ppdiaz .gt. 0)            &
        used = g_send_data(cobalt%id_ppdiaz,  phyto(DIAZO)%jprod_n * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_pppico .gt. 0)            &
        used = g_send_data(cobalt%id_pppico,  phyto(SMALL)%jprod_n * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_ppmisc .gt. 0)            &
        used = g_send_data(cobalt%id_ppmisc,  (phyto(LARGE)%jprod_n * (1 - phyto(LARGE)%silim)) * rho_dzt * cobalt%c_2_n / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CHECK3 _z regridding 
! CAS fixed conversion for all bddt terms 
    if (cobalt%id_bddtdic .gt. 0)            &
        used = g_send_data(cobalt%id_bddtdic,  cobalt%jdic_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bddtdin .gt. 0)            &
        used = g_send_data(cobalt%id_bddtdin,  cobalt%jdin_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bddtdip .gt. 0)            &
        used = g_send_data(cobalt%id_bddtdip,  cobalt%jpo4_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bddtdife .gt. 0)            &
        used = g_send_data(cobalt%id_bddtdife,  cobalt%jfed_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bddtdisi .gt. 0)            &
        used = g_send_data(cobalt%id_bddtdisi,  cobalt%jsio4_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

    if (cobalt%id_bddtalk .gt. 0)            &
        used = g_send_data(cobalt%id_bddtalk,  cobalt%jalk_plus_btm * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: fixed conversion
    if (cobalt%id_fescav .gt. 0)            &
        used = g_send_data(cobalt%id_fescav,  cobalt%jfe_ads * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! CAS: fixed conversion
    if (cobalt%id_fediss .gt. 0)            &
        used = g_send_data(cobalt%id_fediss,  cobalt%jremin_fedet * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

! also in Omon
! CAS: fixed conversion
    if (cobalt%id_graz .gt. 0)            &
        used = g_send_data(cobalt%id_graz,  (phyto(DIAZO)%jzloss_n +  phyto(LARGE)%jzloss_n +  &
        phyto(SMALL)%jzloss_n) * cobalt%c_2_n  * rho_dzt / dzt,  &
        model_time, rmask = grid_tmask,&
        is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Omon: Marine Biogeochemical 2-D Surface Fields 
!  Identical to Oyr 3-D Tracer fields but for surface only

    if (cobalt%id_dissicos .gt. 0)            &
        used = g_send_data(cobalt%id_dissicos,  cobalt%p_dic(:,:,1,tau) * cobalt%Rho_0,           &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_dissicnatos .gt. 0)            &
!        used = g_send_data(cobalt%id_dissicnatos,                                                    &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_dissicabioos .gt. 0)            
!        used = g_send_data(cobalt%id_dissicabioos,                                                    &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_dissi14cabioos .gt. 0)            
!        used = g_send_data(cobalt%id_dissi14cabioos,                                                    &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_dissocos .gt. 0)            &
        used = g_send_data(cobalt%id_dissocos,  cobalt%dissoc(:,:,1) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phycos .gt. 0)            &
        used = g_send_data(cobalt%id_phycos,  (cobalt%p_nlg(:,:,1,tau) + cobalt%p_nsm(:,:,1,tau) +  &
        cobalt%p_ndi(:,:,1,tau)) * cobalt%c_2_n * cobalt%Rho_0, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_zoocos .gt. 0)            &
        used = g_send_data(cobalt%id_zoocos,  (cobalt%p_nlgz(:,:,1,tau) + cobalt%p_nsmz(:,:,1,tau) +  &
        cobalt%p_nmdz(:,:,1,tau)) * cobalt%c_2_n * cobalt%Rho_0, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_baccos .gt. 0)            &
        used = g_send_data(cobalt%id_baccos,  cobalt%p_nbact(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_detocos .gt. 0)            &
        used = g_send_data(cobalt%id_detocos,  cobalt%p_ndet(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_calcos .gt. 0)            &
        used = g_send_data(cobalt%id_calcos,  cobalt%p_cadet_calc(:,:,1,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_aragos .gt. 0)            &
        used = g_send_data(cobalt%id_aragos,  cobalt%p_cadet_arag(:,:,1,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phydiatos.gt. 0)  &
        used = g_send_data(cobalt%id_phydiatos,  cobalt%nlg_diatoms(:,:,1) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phydiazos.gt. 0)  &
        used = g_send_data(cobalt%id_phydiazos,  cobalt%p_ndi(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phypicoos.gt. 0)  &
        used = g_send_data(cobalt%id_phypicoos,  cobalt%p_nsm(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phymiscos.gt. 0)  &
        used = g_send_data(cobalt%id_phymiscos,  (cobalt%p_nlg(:,:,1,tau)-cobalt%nlg_diatoms(:,:,1)) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_zmicroos.gt. 0)  &
        used = g_send_data(cobalt%id_zmicroos,  cobalt%p_nsmz(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_zmesoos.gt. 0)  &
        used = g_send_data(cobalt%id_zmesoos,  (cobalt%p_nlgz(:,:,1,tau)+cobalt%p_nmdz(:,:,1,tau)) * cobalt%c_2_n * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_talkos .gt. 0)            &
        used = g_send_data(cobalt%id_talkos,  cobalt%p_alk(:,:,1,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_talknatos .gt. 0)            &
!        used = g_send_data(cobalt%id_talknatos,                     
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: this is using ntau=1
    if (cobalt%id_phos .gt. 0)            &
        used = g_send_data(cobalt%id_phos,  log10(cobalt%f_htotal(:,:,1)+epsln) * (-1.0),       &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING: 
!    if (cobalt%id_phnatos .gt. 0)            &
!        used = g_send_data(cobalt%id_phnatos,  log10(cobalt%f_htotal(:,:,1)+epsln) * -1.0,       &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING: 
!    if (cobalt%id_phabioos .gt. 0)            &
!        used = g_send_data(cobalt%id_phabioos,  log10(cobalt%f_htotal(:,:,1)+epsln) * -1.0,       &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_o2os .gt. 0)            &
        used = g_send_data(cobalt%id_o2os,  cobalt%p_o2(:,:,1,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_o2satos .gt. 0)            &
        used = g_send_data(cobalt%id_o2satos,  cobalt%o2sat(:,:,1), &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_no3os .gt. 0)            &
        used = g_send_data(cobalt%id_no3os,  cobalt%p_no3(:,:,1,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_nh4os.gt. 0)            &
        used = g_send_data(cobalt%id_nh4os,  cobalt%p_nh4(:,:,1,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_po4os.gt. 0)            &
        used = g_send_data(cobalt%id_po4os,  cobalt%p_po4(:,:,1,tau) * cobalt%Rho_0,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_dfeos .gt. 0)            &
        used = g_send_data(cobalt%id_dfeos,  cobalt%p_fed(:,:,1,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_sios .gt. 0)            &
        used = g_send_data(cobalt%id_sios,  cobalt%p_sio4(:,:,1,tau) * cobalt%Rho_0,       &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_chlos .gt. 0)            &
        used = g_send_data(cobalt%id_chlos,  cobalt%f_chl(:,:,1) * cobalt%Rho_0 / 1e9,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_chldiatos .gt. 0)            &
        used = g_send_data(cobalt%id_chldiatos,  phyto(LARGE)%theta(:,:,1) * cobalt%nlg_diatoms(:,:,1) * cobalt%c_2_n * cobalt%Rho_0 * 12e3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_chldiazos .gt. 0)            &
        used = g_send_data(cobalt%id_chldiazos,  phyto(DIAZO)%theta(:,:,1) * cobalt%p_ndi(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0 * 12e3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_chlpicoos .gt. 0)            &
        used = g_send_data(cobalt%id_chlpicoos,  phyto(SMALL)%theta(:,:,1) * cobalt%p_nsm(:,:,1,tau) * cobalt%c_2_n * cobalt%Rho_0 * 12e3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_chlmiscos .gt. 0)            &
        used = g_send_data(cobalt%id_chlmiscos,  phyto(LARGE)%theta(:,:,1) * (cobalt%p_nlg(:,:,1,tau)-cobalt%nlg_diatoms(:,:,1)) *  &
        cobalt%c_2_n * cobalt%Rho_0 * 12e3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_ponos .gt. 0)            &
        used = g_send_data(cobalt%id_ponos,  (cobalt%p_ndi(:,:,1,tau) + cobalt%p_nlg(:,:,1,tau) + &
        cobalt%p_nsm(:,:,1,tau) + cobalt%p_nbact(:,:,1,tau) +  cobalt%p_ndet(:,:,1,tau) + &
        cobalt%p_nsmz(:,:,1,tau) + cobalt%p_nmdz(:,:,1,tau) + cobalt%p_nlgz(:,:,1,tau)) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_popos .gt. 0)            &
        used = g_send_data(cobalt%id_popos,  ((phyto(DIAZO)%p_2_n_static * cobalt%p_ndi(:,:,1,tau)) + cobalt%p_nlg(:,:,1,tau) + &
        cobalt%p_nsm(:,:,1,tau) + cobalt%p_nbact(:,:,1,tau) +  cobalt%p_pdet(:,:,1,tau) + &
        cobalt%p_nsmz(:,:,1,tau) + cobalt%p_nmdz(:,:,1,tau) + cobalt%p_nlgz(:,:,1,tau)) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_bfeos .gt. 0)            &
        used = g_send_data(cobalt%id_bfeos,  (cobalt%p_fedi(:,:,1,tau) + cobalt%p_felg(:,:,1,tau) + cobalt%p_fesm(:,:,1,tau) + & 
        cobalt%p_fedet(:,:,1,tau))  * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_bsios .gt. 0)            &
        used = g_send_data(cobalt%id_bsios,  (cobalt%p_silg(:,:,1,tau) + cobalt%p_sidet(:,:,1,tau))  * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phynos .gt. 0)            &
        used = g_send_data(cobalt%id_phynos,  (cobalt%p_nlg(:,:,1,tau) + cobalt%p_nsm(:,:,1,tau) +  &
        cobalt%p_ndi(:,:,1,tau)) * cobalt%Rho_0, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phypos .gt. 0)            &
        used = g_send_data(cobalt%id_phypos,  ((phyto(DIAZO)%p_2_n_static * cobalt%p_ndi(:,:,1,tau)) + cobalt%p_nlg(:,:,1,tau) + &
        cobalt%p_nsm(:,:,1,tau)) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_phyfeos .gt. 0)            &
        used = g_send_data(cobalt%id_phyfeos,  (cobalt%p_fedi(:,:,1,tau) + cobalt%p_felg(:,:,1,tau) +  &
        cobalt%p_fesm(:,:,1,tau)) * cobalt%Rho_0, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_physios .gt. 0)            &
        used = g_send_data(cobalt%id_physios,  cobalt%p_silg(:,:,1,tau) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_co3os .gt. 0)            &
        used = g_send_data(cobalt%id_co3os,  cobalt%f_co3_ion(:,:,1) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_co3natos .gt. 0)            &
         used = g_send_data(cobalt%id_co3natos,  cobalt%f_co3_ion(:,:,1) * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (cobalt%id_co3abioos .gt. 0)            &
         used = g_send_data(cobalt%id_co3abioos,  cobalt%f_co3_ion(:,:,1) * cobalt%Rho_0,  &
         model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_co3satcalcos .gt. 0)            &
        used = g_send_data(cobalt%id_co3satcalcos,  cobalt%co3_sol_calc(:,:,1) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_co3sataragos .gt. 0)            &
        used = g_send_data(cobalt%id_co3sataragos,  cobalt%co3_sol_arag(:,:,1) * cobalt%Rho_0,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Omon: 3-D Marine Biogeochemical 3-D Fields  
! Tracers and rates above
! Limitation terms below 
! 2018/07/18 changed to 2d terms
!
    if (cobalt%id_limndiat .gt. 0)            &
        used = g_send_data(cobalt%id_limndiat,  phyto(LARGE)%nlim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS, JGJ: not outputting/serving limndiaz because diazotrophs are not N limited
!    if (cobalt%id_limndiaz .gt. 0)            &
!        used = g_send_data(cobalt%id_limndiaz,  phyto(DIAZO)%nlim_bw_100, &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limnpico .gt. 0)            &
        used = g_send_data(cobalt%id_limnpico,  phyto(SMALL)%nlim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limnmisc .gt. 0)            &
        used = g_send_data(cobalt%id_limnmisc,  phyto(LARGE)%nlim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limirrdiat .gt. 0)            &
        used = g_send_data(cobalt%id_limirrdiat,  phyto(LARGE)%irrlim_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limirrdiaz .gt. 0)            &
        used = g_send_data(cobalt%id_limirrdiaz,  phyto(DIAZO)%irrlim_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limirrpico .gt. 0)            &
        used = g_send_data(cobalt%id_limirrpico,  phyto(SMALL)%irrlim_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limirrmisc .gt. 0)            &
        used = g_send_data(cobalt%id_limirrmisc,  phyto(LARGE)%irrlim_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limfediat .gt. 0)            &
        used = g_send_data(cobalt%id_limfediat,  phyto(LARGE)%def_fe_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limfediaz .gt. 0)            &
        used = g_send_data(cobalt%id_limfediaz,  phyto(DIAZO)%def_fe_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limfepico .gt. 0)            &
        used = g_send_data(cobalt%id_limfepico,  phyto(SMALL)%def_fe_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limfemisc .gt. 0)            &
        used = g_send_data(cobalt%id_limfemisc,  phyto(LARGE)%def_fe_bw_100,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! added by JGJ, not requested for CMIP6
    if (cobalt%id_limpdiat .gt. 0)            &
        used = g_send_data(cobalt%id_limpdiat,  phyto(LARGE)%plim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limpdiaz .gt. 0)            &
        used = g_send_data(cobalt%id_limpdiaz,  phyto(DIAZO)%plim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limppico .gt. 0)            &
        used = g_send_data(cobalt%id_limppico,  phyto(SMALL)%plim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_limpmisc .gt. 0)            &
        used = g_send_data(cobalt%id_limpmisc,  phyto(LARGE)%plim_bw_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem Omon: Marine Biogeochemical 2-D Fields  

    if (cobalt%id_intpp .gt. 0)            &
        used = g_send_data(cobalt%id_intpp,  cobalt%jprod_allphytos_100 * cobalt%c_2_n, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intppnitrate .gt. 0)            &
        used = g_send_data(cobalt%id_intppnitrate,  (phyto(DIAZO)%jprod_n_new_100 +  phyto(LARGE)%jprod_n_new_100 +  &
        phyto(SMALL)%jprod_n_new_100) * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: defined jprod_diat_100, removed extra "p" from individual group production
    if (cobalt%id_intppdiat .gt. 0)            &
        used = g_send_data(cobalt%id_intppdiat,  cobalt%jprod_diat_100 * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intppdiaz .gt. 0)            &
        used = g_send_data(cobalt%id_intppdiaz,  phyto(DIAZO)%jprod_n_100 * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intpppico .gt. 0)            &
        used = g_send_data(cobalt%id_intpppico,  phyto(SMALL)%jprod_n_100 * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: can now use jprod_diat_100 to back out misc production
    if (cobalt%id_intppmisc .gt. 0)            &
        used = g_send_data(cobalt%id_intppmisc,  (phyto(LARGE)%jprod_n_100 - cobalt%jprod_diat_100)  * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: I think this is fine
    if (cobalt%id_intpbn .gt. 0)            &
        used = g_send_data(cobalt%id_intpbn,  cobalt%jprod_allphytos_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: I've added juptake_po4_100 in a manner analagous to iron below
    if (cobalt%id_intpbp .gt. 0)             & 
        used = g_send_data(cobalt%id_intpbp, (phyto(DIAZO)%juptake_po4_100 +  phyto(LARGE)%juptake_po4_100 +  &
        phyto(SMALL)%juptake_po4_100),  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intpbfe .gt. 0)            &
        used = g_send_data(cobalt%id_intpbfe,  (phyto(DIAZO)%juptake_fe_100 +  phyto(LARGE)%juptake_fe_100 +  &
        phyto(SMALL)%juptake_fe_100),  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intpbsi .gt. 0)            &
        used = g_send_data(cobalt%id_intpbsi,  phyto(LARGE)%juptake_sio4_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intpcalcite .gt. 0)            &
        used = g_send_data(cobalt%id_intpcalcite,  cobalt%jprod_cadet_calc_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intparag .gt. 0)            &
        used = g_send_data(cobalt%id_intparag,  cobalt%jprod_cadet_arag_100, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epc100 .gt. 0)            &
        used = g_send_data(cobalt%id_epc100,  cobalt%fndet_100 * cobalt%c_2_n,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epn100 .gt. 0)            &
        used = g_send_data(cobalt%id_epn100,  cobalt%fndet_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epp100 .gt. 0)            &
        used = g_send_data(cobalt%id_epp100,  cobalt%fpdet_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epfe100 .gt. 0)            &
        used = g_send_data(cobalt%id_epfe100,  cobalt%ffedet_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epsi100 .gt. 0)            &
        used = g_send_data(cobalt%id_epsi100,  cobalt%fsidet_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_epcalc100 .gt. 0)            &
        used = g_send_data(cobalt%id_epcalc100,  cobalt%fcadet_calc_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: should be AT 100m 
    if (cobalt%id_eparag100 .gt. 0)            &
        used = g_send_data(cobalt%id_eparag100,  cobalt%fcadet_arag_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: should be wc_vert_int_dic?, *12e-3 to go from moles C m-2 to kg C m-2
    if (cobalt%id_intdic .gt. 0)            &
        used = g_send_data(cobalt%id_intdic,  cobalt%wc_vert_int_dic*12e-3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: added wc_vert_int_doc, *12e-3 to go from moles C m-2 to kg C m-2
    if (cobalt%id_intdoc .gt. 0)            &          
        used = g_send_data(cobalt%id_intdoc,  cobalt%wc_vert_int_doc*12e-3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: added wc_vert_int_poc, *12e-3 to go from moles C m-2 to kg C m-2
    if (cobalt%id_intpoc .gt. 0)            &
        used = g_send_data(cobalt%id_intpoc,  cobalt%wc_vert_int_poc*12e-3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_spco2 .gt. 0)            &
        used = g_send_data(cobalt%id_spco2,  cobalt%pco2_csurf * 0.1013,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_spco2nat .gt. 0)            &
!        used = g_send_data(cobalt%id_spco2nat,  cobalt%pco2_csurf * 0.1013,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_spco2abio .gt. 0)            &
!        used = g_send_data(cobalt%id_spco2abio,  cobalt%pco2_csurf * 0.1013,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3:
    if (cobalt%id_dpco2 .gt. 0)            &
        used = g_send_data(cobalt%id_dpco2,  cobalt%deltap_dic * 0.1013,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_dpco2nat .gt. 0)            &
!        used = g_send_data(cobalt%id_dpco2nat,  cobalt%dic_deltap * 0.1013,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_dpco2abio .gt. 0)            &
!        used = g_send_data(cobalt%id_dpco2abio,  cobalt%dic_deltap * 0.1013,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_dpo2 .gt. 0)            &
        used = g_send_data(cobalt%id_dpo2,  cobalt%deltap_o2 * 0.1013,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_fgco2 .gt. 0)            &
        used = g_send_data(cobalt%id_fgco2,  cobalt%stf_gas_dic * 12e-3,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! PENDING:
!    if (cobalt%id_fgco2nat .gt. 0)            &
!        used = g_send_data(cobalt%id_fgco2nat,  
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fgco2abio .gt. 0)            &
!        used = g_send_data(cobalt%id_fgco2abio,  
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fg14co2abio .gt. 0)            &
!        used = g_send_data(cobalt%id_fg14co2abio,  
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_fgo2 .gt. 0)            &
        used = g_send_data(cobalt%id_fgo2,  cobalt%stf_gas_o2,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: fcased_redis is accounted for elsewhere, just keep runoff
    if (cobalt%id_icfriver .gt. 0)            &
        used = g_send_data(cobalt%id_icfriver,  cobalt%runoff_flux_dic,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: Updated based on 7/19 e-mails with jpd and jgj
    if (cobalt%id_fric .gt. 0)            &
        used = g_send_data(cobalt%id_fric,  cobalt%fcadet_calc_btm -  cobalt%fcased_redis,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: organic nitrogen runoff from rivers*c_2_n ratio
    if (cobalt%id_ocfriver .gt. 0)            &
        used = g_send_data(cobalt%id_ocfriver, cobalt%c_2_n* &
        (cobalt%runoff_flux_ldon+cobalt%runoff_flux_sldon+cobalt%runoff_flux_srdon),&  
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CAS: equal to ndet burial*c_2_n
    if (cobalt%id_froc .gt. 0)            &
        used = g_send_data(cobalt%id_froc,cobalt%c_2_n*cobalt%fndet_burial, & 
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_intpn2 .gt. 0)            &
        used = g_send_data(cobalt%id_intpn2,  cobalt%wc_vert_int_nfix,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: TOPAZ: nongas_source_n=runoff_flux_no3+runoff_flux_nh4 +dry_no3+wet_no3+dry_nh4+wet_nh4+wc_vert_int_nfix
! CAS: should we include 1) don fluxes from rivers? 2) nh4 fluxes from rivers? 3) nh4 deposition?
    if (cobalt%id_fsn .gt. 0)            &
        used = g_send_data(cobalt%id_fsn,  cobalt%runoff_flux_no3 + cobalt%dry_no3 + cobalt%wet_no3 + cobalt%wc_vert_int_nfix,  &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: TOPAZ: fno3denit_tot=fno3denit_sed+wc_vert_int_jno3denit,  where wc_vert_int_jno3denit=jno3denit_wc*rho_dzt (sum over k)
! CAS: added burial
    if (cobalt%id_frn .gt. 0)            &
        used = g_send_data(cobalt%id_frn,  cobalt%fno3denit_sed + cobalt%wc_vert_int_jno3denit + cobalt%fndet_burial, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK3: TOPAZ: nongas_source_fe=runoff_flux_fed+wc_vert_int_jfe_coast+dry_fed+wet_fed+ffe_sed, where wc_vert_int_jfe_coast=jfe_coast*rho_dzt (sum over k)
    if (cobalt%id_fsfe .gt. 0)            &
        used = g_send_data(cobalt%id_fsfe,  cobalt%runoff_flux_fed + cobalt%dry_fed + cobalt%wet_fed + &
        cobalt%ffe_sed+cobalt%ffe_geotherm, &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_frfe .gt. 0)            &
        used = g_send_data(cobalt%id_frfe,  cobalt%ffedet_btm,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! revise calculation if providing to CMIP6
! 2016/08/15 - we will not be providing o2min, zo2min
!    if (cobalt%id_o2min .gt. 0)            &
!        used = g_send_data(cobalt%id_o2min,  cobalt%o2min,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_zo2min .gt. 0)            &
!        used = g_send_data(cobalt%id_zo2min,  cobalt%zo2min,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! 2016/08/15 - we will not be providing zsatcalc, zsatarag
    if (cobalt%id_zsatcalc .gt. 0)            &
        used = g_send_data(cobalt%id_zsatcalc,  cobalt%z_sat_calc,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_zsatarag .gt. 0)            &
        used = g_send_data(cobalt%id_zsatarag,  cobalt%z_sat_arag,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! 2016/08/15 - we will not be providing these fields
! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtdic .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtdic,  cobalt%f_dic_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtdin .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtdin,  cobalt%f_din_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtdip .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtdip,  cobalt%f_po4_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtdife .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtdife,  cobalt%f_fed_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtdisi .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtdisi,  cobalt%f_sio4_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

! CHECK: rate was computed offline for TOPAZ by saving a reference history file, dividing by secs_per_month and differencing monthly averages
! can we compute rates in the code this time?
!    if (cobalt%id_fddtalk .gt. 0)            &
!        used = g_send_data(cobalt%id_fddtalk,  cobalt%f_alk_int_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtdic .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtdic,  cobalt%jdic_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtdin .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtdin,  cobalt%jdin_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtdip .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtdip,  cobalt%jpo4_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtdife .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtdife,  cobalt%jfed_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtdisi .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtdisi,  cobalt%jsio4_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!    if (cobalt%id_fbddtalk .gt. 0)            &
!        used = g_send_data(cobalt%id_fbddtalk,  cobalt%jalk_100,   &
!        model_time, rmask = grid_tmask(:,:,1),&
!        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem day: Marine Biogeochemical daily fields
! chlos and phycos - in Omon and Oday
!==============================================================================================================
! 2016/08/15 JGJ: 100m integrals w/o CMOR conversion 

    if (cobalt%id_f_dic_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_dic_int_100,  cobalt%f_dic_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_f_din_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_din_int_100,  cobalt%f_din_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_f_po4_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_po4_int_100,  cobalt%f_po4_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_f_fed_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_fed_int_100,  cobalt%f_fed_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_f_sio4_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_sio4_int_100,  cobalt%f_sio4_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_f_alk_int_100 .gt. 0)            &
        used = g_send_data(cobalt%id_f_alk_int_100,  cobalt%f_alk_int_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jdic_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jdic_100,  cobalt%jdic_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jdin_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jdin_100,  cobalt%jdin_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jpo4_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jpo4_100,  cobalt%jpo4_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jfed_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jfed_100,  cobalt%jfed_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jsio4_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jsio4_100,  cobalt%jsio4_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

    if (cobalt%id_jalk_100 .gt. 0)            &
        used = g_send_data(cobalt%id_jalk_100,  cobalt%jalk_100,   &
        model_time, rmask = grid_tmask(:,:,1),&
        is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

!==============================================================================================================

    call mpp_clock_end(id_clock_cobalt_send_diagnostics)

  end subroutine generic_COBALT_update_from_source


  ! <SUBROUTINE NAME="generic_COBALT_set_boundary_values">
  !  <OVERVIEW>
  !   Calculate and set coupler values at the surface / bottom
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau,dzt)
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
  subroutine generic_COBALT_set_boundary_values(tracer_list,SST,SSS,rho,ilb,jlb,tau,dzt,model_time)
    type(g_tracer_type),          pointer    :: tracer_list
    real, dimension(ilb:,jlb:),   intent(in)   :: SST, SSS
    real, dimension(ilb:,jlb:,:,:), intent(in) :: rho
    integer,                        intent(in) :: ilb,jlb,tau
    real, dimension(ilb:,jlb:,:), optional, intent(in) :: dzt
    type(time_type),                intent(in) :: model_time
    integer :: isc,iec, jsc,jec,isd,ied,jsd,jed,nk,ntau , i, j
    real    :: sal,ST,o2_saturation
    real    :: tt,tk,ts,ts2,ts3,ts4,ts5
    real, dimension(:,:,:)  ,pointer  :: grid_tmask
    real, dimension(:,:,:,:), pointer :: o2_field,dic_field,po4_field,sio4_field,alk_field,di14c_field,nh4_field
    real, dimension(:,:,:), ALLOCATABLE :: htotal_field,co3_ion_field
    real, dimension(:,:), ALLOCATABLE :: co2_alpha,co2_csurf,co2_sc_no,o2_alpha,o2_csurf,o2_sc_no,nh3_alpha,nh3_csurf,nh3_sc_no,phos_nh3_exchange
    real, dimension(:,:), ALLOCATABLE :: c14o2_alpha,c14o2_csurf
    real :: pka_nh3,tr,ltr

    logical :: phos_nh3_override

    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_set_boundary_values'

    !
    !
    !Get the necessary properties
    !
    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    call g_tracer_get_pointer(tracer_list,'o2' ,'field',  o2_field)

    allocate(co2_alpha(isd:ied, jsd:jed)); co2_alpha=0.0
    allocate(co2_csurf(isd:ied, jsd:jed)); co2_csurf=0.0
    allocate(co2_sc_no(isd:ied, jsd:jed)); co2_sc_no=0.0
    allocate(nh3_alpha(isd:ied, jsd:jed)); nh3_alpha=0.0
    allocate(nh3_csurf(isd:ied, jsd:jed)); nh3_csurf=0.0
    allocate(nh3_sc_no(isd:ied, jsd:jed)); nh3_sc_no=0.0
    !for nh3 ph emission override
    allocate(phos_nh3_exchange(isd:ied, jsd:jed)); phos_nh3_exchange=0.0
    allocate(c14o2_alpha(isd:ied, jsd:jed)); c14o2_alpha=0.0
    allocate(c14o2_csurf(isd:ied, jsd:jed)); c14o2_csurf=0.0
    allocate(o2_alpha(isd:ied, jsd:jed)); o2_alpha=0.0
    allocate(o2_csurf(isd:ied, jsd:jed)); o2_csurf=0.0
    allocate(o2_sc_no(isd:ied, jsd:jed)); o2_sc_no=0.0
    allocate(htotal_field(isd:ied,jsd:jed,nk),co3_ion_field(isd:ied,jsd:jed,nk))
    htotal_field=0.0 ; co3_ion_field=0.0

    !nnz: Since the generic_COBALT_update_from_source() subroutine is called by this time
    !     the following if block is not really necessary (since this calculation is already done in source).
    !     It is only neccessary if source routine is commented out for debugging.
    !Note: In order for this to work we should NOT zero out the coupler values for generic tracers
    !      This zeroing is done for non-generic TOPAZ by calling zero_ocean_sfc.
    !      Since the coupler values here are non-cumulative there is no need to zero them out anyway.
    
    if (cobalt%init .OR. cobalt%force_update_fluxes) then
       !Get necessary fields
       call g_tracer_get_pointer(tracer_list,'dic'   ,'field', dic_field)
       call g_tracer_get_pointer(tracer_list,'po4'   ,'field', po4_field)
       call g_tracer_get_pointer(tracer_list,'sio4'  ,'field', sio4_field)
       call g_tracer_get_pointer(tracer_list,'alk'   ,'field', alk_field)
       call g_tracer_get_pointer(tracer_list,'nh4'   ,'field', nh4_field)

       call g_tracer_get_values(tracer_list,'htotal' ,'field', htotal_field,isd,jsd,ntau=1)
       call g_tracer_get_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       do j = jsc, jec ; do i = isc, iec  !{
          cobalt%htotallo(i,j) = cobalt%htotal_scale_lo * htotal_field(i,j,1)
          cobalt%htotalhi(i,j) = cobalt%htotal_scale_hi * htotal_field(i,j,1)
       enddo; enddo ; !} i, j

       if(present(dzt)) then
          ! 2017/08/11 jgj is cobalt type defined/passed here ?
          !        do j = jsc, jec ; do i = isc, iec  !{
          !         cobalt%zt(i,j,1) = dzt(i,j,1)
          !        enddo; enddo ; !} i, j
       elseif (trim(co2_calc) == 'mocsy') then
          call mpp_error(FATAL,"mocsy method of co2_calc needs dzt to be passed to the FMS_ocmip2_co2calc subroutine.")
       endif

       call FMS_ocmip2_co2calc(CO2_dope_vec,grid_tmask(:,:,1), &
            SST(:,:), SSS(:,:),                            &
            dic_field(:,:,1,tau),                          &
            po4_field(:,:,1,tau),                          &
            sio4_field(:,:,1,tau),                         &
            alk_field(:,:,1,tau),                          &
            cobalt%htotallo, cobalt%htotalhi,              &
                                !InOut
            htotal_field(:,:,1),                           &
                                !Optional In
            co2_calc=trim(co2_calc),                       & 
            !! jgj 2017/08/11
            !!zt=cobalt%zt(:,:,1),                           & 
            zt=dzt(:,:,1),                                 & 
                                !OUT
            co2star=co2_csurf(:,:), alpha=co2_alpha(:,:),  &
            pCO2surf=cobalt%pco2_csurf(:,:), &
            co3_ion=co3_ion_field(:,:,1), &
            omega_arag=cobalt%omegaa(:,:,1), &
            omega_calc=cobalt%omegac(:,:,1))

       !Set fields !nnz: if These are pointers do I need to do this?
       call g_tracer_set_values(tracer_list,'htotal' ,'field',htotal_field ,isd,jsd,ntau=1)
       call g_tracer_set_values(tracer_list,'co3_ion','field',co3_ion_field,isd,jsd,ntau=1)

       call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha    ,isd,jsd)
       call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf    ,isd,jsd)


      if (do_14c) then                                        !<<RADIOCARBON
      
        ! Normally, the alpha would be multiplied by the atmospheric 14C/12C ratio. However,
        ! here that is set to 1, so that alpha_14C = alpha_12C. This needs to be changed!

       call g_tracer_get_pointer(tracer_list,'di14c'   ,'field', di14c_field)
    
        do j = jsc, jec ; do i = isc, iec  !{
          c14o2_csurf(i,j) =  co2_csurf(i,j) *                &
            di14c_field(i,j,1,tau) / (dic_field(i,j,1,tau) + epsln)
          c14o2_alpha(i,j) =  co2_alpha(i,j)
        enddo; enddo ; !} i, j

        call g_tracer_set_values(tracer_list,'di14c','alpha',c14o2_alpha      ,isd,jsd)
        call g_tracer_set_values(tracer_list,'di14c','csurf',c14o2_csurf      ,isd,jsd)

      endif 

       if (do_nh3_atm_ocean_exchange) then
          !          write(*,*) 'min htot ',minval(htotal_field(:,:,1))

          call data_override('OCN', 'phos_nh3_exchange', phos_nh3_exchange(isc:iec,jsc:jec), model_time,override=phos_nh3_override)       

          do j = jsc, jec ; do i = isc, iec  !{
             !nh3
             pka_nh3        = calc_pka_nh3(SST(i,j),SSS(i,j))
             tr             = 298.15/(SST(i,j)+273.15)-1.
             ltr            = -tr+log(298.15/(SST(i,j)+273.15))       
             !mol/kg/atm from Jacobson 2005 (fundamental of atmospheric modeling)
             !997/1035 is to convert pure water to salt water
             nh3_alpha(i,j) = 5.76e1*exp(13.79*tr-5.39*ltr)*997. !in mol/kg(water)/atm -> mol/m3/atm
             nh3_alpha(i,j) = nh3_alpha(i,j)/saltout_correction(101325./(1e-3*rdgas*wtmair*(SST(i,j)+273.15)*nh3_alpha(i,j)),vb_nh3,SSS(i,j)) * 1./cobalt%Rho_0 !mol/kg/atm
             if (phos_nh3_override) then
                nh3_csurf(i,j) = nh4_field(i,j,1,tau)/(1.+10**(pka_nh3-(max(min(11.,phos_nh3_exchange(i,j)),3.))))
             else
                nh3_csurf(i,j) = nh4_field(i,j,1,tau)/(1.+10**(pka_nh3+log10(min(max(1e-11,htotal_field(i,j,1)),1e-3))))
             end if
             cobalt%pnh3_csurf(i,j)  =  cobalt%nh3_csurf(i,j)/nh3_alpha(i,j)*1.e6 !in uatm
          enddo;enddo

          call g_tracer_set_values(tracer_list,'nh4','alpha',nh3_alpha    ,isd,jsd)
          call g_tracer_set_values(tracer_list,'nh4','csurf',nh3_csurf    ,isd,jsd)

       end if   
       !!nnz: If source is called uncomment the following
       cobalt%init = .false. !nnz: This is necessary since the above two calls appear in source subroutine too.
    endif

    call g_tracer_get_values(tracer_list,'dic','alpha', co2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'dic','csurf', co2_csurf ,isd,jsd)

    call g_tracer_get_values(tracer_list,'o2','alpha', o2_alpha ,isd,jsd)
    call g_tracer_get_values(tracer_list,'o2','csurf', o2_csurf ,isd,jsd)


    do j=jsc,jec ; do i=isc,iec
       !This calculation needs an input of SST and SSS
       sal = SSS(i,j) ; ST = SST(i,j)

       !nnz:
       !Note: In the following calculations in order to get results for co2 and o2
       !      identical with cobalt code in MOM cobalt%Rho_0 must be replaced with rho(i,j,1,tau)
       !      This is achieved by uncommenting the following if desired.
       !! cobalt%Rho_0 = rho(i,j,1,tau)
       !      But since %Rho_0 plays the role of a unit conversion factor in this module
       !      it may be safer to keep it as a constant (1035.0) rather than the actual variable
       !      surface density rho(i,j,1,tau)

       !---------------------------------------------------------------------
       !     CO2
       !---------------------------------------------------------------------

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of CO2 in seawater using the
       !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
       !  7373-7382).
       !  2018/01/17 jgj  update Schmidt number for CO2 to use
       !  Wanninkhof, Limnol. Oceanogr: Methods, 12, 2014, 351-362
       !---------------------------------------------------------------------
       if (trim(as_param_cobalt) == 'W92') then
         co2_sc_no(i,j) = cobalt%a1_co2 + ST*(cobalt%a2_co2 + ST*(cobalt%a3_co2 + ST*cobalt%a4_co2)) * & 
                          grid_tmask(i,j,1)
       else if ((trim(as_param_cobalt) == 'W14') .or. (trim(as_param_cobalt) == 'gfdl_cmip6')) then
         co2_sc_no(i,j) = cobalt%a1_co2 + ST*(cobalt%a2_co2 + ST*(cobalt%a3_co2 + & 
                          ST*(cobalt%a4_co2 + ST*cobalt%a5_co2))) * grid_tmask(i,j,1)
       endif
!       sc_no_term = sqrt(660.0 / (sc_co2 + epsln))
!
!       co2_alpha(i,j) = co2_alpha(i,j)* sc_no_term * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
!       co2_csurf(i,j) = co2_csurf(i,j)* sc_no_term * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
!
! in 'ocmip2_new' atmos_ocean_fluxes.F90 coupler formulation, the schmidt number is carried in explicitly
!
       co2_alpha(i,j) = co2_alpha(i,j) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)
       co2_csurf(i,j) = co2_csurf(i,j) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)

       !---------------------------------------------------------------------
       !     O2
       !---------------------------------------------------------------------
       !  Compute the oxygen saturation concentration at 1 atm total
       !  pressure in mol/kg given the temperature (t, in deg C) and
       !  the salinity (s, in permil)
       !
       !  From Garcia and Gordon (1992), Limnology and Oceonography.
       !  The formula used is from page 1310, eq (8).
       !
       !  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
       !  *** It shouldn't be there.                                ***
       !
       !  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
       !                                   0 permil <= S <= 42 permil
       ! 2015/05/15  jgj ESM2.6 has values of 60+ in Red Sea - impose
       ! bounds on salinity to keep it in 0-42 range
       !
       ! check value: T = 10 deg C, S = 35 permil,
       !              o2_saturation = 0.282015 mol m-3
       !---------------------------------------------------------------------
       !
       ! jgj 2015/05/14 impose temperature and salinity bounds for o2sat 
       sal = min(42.0,max(0.0,sal))
       tt = 298.15 - min(40.0,max(0.0,ST))
       tk = 273.15 + min(40.0,max(0.0,ST))
       ts = log(tt / tk)
       ts2 = ts  * ts
       ts3 = ts2 * ts
       ts4 = ts3 * ts
       ts5 = ts4 * ts

       !The atmospheric code needs solubilities in units of mol/m3/atm
       o2_saturation = (1000.0/22391.6) * grid_tmask(i,j,1) *  & !convert from ml/l to mol m-3
            exp( cobalt%a_0 + cobalt%a_1*ts + cobalt%a_2*ts2 + cobalt%a_3*ts3 + cobalt%a_4*ts4 + cobalt%a_5*ts5 + &
            (cobalt%b_0 + cobalt%b_1*ts + cobalt%b_2*ts2 + cobalt%b_3*ts3 + cobalt%c_0*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !  2018/01/17 jgj  update Schmidt number for O2 to use
       !  Wanninkhof, Limnol. Oceanogr: Methods, 12, 2014, 351-362
       !---------------------------------------------------------------------
       !
       ! In 'ocmip2_generic' atmos_ocean_fluxes.F90 coupler formulation,
       ! the schmidt number is carried in explicitly
       !
       if (trim(as_param_cobalt) == 'W92') then
         o2_sc_no(i,j)  = cobalt%a1_o2 + ST * (cobalt%a2_o2 + ST * (cobalt%a3_o2 + ST * cobalt%a4_o2 )) * &
                          grid_tmask(i,j,1)
       else if ((trim(as_param_cobalt) == 'W14') .or. (trim(as_param_cobalt) == 'gfdl_cmip6')) then
         o2_sc_no(i,j) = cobalt%a1_o2 + ST*(cobalt%a2_o2 + ST*(cobalt%a3_o2 + &
                         ST*(cobalt%a4_o2 + ST*cobalt%a5_o2))) * grid_tmask(i,j,1)
       endif
       !
       !      renormalize the alpha value for atm o2
       !      data table override for o2_flux_pcair_atm is now set to 0.21
       !
       o2_alpha(i,j) = (o2_saturation / 0.21)
       o2_csurf(i,j) = o2_field(i,j,1,tau) * cobalt%Rho_0 !nnz: MOM has rho(i,j,1,tau)

    enddo; enddo

    !
    ! Set %csurf, %alpha and %sc_no for these tracers. This will mark them
    ! for sending fluxes to coupler
    !
    call g_tracer_set_values(tracer_list,'dic','alpha',co2_alpha,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','csurf',co2_csurf,isd,jsd)
    call g_tracer_set_values(tracer_list,'dic','sc_no',co2_sc_no,isd,jsd)

    call g_tracer_set_values(tracer_list,'o2', 'alpha',o2_alpha, isd,jsd)
    call g_tracer_set_values(tracer_list,'o2', 'csurf',o2_csurf, isd,jsd)
    call g_tracer_set_values(tracer_list,'o2', 'sc_no',o2_sc_no, isd,jsd)

    if (do_nh3_atm_ocean_exchange) then
       call g_tracer_get_values(tracer_list,'nh4','alpha', nh3_alpha ,isd,jsd)
       call g_tracer_get_values(tracer_list,'nh4','csurf', nh3_csurf ,isd,jsd)

       do j=jsc,jec ; do i=isc,iec
       !nh3
       !f1p
          nh3_sc_no(i,j) = schmidt_w(sst(i,j),sss(i,j),vb_nh3)*grid_tmask(i,j,1)
          nh3_csurf(i,j) = nh3_csurf(i,j)*cobalt%Rho_0
          nh3_alpha(i,j) = nh3_alpha(i,j)*cobalt%Rho_0
       end do;end do

       call g_tracer_set_values(tracer_list,'nh4', 'alpha',nh3_alpha, isd,jsd)
       call g_tracer_set_values(tracer_list,'nh4', 'csurf',nh3_csurf, isd,jsd)
       call g_tracer_set_values(tracer_list,'nh4', 'sc_no',nh3_sc_no, isd,jsd)

    end if

    if (do_14c) then                                      !<<RADIOCARBON
      
       call g_tracer_get_values(tracer_list,'di14c','alpha', c14o2_alpha ,isd,jsd)
       call g_tracer_get_values(tracer_list,'di14c','csurf', c14o2_csurf ,isd,jsd)
       
       do j=jsc,jec ; do i=isc,iec
         !---------------------------------------------------------------------
         !     14CO2 - schmidt number is calculated same as before, as is alpha.
         !   Need to multiply alpha by frac_14catm to get the real alpha.
         !   NOTE: FOR NOW, D14C fixed at 0 permil!! Need to fix this later.
         !---------------------------------------------------------------------

          sal = SSS(i,j) ; ST = SST(i,j)
          co2_sc_no(i,j) = cobalt%a1_co2 + ST * (cobalt%a2_co2 + ST * (cobalt%a3_co2 + ST * cobalt%a4_co2)) * &
               grid_tmask(i,j,1)
       
          c14o2_alpha(i,j) = c14o2_alpha(i,j) * cobalt%Rho_0 
          c14o2_csurf(i,j) = c14o2_csurf(i,j) * cobalt%Rho_0 
          
       enddo; enddo
       
       call g_tracer_set_values(tracer_list,'di14c','alpha',c14o2_alpha,isd,jsd)
       call g_tracer_set_values(tracer_list,'di14c','csurf',c14o2_csurf,isd,jsd)
       call g_tracer_set_values(tracer_list,'di14c','sc_no',co2_sc_no,isd,jsd)

    endif                                                  !RADIOCARBON>>

    deallocate(co2_alpha,co2_csurf,&
         co2_sc_no,o2_alpha,          &
         c14o2_alpha,c14o2_csurf,     &
         o2_csurf,o2_sc_no, nh3_alpha,nh3_csurf,&
         nh3_sc_no,phos_nh3_exchange)

  end subroutine generic_COBALT_set_boundary_values


  ! <SUBROUTINE NAME="generic_COBALT_end">
  !  <OVERVIEW>
  !   End the module.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Deallocate all work arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call generic_COBALT_end
  !  </TEMPLATE>
  ! </SUBROUTINE>


  subroutine generic_COBALT_end
    character(len=fm_string_len), parameter :: sub_name = 'generic_COBALT_end'
    call user_deallocate_arrays
  end subroutine generic_COBALT_end

  !
  !   This is an internal sub, not a public interface.
  !   Allocate all the work arrays to be used in this module.
  !
  subroutine user_allocate_arrays
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,n

    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau) 

    !Allocate all the private arrays.

    !Used in ocmip2_co2calc
    CO2_dope_vec%isc = isc ; CO2_dope_vec%iec = iec 
    CO2_dope_vec%jsc = jsc ; CO2_dope_vec%jec = jec
    CO2_dope_vec%isd = isd ; CO2_dope_vec%ied = ied
    CO2_dope_vec%jsd = jsd ; CO2_dope_vec%jed = jed    

    allocate(cobalt%htotallo(isd:ied,jsd:jed))
    allocate(cobalt%htotalhi(isd:ied,jsd:jed))

    !
    ! allocate and initialize array elements of all phytoplankton groups
    ! CAS: add fluxes for additional explicit phytoplankton loss terms

    do n = 1, NUM_PHYTO
       allocate(phyto(n)%def_fe(isd:ied,jsd:jed,nk))       ; phyto(n)%def_fe         = 0.0
       allocate(phyto(n)%def_p(isd:ied,jsd:jed,nk))        ; phyto(n)%def_p          = 0.0
       allocate(phyto(n)%f_fe(isd:ied,jsd:jed,nk))         ; phyto(n)%f_fe           = 0.0
       allocate(phyto(n)%f_n(isd:ied,jsd:jed,nk))          ; phyto(n)%f_n            = 0.0
       allocate(phyto(n)%felim(isd:ied,jsd:jed,nk))        ; phyto(n)%felim          = 0.0
       allocate(phyto(n)%irrlim(isd:ied,jsd:jed,nk))       ; phyto(n)%irrlim         = 0.0
       allocate(phyto(n)%jzloss_fe(isd:ied,jsd:jed,nk))    ; phyto(n)%jzloss_fe      = 0.0
       allocate(phyto(n)%jzloss_n(isd:ied,jsd:jed,nk))     ; phyto(n)%jzloss_n       = 0.0
       allocate(phyto(n)%jzloss_p(isd:ied,jsd:jed,nk))     ; phyto(n)%jzloss_p       = 0.0
       allocate(phyto(n)%jzloss_sio2(isd:ied,jsd:jed,nk))  ; phyto(n)%jzloss_sio2    = 0.0
       allocate(phyto(n)%jaggloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jaggloss_fe    = 0.0
       allocate(phyto(n)%jaggloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jaggloss_n     = 0.0
       allocate(phyto(n)%jaggloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jaggloss_p     = 0.0
       allocate(phyto(n)%jaggloss_sio2(isd:ied,jsd:jed,nk)); phyto(n)%jaggloss_sio2  = 0.0
       allocate(phyto(n)%jvirloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jvirloss_fe    = 0.0
       allocate(phyto(n)%jvirloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jvirloss_n     = 0.0
       allocate(phyto(n)%jvirloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jvirloss_p     = 0.0
       allocate(phyto(n)%jvirloss_sio2(isd:ied,jsd:jed,nk)); phyto(n)%jvirloss_sio2  = 0.0
       allocate(phyto(n)%jexuloss_fe(isd:ied,jsd:jed,nk))  ; phyto(n)%jexuloss_fe    = 0.0
       allocate(phyto(n)%jexuloss_n(isd:ied,jsd:jed,nk))   ; phyto(n)%jexuloss_n     = 0.0
       allocate(phyto(n)%jexuloss_p(isd:ied,jsd:jed,nk))   ; phyto(n)%jexuloss_p     = 0.0
       allocate(phyto(n)%jhploss_fe(isd:ied,jsd:jed,nk))   ; phyto(n)%jhploss_fe     = 0.0
       allocate(phyto(n)%jhploss_n(isd:ied,jsd:jed,nk))    ; phyto(n)%jhploss_n      = 0.0
       allocate(phyto(n)%jhploss_p(isd:ied,jsd:jed,nk))    ; phyto(n)%jhploss_p      = 0.0
       allocate(phyto(n)%jhploss_sio2(isd:ied,jsd:jed,nk)) ; phyto(n)%jhploss_sio2   = 0.0
       allocate(phyto(n)%juptake_fe(isd:ied,jsd:jed,nk))   ; phyto(n)%juptake_fe     = 0.0
       allocate(phyto(n)%juptake_nh4(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_nh4    = 0.0
       allocate(phyto(n)%juptake_no3(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_no3    = 0.0
       allocate(phyto(n)%juptake_po4(isd:ied,jsd:jed,nk))  ; phyto(n)%juptake_po4    = 0.0
       allocate(phyto(n)%jprod_n(isd:ied,jsd:jed,nk))      ; phyto(n)%jprod_n        = 0.0
       allocate(phyto(n)%liebig_lim(isd:ied,jsd:jed,nk))   ; phyto(n)%liebig_lim     = 0.0
       allocate(phyto(n)%mu(isd:ied,jsd:jed,nk))           ; phyto(n)%mu             = 0.0
       allocate(phyto(n)%po4lim(isd:ied,jsd:jed,nk))       ; phyto(n)%po4lim         = 0.0
       allocate(phyto(n)%q_fe_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_fe_2_n       = 0.0
       allocate(phyto(n)%q_p_2_n(isd:ied,jsd:jed,nk))      ; phyto(n)%q_p_2_n        = 0.0
       allocate(phyto(n)%q_si_2_n(isd:ied,jsd:jed,nk))     ; phyto(n)%q_si_2_n       = 0.0
       allocate(phyto(n)%theta(isd:ied,jsd:jed,nk))        ; phyto(n)%theta          = 0.0
       allocate(phyto(n)%f_mu_mem(isd:ied,jsd:jed,nk))     ; phyto(n)%f_mu_mem       = 0.0
       allocate(phyto(n)%mu_mix(isd:ied,jsd:jed,nk))       ; phyto(n)%mu_mix         = 0.0
       allocate(phyto(n)%agg_lim(isd:ied,jsd:jed,nk))      ; phyto(n)%agg_lim        = 0.0
       allocate(phyto(n)%nh4lim(isd:ied,jsd:jed,nk))       ; phyto(n)%nh4lim         = 0.0
       allocate(phyto(n)%no3lim(isd:ied,jsd:jed,nk))       ; phyto(n)%no3lim         = 0.0
    enddo
    !
    ! allocate and initialize array elements of only one phytoplankton group
    !
    allocate(phyto(DIAZO)%juptake_n2(isd:ied,jsd:jed,nk))   ; phyto(DIAZO)%juptake_n2   = 0.0
    allocate(phyto(DIAZO)%o2lim(isd:ied,jsd:jed,nk))        ; phyto(DIAZO)%o2lim        = 0.0
    allocate(phyto(LARGE)%juptake_sio4(isd:ied,jsd:jed,nk)) ; phyto(LARGE)%juptake_sio4 = 0.0
    allocate(phyto(LARGE)%silim(isd:ied,jsd:jed,nk))        ; phyto(LARGE)%silim      = 0.0
    !
    ! allocate and initialize arrays for bacteria
    !
    allocate(bact(1)%f_n(isd:ied,jsd:jed,nk))              ; bact(1)%f_n             = 0.0
    allocate(bact(1)%jzloss_n(isd:ied,jsd:jed,nk))         ; bact(1)%jzloss_n        = 0.0
    allocate(bact(1)%jzloss_p(isd:ied,jsd:jed,nk))         ; bact(1)%jzloss_p        = 0.0
    allocate(bact(1)%jvirloss_n(isd:ied,jsd:jed,nk))       ; bact(1)%jvirloss_n      = 0.0
    allocate(bact(1)%jvirloss_p(isd:ied,jsd:jed,nk))       ; bact(1)%jvirloss_p      = 0.0
    allocate(bact(1)%jhploss_n(isd:ied,jsd:jed,nk))        ; bact(1)%jhploss_n       = 0.0
    allocate(bact(1)%jhploss_p(isd:ied,jsd:jed,nk))        ; bact(1)%jhploss_p       = 0.0
    allocate(bact(1)%juptake_ldon(isd:ied,jsd:jed,nk))     ; bact(1)%juptake_ldon    = 0.0
    allocate(bact(1)%juptake_ldop(isd:ied,jsd:jed,nk))     ; bact(1)%juptake_ldop    = 0.0
    allocate(bact(1)%jprod_nh4(isd:ied,jsd:jed,nk))        ; bact(1)%jprod_nh4       = 0.0
    allocate(bact(1)%jprod_po4(isd:ied,jsd:jed,nk))        ; bact(1)%jprod_po4       = 0.0
    allocate(bact(1)%jprod_n(isd:ied,jsd:jed,nk))          ; bact(1)%jprod_n         = 0.0
    allocate(bact(1)%o2lim(isd:ied,jsd:jed,nk))            ; bact(1)%o2lim           = 0.0
    allocate(bact(1)%ldonlim(isd:ied,jsd:jed,nk))          ; bact(1)%ldonlim         = 0.0
    allocate(bact(1)%temp_lim(isd:ied,jsd:jed,nk))         ; bact(1)%temp_lim        = 0.0
    !
    ! CAS: allocate and initialize array elements for all zooplankton groups
    !
    do n = 1, NUM_ZOO
       allocate(zoo(n)%f_n(isd:ied,jsd:jed,nk))           ; zoo(n)%f_n            = 0.0
       allocate(zoo(n)%jzloss_n(isd:ied,jsd:jed,nk))      ; zoo(n)%jzloss_n       = 0.0
       allocate(zoo(n)%jzloss_p(isd:ied,jsd:jed,nk))      ; zoo(n)%jzloss_p       = 0.0
       allocate(zoo(n)%jhploss_n(isd:ied,jsd:jed,nk))     ; zoo(n)%jhploss_n      = 0.0
       allocate(zoo(n)%jhploss_p(isd:ied,jsd:jed,nk))     ; zoo(n)%jhploss_p      = 0.0
       allocate(zoo(n)%jingest_n(isd:ied,jsd:jed,nk))     ; zoo(n)%jingest_n      = 0.0
       allocate(zoo(n)%jingest_p(isd:ied,jsd:jed,nk))     ; zoo(n)%jingest_p      = 0.0
       allocate(zoo(n)%jingest_sio2(isd:ied,jsd:jed,nk))  ; zoo(n)%jingest_sio2   = 0.0
       allocate(zoo(n)%jingest_fe(isd:ied,jsd:jed,nk))    ; zoo(n)%jingest_fe     = 0.0
       allocate(zoo(n)%jprod_fed(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_fed      = 0.0
       allocate(zoo(n)%jprod_fedet(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_fedet    = 0.0
       allocate(zoo(n)%jprod_ndet(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ndet     = 0.0
       allocate(zoo(n)%jprod_pdet(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_pdet     = 0.0
       allocate(zoo(n)%jprod_ldon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ldon     = 0.0
       allocate(zoo(n)%jprod_ldop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_ldop     = 0.0
       allocate(zoo(n)%jprod_srdon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_srdon   = 0.0
       allocate(zoo(n)%jprod_srdop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_srdop   = 0.0
       allocate(zoo(n)%jprod_sldon(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_sldon   = 0.0
       allocate(zoo(n)%jprod_sldop(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_sldop   = 0.0
       allocate(zoo(n)%jprod_sidet(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_sidet    = 0.0
       allocate(zoo(n)%jprod_sio4(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_sio4      = 0.0
       allocate(zoo(n)%jprod_po4(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_po4      = 0.0
       allocate(zoo(n)%jprod_nh4(isd:ied,jsd:jed,nk))     ; zoo(n)%jprod_nh4      = 0.0
       allocate(zoo(n)%jprod_n(isd:ied,jsd:jed,nk))      ; zoo(n)%jprod_n         = 0.0
       allocate(zoo(n)%o2lim(isd:ied,jsd:jed,nk))        ; zoo(n)%o2lim           = 0.0
       allocate(zoo(n)%temp_lim(isd:ied,jsd:jed,nk))      ; zoo(n)%temp_lim       = 0.0
    enddo

    ! higher predator ingestion
    allocate(cobalt%hp_jingest_n(isd:ied,jsd:jed,nk))     ; cobalt%hp_jingest_n      = 0.0
    allocate(cobalt%hp_jingest_p(isd:ied,jsd:jed,nk))     ; cobalt%hp_jingest_p      = 0.0
    allocate(cobalt%hp_jingest_sio2(isd:ied,jsd:jed,nk))  ; cobalt%hp_jingest_sio2   = 0.0
    allocate(cobalt%hp_jingest_fe(isd:ied,jsd:jed,nk))    ; cobalt%hp_jingest_fe     = 0.0 

    allocate(cobalt%f_alk(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_alk=0.0
    allocate(cobalt%f_cadet_arag(isd:ied, jsd:jed, 1:nk)) ; cobalt%f_cadet_arag=0.0
    allocate(cobalt%f_cadet_calc(isd:ied, jsd:jed, 1:nk)) ; cobalt%f_cadet_calc=0.0
    allocate(cobalt%f_dic(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_dic=0.0
    allocate(cobalt%f_fed(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_fed=0.0
    allocate(cobalt%f_fedet(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_fedet=0.0
    allocate(cobalt%f_ldon(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ldon=0.0
    allocate(cobalt%f_ldop(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ldop=0.0
    allocate(cobalt%f_lith(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_lith=0.0
    allocate(cobalt%f_lithdet(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_lithdet=0.0
    allocate(cobalt%f_ndet(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_ndet=0.0
    allocate(cobalt%f_nh4(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_nh4=0.0
    allocate(cobalt%f_no3(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_no3=0.0
    allocate(cobalt%f_o2(isd:ied, jsd:jed, 1:nk))         ; cobalt%f_o2=0.0
    allocate(cobalt%f_pdet(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_pdet=0.0
    allocate(cobalt%f_po4(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_po4=0.0
    allocate(cobalt%f_srdon(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_srdon=0.0
    allocate(cobalt%f_srdop(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_srdop=0.0
    allocate(cobalt%f_sldon(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_sldon=0.0
    allocate(cobalt%f_sldop(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_sldop=0.0
    allocate(cobalt%f_sidet(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_sidet=0.0
    allocate(cobalt%f_silg(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_silg=0.0
    allocate(cobalt%f_sio4(isd:ied, jsd:jed, 1:nk))       ; cobalt%f_sio4=0.0
    allocate(cobalt%co3_sol_arag(isd:ied, jsd:jed, 1:nk)) ; cobalt%co3_sol_arag=0.0
    allocate(cobalt%co3_sol_calc(isd:ied, jsd:jed, 1:nk)) ; cobalt%co3_sol_calc=0.0
    allocate(cobalt%f_chl(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_chl=0.0
    if (do_nh3_diag) allocate(cobalt%f_nh3(isd:ied, jsd:jed, 1:nk))        ; cobalt%f_nh3=0.0
    allocate(cobalt%f_co3_ion(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_co3_ion=0.0
    allocate(cobalt%f_htotal(isd:ied, jsd:jed, 1:nk))     ; cobalt%f_htotal=0.0
    allocate(cobalt%f_irr_mem(isd:ied, jsd:jed, 1:nk))    ; cobalt%f_irr_mem=0.0
    allocate(cobalt%f_cased(isd:ied, jsd:jed, 1:nk))      ; cobalt%f_cased=0.0
    allocate(cobalt%f_cadet_arag_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_cadet_arag_btf=0.0
    allocate(cobalt%f_cadet_calc_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_cadet_calc_btf=0.0
    allocate(cobalt%f_fedet_btf(isd:ied, jsd:jed, 1:nk))  ; cobalt%f_fedet_btf=0.0
    allocate(cobalt%f_lithdet_btf(isd:ied, jsd:jed, 1:nk)); cobalt%f_lithdet_btf=0.0
    allocate(cobalt%f_ndet_btf(isd:ied, jsd:jed, 1:nk))   ; cobalt%f_ndet_btf=0.0
    allocate(cobalt%f_pdet_btf(isd:ied, jsd:jed, 1:nk))   ; cobalt%f_pdet_btf=0.0
    allocate(cobalt%f_sidet_btf(isd:ied, jsd:jed, 1:nk))  ; cobalt%f_sidet_btf=0.0

    allocate(cobalt%jnbact(isd:ied, jsd:jed, 1:nk))       ; cobalt%jnbact=0.0
    allocate(cobalt%jndi(isd:ied, jsd:jed, 1:nk))         ; cobalt%jndi=0.0
    allocate(cobalt%jnsm(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnsm=0.0
    allocate(cobalt%jnlg(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnlg=0.0
    allocate(cobalt%jnsmz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnsmz=0.0
    allocate(cobalt%jnmdz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnmdz=0.0
    allocate(cobalt%jnlgz(isd:ied, jsd:jed, 1:nk))        ; cobalt%jnlgz=0.0
    allocate(cobalt%jalk(isd:ied, jsd:jed, 1:nk))         ; cobalt%jalk=0.0
    allocate(cobalt%jalk_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jalk_plus_btm=0.0
    allocate(cobalt%jcadet_arag(isd:ied, jsd:jed, 1:nk))  ; cobalt%jcadet_arag=0.0
    allocate(cobalt%jcadet_calc(isd:ied, jsd:jed, 1:nk))  ; cobalt%jcadet_calc=0.0
    allocate(cobalt%jdic(isd:ied, jsd:jed, 1:nk))         ; cobalt%jdic=0.0
    allocate(cobalt%jdic_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jdic_plus_btm=0.0
    allocate(cobalt%jdin_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jdin_plus_btm=0.0
    allocate(cobalt%jfed(isd:ied, jsd:jed, 1:nk))         ; cobalt%jfed=0.0
    allocate(cobalt%jfed_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jfed_plus_btm=0.0
    allocate(cobalt%jfedi(isd:ied, jsd:jed, 1:nk))        ; cobalt%jfedi=0.0
    allocate(cobalt%jfelg(isd:ied, jsd:jed, 1:nk))        ; cobalt%jfelg=0.0
    allocate(cobalt%jfesm(isd:ied, jsd:jed, 1:nk))        ; cobalt%jfesm=0.0
    allocate(cobalt%jfedet(isd:ied, jsd:jed, 1:nk))       ; cobalt%jfedet=0.0
    allocate(cobalt%jldon(isd:ied, jsd:jed, 1:nk))        ; cobalt%jldon=0.0
    allocate(cobalt%jldop(isd:ied, jsd:jed, 1:nk))        ; cobalt%jldop=0.0
    allocate(cobalt%jlith(isd:ied, jsd:jed, 1:nk))        ; cobalt%jlith=0.0
    allocate(cobalt%jlithdet(isd:ied, jsd:jed, 1:nk))     ; cobalt%jlithdet=0.0
    allocate(cobalt%jndet(isd:ied, jsd:jed, 1:nk))        ; cobalt%jndet=0.0
    allocate(cobalt%jnh4(isd:ied, jsd:jed, 1:nk))         ; cobalt%jnh4=0.0
    allocate(cobalt%jnh4_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jnh4_plus_btm=0.0
    allocate(cobalt%jno3(isd:ied, jsd:jed, 1:nk))         ; cobalt%jno3=0.0
    allocate(cobalt%jno3_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jno3_plus_btm=0.0
    allocate(cobalt%jo2(isd:ied, jsd:jed, 1:nk))          ; cobalt%jo2=0.0
    allocate(cobalt%jo2_plus_btm(isd:ied, jsd:jed, 1:nk)) ; cobalt%jo2_plus_btm=0.0
    allocate(cobalt%jpdet(isd:ied, jsd:jed, 1:nk))        ; cobalt%jpdet=0.0
    allocate(cobalt%jpo4(isd:ied, jsd:jed, 1:nk))         ; cobalt%jpo4=0.0
    allocate(cobalt%jpo4_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jpo4_plus_btm=0.0
    allocate(cobalt%jsrdon(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsrdon=0.0
    allocate(cobalt%jsrdop(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsrdop=0.0
    allocate(cobalt%jsldon(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsldon=0.0
    allocate(cobalt%jsldop(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsldop=0.0
    allocate(cobalt%jsidet(isd:ied, jsd:jed, 1:nk))       ; cobalt%jsidet=0.0
    allocate(cobalt%jsilg(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsilg=0.0
    allocate(cobalt%jsio4(isd:ied, jsd:jed, 1:nk))        ; cobalt%jsio4=0.0
    allocate(cobalt%jsio4_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jsio4_plus_btm=0.0
    allocate(cobalt%jprod_fed(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_fed=0.0
    allocate(cobalt%jprod_fedet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_fedet=0.0
    allocate(cobalt%jprod_ndet(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ndet=0.0
    allocate(cobalt%jprod_pdet(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_pdet=0.0
    allocate(cobalt%jprod_ldon(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ldon=0.0
    allocate(cobalt%jprod_ldop(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_ldop=0.0
    allocate(cobalt%jprod_sldon(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_sldon=0.0
    allocate(cobalt%jprod_sldop(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_sldop=0.0
    allocate(cobalt%jprod_srdon(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_srdon=0.0
    allocate(cobalt%jprod_srdop(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_srdop=0.0
    allocate(cobalt%jprod_sidet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jprod_sidet=0.0
    allocate(cobalt%jprod_sio4(isd:ied, jsd:jed, 1:nk))   ; cobalt%jprod_sio4=0.0
    allocate(cobalt%jprod_lithdet(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_lithdet=0.0
    allocate(cobalt%jprod_cadet_arag(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_cadet_arag=0.0
    allocate(cobalt%jprod_cadet_calc(isd:ied, jsd:jed, 1:nk)); cobalt%jprod_cadet_calc=0.0
    allocate(cobalt%jprod_nh4(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_nh4=0.0
    allocate(cobalt%jprod_nh4_plus_btm(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_nh4_plus_btm=0.0
    allocate(cobalt%jprod_po4(isd:ied, jsd:jed, 1:nk))    ; cobalt%jprod_po4=0.0
    allocate(cobalt%det_jzloss_n(isd:ied, jsd:jed, 1:nk)) ; cobalt%det_jzloss_n=0.0
    allocate(cobalt%det_jzloss_p(isd:ied, jsd:jed, 1:nk)) ; cobalt%det_jzloss_p=0.0
    allocate(cobalt%det_jzloss_fe(isd:ied, jsd:jed, 1:nk)); cobalt%det_jzloss_fe=0.0
    allocate(cobalt%det_jzloss_si(isd:ied, jsd:jed, 1:nk)); cobalt%det_jzloss_si=0.0
    allocate(cobalt%det_jhploss_n(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_n=0.0
    allocate(cobalt%det_jhploss_p(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_p=0.0
    allocate(cobalt%det_jhploss_fe(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_fe=0.0
    allocate(cobalt%det_jhploss_si(isd:ied, jsd:jed, 1:nk)); cobalt%det_jhploss_si=0.0
    allocate(cobalt%jdiss_cadet_arag(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_arag=0.0
    allocate(cobalt%jdiss_cadet_arag_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_arag_plus_btm=0.0
    allocate(cobalt%jdiss_cadet_calc(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_calc=0.0
    allocate(cobalt%jdiss_cadet_calc_plus_btm(isd:ied, jsd:jed, 1:nk)); cobalt%jdiss_cadet_calc_plus_btm=0.0
    allocate(cobalt%jdiss_sidet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jdiss_sidet=0.0
    allocate(cobalt%jremin_ndet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jremin_ndet=0.0
    allocate(cobalt%jremin_pdet(isd:ied, jsd:jed, 1:nk))  ; cobalt%jremin_pdet=0.0
    allocate(cobalt%jremin_fedet(isd:ied, jsd:jed, 1:nk)) ; cobalt%jremin_fedet=0.0
    allocate(cobalt%jfe_ads(isd:ied, jsd:jed, 1:nk))      ; cobalt%jfe_ads=0.0
    allocate(cobalt%jfe_coast(isd:ied, jsd:jed, 1:nk))    ; cobalt%jfe_coast=0.0
    allocate(cobalt%kfe_eq_lig(isd:ied, jsd:jed, 1:nk))   ; cobalt%kfe_eq_lig=0.0
    allocate(cobalt%feprime(isd:ied, jsd:jed, 1:nk))      ; cobalt%feprime=0.0
    allocate(cobalt%ligand(isd:ied, jsd:jed, 1:nk))       ; cobalt%ligand=0.0
    allocate(cobalt%fe_sol(isd:ied, jsd:jed, 1:nk))       ; cobalt%fe_sol=0.0
    allocate(cobalt%expkT(isd:ied, jsd:jed, 1:nk))        ; cobalt%expkT=0.0
    allocate(cobalt%expkreminT(isd:ied, jsd:jed, 1:nk))   ; cobalt%expkreminT=0.0
    allocate(cobalt%hp_o2lim(isd:ied, jsd:jed, 1:nk))     ; cobalt%hp_o2lim=0.0
    allocate(cobalt%hp_temp_lim(isd:ied, jsd:jed, 1:nk))  ; cobalt%hp_temp_lim=0.0
    allocate(cobalt%irr_inst(isd:ied, jsd:jed, 1:nk))     ; cobalt%irr_inst=0.0
    allocate(cobalt%irr_mix(isd:ied, jsd:jed, 1:nk))      ; cobalt%irr_mix=0.0
    allocate(cobalt%jno3denit_wc(isd:ied, jsd:jed, 1:nk)) ; cobalt%jno3denit_wc=0.0
    allocate(cobalt%jo2resp_wc(isd:ied, jsd:jed, 1:nk))   ; cobalt%jo2resp_wc=0.0
    allocate(cobalt%jnitrif(isd:ied, jsd:jed, 1:nk))      ; cobalt%jnitrif=0.0
    allocate(cobalt%omega_arag(isd:ied, jsd:jed, 1:nk))   ; cobalt%omega_arag=0.0
    allocate(cobalt%omega_calc(isd:ied, jsd:jed, 1:nk))   ; cobalt%omega_calc=0.0
    allocate(cobalt%omegaa(isd:ied, jsd:jed, 1:nk))       ; cobalt%omegaa=0.0
    allocate(cobalt%omegac(isd:ied, jsd:jed, 1:nk))       ; cobalt%omegac=0.0
    allocate(cobalt%tot_layer_int_c(isd:ied, jsd:jed,1:nk))  ; cobalt%tot_layer_int_c=0.0
    allocate(cobalt%tot_layer_int_fe(isd:ied, jsd:jed,1:nk)) ; cobalt%tot_layer_int_fe=0.0
    allocate(cobalt%tot_layer_int_n(isd:ied, jsd:jed, 1:nk)) ; cobalt%tot_layer_int_n=0.0
    allocate(cobalt%tot_layer_int_p(isd:ied, jsd:jed, 1:nk)) ; cobalt%tot_layer_int_p=0.0
    allocate(cobalt%tot_layer_int_si(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_si=0.0
    allocate(cobalt%tot_layer_int_o2(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_o2=0.0
    allocate(cobalt%tot_layer_int_alk(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_alk=0.0
    allocate(cobalt%total_filter_feeding(isd:ied,jsd:jed,1:nk)); cobalt%total_filter_feeding=0.0
    allocate(cobalt%nlg_diatoms(isd:ied, jsd:jed, 1:nk)); cobalt%nlg_diatoms=0.0
    allocate(cobalt%q_si_2_n_lg_diatoms(isd:ied, jsd:jed, 1:nk)); cobalt%q_si_2_n_lg_diatoms=0.0
    allocate(cobalt%zt(isd:ied, jsd:jed, 1:nk))           ; cobalt%zt=0.0
    allocate(cobalt%zm(isd:ied, jsd:jed, 1:nk))           ; cobalt%zm=0.0
    allocate(cobalt%b_alk(isd:ied, jsd:jed))              ; cobalt%b_alk=0.0
    allocate(cobalt%b_dic(isd:ied, jsd:jed))              ; cobalt%b_dic=0.0
    allocate(cobalt%b_fed(isd:ied, jsd:jed))              ; cobalt%b_fed=0.0
    allocate(cobalt%b_nh4(isd:ied, jsd:jed))              ; cobalt%b_nh4=0.0
    allocate(cobalt%b_no3(isd:ied, jsd:jed))              ; cobalt%b_no3=0.0
    allocate(cobalt%b_o2(isd:ied, jsd:jed))               ; cobalt%b_o2=0.0
    allocate(cobalt%b_po4(isd:ied, jsd:jed))              ; cobalt%b_po4=0.0
    allocate(cobalt%b_sio4(isd:ied, jsd:jed))             ; cobalt%b_sio4=0.0
    allocate(cobalt%pco2_csurf(isd:ied, jsd:jed))         ; cobalt%pco2_csurf=0.0
    allocate(cobalt%co2_csurf(isd:ied, jsd:jed))          ; cobalt%co2_csurf=0.0
    allocate(cobalt%co2_alpha(isd:ied, jsd:jed))          ; cobalt%co2_alpha=0.0
    allocate(cobalt%nh3_csurf(isd:ied, jsd:jed))          ; cobalt%nh3_csurf=0.0
    allocate(cobalt%nh3_alpha(isd:ied, jsd:jed))          ; cobalt%nh3_alpha=0.0
    allocate(cobalt%pnh3_csurf(isd:ied, jsd:jed))         ; cobalt%pnh3_csurf=0.0
    allocate(cobalt%fcadet_arag_btm(isd:ied, jsd:jed))    ; cobalt%fcadet_arag_btm=0.0
    allocate(cobalt%fcadet_calc_btm(isd:ied, jsd:jed))    ; cobalt%fcadet_calc_btm=0.0
    allocate(cobalt%ffedet_btm(isd:ied, jsd:jed))         ; cobalt%ffedet_btm=0.0
    allocate(cobalt%flithdet_btm(isd:ied, jsd:jed))       ; cobalt%flithdet_btm=0.0
    allocate(cobalt%fpdet_btm(isd:ied, jsd:jed))          ; cobalt%fpdet_btm=0.0
    allocate(cobalt%fndet_btm(isd:ied, jsd:jed))          ; cobalt%fndet_btm=0.0
    allocate(cobalt%fsidet_btm(isd:ied, jsd:jed))         ; cobalt%fsidet_btm=0.0
    allocate(cobalt%fcased_burial(isd:ied, jsd:jed))      ; cobalt%fcased_burial=0.0
    allocate(cobalt%fcased_redis(isd:ied, jsd:jed))       ; cobalt%fcased_redis=0.0  
    allocate(cobalt%fcased_redis_surfresp(isd:ied, jsd:jed)) ; cobalt%fcased_redis_surfresp=0.0
    allocate(cobalt%cased_redis_coef(isd:ied, jsd:jed))   ; cobalt%cased_redis_coef=0.0
    allocate(cobalt%cased_redis_delz(isd:ied, jsd:jed))   ; cobalt%cased_redis_delz=0.0
    allocate(cobalt%ffe_sed(isd:ied, jsd:jed))            ; cobalt%ffe_sed=0.0
    allocate(cobalt%ffe_geotherm(isd:ied, jsd:jed))       ; cobalt%ffe_geotherm=0.0
    allocate(cobalt%ffe_iceberg(isd:ied, jsd:jed))        ; cobalt%ffe_iceberg=0.0
    allocate(cobalt%fnfeso4red_sed(isd:ied, jsd:jed))     ; cobalt%fnfeso4red_sed=0.0
    allocate(cobalt%fno3denit_sed(isd:ied, jsd:jed))      ; cobalt%fno3denit_sed=0.0
    allocate(cobalt%fnoxic_sed(isd:ied, jsd:jed))         ; cobalt%fnoxic_sed=0.0
    allocate(cobalt%frac_burial(isd:ied, jsd:jed))        ; cobalt%frac_burial=0.0
    allocate(cobalt%fndet_burial(isd:ied, jsd:jed))       ; cobalt%fndet_burial=0.0
    allocate(cobalt%fpdet_burial(isd:ied, jsd:jed))       ; cobalt%fpdet_burial=0.0
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem 
    allocate(cobalt%dissoc(isd:ied, jsd:jed, 1:nk))        ; cobalt%dissoc=0.0
    allocate(cobalt%o2sat(isd:ied, jsd:jed, 1:nk))         ; cobalt%o2sat=0.0
    allocate(cobalt%remoc(isd:ied, jsd:jed, 1:nk))         ; cobalt%remoc=0.0
    allocate(cobalt%tot_layer_int_doc(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_doc=0.0
    allocate(cobalt%tot_layer_int_poc(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_poc=0.0
    allocate(cobalt%tot_layer_int_dic(isd:ied, jsd:jed, 1:nk)); cobalt%tot_layer_int_dic=0.0
    allocate(cobalt%f_alk_int_100(isd:ied, jsd:jed))       ; cobalt%f_alk_int_100=0.0
    allocate(cobalt%f_dic_int_100(isd:ied, jsd:jed))       ; cobalt%f_dic_int_100=0.0
    allocate(cobalt%f_din_int_100(isd:ied, jsd:jed))       ; cobalt%f_din_int_100=0.0
    allocate(cobalt%f_fed_int_100(isd:ied, jsd:jed))       ; cobalt%f_fed_int_100=0.0
    allocate(cobalt%f_po4_int_100(isd:ied, jsd:jed))       ; cobalt%f_po4_int_100=0.0
    allocate(cobalt%f_sio4_int_100(isd:ied, jsd:jed))      ; cobalt%f_sio4_int_100=0.0
    allocate(cobalt%jalk_100(isd:ied, jsd:jed))            ; cobalt%jalk_100=0.0
    allocate(cobalt%jdic_100(isd:ied, jsd:jed))            ; cobalt%jdic_100=0.0
    allocate(cobalt%jdin_100(isd:ied, jsd:jed))            ; cobalt%jdin_100=0.0
    allocate(cobalt%jfed_100(isd:ied, jsd:jed))            ; cobalt%jfed_100=0.0
    allocate(cobalt%jpo4_100(isd:ied, jsd:jed))            ; cobalt%jpo4_100=0.0
    allocate(cobalt%jsio4_100(isd:ied, jsd:jed))           ; cobalt%jsio4_100=0.0
    allocate(cobalt%jprod_ptot_100(isd:ied, jsd:jed))      ; cobalt%jprod_ptot_100=0.0
    allocate(cobalt%wc_vert_int_c(isd:ied, jsd:jed))       ; cobalt%wc_vert_int_c=0.0
    allocate(cobalt%wc_vert_int_dic(isd:ied, jsd:jed))        ; cobalt%wc_vert_int_dic=0.0
    allocate(cobalt%wc_vert_int_doc(isd:ied, jsd:jed))        ; cobalt%wc_vert_int_doc=0.0
    allocate(cobalt%wc_vert_int_poc(isd:ied, jsd:jed))        ; cobalt%wc_vert_int_poc=0.0
    allocate(cobalt%wc_vert_int_n(isd:ied, jsd:jed))          ; cobalt%wc_vert_int_n=0.0
    allocate(cobalt%wc_vert_int_p(isd:ied, jsd:jed))          ; cobalt%wc_vert_int_p=0.0
    allocate(cobalt%wc_vert_int_fe(isd:ied, jsd:jed))         ; cobalt%wc_vert_int_fe=0.0
    allocate(cobalt%wc_vert_int_si(isd:ied, jsd:jed))         ; cobalt%wc_vert_int_si=0.0
    allocate(cobalt%wc_vert_int_o2(isd:ied, jsd:jed))         ; cobalt%wc_vert_int_o2=0.0
    allocate(cobalt%wc_vert_int_alk(isd:ied, jsd:jed))        ; cobalt%wc_vert_int_alk=0.0
    allocate(cobalt%wc_vert_int_jdiss_sidet(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_jdiss_sidet=0.0
    allocate(cobalt%wc_vert_int_jdiss_cadet(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_jdiss_cadet=0.0
    allocate(cobalt%wc_vert_int_jo2resp(isd:ied, jsd:jed))      ; cobalt%wc_vert_int_jo2resp=0.0
    allocate(cobalt%wc_vert_int_jprod_cadet(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_jprod_cadet=0.0
    allocate(cobalt%wc_vert_int_jno3denit(isd:ied, jsd:jed))    ; cobalt%wc_vert_int_jno3denit=0.0
    allocate(cobalt%wc_vert_int_jnitrif(isd:ied, jsd:jed))      ; cobalt%wc_vert_int_jnitrif=0.0
    allocate(cobalt%wc_vert_int_juptake_nh4(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_juptake_nh4=0.0
    allocate(cobalt%wc_vert_int_jprod_nh4(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_jprod_nh4=0.0
    allocate(cobalt%wc_vert_int_juptake_no3(isd:ied, jsd:jed))  ; cobalt%wc_vert_int_juptake_no3=0.0
    allocate(cobalt%wc_vert_int_nfix(isd:ied, jsd:jed))         ; cobalt%wc_vert_int_nfix=0.0
!==============================================================================================================
    !
    ! allocate 100m integrated quantities
    !
    do n = 1, NUM_PHYTO
       allocate(phyto(n)%jprod_n_100(isd:ied,jsd:jed))      ; phyto(n)%jprod_n_100      = 0.0
       allocate(phyto(n)%jprod_n_new_100(isd:ied,jsd:jed))  ; phyto(n)%jprod_n_new_100  = 0.0
       allocate(phyto(n)%jzloss_n_100(isd:ied,jsd:jed))     ; phyto(n)%jzloss_n_100  = 0.0
       allocate(phyto(n)%jexuloss_n_100(isd:ied,jsd:jed))   ; phyto(n)%jexuloss_n_100  = 0.0
       allocate(phyto(n)%f_n_100(isd:ied,jsd:jed))          ; phyto(n)%f_n_100  = 0.0
       allocate(phyto(n)%juptake_fe_100(isd:ied,jsd:jed))   ; phyto(n)%juptake_fe_100  = 0.0
       allocate(phyto(n)%juptake_po4_100(isd:ied,jsd:jed))  ; phyto(n)%juptake_po4_100  = 0.0
       ! Biomass-weighted limitation terms (for cmip)
       allocate(phyto(n)%nlim_bw_100(isd:ied,jsd:jed)) ; phyto(n)%nlim_bw_100 = 0.0
       allocate(phyto(n)%plim_bw_100(isd:ied,jsd:jed)) ; phyto(n)%plim_bw_100 = 0.0
       allocate(phyto(n)%irrlim_bw_100(isd:ied,jsd:jed)) ; phyto(n)%irrlim_bw_100 = 0.0
       allocate(phyto(n)%def_fe_bw_100(isd:ied,jsd:jed)) ; phyto(n)%def_fe_bw_100 = 0.0
    enddo
    allocate(phyto(DIAZO)%jprod_n_n2_100(isd:ied,jsd:jed)); phyto(DIAZO)%jprod_n_n2_100 = 0.0
    allocate(phyto(SMALL)%jvirloss_n_100(isd:ied,jsd:jed))  ; phyto(SMALL)%jvirloss_n_100 = 0.0
    allocate(phyto(SMALL)%jaggloss_n_100(isd:ied,jsd:jed))  ; phyto(SMALL)%jaggloss_n_100 = 0.0
    allocate(phyto(LARGE)%jaggloss_n_100(isd:ied,jsd:jed))  ; phyto(LARGE)%jaggloss_n_100 = 0.0
    allocate(cobalt%jprod_allphytos_100(isd:ied,jsd:jed))   ; cobalt%jprod_allphytos_100 = 0.0
    allocate(cobalt%jprod_diat_100(isd:ied,jsd:jed))   ; cobalt%jprod_diat_100 = 0.0
    allocate(phyto(LARGE)%juptake_sio4_100(isd:ied,jsd:jed)) ; phyto(LARGE)%juptake_sio4_100 = 0.0

    do n = 1, NUM_ZOO
       allocate(zoo(n)%jprod_n_100(isd:ied,jsd:jed))      ; zoo(n)%jprod_n_100      = 0.0
       allocate(zoo(n)%jingest_n_100(isd:ied,jsd:jed))    ; zoo(n)%jingest_n_100    = 0.0
       allocate(zoo(n)%jremin_n_100(isd:ied,jsd:jed))     ; zoo(n)%jremin_n_100     = 0.0
       allocate(zoo(n)%f_n_100(isd:ied,jsd:jed))          ; zoo(n)%f_n_100          = 0.0
    enddo

   do n = 1,2
       allocate(zoo(n)%jzloss_n_100(isd:ied,jsd:jed))     ; zoo(n)%jzloss_n_100     = 0.0
       allocate(zoo(n)%jprod_don_100(isd:ied,jsd:jed))    ; zoo(n)%jprod_don_100    = 0.0
   enddo

   do n = 2,3
       allocate(zoo(n)%jhploss_n_100(isd:ied,jsd:jed))    ; zoo(n)%jhploss_n_100     = 0.0
       allocate(zoo(n)%jprod_ndet_100(isd:ied,jsd:jed))   ; zoo(n)%jprod_ndet_100    = 0.0
   enddo

   allocate(cobalt%hp_jingest_n_100(isd:ied,jsd:jed))    ; cobalt%hp_jingest_n_100    = 0.0
   allocate(cobalt%hp_jremin_n_100(isd:ied,jsd:jed))     ; cobalt%hp_jremin_n_100     = 0.0
   allocate(cobalt%hp_jprod_ndet_100(isd:ied,jsd:jed))   ; cobalt%hp_jprod_ndet_100   = 0.0

   allocate(bact(1)%jprod_n_100(isd:ied,jsd:jed))   ; bact(1)%jprod_n_100 = 0.0
   allocate(bact(1)%jzloss_n_100(isd:ied,jsd:jed))  ; bact(1)%jzloss_n_100 = 0.0
   allocate(bact(1)%jvirloss_n_100(isd:ied,jsd:jed)); bact(1)%jvirloss_n_100 = 0.0
   allocate(bact(1)%jremin_n_100(isd:ied,jsd:jed))  ; bact(1)%jremin_n_100 = 0.0
   allocate(bact(1)%juptake_ldon_100(isd:ied,jsd:jed)) ; bact(1)%juptake_ldon_100 = 0.0
   allocate(bact(1)%f_n_100(isd:ied,jsd:jed))       ; bact(1)%f_n_100 = 0.0

   allocate(cobalt%jprod_lithdet_100(isd:ied,jsd:jed))      ; cobalt%jprod_lithdet_100 = 0.0
   allocate(cobalt%jprod_sidet_100(isd:ied,jsd:jed))        ; cobalt%jprod_sidet_100 = 0.0
   allocate(cobalt%jprod_cadet_calc_100(isd:ied,jsd:jed))   ; cobalt%jprod_cadet_calc_100 = 0.0
   allocate(cobalt%jprod_cadet_arag_100(isd:ied,jsd:jed))   ; cobalt%jprod_cadet_arag_100 = 0.0
   allocate(cobalt%jremin_ndet_100(isd:ied,jsd:jed))        ; cobalt%jremin_ndet_100 = 0.0
   allocate(cobalt%jprod_mesozoo_200(isd:ied,jsd:jed))      ; cobalt%jprod_mesozoo_200 = 0.0

   allocate(cobalt%f_ndet_100(isd:ied,jsd:jed))             ; cobalt%f_ndet_100 = 0.0
   allocate(cobalt%f_don_100(isd:ied,jsd:jed))              ; cobalt%f_don_100  = 0.0
   allocate(cobalt%f_silg_100(isd:ied,jsd:jed))             ; cobalt%f_silg_100 = 0.0
   allocate(cobalt%f_mesozoo_200(isd:ied,jsd:jed))          ; cobalt%f_mesozoo_200 = 0.0

   allocate(cobalt%fndet_100(isd:ied,jsd:jed))             ; cobalt%fndet_100 = 0.0
   allocate(cobalt%fpdet_100(isd:ied,jsd:jed))             ; cobalt%fpdet_100 = 0.0
   allocate(cobalt%fsidet_100(isd:ied,jsd:jed))            ; cobalt%fsidet_100 = 0.0
   allocate(cobalt%flithdet_100(isd:ied,jsd:jed))          ; cobalt%flithdet_100 = 0.0
   allocate(cobalt%fcadet_calc_100(isd:ied,jsd:jed))       ; cobalt%fcadet_calc_100 = 0.0
   allocate(cobalt%fcadet_arag_100(isd:ied,jsd:jed))       ; cobalt%fcadet_arag_100 = 0.0
   allocate(cobalt%ffedet_100(isd:ied,jsd:jed))            ; cobalt%ffedet_100 = 0.0

   allocate(cobalt%btm_temp(isd:ied,jsd:jed))              ; cobalt%btm_temp = 0.0
   allocate(cobalt%btm_o2(isd:ied,jsd:jed))                ; cobalt%btm_o2 = 0.0
   allocate(cobalt%btm_htotal(isd:ied,jsd:jed))                ; cobalt%btm_htotal = 0.0
   allocate(cobalt%btm_co3_ion(isd:ied,jsd:jed))           ; cobalt%btm_co3_ion = 0.0
   allocate(cobalt%btm_co3_sol_arag(isd:ied,jsd:jed))      ; cobalt%btm_co3_sol_arag = 0.0
   allocate(cobalt%btm_co3_sol_calc(isd:ied,jsd:jed))      ; cobalt%btm_co3_sol_calc = 0.0
   allocate(cobalt%cased_2d(isd:ied,jsd:jed))              ; cobalt%cased_2d = 0.0

   allocate(cobalt%o2min(isd:ied, jsd:jed))                ; cobalt%o2min=0.0
   allocate(cobalt%z_o2min(isd:ied, jsd:jed))              ; cobalt%z_o2min=0.0
   allocate(cobalt%z_sat_arag(isd:ied, jsd:jed))           ; cobalt%z_sat_arag=0.0
   allocate(cobalt%z_sat_calc(isd:ied, jsd:jed))           ; cobalt%z_sat_calc=0.0
   allocate(cobalt%mask_z_sat_arag(isd:ied, jsd:jed))      ; cobalt%mask_z_sat_arag = .FALSE.
   allocate(cobalt%mask_z_sat_calc(isd:ied, jsd:jed))      ; cobalt%mask_z_sat_calc = .FALSE.
   if (do_14c) then                                        !<<RADIOCARBON
      allocate(cobalt%c14_2_n(isd:ied, jsd:jed, 1:nk));        cobalt%c14_2_n=0.0
      allocate(cobalt%f_di14c(isd:ied, jsd:jed, 1:nk));        cobalt%f_di14c=0.0
      allocate(cobalt%f_do14c(isd:ied, jsd:jed, 1:nk));        cobalt%f_do14c=0.0
      allocate(cobalt%fpo14c(isd:ied, jsd:jed, 1:nk));         cobalt%fpo14c=0.0
      allocate(cobalt%j14c_decay_dic(isd:ied, jsd:jed, 1:nk)); cobalt%j14c_decay_dic=0.0
      allocate(cobalt%j14c_decay_doc(isd:ied, jsd:jed, 1:nk)); cobalt%j14c_decay_doc=0.0
      allocate(cobalt%j14c_reminp(isd:ied, jsd:jed, 1:nk));    cobalt%j14c_reminp=0.0
      allocate(cobalt%jdi14c(isd:ied, jsd:jed, 1:nk));         cobalt%jdi14c=0.0
      allocate(cobalt%jdo14c(isd:ied, jsd:jed, 1:nk));         cobalt%jdo14c=0.0
      allocate(cobalt%c14o2_csurf  (isd:ied, jsd:jed));        cobalt%c14o2_csurf=0.0
      allocate(cobalt%c14o2_alpha  (isd:ied, jsd:jed));        cobalt%c14o2_alpha=0.0
      allocate(cobalt%b_di14c      (isd:ied, jsd:jed));        cobalt%b_di14c=0.0
   endif                                                   !RADIOCARBON>>
      allocate(cobalt%runoff_flux_alk(isd:ied, jsd:jed));      cobalt%runoff_flux_alk=0.0
      allocate(cobalt%runoff_flux_dic(isd:ied, jsd:jed));      cobalt%runoff_flux_dic=0.0
      allocate(cobalt%runoff_flux_di14c(isd:ied, jsd:jed));    cobalt%runoff_flux_di14c=0.0
      allocate(cobalt%runoff_flux_lith(isd:ied, jsd:jed));     cobalt%runoff_flux_lith=0.0
      allocate(cobalt%runoff_flux_fed(isd:ied, jsd:jed));      cobalt%runoff_flux_fed=0.0
      allocate(cobalt%runoff_flux_no3(isd:ied, jsd:jed));      cobalt%runoff_flux_no3=0.0
      allocate(cobalt%runoff_flux_ldon(isd:ied, jsd:jed));     cobalt%runoff_flux_ldon=0.0
      allocate(cobalt%runoff_flux_sldon(isd:ied, jsd:jed));    cobalt%runoff_flux_sldon=0.0
      allocate(cobalt%runoff_flux_srdon(isd:ied, jsd:jed));    cobalt%runoff_flux_srdon=0.0
      allocate(cobalt%runoff_flux_ndet(isd:ied, jsd:jed));     cobalt%runoff_flux_ndet=0.0
      allocate(cobalt%runoff_flux_po4(isd:ied, jsd:jed));      cobalt%runoff_flux_po4=0.0
      allocate(cobalt%runoff_flux_ldop(isd:ied, jsd:jed));     cobalt%runoff_flux_ldop=0.0
      allocate(cobalt%runoff_flux_sldop(isd:ied, jsd:jed));    cobalt%runoff_flux_sldop=0.0
      allocate(cobalt%runoff_flux_srdop(isd:ied, jsd:jed));    cobalt%runoff_flux_srdop=0.0
      allocate(cobalt%dry_fed(isd:ied, jsd:jed));              cobalt%dry_fed=0.0
      allocate(cobalt%wet_fed(isd:ied, jsd:jed));              cobalt%wet_fed=0.0
      allocate(cobalt%dry_lith(isd:ied, jsd:jed));             cobalt%dry_lith=0.0
      allocate(cobalt%wet_lith(isd:ied, jsd:jed));             cobalt%wet_lith=0.0
      allocate(cobalt%dry_no3(isd:ied, jsd:jed));              cobalt%dry_no3=0.0
      allocate(cobalt%wet_no3(isd:ied, jsd:jed));              cobalt%wet_no3=0.0
      allocate(cobalt%dry_nh4(isd:ied, jsd:jed));              cobalt%dry_nh4=0.0
      allocate(cobalt%wet_nh4(isd:ied, jsd:jed));              cobalt%wet_nh4=0.0
      allocate(cobalt%dry_po4(isd:ied, jsd:jed));              cobalt%dry_po4=0.0
      allocate(cobalt%wet_po4(isd:ied, jsd:jed));              cobalt%wet_po4=0.0
      allocate(cobalt%stf_gas_dic(isd:ied, jsd:jed));          cobalt%stf_gas_dic=0.0
      allocate(cobalt%stf_gas_o2(isd:ied, jsd:jed));           cobalt%stf_gas_o2=0.0
      allocate(cobalt%deltap_dic(isd:ied, jsd:jed));           cobalt%deltap_dic=0.0
      allocate(cobalt%deltap_o2(isd:ied, jsd:jed));            cobalt%deltap_o2=0.0


  end subroutine user_allocate_arrays

  !
  !   This is an internal sub, not a public interface.
  !   Deallocate all the work arrays allocated by user_allocate_arrays.
  !
  subroutine user_deallocate_arrays
    integer n

    deallocate(cobalt%htotalhi,cobalt%htotallo)

    do n = 1, NUM_PHYTO
       deallocate(phyto(n)%def_fe)
       deallocate(phyto(n)%def_p)
       deallocate(phyto(n)%f_fe)
       deallocate(phyto(n)%f_n)
       deallocate(phyto(n)%felim)
       deallocate(phyto(n)%irrlim)
       deallocate(phyto(n)%jzloss_fe)
       deallocate(phyto(n)%jzloss_n)
       deallocate(phyto(n)%jzloss_p)
       deallocate(phyto(n)%jzloss_sio2)
       deallocate(phyto(n)%jaggloss_n) 
       deallocate(phyto(n)%jaggloss_p)
       deallocate(phyto(n)%jaggloss_fe)
       deallocate(phyto(n)%jaggloss_sio2)
       deallocate(phyto(n)%jvirloss_n)
       deallocate(phyto(n)%jvirloss_p)
       deallocate(phyto(n)%jvirloss_fe)
       deallocate(phyto(n)%jvirloss_sio2)
       deallocate(phyto(n)%jexuloss_n)
       deallocate(phyto(n)%jexuloss_p)
       deallocate(phyto(n)%jexuloss_fe)
       deallocate(phyto(n)%jhploss_fe)
       deallocate(phyto(n)%jhploss_n)
       deallocate(phyto(n)%jhploss_p)
       deallocate(phyto(n)%juptake_fe)
       deallocate(phyto(n)%juptake_nh4)
       deallocate(phyto(n)%juptake_no3)
       deallocate(phyto(n)%juptake_po4)
       deallocate(phyto(n)%jprod_n)
       deallocate(phyto(n)%liebig_lim)
       deallocate(phyto(n)%mu)
       deallocate(phyto(n)%po4lim)
       deallocate(phyto(n)%q_fe_2_n)
       deallocate(phyto(n)%q_p_2_n)
       deallocate(phyto(n)%q_si_2_n)
       deallocate(phyto(n)%theta)
       deallocate(phyto(n)%f_mu_mem)
       deallocate(phyto(n)%mu_mix)
       deallocate(phyto(n)%agg_lim)
       deallocate(phyto(n)%juptake_fe_100)
       deallocate(phyto(n)%juptake_po4_100)
       deallocate(phyto(n)%nh4lim)
       deallocate(phyto(n)%no3lim)
       deallocate(phyto(n)%nlim_bw_100)
       deallocate(phyto(n)%plim_bw_100)
       deallocate(phyto(n)%irrlim_bw_100)
       deallocate(phyto(n)%def_fe_bw_100)
    enddo
    deallocate(phyto(DIAZO)%juptake_n2)
    deallocate(phyto(DIAZO)%o2lim)
    deallocate(phyto(LARGE)%juptake_sio4)
    deallocate(phyto(LARGE)%juptake_sio4_100)
    deallocate(phyto(LARGE)%silim)

    ! bacteria
    deallocate(bact(1)%f_n)
    deallocate(bact(1)%jzloss_n)
    deallocate(bact(1)%jzloss_p)
    deallocate(bact(1)%jvirloss_n) 
    deallocate(bact(1)%jvirloss_p)   
    deallocate(bact(1)%jhploss_n)  
    deallocate(bact(1)%jhploss_p)
    deallocate(bact(1)%juptake_ldon)
    deallocate(bact(1)%juptake_ldop)
    deallocate(bact(1)%jprod_nh4)
    deallocate(bact(1)%jprod_po4)
    deallocate(bact(1)%jprod_n)
    deallocate(bact(1)%o2lim)
    deallocate(bact(1)%ldonlim)
    deallocate(bact(1)%temp_lim)

    ! zooplankton
    do n = 1, NUM_ZOO
       deallocate(zoo(n)%f_n)
       deallocate(zoo(n)%jzloss_n)
       deallocate(zoo(n)%jzloss_p)
       deallocate(zoo(n)%jhploss_n)
       deallocate(zoo(n)%jhploss_p)
       deallocate(zoo(n)%jingest_n)
       deallocate(zoo(n)%jingest_p)
       deallocate(zoo(n)%jingest_sio2)
       deallocate(zoo(n)%jingest_fe)
       deallocate(zoo(n)%jprod_ndet)
       deallocate(zoo(n)%jprod_pdet)
       deallocate(zoo(n)%jprod_ldon)
       deallocate(zoo(n)%jprod_ldop)
       deallocate(zoo(n)%jprod_srdon)
       deallocate(zoo(n)%jprod_srdop)
       deallocate(zoo(n)%jprod_sldon)
       deallocate(zoo(n)%jprod_sldop)
       deallocate(zoo(n)%jprod_fedet)
       deallocate(zoo(n)%jprod_fed)
       deallocate(zoo(n)%jprod_sidet)
       deallocate(zoo(n)%jprod_sio4)
       deallocate(zoo(n)%jprod_po4)
       deallocate(zoo(n)%jprod_nh4)
       deallocate(zoo(n)%jprod_n)
       deallocate(zoo(n)%o2lim)
       deallocate(zoo(n)%temp_lim)
    enddo

    deallocate(cobalt%f_alk)
    deallocate(cobalt%f_cadet_arag)
    deallocate(cobalt%f_cadet_calc)  
    deallocate(cobalt%f_dic)  
    deallocate(cobalt%f_fed)  
    deallocate(cobalt%f_fedet)  
    deallocate(cobalt%f_ldon)  
    deallocate(cobalt%f_ldop)  
    deallocate(cobalt%f_lith)  
    deallocate(cobalt%f_lithdet)  
    deallocate(cobalt%f_ndet)  
    deallocate(cobalt%f_nh4)  
    deallocate(cobalt%f_no3)  
    deallocate(cobalt%f_o2)  
    deallocate(cobalt%f_pdet)  
    deallocate(cobalt%f_po4)  
    deallocate(cobalt%f_srdon)  
    deallocate(cobalt%f_srdop)  
    deallocate(cobalt%f_sldon)  
    deallocate(cobalt%f_sldop)  
    deallocate(cobalt%f_sidet)  
    deallocate(cobalt%f_silg)  
    deallocate(cobalt%f_sio4)  
    deallocate(cobalt%co3_sol_arag)  
    deallocate(cobalt%co3_sol_calc)  
    deallocate(cobalt%f_chl)  
    if (allocated(cobalt%f_nh3)) deallocate(cobalt%f_nh3)
    deallocate(cobalt%f_co3_ion)  
    deallocate(cobalt%f_htotal)  
    deallocate(cobalt%f_irr_mem)  
    deallocate(cobalt%f_cased)  
    deallocate(cobalt%f_cadet_arag_btf)  
    deallocate(cobalt%f_cadet_calc_btf)  
    deallocate(cobalt%f_fedet_btf)  
    deallocate(cobalt%f_lithdet_btf)  
    deallocate(cobalt%f_ndet_btf)  
    deallocate(cobalt%f_pdet_btf)  
    deallocate(cobalt%f_sidet_btf)  
    deallocate(cobalt%jnbact)  
    deallocate(cobalt%jndi)  
    deallocate(cobalt%jnsm)  
    deallocate(cobalt%jnlg)  
    deallocate(cobalt%jnsmz)  
    deallocate(cobalt%jnmdz)  
    deallocate(cobalt%jnlgz)  
    deallocate(cobalt%jalk)  
    deallocate(cobalt%jalk_plus_btm)  
    deallocate(cobalt%jcadet_arag)  
    deallocate(cobalt%jcadet_calc)  
    deallocate(cobalt%jdic)  
    deallocate(cobalt%jdic_plus_btm)  
    deallocate(cobalt%jdin_plus_btm)  
    deallocate(cobalt%jfed)  
    deallocate(cobalt%jfed_plus_btm)  
    deallocate(cobalt%jfedi)  
    deallocate(cobalt%jfelg)  
    deallocate(cobalt%jfesm)  
    deallocate(cobalt%jfedet)  
    deallocate(cobalt%jldon)  
    deallocate(cobalt%jldop)  
    deallocate(cobalt%jlith)  
    deallocate(cobalt%jlithdet)  
    deallocate(cobalt%jndet)  
    deallocate(cobalt%jnh4)  
    deallocate(cobalt%jnh4_plus_btm)  
    deallocate(cobalt%jno3)  
    deallocate(cobalt%jno3_plus_btm)  
    deallocate(cobalt%jo2)  
    deallocate(cobalt%jo2_plus_btm)  
    deallocate(cobalt%jpdet)  
    deallocate(cobalt%jpo4)  
    deallocate(cobalt%jpo4_plus_btm)  
    deallocate(cobalt%jsrdon)  
    deallocate(cobalt%jsrdop)  
    deallocate(cobalt%jsldon)  
    deallocate(cobalt%jsldop)  
    deallocate(cobalt%jsidet)  
    deallocate(cobalt%jsilg)  
    deallocate(cobalt%jsio4)  
    deallocate(cobalt%jsio4_plus_btm)  
    deallocate(cobalt%jprod_ndet)  
    deallocate(cobalt%jprod_pdet)  
    deallocate(cobalt%jprod_ldon)  
    deallocate(cobalt%jprod_ldop)  
    deallocate(cobalt%jprod_sldop)  
    deallocate(cobalt%jprod_sldon)  
    deallocate(cobalt%jprod_srdon)  
    deallocate(cobalt%jprod_srdop)  
    deallocate(cobalt%jprod_fed)  
    deallocate(cobalt%jprod_fedet)  
    deallocate(cobalt%jprod_sidet)  
    deallocate(cobalt%jprod_sio4)  
    deallocate(cobalt%jprod_lithdet)  
    deallocate(cobalt%jprod_cadet_arag)  
    deallocate(cobalt%jprod_cadet_calc)  
    deallocate(cobalt%jprod_nh4)  
    deallocate(cobalt%jprod_nh4_plus_btm)  
    deallocate(cobalt%jprod_po4)  
    deallocate(cobalt%det_jzloss_n)  
    deallocate(cobalt%det_jzloss_p)  
    deallocate(cobalt%det_jzloss_fe)  
    deallocate(cobalt%det_jzloss_si)  
    deallocate(cobalt%det_jhploss_n)  
    deallocate(cobalt%det_jhploss_p)  
    deallocate(cobalt%det_jhploss_fe)  
    deallocate(cobalt%det_jhploss_si)  
    deallocate(cobalt%jdiss_cadet_arag)  
    deallocate(cobalt%jdiss_cadet_arag_plus_btm)  
    deallocate(cobalt%jdiss_cadet_calc)  
    deallocate(cobalt%jdiss_cadet_calc_plus_btm)  
    deallocate(cobalt%jdiss_sidet)  
    deallocate(cobalt%jremin_ndet)  
    deallocate(cobalt%jremin_pdet)  
    deallocate(cobalt%jremin_fedet)  
    deallocate(cobalt%jfe_ads)  
    deallocate(cobalt%jfe_coast)  
    deallocate(cobalt%kfe_eq_lig)
    deallocate(cobalt%feprime)
    deallocate(cobalt%ligand)
    deallocate(cobalt%fe_sol)  
    deallocate(cobalt%expkT)  
    deallocate(cobalt%expkreminT) 
    deallocate(cobalt%hp_o2lim)
    deallocate(cobalt%hp_temp_lim)  
    deallocate(cobalt%hp_jingest_n)
    deallocate(cobalt%hp_jingest_p)
    deallocate(cobalt%hp_jingest_sio2)
    deallocate(cobalt%hp_jingest_fe)
    deallocate(cobalt%irr_inst)  
    deallocate(cobalt%irr_mix)  
    deallocate(cobalt%jno3denit_wc)  
    deallocate(cobalt%jo2resp_wc)  
    deallocate(cobalt%jnitrif)  
    deallocate(cobalt%omega_arag)  
    deallocate(cobalt%omega_calc)  
    deallocate(cobalt%omegaa)  
    deallocate(cobalt%omegac)  
    deallocate(cobalt%tot_layer_int_c)  
    deallocate(cobalt%tot_layer_int_fe)  
    deallocate(cobalt%tot_layer_int_n)  
    deallocate(cobalt%tot_layer_int_p)  
    deallocate(cobalt%tot_layer_int_si)
    deallocate(cobalt%tot_layer_int_alk)
    deallocate(cobalt%tot_layer_int_o2)  
    deallocate(cobalt%total_filter_feeding)  
    deallocate(cobalt%nlg_diatoms)  
    deallocate(cobalt%q_si_2_n_lg_diatoms)  
    deallocate(cobalt%zt)  
    deallocate(cobalt%zm)  
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem 
    deallocate(cobalt%dissoc)  
    deallocate(cobalt%o2sat)  
    deallocate(cobalt%remoc)  
    deallocate(cobalt%tot_layer_int_doc)  
    deallocate(cobalt%tot_layer_int_poc)  
    deallocate(cobalt%tot_layer_int_dic)

!==============================================================================================================

    deallocate(cobalt%b_alk)  
    deallocate(cobalt%b_dic)  
    deallocate(cobalt%b_fed)  
    deallocate(cobalt%b_nh4)  
    deallocate(cobalt%b_no3)  
    deallocate(cobalt%b_o2)   
    deallocate(cobalt%b_po4)  
    deallocate(cobalt%b_sio4)  
    deallocate(cobalt%pco2_csurf)  
    deallocate(cobalt%co2_csurf)  
    deallocate(cobalt%co2_alpha)  
    deallocate(cobalt%nh3_csurf)  
    deallocate(cobalt%pnh3_csurf)  
    deallocate(cobalt%nh3_alpha)  
    deallocate(cobalt%fcadet_arag_btm)  
    deallocate(cobalt%fcadet_calc_btm)  
    deallocate(cobalt%ffedet_btm)  
    deallocate(cobalt%flithdet_btm)  
    deallocate(cobalt%fpdet_btm)  
    deallocate(cobalt%fndet_btm)  
    deallocate(cobalt%fsidet_btm)  
    deallocate(cobalt%fcased_burial)  
    deallocate(cobalt%fcased_redis)
    deallocate(cobalt%fcased_redis_surfresp)
    deallocate(cobalt%cased_redis_coef)
    deallocate(cobalt%cased_redis_delz)  
    deallocate(cobalt%ffe_sed)
    deallocate(cobalt%ffe_geotherm)
    deallocate(cobalt%ffe_iceberg)  
    deallocate(cobalt%fnfeso4red_sed)  
    deallocate(cobalt%fno3denit_sed)  
    deallocate(cobalt%fnoxic_sed)  
    deallocate(cobalt%frac_burial)  
    deallocate(cobalt%fndet_burial)  
    deallocate(cobalt%fpdet_burial)
    deallocate(cobalt%jprod_allphytos_100)
    deallocate(cobalt%jprod_diat_100)
    deallocate(cobalt%hp_jingest_n_100) 
    deallocate(cobalt%hp_jremin_n_100)  
    deallocate(cobalt%hp_jprod_ndet_100)  
    deallocate(cobalt%jprod_lithdet_100)
    deallocate(cobalt%jprod_sidet_100)
    deallocate(cobalt%jprod_cadet_arag_100)
    deallocate(cobalt%jprod_cadet_calc_100)
    deallocate(cobalt%jprod_mesozoo_200)
    deallocate(cobalt%jremin_ndet_100)
    deallocate(cobalt%f_ndet_100)
    deallocate(cobalt%f_don_100)
    deallocate(cobalt%f_silg_100)
    deallocate(cobalt%f_mesozoo_200)
    deallocate(cobalt%fndet_100)
    deallocate(cobalt%fpdet_100)
    deallocate(cobalt%fsidet_100)
    deallocate(cobalt%fcadet_calc_100)
    deallocate(cobalt%fcadet_arag_100)
    deallocate(cobalt%ffedet_100)
    deallocate(cobalt%flithdet_100)
    deallocate(cobalt%btm_temp)
    deallocate(cobalt%btm_o2)
    deallocate(cobalt%btm_htotal)
    deallocate(cobalt%btm_co3_ion)
    deallocate(cobalt%btm_co3_sol_arag)
    deallocate(cobalt%btm_co3_sol_calc)
    deallocate(cobalt%cased_2d)
    deallocate(cobalt%o2min)
    deallocate(cobalt%z_o2min)
    deallocate(cobalt%z_sat_arag)
    deallocate(cobalt%z_sat_calc)
    deallocate(cobalt%mask_z_sat_arag)
    deallocate(cobalt%mask_z_sat_calc)
!==============================================================================================================
! JGJ 2016/08/08 CMIP6 OcnBgchem 
    deallocate(cobalt%f_alk_int_100)  
    deallocate(cobalt%f_dic_int_100)  
    deallocate(cobalt%f_din_int_100)  
    deallocate(cobalt%f_fed_int_100)  
    deallocate(cobalt%f_po4_int_100)  
    deallocate(cobalt%f_sio4_int_100)  
    deallocate(cobalt%jalk_100)  
    deallocate(cobalt%jdic_100)  
    deallocate(cobalt%jdin_100)  
    deallocate(cobalt%jfed_100)  
    deallocate(cobalt%jpo4_100)  
    deallocate(cobalt%jsio4_100)  
    deallocate(cobalt%jprod_ptot_100)
    deallocate(cobalt%wc_vert_int_c)  
    deallocate(cobalt%wc_vert_int_dic)  
    deallocate(cobalt%wc_vert_int_doc)
    deallocate(cobalt%wc_vert_int_poc)  
    deallocate(cobalt%wc_vert_int_n)
    deallocate(cobalt%wc_vert_int_p)
    deallocate(cobalt%wc_vert_int_fe)
    deallocate(cobalt%wc_vert_int_si)
    deallocate(cobalt%wc_vert_int_o2)
    deallocate(cobalt%wc_vert_int_alk)
    deallocate(cobalt%wc_vert_int_jdiss_sidet)
    deallocate(cobalt%wc_vert_int_jdiss_cadet)
    deallocate(cobalt%wc_vert_int_jo2resp)
    deallocate(cobalt%wc_vert_int_jprod_cadet)
    deallocate(cobalt%wc_vert_int_jno3denit)
    deallocate(cobalt%wc_vert_int_jnitrif)
    deallocate(cobalt%wc_vert_int_juptake_nh4)
    deallocate(cobalt%wc_vert_int_jprod_nh4)
    deallocate(cobalt%wc_vert_int_juptake_no3)
    deallocate(cobalt%wc_vert_int_nfix)  
!==============================================================================================================

    do n = 1, NUM_PHYTO
       deallocate(phyto(n)%jprod_n_100)
       deallocate(phyto(n)%jprod_n_new_100)
       deallocate(phyto(n)%jzloss_n_100)
       deallocate(phyto(n)%jexuloss_n_100)
       deallocate(phyto(n)%f_n_100)
    enddo
    deallocate(phyto(DIAZO)%jprod_n_n2_100)
    deallocate(phyto(SMALL)%jvirloss_n_100)
    deallocate(phyto(SMALL)%jaggloss_n_100)
    deallocate(phyto(LARGE)%jaggloss_n_100)

    do n = 1, NUM_ZOO
       deallocate(zoo(n)%jprod_n_100)
       deallocate(zoo(n)%jingest_n_100)
       deallocate(zoo(n)%jremin_n_100)
       deallocate(zoo(n)%f_n_100)
    enddo

    do n = 1,2
       deallocate(zoo(n)%jzloss_n_100)
       deallocate(zoo(n)%jprod_don_100)
    enddo

    do n = 2,3
       deallocate(zoo(n)%jhploss_n_100)
       deallocate(zoo(n)%jprod_ndet_100)
    enddo

    deallocate(bact(1)%jprod_n_100)
    deallocate(bact(1)%jzloss_n_100)
    deallocate(bact(1)%jvirloss_n_100)
    deallocate(bact(1)%jremin_n_100)
    deallocate(bact(1)%juptake_ldon_100)
    deallocate(bact(1)%f_n_100)

    if (do_14c) then                                        !<<RADIOCARBON
      deallocate(cobalt%c14_2_n)  
      deallocate(cobalt%f_di14c)  
      deallocate(cobalt%f_do14c)  
      deallocate(cobalt%fpo14c)  
      deallocate(cobalt%j14c_decay_dic)  
      deallocate(cobalt%j14c_decay_doc)  
      deallocate(cobalt%j14c_reminp)  
      deallocate(cobalt%jdi14c)  
      deallocate(cobalt%jdo14c)  
      deallocate(cobalt%c14o2_alpha)  
      deallocate(cobalt%c14o2_csurf)  
      deallocate(cobalt%b_di14c )
    endif                                                   !RADIOCARBON>>
      deallocate(cobalt%runoff_flux_alk)
      deallocate(cobalt%runoff_flux_dic)
      deallocate(cobalt%runoff_flux_di14c)
      deallocate(cobalt%runoff_flux_lith)
      deallocate(cobalt%runoff_flux_fed)
      deallocate(cobalt%runoff_flux_no3)
      deallocate(cobalt%runoff_flux_ldon)
      deallocate(cobalt%runoff_flux_sldon)
      deallocate(cobalt%runoff_flux_srdon)
      deallocate(cobalt%runoff_flux_ndet)
      deallocate(cobalt%runoff_flux_po4)
      deallocate(cobalt%runoff_flux_ldop)
      deallocate(cobalt%runoff_flux_sldop)
      deallocate(cobalt%runoff_flux_srdop)
      deallocate(cobalt%dry_fed)
      deallocate(cobalt%wet_fed)
      deallocate(cobalt%dry_lith)
      deallocate(cobalt%wet_lith)
      deallocate(cobalt%dry_no3)
      deallocate(cobalt%wet_no3)
      deallocate(cobalt%dry_nh4)
      deallocate(cobalt%wet_nh4)
      deallocate(cobalt%dry_po4)
      deallocate(cobalt%wet_po4)
      deallocate(cobalt%stf_gas_dic)
      deallocate(cobalt%stf_gas_o2)
      deallocate(cobalt%deltap_dic)
      deallocate(cobalt%deltap_o2)

  end subroutine user_deallocate_arrays


!f1p
 function calc_pka_nh3(tc,salt) result(pka)
    !temperature, salinity
    real, intent(in) :: tc,salt
    real :: pka,tk
!Bell 2007
!    pka = 10.0423-0.0315536*tc+0.003071*salt

!Clegg 1995
    real, parameter :: a1=0.0500616
    real, parameter :: a2=-9.412696
    real, parameter :: a3=-2.029559e-7
    real, parameter :: a4=-0.0142372
    real, parameter :: a5=1.46041e-5
    real, parameter :: a6=3.730005
    real, parameter :: a7=7.14045e-5
    real, parameter :: a8=-0.0229021
    real, parameter :: a9=-5.521278e-7
    real, parameter :: a10=1.95413e-4
    
    tk=tc+273.15;
    pka     = 9.244605-2729.33*(1/298.15-1./tk)  &
            + (a1+a2/tk+a3*tk**2.)*salt**0.5       &
            + (a4+a5*tk+a6/tk)*salt              &
            + (a7+a8/tk)*salt**2.                 &
            + (a9+a10/tk)*salt**3.;

  end function calc_pka_nh3

!salting out correction for solubility (Johnson 2010, Ocean Science)
  function saltout_correction(kh,vb,salt) result(C)
    real, intent(in) :: Kh,vb,salt
    real*8 :: T,log_kh,theta2    
    real :: theta,C
    log_kh = log(kh)
    theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log_kh) + (-2.4088830102075734e-06*(log_kh)**2) + (1.5711393120941302e-07*(log_kh)**3))*log(vb)
    C = 10**(theta*salt)    
  end function saltout_correction

  !schmidt number in water
  function schmidt_w(t,s,vb,rho) result(sc)
    !schmidt number of the gas in the water                                                    
    real, intent(in) :: t,s,vb
    real, intent(in), optional :: rho
    real             :: sc

    sc=2.*v_sw(t,s,rho)/(d_hm(t,s,vb)+d_wc(t,s,vb))    
  end function schmidt_w

  function v_sw(t,s,rho) result(v)
    real, intent(in) :: t,s
    real, intent(in), optional :: rho
    real             :: n,p,v
    n=n_sw(t,s)*1e-3
    if (present(rho)) then
       p=rho
    else
       p=p_sw(t,s)
    end if
    v = 1e4*n/p
  end function  v_sw

  function p_sw(t,s) result(p)
    !density of sea water 
    !millero and poisson (1981)                                                               
    real, intent(in) :: t,s
    real             :: p, a, b, c
    a = 0.824493-(4.0899e-3*t)+(7.6438e-5*(t**2))-(8.2467e-7*(t**3))+(5.3875e-9*(t**4))
    b = -5.72466e-3+(1.0277e-4*t)-(1.6546e-6*(t**2))
    c = 4.8314e-4
    ! density of pure water
    p = 999.842594+(6.793952e-2*t)-(9.09529e-3*(t**2))+(1.001685e-4*(t**3))-(1.120083e-6*(t**4))+(6.536332e-9*(t**5))
    !salinity correction
    p = (p+(a*s)+(b*(s**(1.5)))+(c*s))
  end function p_sw

  function d_wc(t,s,vb) result(d)
    real, intent(in) :: t,s,vb
    real             :: d
    real, parameter  :: phi = 2.6
    !wilkie and chang 1955
    d = ((t+273.15)*7.4e-8*(phi*18.01)**0.5)/((n_sw(t,s))*(vb**0.6))
  end function d_wc

  function d_hm(t,s,vb) result(d)
    real, intent(in) :: t,s,vb
    real             :: d, epsilonstar
    ! hayduk 1982                                                                              
    epsilonstar = (9.58/vb)-1.12
    d=1.25e-8*(vb**(-0.19)-0.292)*((t+273.15)**(1.52))*((n_sw(t,s))**epsilonstar)
  end function d_hm

  function n_sw(t,s) result(n)
    !dynamic viscosity                                   
    !laliberte 2007                                      
    real :: n
    real, intent(in) :: t,s !temperature (c) and salinity
    !salt in the order nacl,kcl,cacl2,mgcl2,mgso4      
    real, parameter :: mass_fraction(5) = (/ 0.798,0.022,0.033,0.047,0.1 /)
    real, parameter :: v1(5) = (/ 16.22 , 6.4883, 32.028, 24.032, 72.269/)
    real, parameter :: v2(5) = (/ 1.3229 , 1.3175, 0.78792, 2.2694, 2.2238/)
    real, parameter :: v3(5) = (/ 1.4849 , -0.7785, -1.1495,  3.7108, 6.6037/)
    real, parameter :: v4(5) = (/ 0.0074691 , 0.09272, 0.0026995,  0.021853, 0.0079004/)
    real, parameter :: v5(5) = (/ 30.78 , -1.3, 780860., -1.1236, 3340.1/)
    real, parameter :: v6(5) = (/ 2.0583 , 2.0811, 5.8442,0.14474, 6.1304/)

    real :: n_0,ln_n_m,w_i_ln_n_i_tot,ni,w_i_tot,w_i
    integer :: i
    w_i_tot=0
    w_i_ln_n_i_tot=0
    do i=1,5
       w_i = mass_fraction(i)*s/1000
       w_i_tot = w_i+w_i_tot
    enddo
    do i=1,5
       w_i = mass_fraction(i)*s/1000
       ni = (exp(((v1(i)*w_i_tot**v2(i))+v3(i))/((v4(i)*t) + 1)))/((v5(i)*(w_i_tot**v6(i)))+1)
       w_i_ln_n_i_tot = w_i_ln_n_i_tot + (w_i*log(ni))
    enddo
    n_0 = (t+246)/(137.37+(5.2842*t)+(0.05594*(t**2)))
    ln_n_m = (1-w_i_tot)*log(n_0)+w_i_ln_n_i_tot
    n = exp(ln_n_m)
  end function n_sw


 
end module generic_COBALT
