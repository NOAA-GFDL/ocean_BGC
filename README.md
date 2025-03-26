# NOAA-GFDL ocean_BGC

This repository contains several generic tracer packages that can interface with MOM5 and MOM6. The generic tracers reley on the FMS coupler types and are generally used with the FMS coupler. 
Details of the different generic tracer packages are given below.

## Interfacing with MOM5 and MOM6

Prior to 2025, `ocean_BGC` was compatible with MOM5 and MOM6 by default. However, a number of changes to the MOM6-ocean_BG interface needed for the updated version of COBALT motivated moving the interface. 
This means that the file `MOM6_generic_tracers.F90` which perviously was included in MOM6 has been copied here. 

To use `ocean_BGC` with MOM5, use version 2024 or earlier.

To use `ocean_BGC` with older versions of MOM6, before [MOM6 commit fcf5fff](https://github.com/NOAA-GFDL/MOM6/pull/790), use version 2024 or earlier.

To use `ocean_BGC` with newer versions of MOM6, **after** [MOM6 commit fcf5fff](https://github.com/NOAA-GFDL/MOM6/pull/790), use version 2025.

## Generic Tracer packages

### generic_COBALT.F90
COBALT simulates the biogeochemical cycling of carbon, nitrogen,
phosphorous, iron, silica, calcium carbonate, and lithogenic
material in the ocean.  The code is built upon the TOPAZ code
developed by John Dunne.  The primary changes to TOPAZ are:
1) the addition of three zooplankton groups
2) The addition of bacteria
3) The expansion of the dissolved organic nitrogen and 
phosphorous groups to include three types each: labile,
semi-labile, and refractory
4) Constant Stoichiometry by plankton functional type 

The primary COBALT reference is:

Stock, CA, Dunne, JP, John, JG. 2014. Global-scale carbon and energy
flows through the marine planktonic food web: An analysis with a'
coupled physical-biological model.  Progress in Oceanography 120, 1-18.  

Version 2.0 has a number of refinements:

1) Ammonia uptake parameters are now based on the "high-affinity" 
settings from Paulot et al., 2015; GBC; 29(8)
2) Phytoplankton aggregation is initiated only when growth rates
fall below 1/4 maximum values
3) The default parameterization has elevated N:P ratios for both
diazotrophs and small phytoplankton
4) Remineralization of sinking detritus is now based on the temperature
and oxygen dependences described in Laufkotter et al., 2017; O2
dependence of other aerobic processes have also been adjusted
for consistency.
5) The default carbon chemistry routine is now MOCSY
6) The iron scavenging has been re-tuned to new atmospheric (Ginoux-AM4),
sediment, river and hydrothermal vent (Tagliabue) sources.
       
### generic_BLING.F90
Biogeochemistry with Light, Iron, Nutrient and Gas (BLING) includes an implicit ecological model of growth limitation by light, temperature, phosphate and iron, along with dissolved organic
phosphorus and O2 pools. Food web processing in the euphotic zone and remineralization/dissolution through the ocean interior are handled as in Dunne et al. (2005).  
O2 equilibria and gas exchange follow OCMIP2 protocols. Additional functionality comes from an optional carbon cycle that is  non-interactive, i.e. does not change the core BLING behaviour, as
well as tracers for radiocarbon (14c), a decomposition of carbon components by gas exchange and remineralization (carbon_pre), and a decomposition of phosphate as preformed and remineralized (po4_pre).

### generic_miniBLING.F90

### generic_TOPAZ.F90
Phytoplankton Biogeochemistry: Includes an explicit ecological model
including three phytoplankton groups (small, large/diatoms and
diazotrophs), growth limitation by light, temperature and a suite of
nutrients including nitrate, ammonia, phosphate, iron and silicate,
dissolved inorganic carbon, alkalinity, two kinds of dissolved organic
material, O2, nitrogen fixation and denitrification. CO2 gas exchange
is function of the biologically and physically forced solubility.
Additionally, changes in the vertical distribution of phytoplankton
affect heat absorption with climate feedbacks. Food web processing in
the euphotic zone and remineralization/dissolution through the ocean
interior are handled as in Dunne et al. (in prep).  CO2 and O2
equilibria and gas exchange follow OCMIP2 protocols.

### generic_CFC.F90

### generic_ERGOM.F90

### generic_SF6.F90

### generic_abiotic.F90

### generic_age.F90

### generic_argon.F90

### generic_blres.F90


