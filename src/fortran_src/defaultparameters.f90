!!This file is used to autogenerate the docs. So please ignore the mess!
!!Double !! lines do not show up in docs, single ones do. 
!!If you add a parameter, please take the time to add a useful descriptor comment on the same line
!!and then re-run utils/generate_param_docs.py to update the docs.
!!note the resuting md file needs manually adding to the website.
MODULE DEFAULTPARAMETERS
USE constants
!---  
!id: parameters
!title: Model Parameters
!---
!UCLCHEM will default to these values unless they are overridden by user. Users can override these by adding the variable name as written here in the param_dict argument of any UCLCHEM model function. param_dict is not case sensitive.
!
!## Physical Variables
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
IMPLICIT NONE
LOGICAL :: useAVDirectly=.True. ! Determins whether we do av=baseAv or av = baseAv + coldens/1.6d21
LOGICAL :: evolving_physical_params = .False.

REAL(dp) :: initialTemp=10.0 !Initial gas temperature in Kelvin for all gas parcels in model.
REAL(dp) :: initialDens=1.00d2 !Initial gas density in H nuclei per cm$^{-3}$ for all gas parcels in model.
REAL(dp) :: initialBaseAv=2.0 !Extinction at cloud edge, Av of a parcel at rout.
REAL(dp) :: initialRadfield=1.0 !Interstellar radiation field in Habing

REAL(dp) :: tempRate=0.0 !The rate at which the temperature changes.
REAL(dp) :: densRate=0.0 !The rate at which the density changes.
REAL(dp) :: baseAvRate=0.0 !The rate at which the extinction changes.
REAL(dp) :: radfieldRate=0.0 !The rate at which the radiation field changes.

REAL(dp) :: radfield=1.0 !Interstellar radiation field in Habing
REAL(dp) :: baseAv=2.0 !Extinction at cloud edge, Av of a parcel at rout.
REAL(dp) :: finalDens=1.00d5 !Final gas density achieved via freefall. 
REAL(dp) :: currentTime=0.0 !Time at start of model in years.
REAL(dp) :: finalTime=5.0d6 !Time to stop model in years, if not using `endAtFinalDensity` below.
REAL(dp) :: zeta=1.0 !Cosmic ray ionisation rate as multiple of $1.3 10^{-17} s^{-1}$
REAL(dp) :: rout=0.05 !Outer radius of cloud being modelled in pc.
REAL(dp) :: rin=0.0 !Minimum radial distance from cloud centre to consider.
INTEGER(dp) :: points=1 !Number of gas parcels equally spaced between rin to rout to consider
REAL(dp) :: bm0=1.0 !magnetic parameter [microgauss]: B0 = bm0*sqrt(initialDens)

!## Behavioural Controls
!*The following parameters generally turn on or off features of the model. If a parameter is set to `True`, then it is turned on. If it is set to `False`, then it is turned off.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: freezeFactor=1.0 !Modify freeze out rate of gas parcels by this factor.
LOGICAL :: endAtFinalDensity=.False. !Choose to end model at final density, otherwise end at final time.
LOGICAL :: freefall=.False. !Controls whether models density increaes following freefall equation.
REAL(dp) :: freefallFactor=1.0 !Modify freefall rate by factor, usually to slow it.
LOGICAL :: desorb=.True. !Toggles all non-thermal desoprtion processes on or off.
LOGICAL :: h2desorb=.True. !Individually toggle non-thermal desorption due to H2 formation.
LOGICAL :: crdesorb=.True. !Individually toggle non-thermal desorption due to cosmic rays.
LOGICAL :: uvdesorb=.True. !Individually toggle non-thermal desorption due to uv photons.
LOGICAL :: thermdesorb=.True. !Toggle continuous thermal desorption.
LOGICAL :: instantSublimation=.False. !Toggle instantaneous sublimation of the ices at t=0
LOGICAL :: cosmicRayAttenuation=.False. !Use column density to attenuate cosmic ray ionisation rate following [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
CHARACTER :: ionModel='L' !L/H model for cosmic ray attenuation [Padovani et al. 2018](https://arxiv.org/abs/1803.09348).
LOGICAL :: improvedH2CRPDissociation=.False. !Use H2 CRP dissociation rate from [Padovani et al. 2018b](https://arxiv.org/abs/1809.04168).
!
!## Input and Output
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
CHARACTER(256) :: outputFile="" !File to write full output of UCLCHEM. This includes physical parameter values and all abundances at every time step.
CHARACTER(256) :: columnFile="" !File to write specific species abundances, see outSpecies.
CHARACTER(256) :: rateFile="" !File to write rate 'constants' at each timestep. This includes physical parameter values.
CHARACTER(256) :: fluxFile="" !File to write reaction rates (flux) at each timestep. This includes physical parameter values.
INTEGER :: writeStep=1 !Writing to columnFile only happens every writeStep timesteps.
CHARACTER(256) :: abundSaveFile="" ! The file to save the abundances to at the end of the model.
CHARACTER(256) :: abundLoadFile="" ! The file to load the abundances from at the start of the model.
!|abundSaveFile |None| File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced.
!|abundLoadFile |None| File from which to load initial abundances for the model, created through `abundSaveFile`. If not provided, the model starts from elemental gas.
!|outSpecies|None| A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions, see python API docs.
!
!## Initial Abundances
!*Unless otherwise specified, we take all abundances from Jenkins et al. 2009, using the heavily depleted case from Table 4.*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: metallicity=1.0 !Scale the abundances of all elements heavier than He by this factor.
INTEGER(dp) :: ion=2 !Sets how much elemental C is initially atomic (0= all atomic/1=50:50/2=fully ionized).
REAL(dp) :: fh=0.5 !Total elemental abundance of H is always 1 by definition because abundances are relative to number of H nuclei. Use fh to set how much to initially put in atomic H, the rest goes to H2.
REAL(dp) :: fhe = 0.1 !Total elemental abundance of He.
REAL(dp) :: fc=1.77d-04 !Total elemental abundance of C.
REAL(dp) :: fo  = 3.34d-04 !Total elemental abundance of O.
REAL(dp) :: fn  = 6.18d-05 !Total elemental abundance of N.
REAL(dp) :: fs  = 3.51d-6 !Total elemental abundance of S.
REAL(dp) :: fmg = 2.256d-06 !Total elemental abundance of Mg.
REAL(dp) :: fsi = 1.78d-06 !Total elemental abundance of Si.
REAL(dp) :: fcl = 3.39d-08 !Total elemental abundance of Cl.
REAL(dp) :: fp =7.78d-08 !Total elemental abundance of P.
REAL(dp) :: ffe =2.01d-7!Total elemental abundance of Fe.
REAL(dp) :: ff = 3.6d-08 !fp depleted 1/100 of solar from Asplund 2009.
REAL(dp) :: fd=0.0 ! The following elements are not typically used. We do not recommend any particular value.
REAL(dp) :: fli=0.0 !Total elemental abundance of Li.
REAL(dp) :: fna=0.0 !Total elemental abundance of Na.
REAL(dp) :: fpah=0.0 !Total initial abundance of PAHs.
REAL(dp) :: f15n=0.0 !Total initial abundance of 15N.
REAL(dp) :: f13c=0.0 !Total initial abundance of 13C.
REAL(dp) :: f18O=0.0 !Total initial abundance of 18O.
!!
!! We used to use Asplund et al. 2009,kept here for reference
!! !initial fractional abundances of elements(from Asplund et al. 2009 ARAA table 1 -SOLAR)
!! !note fh is fraction of H initially in H atoms. Total H is always 1.
!! !fh=0.5;fhe = 0.1;fc  = 2.6d-04;fo  = 4.6d-04;fn  = 6.1d-05
!! fs  = 1.318d-05;fmg = 3.981d-05;fsi = 1.0d-07;fcl = 3.162d-07;
!! fp=2.57d-09 ; ff = 3.6d-08 !fp depleted 1/100 of solar
!
!## Integration Controls
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: reltol=1d-8 !Relative tolerance for integration, see [integration docs](/docs/trouble-integration) for advice.
REAL(dp) :: abstol_factor=1.0d-14 !Absolute tolerance for integration is calculated by multiplying species abundance by this factor.
REAL(dp) :: abstol_min=1.0d-25 !Minimum value absolute tolerances can take.
INTEGER :: MXSTEP=10000 !Maximum steps allowed in integration before warning is thrown. ! HAS TO BE INT4 instead of INT8
!
!## Here be Dragons
!*These are not recommended to be changed unless you know what you are doing*
!
!|Parameter|Default Value |Description|
!| ----- | ------| ------ |
REAL(dp) :: ebmaxh2=1.21d3 ! Maximum binding energy of species desorbed by H2 formation.
REAL(dp) :: ebmaxcr=1.21d3 ! Maximum binding energy of species desorbed by cosmic ray ionisation.
REAL(dp) :: ebmaxuvcr=1.0d4 ! Maximum binding energy of species desorbed by UV photons.
REAL(dp) :: epsilon=0.01 !Number of molecules desorbed per H2 formation.
REAL(dp) :: uv_yield=0.03 !Number of molecules desorbed per UV photon. The yield is extrapolated from Oberg et al. 2009
REAL(dp) :: phi=1.0d5 !Number of molecules desorbed per cosmic ray ionisation.
REAL(dp) :: uvcreff=1.0d-3 !Ratio of CR induced UV photons to ISRF UV photons.
REAL(dp) :: omega=0.5 !Dust grain albedo.
!|alpha|{1:0.0,2:0.0}| Set alpha coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how alpha is used for each reaction type.|
!|beta|{1:0.0,2:0.0}| Set beta coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how beta is used for each reaction type.|
!|gama|{1:0.0,2:0.0}| Set gama coeffecients of reactions using a python dictionary where keys are reaction numbers and values are the coefficients. Once you do this, you cannot return to the default value in the same python script or without restarting the kernel in iPython. See the chemistry docs for how gama is used for each reaction type.|
!
!
!
! ## Moved the rates.f90 constants to be here. 
!Variables controlling chemistry:
LOGICAL :: PARAMETERIZE_H2FORM=.True.
REAL(dp) :: grainArea,cion,h2dis,lastTemp=0.0

! Controlling ice chemistry
REAL(dp), PARAMETER :: h2StickingZero=0.87d0,hStickingZero=1.0d0, h2StickingTemp=87.0d0,hStickingTemp=52.0d0
!Flags to control desorption processes
REAL(dp) :: turbVel=1.0
!
!
!
!
! ## Moved the photoreactions.f90 constants to be here.
REAL(dp) :: UV_FAC=3.02

!Below are arrays for self-shielding of CO and H2
LOGICAL :: start=.True.
INTEGER :: NUM_LAMBDA=30
REAL(dp), DIMENSION(30) :: LAMBDA_GRID=(/ &
    &  910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
    & 1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
    & 1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
    & 2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
    & 3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
    & 9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)
REAL(dp), DIMENSION(30) :: XLAMBDA_GRID=(/ &
    & 5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
    & 3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
    & 2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
    & 3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
    & 1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
    & 0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)
REAL(dp), DIMENSION(30) :: XLAMBDA_DERIV

logical :: startr=.True.

!  12CO line shielding data from van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
INTEGER, PARAMETER ::  DIMCO=7, DIMH2=6
REAL(KIND=DP), DIMENSION(8) :: NCO_GRID=(/12.0D0,13.0D0,14.0D0,15.0D0,16.0D0,17.0D0,18.0D0,19.0D0/)
REAL(KIND=DP), DIMENSION(6) :: NH2_GRID=(/18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
REAL(KIND=DP), DIMENSION(8,6) :: SCO_GRID=RESHAPE((/ &
    &  0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00, &
    & -8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00, &
    & -1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00, &
    & -4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00, &
    & -1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00, &
    & -3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/), (/8,6/))
REAL(KIND=DP), DIMENSION(8,6) :: SCO_DERIV

REAL(dp) :: ICE_GAS_PHOTO_CROSSSECTION_RATIO = 0.3 ! Kalvans 2018
!
!
!
!
!
! ## Moved the surfacereactions.f90 constants to be here.
!Silicate grain properties for H2 Formation
REAL(dp),PARAMETER :: SILICATE_MU=0.005D0 ! Fraction of newly formed H2 that stays on the grain surface
REAL(dp),PARAMETER :: SILICATE_E_S=110.0D0 ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
REAL(dp),PARAMETER :: SILICATE_E_H2=320.0D0 ! Desorption energy of H2 molecules (K)
REAL(dp),PARAMETER :: SILICATE_E_HP=450.0D0 ! Desorption energy of physisorbed H atoms (K)
REAL(dp),PARAMETER :: SILICATE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
REAL(dp),PARAMETER :: SILICATE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
REAL(dp),PARAMETER :: SILICATE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)
REAL(dp),PARAMETER :: SILICATE_CROSS_SECTION=8.473D-22!*CROSS_SECTION_SCALE ! Silicate grain cross section per H nucleus (cm^-2/nucleus)

!Graphite grain properties for H2 Formation
REAL(dp),PARAMETER :: GRAPHITE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
REAL(dp),PARAMETER :: GRAPHITE_E_S=260.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
REAL(dp),PARAMETER :: GRAPHITE_E_H2=520.0D0 ! Desorption energy of H2 molecules (K)
REAL(dp),PARAMETER :: GRAPHITE_E_HP=800.0D0 ! Desorption energy of physisorbed H atoms (K)
REAL(dp),PARAMETER :: GRAPHITE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
REAL(dp),PARAMETER :: GRAPHITE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
REAL(dp),PARAMETER :: GRAPHITE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)
REAL(dp),PARAMETER :: GRAPHITE_CROSS_SECTION=7.908D-22!*CROSS_SECTION_SCALE ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Grain surface parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp), PARAMETER :: GAS_DUST_MASS_RATIO=100.0,GRAIN_RADIUS=1.d-5, GRAIN_DENSITY = 3.0 ! Mass density of a dust grain
REAL(dp), PARAMETER :: THERMAL_VEL= SQRT(8.0d0*K_BOLTZ/(PI*AMU)) !Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

!reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
REAL(dp), PARAMETER :: GAS_DUST_DENSITY_RATIO = (4.0*PI*(GRAIN_RADIUS**3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0 * AMU)
!Grain area per h nuclei, values taken from Cazaux & Tielens 2004 via UCL-PDR to match H2 formation rate
REAL(dp), PARAMETER :: GRAIN_CROSSSECTION_PER_H=0.5*(7.908D-22+8.473D-22)
REAL(dp), PARAMETER :: GRAIN_SURFACEAREA_PER_H=4.0*GRAIN_CROSSSECTION_PER_H!2.0*4.0*PI*GRAIN_RADIUS*GRAIN_RADIUS/GAS_DUST_DENSITY_RATIO

!Below are values for grain surface reactions
LOGICAL, PARAMETER :: DIFFUSE_REACT_COMPETITION=.True., GRAINS_HAVE_ICE=.True.
REAL(dp), PARAMETER :: DIFFUSION_BIND_RATIO=0.5 ! Ratio between diffusion barrier and binding energy of a species
REAL(dp), PARAMETER :: CHEMICAL_BARRIER_THICKNESS = 1.40d-8! Parameter used to compute the probability for a surface reaction with 
!! activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
REAL(dp), PARAMETER :: SURFACE_SITE_DENSITY = 1.5d15 ! site density on one grain [cm-2]
REAL(dp), PARAMETER :: VDIFF_PREFACTOR=2.0*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
REAL(dp), PARAMETER :: NUM_SITES_PER_GRAIN = GRAIN_RADIUS*GRAIN_RADIUS*SURFACE_SITE_DENSITY*4.0*PI

REAL(dp), PARAMETER :: MAX_GRAIN_TEMP=150.0, MIN_SURFACE_ABUND=1.0d-20

!
!
CONTAINS
! Add a dummy subroutine to help f2py compile: https://github.com/numpy/numpy/issues/27167
SUBROUTINE DUMMY_TWO(dummy_two_output)
        integer, intent(out) :: dummy_two_output
        dummy_two_output = 2
    END SUBROUTINE DUMMY_TWO
END MODULE DEFAULTPARAMETERS