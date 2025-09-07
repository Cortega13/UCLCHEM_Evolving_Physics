MODULE SurfaceReactions
  USE constants
  USE DEFAULTPARAMETERS
  !f2py INTEGER, parameter :: dp
  USE f2py_constants
  USE network
  IMPLICIT NONE
  REAL(dp) :: surfaceCoverage,totalSwap,bulkLayersReciprocal
  REAL(dp) :: safeMantle,safeBulk
  
  REAL(dp), ALLOCATABLE ::vdiff(:)
CONTAINS
  !=======================================================================
  !
  !  Calculate the rate of molecular hydrogen (H2) formation on grains
  !  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
  !  Cazaux & Tielens (2004, ApJ, 604, 222).
  !
  !-----------------------------------------------------------------------
  FUNCTION h2FormEfficiency(gasTemp,dustTemp) RESULT(rate)
    REAL(dp) :: rate
    REAL(dp), INTENT(IN) :: gasTemp,dustTemp

    REAL(dp) :: THERMAL_VELOCITY,STICKING_COEFFICIENT,CROSS_SECTION_SCALE
    REAL(dp) :: HFLUX,FACTOR1,FACTOR2,EPSILON
    REAL(dp) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
    !  Mean thermal velocity of hydrogen atoms (cm s^-1)
    THERMAL_VELOCITY=1.45D5*SQRT(gasTemp/1.0D2)

    !  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
    !  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
    STICKING_COEFFICIENT=1.0D0/(1.0D0+0.04D0*SQRT(gasTemp+dustTemp) &
                    & + 0.2D0*(gasTemp/1.0D2)+0.08D0*(gasTemp/1.0D2)**2)

    HFLUX=1.0D-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

    FACTOR1=SILICATE_MU*HFLUX/(2*SILICATE_NU_H2*EXP(-SILICATE_E_H2/dustTemp))

   FACTOR2=1.0D0*(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0D0*EXP(-SILICATE_E_S/dustTemp)

   EPSILON=1.0D0/(1.0D0+SILICATE_NU_HC/(2*HFLUX)*EXP(-1.5*SILICATE_E_HC/dustTemp) &
              & *(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON


   FACTOR1=GRAPHITE_MU*HFLUX/(2*GRAPHITE_NU_H2*EXP(-GRAPHITE_E_H2/dustTemp))

   FACTOR2=1.0D0*(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0D0*EXP(-GRAPHITE_E_S/dustTemp)

   EPSILON=1.0D0/(1.0D0+GRAPHITE_NU_HC/(2*HFLUX)*EXP(-1.5*GRAPHITE_E_HC/dustTemp) &
              & *(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   rate=0.5D0*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT

     RETURN
  END FUNCTION h2FormEfficiency

  !surface abundance multiplied by this value gives fraction of surface covered by material
  FUNCTION bulkGainFromMantleBuildUp() RESULT(rate)
    REAL(dp) :: rate
    rate=0.5*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN
  END FUNCTION bulkGainFromMantleBuildUp

  !----------------------------------------------------------------------------------------------------
!Reactions on the surface treated by evaluating diffusion rates across grains and accounting
!For competition with chemical desorption. Products remain bound ('DIFF') or are desorbed ('CHEMDES')
!Units of s-1. 
!David Quenard 2017 Arxiv:1711.05184
!----------------------------------------------------------------------------------------------------
double precision FUNCTION diffusionReactionRate(reacIndx,dustTemperature)
    double precision :: reducedMass,tunnelProb,dustTemperature
    double precision :: diffuseProb,desorbProb,reacProb,n_dust
    integer(dp) :: index1,index2,reacIndx,i

    !want position of species in the grain array but gas phase species aren't in there
    !so store species index
    index1=re1(reacIndx)
    index2=re2(reacIndx)

    !then try to overwrite with position in grain array
    DO i=lbound(iceList,1),ubound(iceList,1)
        IF (iceList(i) .eq. index1) index1 = i
        IF (iceList(i) .eq. index2) index2 = i
    END DO

    !Hasegawa 1992 diffusion rate. Rate that two species diffuse and meet on grain surface
    diffuseProb = vdiff(index1)*dexp(-DIFFUSION_BIND_RATIO*bindingEnergy(index1)/dustTemperature)
    diffuseProb = diffuseProb+ (vdiff(index2)*dexp(-DIFFUSION_BIND_RATIO*bindingEnergy(index2)/dustTemperature))

    !probability a reactant will just desorb
    desorbProb = vdiff(index1)*dexp(-bindingEnergy(index1)/dustTemperature)
    desorbProb = desorbProb + vdiff(index2)*dexp(-bindingEnergy(index2)/dustTemperature) 

    !Calculate classical activation energy barrier exponent
    reacProb = gama(reacIndx)/dustTemperature
    !Calculate quantum activation energy barrier exponent
    reducedMass = reducedMasses(reacIndx)
    IF (reducedMass .eq. 0.0) THEN 
        ! Should never happen, just as a backup
        ! If no reducedMass was supplied in the array, calculate it from the two reacting species
        reducedMass = mass(icelist(index1)) * mass(icelist(index2)) / (mass(icelist(index1)) + mass(icelist(index2)))
    END IF
    tunnelProb = 2.0d0 *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * dsqrt(2.0d0*AMU*reducedMass*K_BOLTZ*gama(reacIndx))

    !Choose fastest between classical and tunnelling
    IF (reacProb.GT.tunnelProb) reacProb=tunnelProb

    !Overall reaction probability is chance of reaction occuring on meeting * diffusion rate
    reacProb = max(vdiff(index1),vdiff(index2)) * dexp(-reacProb)       


    ! Keff from Garrod & Pauly 2011 and Ruaud+2016
    ! Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    IF(DIFFUSE_REACT_COMPETITION) THEN
       reacProb = reacProb/(reacProb + desorbProb + diffuseProb)
    END IF
    
    !see Eq A1 of Quenard et al. 2018
    !NUM_SITES_PER_GRAIN should be multiplied by n_dust as in A1
    !n_dust=density/GAS_DUST_DENSITY_RATIO so we use the 1/density to cancel the density in odes.f90 and drop it here
    diffusionReactionRate=alpha(reacIndx) *reacProb* diffuseProb*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN

END FUNCTION diffusionReactionRate

! ---------------------------------------------------------------------
!  Chemical Reactive Desorption (CRD)
! David Quenard 2017 Arxiv:1711.05184
! From Minissalle+ 2016 and Vasyunin+ 2016
! ---------------------------------------------------------------------
double precision FUNCTION desorptionFraction(reacIndx)
    integer(dp) :: reacIndx,reactIndex1,reactIndex2,degreesOfFreedom,i
    integer(dp) :: productIndex(4)

    double precision :: deltaEnthalpy,maxBindingEnergy,epsilonCd,productEnthalpy
    double precision, parameter :: EFFECTIVE_SURFACE_MASS = 120.0

    
    !Get indices of grain surface version of products products 
    productIndex=0
    !Arrays like binding energy and formation enthalpy are indexed by position in iceList
    !rather than species list. Need to find position in grain list where every reactant and product appears
    !bearing in mind that Eley-Rideal reactions can have reactants in gas phase and CHEMDES has products in gas
    DO i=lbound(iceList,1),ubound(iceList,1)
        !check grain lists for reactants
        IF (iceList(i) .eq. re1(reacIndx)) reactIndex1 = i
        IF (gasiceList(i) .eq. re1(reacIndx)) reactIndex1 = i
        !check equivalent gas list in case of ER reaction.
        IF (iceList(i) .eq. re2(reacIndx)) reactIndex2 = i
        IF (gasiceList(i) .eq. re2(reacIndx)) reactIndex2 = i

        IF (iceList(i) .eq. p1(reacIndx)) productIndex(1) = i
        IF (iceList(i) .eq. p2(reacIndx)) productIndex(2) = i
        IF (iceList(i) .eq. p3(reacIndx)) productIndex(3) = i
        IF (iceList(i) .eq. p4(reacIndx)) productIndex(4) = i

        IF (gasiceList(i) .eq. p1(reacIndx)) productIndex(1) = i
        IF (gasiceList(i) .eq. p2(reacIndx)) productIndex(2) = i
        IF (gasiceList(i) .eq. p3(reacIndx)) productIndex(3) = i
        IF (gasiceList(i) .eq. p4(reacIndx)) productIndex(4) = i
    END DO

    maxBindingEnergy=0.0
    productEnthalpy=0.0
    epsilonCd=0.0
    DO i=1,4
        IF (productIndex(i) .ne. 0) THEN
            maxBindingEnergy=MAX(maxBindingEnergy,bindingEnergy(productIndex(i)))
            productEnthalpy=productEnthalpy+formationEnthalpy(productIndex(i))
            epsilonCd=epsilonCd + mass(productIndex(i))
        END IF
    END DO

    !epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))**2
    
    !Now calculate the change in enthalpy of the reaction.
    deltaEnthalpy= formationEnthalpy(reactIndex1)+formationEnthalpy(reactIndex2)-productEnthalpy
    
    !Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*4.184d03/(1.38054D-23*6.02214129d23)
    ! Total energy change includes activation energy of the reaction !
    deltaEnthalpy = deltaEnthalpy + gama(reacIndx)


    IF (deltaEnthalpy.eq.0.00) deltaEnthalpy = 1e-30 

    !Degrees of freedom = 3 * number of atoms in the molecule
    degreesOfFreedom = atomCounts(productIndex(1))
    if (productIndex(2).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(2)))
    if (productIndex(3).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(3)))
    if (productIndex(4).NE.0) degreesOfFreedom = max(degreesOfFreedom,atomCounts(productIndex(4)))                    
    degreesOfFreedom = 3 * degreesOfFreedom
        
    desorptionFraction = dexp((-maxBindingEnergy*real(degreesOfFreedom)) / (epsilonCd * deltaEnthalpy))
    
   IF (deltaEnthalpy.lt.0.d0) THEN        !< If reaction is endothermic, no CRD
        desorptionFraction = 0.d0
    END IF
    
    IF (GRAINS_HAVE_ICE) THEN
        desorptionFraction = desorptionFraction/10    !< See Minisalle et al. 2016 for icy grain surface.
        ! Special case of OH+H, O+H, N+N on ices, see same paper
        if (re1(reacIndx).eq.ngn.and.re2(reacIndx).eq.ngn) desorptionFraction = 0.5
        if ((re1(reacIndx).eq.ngo.and.re2(reacIndx).eq.nh) &
            &.or. (re1(reacIndx).eq. nh.and.re2(reacIndx).eq.ngo)) desorptionFraction = 0.3
        if ((re1(reacIndx).eq.ngoh.and.re2(reacIndx).eq.nh) &
            &.or. (re1(reacIndx).eq.nh.and.re2(reacIndx).eq.ngoh)) desorptionFraction = 0.25
    ENDIF
END FUNCTION desorptionFraction
END MODULE SurfaceReactions
