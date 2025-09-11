MODULE RATES
    USE constants
    USE DEFAULTPARAMETERS
    !f2py INTEGER, parameter :: dp
    USE network
    USE physicscore
    USE SurfaceReactions
    use photoreactions, only: H2PhotoDissRate, COPhotoDissRate, cIonizationRate, ICE_GAS_PHOTO_CROSSSECTION_RATIO
    IMPLICIT NONE
CONTAINS
    SUBROUTINE calculateConditionalRates(Y, safemantle, rate, Td, Tg, D)
        REAL(dp), INTENT(IN) :: Y(:), safemantle, Td, Tg, D
        REAL(dp), INTENT(INOUT) :: rate(:)
        INTEGER:: idx1,idx2
        INTEGER(dp) :: i,j

        ! lhdes

        idx1=descrReacs(1)
        idx2=descrReacs(2)
        IF (safeMantle .gt. MIN_SURFACE_ABUND) THEN
            rate(idx1:idx2) = 4.0*pi*zeta*1.64d-4*(GRAIN_SURFACEAREA_PER_H)*phi
            WHERE(gama(idx1:idx2) .gt. ebmaxcr) rate(idx1:idx2)=0.0 
        ELSE
            rate(idx1:idx2) = 0.0
        ENDIF
        WHERE((rate(freezePartners)*Y(re1(freezePartners))*D)&
        <MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0

        idx1 = bulkswapReacs(1)
        idx2 = bulkswapReacs(2)
        IF ((Td .gt. MAX_GRAIN_TEMP) .or. (safeMantle .lt. MIN_SURFACE_ABUND)) THEN
            rate(idx1:idx2) = 0.0
        ELSE
            DO i=idx1,idx2
                DO j=lbound(iceList,1),ubound(iceList,1)
                    IF (iceList(j) .eq. re1(i)) THEN
                    rate(i)= 1 
                    END IF
                END DO
            END DO
        END IF

        idx1=deuvcrReacs(1)
        idx2=deuvcrReacs(2)
        IF (safeMantle .gt. MIN_SURFACE_ABUND) THEN
                rate(idx1:idx2) = 1
            WHERE(gama(idx1:idx2) .gt. ebmaxuvcr) rate(idx1:idx2)=0.0 
        ELSE
            rate(idx1:idx2) = 0.0
        ENDIF
        WHERE((rate(freezePartners)*Y(re1(freezePartners))*D)&
        &<MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0

        idx1=thermReacs(1)
        idx2=thermReacs(2)
        DO j=idx1,idx2
            DO i=lbound(iceList,1),ubound(iceList,1)
                IF (iceList(i) .eq. re1(j)) THEN
                    rate(j)= 1
                END IF
            END DO
        END DO
        WHERE(rate(freezePartners)*Y(re1(freezePartners))*D&
        &<MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        IF (safeMantle .lt. MIN_SURFACE_ABUND) rate(idx1:idx2)=0.0

        idx1=desoh2Reacs(1)
        idx2=desoh2Reacs(2)
        IF (safeMantle .gt. MIN_SURFACE_ABUND) THEN
            rate(idx1:idx2) = 1
            WHERE(gama(idx1:idx2) .gt. ebmaxh2) rate(idx1:idx2)=0.0 
        ELSE
            rate(idx1:idx2) = 0.0
        ENDIF
        WHERE((rate(freezePartners)*Y(re1(freezePartners)))<&
        &MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0

        idx1 = surfSwapReacs(1)
        idx2 = surfSwapReacs(2)
        IF ((Td .gt. MAX_GRAIN_TEMP) .or. (safeMantle .lt. MIN_SURFACE_ABUND)) THEN
            rate(idx1:idx2) = 0.0
        ELSE
            rate(idx1:idx2) = 1.0
        END IF

    END SUBROUTINE calculateConditionalrates

    SUBROUTINE calculateReactionRates(abund, safemantle, rate)
        REAL(dp), INTENT(IN) :: abund(:, :), safemantle
        REAL(dp), INTENT(INOUT) :: rate(:)
        INTEGER:: idx1,idx2
        REAL(dp) :: vA,vB
        INTEGER(dp) :: i,j
        ! REAL(dp) :: vdiff(:)

        ! Up here we will include all the reactions which are included in the odes.f90 already by carlos. Note: We have to manually set rate=1 to initalize them.
        !Initialize twobody reactions. 
        idx1=twobodyReacs(1)
        idx2=twobodyReacs(2)
        rate(idx1:idx2) = 1 ! alpha(idx1:idx2)*((gasTemp(dstep)/300.)**beta(idx1:idx2))*dexp(-gama(idx1:idx2)/gasTemp(dstep)) 

        !UV photons, radfield has (factor of 1.7 conversion from habing to Draine)
        idx1=photonReacs(1)
        idx2=photonReacs(2)
        IF (idx1 .ne. idx2) THEN
            rate(idx1:idx2) = 1 !alpha(idx1:idx2)*dexp(-gama(idx1:idx2)*av(dstep))*radfield/1.7
            ! For all solid species, decrease rate by 0.3 (Kalvans 2018)
            ! For bulk species, also decrease rate by (1-Pabs)**(Bs+0.5*Bb) (Kalvans 2014)
            DO j=idx1,idx2
                IF (ANY(bulkList==re1(j))) THEN
                    rate(j) = 1 !rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO * (1.0-0.007)**(1.0+0.5/bulkLayersReciprocal)
                ELSE IF (ANY(surfaceList==re1(j))) THEN
                    rate(j) = 1 !rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO
                END IF
            END DO 
        END IF

        !Reactions involving cosmic ray induced photon
        idx1=crphotReacs(1)
        idx2=crphotReacs(2)
        IF (idx1 .ne. idx2) THEN
            rate(idx1:idx2)= 1 !alpha(idx1:idx2)*gama(idx1:idx2)*1.0/(1.0-omega)*zeta*(gasTemp(dstep)/300)**beta(idx1:idx2)
            ! For all solid species, decrease rate by 0.3 (Kalvans 2018)
            ! For bulk species, also decrease rate by (1-Pabs)**(Bs+0.5*Bb) (Kalvans 2014)
            DO j=idx1,idx2
                IF (ANY(bulkList==re1(j))) THEN
                    rate(j) = 1 !rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO * (1-0.007)**(1+0.5/bulkLayersReciprocal)
                ELSE IF (ANY(surfaceList==re1(j))) THEN
                    rate(j) = 1 !rate(j) * ICE_GAS_PHOTO_CROSSSECTION_RATIO
                END IF
            END DO 
        END IF

        !freeze out only happens if freezeFactor>0 and depending on evap choice
        idx1=freezeReacs(1)
        idx2=freezeReacs(2)
        IF (idx1 .ne. idx2) THEN
            rate(idx1:idx2)= 1 !freezeOutRate(idx1,idx2)
            !freeze out rate uses thermal velocity but mass of E is 0 giving us infinite rates
            !just assume it's same as H
            rate(nR_EFreeze)= 1 !rate(nR_HFreeze)
            rate(nR_H2Freeze)= 1 !stickingCoefficient(h2StickingZero,h2StickingTemp,gasTemp(dstep))*rate(nR_H2Freeze)
            rate(nR_HFreeze)= 1 ! stickingCoefficient(hStickingZero,hStickingTemp,gasTemp(dstep))*rate(nR_HFreeze)
        END IF

        !Desorption due to energy released by H2 Formations. !The below desorption mechanisms are from Roberts et al. 2007 MNRAS with the addition of direct UV photodesorption. DESOH2,DESCR1,DEUVCR
        idx1=desoh2Reacs(1)
        idx2=desoh2Reacs(2)
        IF (idx1 .ne. idx2) THEN
            IF ((desorb) .and. (h2desorb) .and. (safeMantle .gt. MIN_SURFACE_ABUND)) THEN
                !Epsilon is efficieny of this process, number of molecules removed per event
                !h2form is formation rate of h2, dependent on hydrogen abundance. 
                rate(idx1:idx2) = 1 !epsilon*h2FormEfficiency(gasTemp(dstep),dustTemp(dstep))

                !Don't remove species with binding energy > max BE removed by this process
                WHERE(gama(idx1:idx2) .gt. ebmaxh2) rate(idx1:idx2)=0.0 
            ELSE
                rate(idx1:idx2) = 0.0
            ENDIF
            !turn off freeze out if desorption due to H2 formation is much faster
            !both rates combine with density to get rate of change so drop that factor
            WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep))<&
            &MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        END IF

        ! Parameterize H2 Formation
        IF (PARAMETERIZE_H2FORM) THEN
            rate(nR_H2Form_CT) = 1 !h2FormEfficiency(dustTemp(dstep),dustTemp(dstep)) !!
            !rate(nR_H2Form_LH)=0.0
            rate(nR_H2Form_ER)=0.0
            !rate(nR_H2Form_LHDes)=0.0
            rate(nR_H2Form_ERDes)=0.0
        ELSE
            rate(nR_H2Form_CT)= 0.0
        END IF

        !Continuous Thermal Desorption. Reactions can be generated through a flag in Makerates
        idx1=thermReacs(1)
        idx2=thermReacs(2)
        IF (idx1 .ne. idx2) THEN
            IF (thermdesorb) THEN
                DO j=idx1,idx2
                    !then try to overwrite with position in grain array
                    DO i=lbound(iceList,1),ubound(iceList,1)
                        !See Cuppen, Walsh et al. 2017 review (section 4.1)
                        IF (iceList(i) .eq. re1(j)) THEN
                            !Basic rate at which thermal desorption occurs

                            rate(j)= 1 !dsqrt(VDIFF_PREFACTOR*bindingEnergy(i)/mass(iceList(i))) * exp(-gama(j)/dustTemp(dstep)) !!

                            !factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                            !becayse GRAIN_SURFACEAREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                            !that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                            !area seems more correct for this process.
                            IF (.NOT. THREE_PHASE) rate(j)=rate(j)*2.0*SURFACE_SITE_DENSITY*GRAIN_SURFACEAREA_PER_H
                        END IF
                    END DO
                END DO
                !At some point, rate is so fast that there's no point freezing out any more
                !Save the integrator some trouble and turn freeze out off
                WHERE(rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep)&
                &<MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
                IF (safeMantle .lt. MIN_SURFACE_ABUND) rate(idx1:idx2)=0.0
            ELSE
                rate(idx1:idx2)=0.0
            END IF
        END IF

        !Desorption due to UV, partially from ISRF and partially from CR creating photons
        idx1=deuvcrReacs(1)
        idx2=deuvcrReacs(2)
        IF (idx1 .ne. idx2) THEN
            IF ((desorb) .and. (uvdesorb) .and. (safeMantle .gt. MIN_SURFACE_ABUND)&
                    &.and.(zeta .gt. 0)) THEN
                !4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
                !UVY is yield per photon.
                rate(idx1:idx2) = 1 !GRAIN_CROSSSECTION_PER_H*uv_yield*4.875d3*zeta
                !additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of 
                !CR induced UV to ISRF UV.
                rate(idx1:idx2) = 1 !rate(idx1:idx2) * (1+(radfield/uvcreff)*(1.0/zeta)*dexp(-1.8*av(dstep)))

                !Don't remove species with binding energy > max BE removed by this process
                WHERE(gama(idx1:idx2) .gt. ebmaxuvcr) rate(idx1:idx2)=0.0 
            ELSE
                rate(idx1:idx2) = 0.0
            ENDIF
            !turn off freeze out if desorption due to UV is much faster
            WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep))&
            &<MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        END IF

        ! bulkSurfaceExchangeReactions
        idx1 = bulkswapReacs(1)
        idx2 = bulkswapReacs(2)
        IF (THREE_PHASE) THEN
            surfaceCoverage=0.5*GAS_DUST_DENSITY_RATIO/NUM_SITES_PER_GRAIN

            ! bulkToSurfaceSwappingRates
            IF ((dustTemp(dstep) .gt. MAX_GRAIN_TEMP) .or. (safeMantle .lt. MIN_SURFACE_ABUND)) THEN
                rate(idx1:idx2) = 0.0
            ELSE
                DO i=idx1,idx2
                    DO j=lbound(iceList,1),ubound(iceList,1)
                        IF (iceList(j) .eq. re1(i)) THEN
                        rate(i)= 1 !dsqrt(VDIFF_PREFACTOR*bindingEnergy(j)/mass(iceList(j)))*DEXP(-bindingEnergy(j)/dustTemp(dstep))
                        END IF
                    END DO
                END DO
            END IF

        END IF




        ! These are reactions which are untouched by carlos.

        ! CRP
        idx1=crpReacs(1)
        idx2=crpReacs(2)
        IF (idx1 .ne. idx2) THEN 
            rate(idx1:idx2) = alpha(idx1:idx2)*zeta
        END IF
        IF (improvedH2CRPDissociation) THEN
            rate(nR_H2_CRP)=h2CRPRate
        END IF

        !Account for Eley-Rideal reactions in a similar way.
        !First calculate overall rate and then split between desorption and sticking
        idx1=erReacs(1)
        idx2=erReacs(2)
        if (idx1 .ne. idx2) THEN
            rate(idx1:idx2)=freezeOutRate(idx1,idx2) *dexp(-gama(idx1:idx2)/dustTemp(dstep))
            rate(erdesReacs(1):erdesReacs(2))=rate(idx1:idx2)
            !calculate fraction of reaction that goes down desorption route
            idx1=erdesReacs(1)
            idx2=erdesReacs(2)
            DO j=idx1,idx2
                rate(j)=desorptionFraction(j)*rate(j)
                IF (ANY(bulkList==re1(j))) rate(j)=0.0 ! Bulk species are not able to chemically desorb
            END DO
            !remove that fraction from total rate of the diffusion route
            rate(erReacs(1):erReacs(2))=rate(erReacs(1):erReacs(2))-rate(idx1:idx2)
        END IF
        
        !Desorption due to energy from cosmic rays
        idx1=descrReacs(1)
        idx2=descrReacs(2)
        IF (idx1 .ne. idx2) THEN
            IF ((desorb) .and. (crdesorb) .and. (safeMantle .gt. MIN_SURFACE_ABUND)) THEN
                !4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
                !as iron nuclei are main cause of CR heating.
                !GRAIN_SURFACEAREA_PER_H is the total surface area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
                !phi is efficieny of this reaction, number of molecules removed per event.
                rate(idx1:idx2) = 4.0*pi*zeta*1.64d-4*(GRAIN_SURFACEAREA_PER_H)*phi

                !Don't remove species with binding energy > max BE removed by this process
                WHERE(gama(idx1:idx2) .gt. ebmaxcr) rate(idx1:idx2)=0.0 

            ELSE
                rate(idx1:idx2) = 0.0
            ENDIF
            !turn off freeze out if desorption due to CR formation is much faster
            WHERE((rate(freezePartners)*abund(re1(freezePartners),dstep)*density(dstep))&
            <MIN_SURFACE_ABUND*rate(idx1:idx2)) rate(freezePartners)=0.0
        END IF
    

        !CRS reactions represent the production of excited species from cosmic ray bombardment
        !rate equations from Shingledecker et. al. 2018
        idx1=crsReacs(1)
        idx2=crsReacs(2)
        IF (idx1 .ne. idx2) THEN
            !8.6 is the Spitzer-Tomasko cosmic ray flux in cm^-2 s^-1
            !1.3 converts to: ionisation rate/10^-17
            rate(idx1:idx2)=alpha(idx1:idx2)*(beta(idx1:idx2)*(gama(idx1:idx2)/100)*(8.6*zeta*1.3))
        END IF

        !EXRELAX, relaxation reactions for each excited species
        idx1=exrelaxReacs(1)
        idx2=exrelaxReacs(2)
        IF (idx1 .ne. idx2) THEN
            DO j=idx1,idx2
                DO i=lbound(iceList,1),ubound(iceList,1)
                    IF (iceList(i) .eq. re1(j)) THEN
                        vA=vdiff(i)
                    END IF
                END DO 
                rate(j)=vA
            END DO 
        END IF  

        !EXSOLID reactions represent the reactions of excited species on the grain
        idx1=exsolidReacs(1)
        idx2=exsolidReacs(2)
        IF (idx1 .ne. idx2) THEN
            !reaction rates calculated outside of UCLCHEM as per Shingledecker et al. 2018 and included in grain network
            !alpha are branching ratios and beta is reaction rate
            DO j=idx1,idx2
                DO i=lbound(iceList,1),ubound(iceList,1)
                    IF (iceList(i) .eq. re1(j)) THEN
                        vA = vdiff(i)
                    END IF
                    IF (iceList(i) .eq. re2(j)) THEN
                        vB = vdiff(i)
                    END IF
                END DO 
                rate(j) = (vB + vA)/(SURFACE_SITE_DENSITY*1.8d-8)
                rate(j) = alpha(j) * rate(j)
            END DO 
        END IF  


        !Reactions on surface can be treated considering diffusion of reactants
        !as in Langmuir-Hinshelwood mechanism
        !See work of David Quenard 2017 Arxiv:1711.05184
        !First calculate rate of the diffusion reaction
        idx1=lhReacs(1)
        idx2=lhReacs(2)
        if (idx1 .ne. idx2) THEN
            if ((dustTemp(dstep) .lt. MAX_GRAIN_TEMP) .and. (safeMantle .gt. MIN_SURFACE_ABUND)) THEN
                DO j=idx1,idx2
                    rate(j)=diffusionReactionRate(j,dustTemp(dstep))
                END DO
                !two routes for every diffusion reaction: products to gas or products remain on surface
                rate(lhdesReacs(1):lhdesReacs(2))=rate(idx1:idx2)

                !calculate fraction of reaction that goes down desorption route
                idx1=lhdesReacs(1)
                idx2=lhdesReacs(2)
                DO j=idx1,idx2
                    rate(j)=desorptionFraction(j)*rate(j)
                    IF (ANY(bulkList==re1(j))) rate(j)=0.0 ! Bulk species are not able to chemically desorb
                    IF (ANY(bulkList==re1(j))) rate(j)=0.0 ! Bulk species are not able to chemically desorb
                END DO
                !remove that fraction from total rate of the diffusion route
                rate(lhReacs(1):lhReacs(2))=rate(lhReacs(1):lhReacs(2))-rate(idx1:idx2)
            ELSE
                rate(idx1:idx2)=0.0
                rate(lhdesReacs(1):lhdesReacs(2))=0.0
            END IF
        END IF

        idx1 = surfSwapReacs(1)
        idx2 = surfSwapReacs(2)
        IF ((dustTemp(dstep) .gt. MAX_GRAIN_TEMP) .or. (safeMantle .lt. MIN_SURFACE_ABUND)) THEN
            rate(idx1:idx2) = 0.0
        ELSE
            rate(idx1:idx2) = 1.0
        END IF

        idx1=ionopol1Reacs(1)
        idx2=ionopol1Reacs(2)
        IF (idx1 .ne. idx2)&
        !This formula including the magic numbers come from KIDA help page.
        &rate(idx1:idx2)=alpha(idx1:idx2)*beta(idx1:idx2)*(0.62d0+0.4767d0*gama(idx1:idx2)*dsqrt(300.0d0/gasTemp(dstep)))

        idx1=ionopol2Reacs(1)
        idx2=ionopol2Reacs(2)
        IF (idx1 .ne. idx2) THEN
            !This formula including the magic numbers come from KIDA help page.
            rate(idx1:idx2)=alpha(idx1:idx2)*beta(idx1:idx2)*(1.0d0+0.0967d0*gama(idx1:idx2)&
            &*dsqrt(300.0d0/gasTemp(dstep))+gama(idx1:idx2)*gama(idx1:idx2)*300.0/(10.526*gasTemp(dstep)))
        END IF

        !turn off reactions outside their temperature range
        WHERE(.not. ExtrapolateRates .and. (gasTemp(dstep) .lt. minTemps)) rate=0.0

        WHERE(.not. ExtrapolateRates .and. (gasTemp(dstep) .gt. maxTemps)) rate=0.0
    END SUBROUTINE calculateReactionRates



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !Freeze out determined by rate of collisions with grain
    !No sticking coefficient is used because typical values are >0.95 below 150 K
    ! eg Le Bourlot et al. 2013, Molpeceres et al. 2020
    !Above 150 K, thermal desorption will completely remove grain species
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION freezeOutRate(idx1,idx2) RESULT(freezeRates)
        REAL(dp) :: freezeRates(idx2-idx1+1)
        INTEGER :: idx1,idx2
        
        !additional factor for ions (beta=0 for neutrals)
        freezeRates=1.0+beta(idx1:idx2)*16.71d-4/(GRAIN_RADIUS*gasTemp(dstep))
        IF ((freezeFactor .eq. 0.0) .or. (dustTemp(dstep) .gt. MAX_GRAIN_TEMP)) then
            freezeRates=0.0
        ELSE
            freezeRates=freezeRates*freezeFactor*alpha(idx1:idx2)*THERMAL_VEL*dsqrt(gasTemp(dstep)/mass(re1(idx1:idx2)))*GRAIN_CROSSSECTION_PER_H
        END IF

        END FUNCTION freezeOutRate


        FUNCTION stickingCoefficient(stickingZero,criticalTemp,gasTemp) RESULT(stickingCoeff)
            !Sticking coefficient for freeze out taken from Chaabouni et al. 2012 A&A 538 Equation 1
            REAL(dp) :: stickingCoeff
            REAL(dp) :: stickingZero,criticalTemp,gasTemp,tempRatio
            
            stickingCoeff=stickingZero*(1.0d0+(2.5d0)*(gasTemp/criticalTemp))/((1.0d0+(gasTemp/criticalTemp))**(2.5d0))
        END FUNCTION stickingCoefficient

END MODULE RATES