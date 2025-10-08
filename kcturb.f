c***********************************************************************
c  KCTURB calculates the turbulent fluxes of heat, momentum, and 
c  humidity and in the atmospheric surface layer (lowest 10% of the 
c  of the planetary boundary layer on Earth, Mars, and Titan.
c  It includes: 
c (a) Monin-Obukhov similarity functions that encompass the entire range 
c     of atmospheric stability expected on terrestrial planetary bodies, 
c (b) the additional shear associated with buoyant plumes under unstable
c     conditions,  
c (c) surface renewal theory to calculate transfer rates within the 
c     interfacial layer adjacent to the surface, and 
c (d) key humidity effects that become especially important when a 
c     volatile is more buoyant than the ambient gas 
c     (e.g., on Mars where H2O is lighter than CO2). 
c  REFERENCES:
c  Brutsaert, W. (1982). Evaporation into the atmosphere: Theory and 
c  applications. Springer.
c  Andreas, E. L., Horst, T. W., Grachev, A. A., Persson, P. O. G., 
c  Fairall, C. W., Guest, P. S., & Jordan, R. E. (2010). Parametrizing 
c  turbulent exchange over summer sea ice and the marginal ice zone. 
c  Quarterly Journal of the Royal Meteorological Society, 136(649), 
c  927–943.
c  Andreas, E. L. (1987). A theory for the scalar roughness and the 
c  scalar transfer coefficients over snow and sea ice. 
c  Boundary-Layer Meteorology, 38(1), 159–184. 
c  Khuller, A. R., & Clow, G. D. (2024). Turbulent fluxes and 
c  evaporation/sublimation rates on earth, mars, titan, and exoplanets. 
c  Journal of Geophysical Research: Planets, 129(4), e2023JE008114.
cc***********************************************************************     
      subroutine kcturb(Cpd0,Fh02,Fw02,H0sen,lambda0,LE,Ls,qlat,qsen,q2,
     & rho_rep,rhocp_rep,Fm,
     & planet,z0,Zmax,H_PBL,P_in,beta,beta2,
     & Fx_opt,z1,z2_in,z3_in,Theta1,Theta2,Rh1,Rh2,wsp,frh,ltype,
     & snowdepth,dliqvol,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

      implicit none

      double precision snowdepth,z1,z2_in,z3_in,P_in
      double precision vk,h,H_PBL,P,z0,gacc,Mw,Ma,b,Zmax,beta,beta2
      double precision Theta1,Theta2,Rh1,Rh2,wsp,U3,frh,dliqvol
      double precision q1,q2,Cp1,Cp2,rho1,rho2,rhoa,rhow,dummy1,dummy2
      double precision rhocp1,rhocp2,Dtheta,Dq,Pr,Sc,nu,h0,hd
      double precision viscLs,viscLr,viscL,ustar,Re,U3_og
      double precision Cpd0,Cpw0,Fh02,Fw02,H0sen,Hq,lambda0
      double precision LE,Ls,qlat,qsen,qlat_a1,qlat_a2
      double precision theta_rep,rh_rep,Cp_rep,rho_rep,q_rep,rhocp_rep
      double precision rho(2),dq2,dtheta2,E0,es
      double precision z0m,minimum_value,dummy3
      double precision md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b
      double precision md_a0,md_a1,md_a2,md_a3,md_a4,md_b0,md_b1,md_b2
      double precision md_b3,md_b4,md_beta_a,md_beta_b,md_gamma_a
      double precision md_gamma_b,md_kappa_a,md_kappa_b, md_m_a,md_m_ab
      double precision md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b
      double precision md_sum_v_a,md_sum_v_b,d0,Fm,fac1,fac2,eb
      double precision Fhs,Fws,Fh,Fw,Fh12,Fw12,Fhi,Fwi,theta_star,q_star
      double precision zeta2,zeta3,L,theta0,rh0,q0,rho0,rhocp0,wstar
      double precision resid,rel,Conv,TESTQ,TESTT,TESTU,Rib
      double precision theta_star_old,q_star_old,ustar_old,Htheta,Eflux
      double precision thetad,rhd,qd,rhod,rhocpd,theta0l,resid0,Cp0 
      double precision Z(200),theta(200),q(200),rhocp(200),rh(200)
      double precision zeta(200),junk(200),z2,z3

      integer boil,imode,ltype,Fx_opt,planet,dryair,volatile,ic,icmax
      integer i,N,p1,p2,imid


c     SNTHERM assumes that the met data heights are measured relative to the ground underlying the snow/ice
      if (ltype .eq. 1) then ! might have to update to account for firn/glacier ice
            z2 = z2_in - snowdepth
            z3 = z3_in - snowdepth
      endif

c     general parameters
c      von karman's constant
      vk = 0.40d0 

c     Height of surface layer (taken to be 1/10th of H, the height of the mixed layer)

      h = 0.1d0*H_PBL

c     Perform input conversions

      P  = 100.0d0 * P_in ! Atmospheric surface pressure (mb -> Pa)
      

c     Set planet specific parameters
c     Planet:   1 = Earth, 2 = Mars, 3 = Titan
c     dryair:   1 = air, 2 = CO2, 3 = N2
c     volatile: 1 = H2O, 2 = CH4
      
      if (planet .eq. 1) then                ! Earth
            gacc     = 9.78d0  ! gravitational acceleration (m/s^2)
            dryair   = 1       ! predominant dry air component (air)
            volatile = 1       ! predominant volatile species  (H2O)
      else if (planet .eq. 2) then           ! Mars
            gacc     = 3.72d0  ! gravitational acceleration (m/s^2)
            dryair   = 2       ! predominant dry air component (CO2)
            volatile = 1       ! predominant volatile species  (H2O)
      else if (planet .eq. 3) then           ! Titan
            gacc     = 1.352d0 ! gravitational acceleration (m/s^2)
            dryair   = 3       ! predominant dry air component (N2)
            volatile = 2       ! predominant volatile species  (CH4)
      endif

c     Set molecular weight of volatile species
      
      ! set molecular weight of dry air
      if (dryair .eq. 1) then
            Ma = 0.028966d0    ! (modern earth)
      else if (dryair .eq. 2) then
            Ma = 0.0440098d0   ! CO2
      else if (dryair .eq. 3) then
            Ma = 0.0280d0      ! N2
      endif
      
      ! set molecular weight of volatile
      if (volatile .eq. 1) then
            Mw = 0.018016d0
      else if (volatile .eq. 2) then
            Mw = 0.016043d0
      endif

      b = Ma/Mw - 1d0


c     Determine whether height z1 is located at the surface (z1 = 0)

      if (z1 < z0) then
          imode = 1 ! z1 is likely at the surface
      else
          imode = 2 ! z1 is likely within the surface sublayer (z1 >= hd)
      endif


c     Estimate interfacial sublayer height hd  (so we can setup vertical grid)

c     Find (q,rho,rhocp) at met heights (z1,z2)
      
      ! Find dry air and water vapor densities at z1
      CALL rho_sub_new(rhoa,rhow,rho1,q1,boil,
     &                       dryair,volatile,P,Theta1,Rh1)

      ! Find dry air and water vapor densities at z1
      CALL rho_sub_new(dummy1,dummy2,rho2,q2,boil,
     &                       dryair,volatile,P,Theta2,Rh2)

      ! Find specific heat of air at z1
      CALL Cp_sub(dummy1,dummy2,Cp1,
     &                       dryair,volatile,Theta1,q1)

      ! Find specific heat of air at z2
      CALL Cp_sub(dummy1,dummy2,Cp2,
     &                       dryair,volatile,Theta2,q2)

      
      rhocp1              = rho1 * Cp1
      rhocp2              = rho2 * Cp2

      
c     Find met data temperature and humidity differences between heights (z1,z2)

      Dtheta = Theta2 - Theta1
      Dq     = q2 - frh*q1
      
c     Estimate (Pr,Sc,nu) from met data at height z1

      
      rho(1) = rhoa ! density of dry air (kg/m^3) 
      rho(2) = rhow ! density of water vapor (kg/m^3)
      
      
      ! Calculate the Prandtl & Schmidt numbers (dimensionless), 
      ! along with the kinematic viscosity (m^2/s)
      
      CALL mol_diffus_v4_sntherm(Pr,Sc,nu,
     & planet,Theta1,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

c     Estimate (viscL,Re) assuming a neutral atmosphere

      h0     = 7.35d0 * z0     ! height of roughness elements for aerodynamically rough case (m)
      d0     = (2d0/3d0) * h0  ! zero-plane displacement height (m)
      
c     Calculate the primitive function for momentum
      
      CALL Fm_sub_update(Fm,
     & Fx_opt,z0,d0,z3,0d0)

      U3    = wsp ! save SNTHERM's original wind speed
      
      U3_og = U3 ! save original wind speed (can be modified inside kcturb)
      
      if (U3 .lt. 0.1d0) then
          U3 = 0.1d0  ! to prevent the formation of a singularity and a wormhole
      endif


      ustar  = U3 / Fm           ! friction velocity (m/s)
      viscLr = nu / ustar        ! viscous length scale (m)
      z0m    = 0.135d0 * viscLr  ! surface roughness length for momentum (aero smooth case)
      d0     = 0d0               ! zero-plane displacement height (m)
      
c     Calculate the primitive function for momentum again, using z0m and a neutral atmosphere
      
      CALL Fm_sub_update(Fm,
     & Fx_opt,z0m,d0,z3,0d0)

      ustar  = u3 / Fm             ! friction velocity (m/s)
      viscLs = nu / ustar          ! viscous length scale (m)
      viscL  = (viscLr+viscLs)/2d0 ! mean of smooth and rough values for viscous length scale (m)
      
      Re     = z0 / viscL          ! Roughness Reynolds number (dimensionless)

      
c     Estimate (z0m,d0,hd) from neutral estimate of (viscL,Re)

      if (Re .le. 0.135d0) then    ! aero smooth
            z0m = 0.135d0 * viscL  ! surface roughness length for momentum (m)
            d0  = 0d0              ! zero-plane displacement height (m)
            hd  = 30d0 * viscL     ! interfacial layer height (m)
      else if (Re .gt. 0.135d0 .and. Re .lt. 2d0) then  ! aero transition
            z0m = z0               ! surface roughness length for momentum (m)
            h0  = 7.35d0 * z0      ! height of roughness elements (m)
            d0  = (2d0/3d0) * h0   ! zero-plane displacement height (m)
            hd  =  (1d0/(1d0+Re**2d0))*(30d0 * viscL 
     &             + Re**2d0*(d0 + 7.4d0 * z0)) ! interfacial layer height (m)
      else if (Re .ge. 2d0) then   ! aero rough
            z0m = z0               ! surface roughness length for momentum (m)
            h0  = 7.35d0 * z0      ! height of roughness elements (m)
            d0  = (2d0/3d0) * h0   ! zero-plane displacement height (m)
            hd  = d0 + 7.4d0 * z0  ! interfacial layer height (m)
      end if

c     Estimate ustar assuming a neutral atmosphere

      ! Calculate the primitive function for momentum again, using z0 and a neutral atmosphere
      
      CALL Fm_sub_update(Fm,
     & Fx_opt,z0,d0,z3,0d0)
      
      ustar  = u3 / Fm             ! friction velocity (m/s)


c set seed values for representative values within the surface sublayer

      theta_rep = Theta2
      q_rep     = q2
      rho_rep   = rho2
      rhocp_rep = rhocp2

c estimate scaling temperature and humidity (theta_star,q_star)
      
      if (imode .eq. 1) then ! z1 is likely at the surface
    
      ! Calculate the primitive function for heat and humidity
      
            CALL Fhw_subi_v2(Fhi,Fwi,
     &           Re,Pr,Sc)                  ! interfacial sublayer
    
            CALL Fhw_sub_update(Fhs,Fws,
     &             Fx_opt,d0,hd,z2,0d0)     ! surface sublayer
        
        Fh12      = (rhocp_rep / rhocp1) * Fhi + Fhs ! Eq. 52
        Fw12      = (rho_rep / rho1)     * Fwi + Fws ! Eq. 51
      
      else if (imode .eq. 2) then ! z1 is likely within the surface sublayer (z1 >= hd)
      
      ! Calculate the primitive function for heat and humidity
      
            CALL Fhw_sub_update(Fhs,Fws,
     &             Fx_opt,d0,z1,z2,0d0)     ! surface sublayer
      
        Fh12      = Fhs
        Fw12      = Fws
      
      end if

      theta_star = Dtheta / Fh12  ! turbulent temperature scale (K) <theta_star>_12
      q_star     = Dq / Fw12      ! turbulent humidity scale (dimensionless) <q_star>_12

c     Estimate turbulent kinetic energy (TKE) buoyancy production eb (m^2/s^3)

      fac1 = (1d0 / theta_rep)       * ustar * theta_star
      fac2 = (b / (1d0 + b * q_rep)) * ustar * q_star
      eb   = -gacc * (fac1 + fac2) ! buoyant production of TKE (m^2/s^3)

c     Estimate Obukhov length L

      L = -ustar**3d0/ (vk * eb) ! Obukhov length (m)
      
      
c     Get improved estimate of ustar

      zeta3 = (z3 - d0) / L ! dimensionless stability parameter, z/L
      
      ! Calculate the primitive function for momentum again, using z0m and zeta3
      CALL Fm_sub_update(Fm,
     &       Fx_opt,z0m,d0,z3,zeta3)
      
      ustar = u3 / Fm ! friction velocity (m/s)
      
c     Get improved value of (viscL,Re)

      viscL = nu / ustar ! viscous length scale (m)
      Re    = z0 / viscL ! Roughness Reynolds number (dimensionless)

c     Get improved value of the interfacial layer height, hd

      if (Re .le. 0.135d0) then                 ! aero smooth
          z0m = 0.135d0* viscL ! surface roughness length for momentum (m)
          d0  = 0d0            ! zero-plane displacement height (m)
          hd  = 30d0 * viscL   ! interfacial layer height (m)
      else if (Re .gt. 0.135d0 .and. Re .lt. 2d0) then  ! aero transition
          z0m = z0             ! surface roughness length for momentum (m)
          h0  = 7.35d0 * z0    ! height of roughness elements (m)
          d0  = (2d0/3d0) * h0 ! zero-plane displacement height (m)
          hd  =  (1d0/(1d0+Re**2d0))*(30d0 * viscL 
     &       + Re**2d0*(d0 + 7.4d0 * z0)) ! interfacial layer height (m)
      elseif (Re .ge. 2d0) then                          ! aero rough
          z0m = z0              ! surface roughness length for momentum (m)
          h0  = 7.35d0 * z0     ! height of roughness elements (m)
          d0  = (2d0/3d0) * h0  ! zero-plane displacement height (m)
          hd  = d0 + 7.4d0 * z0 ! interfacial layer height (m)
      end if


c  Setup vertical grid Z
      
      CALL SETUP_GRID(hd, Zmax, Z, N) ! Z(1:N) contains grid points
      
!***********************************************************************
C     surf_layer_advan
      
       U3 = U3_og  !reset U3 to input value to follow surf_layer_advan method

C     Find met height indices

      do 123 i = 1,N
            if (Z(i) .ge. z1) then
            p1 = i
            exit
            end if
123   continue

      do 201 i = 1,N
            if (Z(i) .ge. z2) then
            p2 = i
            exit
            end if
201   continue


c     Set surface values if imode = 1 (z1 = 0)

      if (imode .eq. 1) then
    
            theta0          = Theta1
            rh0             = rh1
            q0              = q1
            rho0            = rho1
            rhocp0          = rhocp1
        
      ! Find dry air and water vapor densities at z1
            CALL rho_sub_new(rhoa,rhow,dummy1,dummy2,boil,
     &                       dryair,volatile,P,theta0,rh0)
            
            rho(1) = rhoa
            rho(2) = rhow

            CALL mol_diffus_v4_sntherm(Pr,Sc,nu,
     & planet,theta0,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

  
      else if (imode .eq. 2) then

            CALL rho_sub_new(rhoa,rhow,dummy1,dummy2,boil,
     &                       dryair,volatile,P,Theta1,Rh1)

            rho(1) = rhoa
            rho(2) = rhow

             CALL mol_diffus_v4_sntherm(dummy1,dummy2,nu,
     & planet,Theta1,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

      end if


c     Estimate Obukhov length L assuming a neutral atmosphere

c     Estimate the viscous length scale using Re

      viscL = z0 / Re
      

c     Estimate (z0m,d0,hd)

      if (Re .le. 0.135d0) then                 ! aero smooth
          z0m = 0.135d0* viscL
          d0  = 0d0
          hd  = 30d0 * viscL
      else if (Re .gt. 0.135d0 .and. Re .lt. 2d0) then  ! aero transition
          z0m = z0
          h0  = 7.35d0 * z0
          d0  = (2d0/3d0) * h0
          hd  =  (1d0/(1d0+Re**2d0))*(30d0 * viscL 
     &       + Re**2d0*(d0 + 7.4d0 * z0))
      elseif (Re .ge. 2d0) then                          ! aero rough
          z0m = z0
          h0  = 7.35d0 * z0
          d0  = (2d0/3d0) * h0
          hd  = d0 + 7.4d0 * z0
      end if


c     Estimate ustar assuming a neutral atmosphere

      CALL Fm_sub_update(Fm,
     & Fx_opt,z0m,d0,z3,0d0)

      if (U3 .lt. 0.1d0) then
          U3 = 0.1d0  ! to prevent the formation of a singularity and a wormhole
      endif


      ustar  = U3 / Fm ! friction velocity (m/s)

c     Set seed values for representative values within the surface sublayer

      theta_rep = theta2
      q_rep     = q2
      rho_rep   = rho2
      rhocp_rep = rhocp2


c     Estimate scaling temperature and humidity (theta_star,q_star)

      if (imode .eq. 1) then
    

            CALL Fhw_subi_v2(Fhi,Fwi,
     &           Re,Pr,Sc)                  ! interfacial sublayer
    
            CALL Fhw_sub_update(Fhs,Fws,
     &             Fx_opt,d0,hd,z2,0d0)     ! surface sublayer
      
            Fh02       = (rhocp_rep / rhocp0) * Fhi + Fhs
            Fw02       =   (rho_rep / rho0)   * Fwi + Fws
            theta_star = Dtheta / Fh02                    ! <theta_star>_02
            q_star     =     Dq / Fw02                    ! <q_star>_02
            
      
      else if (imode .eq. 2) then

            CALL Fhw_sub_update(Fhs,Fws,
     &             Fx_opt,d0,z1,z2,0d0)     ! surface sublayer
      
            Fh12       = Fhs
            Fw12       = Fws
            theta_star = Dtheta / Fh12                     ! <theta_star>_12
            q_star     =     Dq / Fw12                    ! <q_star>_12
      
      end if


c     Estimate TKE buoyancy production eb

      fac1 = (1d0 / theta_rep)       * ustar * theta_star
      fac2 = (b / (1d0 + b * q_rep)) * ustar * q_star
      eb   = -gacc * (fac1 + fac2) ! buoyant production of TKE (m^2/s^3)

c     Estimate Obukhov length L

      L = -ustar**3d0 / (vk * eb) ! Obukhov length (m)

      

c     Estimate the bulk Richardson number Rib

      fac1 = gacc / (theta_rep * (z2 - z1)) * (z3 / u3)**2d0
      fac2 = Dtheta  + b * theta_rep / (1d0 + b * q_rep) * Dq
      Rib  = fac1 * fac2

      
      
c     Estimate properties within the surface sublayer
      Conv   = 1d-2             ! star convergence criteria
      resid  = 9999d0           ! absolute convergence criteria
      rel    = 9999d0           ! relative convergence criteria
      ic     = 1                ! iteration counter
      icmax  = 400              ! max number of iterations
      TESTU  = 9999d0           ! ustar convergence criteria
      TESTT  = 9999d0           ! thetastar convergence criteria
      TESTQ  = 9999d0           ! qstar convergence criteria


c ------- begin iteration loop ------
      do 223 while (ic .le. icmax .and. (TESTT .gt. Conv .or. 
     &       TESTU .gt. Conv .or. TESTQ .gt. Conv))   

            
            ustar_old      = ustar
            theta_star_old = theta_star
            q_star_old     = q_star

c     Find Monin-Obukhov stability parameter zeta on the Z-grid

            do 245 i = 1,N
                  zeta(i)    = (Z(i) - d0)  / L
245         continue        
        
            zeta2   = (z2 - d0) / L          ! zeta(z2)
            zeta3   = (z3  - d0) / L         ! zeta(h)

            zeta(1) = 0d0

c     Find ustar

            CALL Fm_sub_update(Fm,
     &       Fx_opt,z0m,d0,z3,zeta3)
      
            ustar = u3 / Fm

c     Find viscous length-scale and roughness Reynolds number

            viscL = nu / ustar
            Re    = z0 / viscL


c     Estimate (z0m,d0,hd)

            if (Re .le. 0.135d0) then                 ! aero smooth
                  z0m = 0.135d0* viscL
                  d0  = 0d0
                  hd  = 30d0 * viscL
            else if (Re .gt. 0.135d0 .and. Re .lt. 2d0) then  ! aero transition
                  z0m = z0
                  h0  = 7.35d0 * z0
                  d0  = (2d0/3d0) * h0
                  hd  =  (1d0/(1d0+Re**2d0))*(30d0 * viscL 
     &                  + Re**2d0*(d0 + 7.4d0 * z0))
            elseif (Re .ge. 2d0) then                 ! aero rough
                  z0m = z0
                  h0  = 7.35d0 * z0
                  d0  = (2d0/3d0) * h0
                  hd  = d0 + 7.4d0 * z0
            end if

c     Find point at the middle of the surface layer

            do 250 i = 1,N
                  junk(i) = abs(Z(i) - (abs(L) - hd)/2d0)
250         continue
        
            
            minimum_value = junk(1)
            imid = 1
            
            do 260 i = 2, N
                  if (junk(i) .lt. minimum_value) then
                        minimum_value = junk(i)
                        imid = i
                  end if
260         continue
            
c     Use values midway through the surface sublayer as representative values
c     for (theta,q,rho,rhocp)

            if (imid .ne. 1) then

                  CALL Fhw_sub_update(Fh,Fw,
     &            Fx_opt,d0,z2,Z(imid),zeta(imid))     ! surface sublayer
            
            
                  dtheta2   = theta_star * Fh
                  dq2      =     q_star * Fw
                  theta_rep = theta2 + dtheta2
                  q_rep     =     q2 + dq2
                  
                  CALL Rh_sub(dummy1,rh_rep,
     &                 dryair,volatile,P,theta_rep,q_rep)

                  CALL rho_sub_new(dummy1,dummy2,rho_rep,dummy3,boil,
     &                 dryair,volatile,P,theta_rep,rh_rep)
            
                  CALL Cp_sub(dummy1,dummy2,Cp_rep,
     &                 dryair,volatile,theta_rep,q_rep)
            
                  rhocp_rep       = rho_rep * Cp_rep

            end if

c     Restore array values at z1,z2 in case they're a little of
            
            theta(p1) = theta1                                         
            q(p1)     = q1                                         
            theta(p2) = theta2
            q(p2)     = q2


c     Find (theta_star,q_star) from met data using improved values for
c     (Re,d0,hd,zeta2,rho_rep,rhocp_rep)

            if (imode .eq. 1) then
    

                  CALL Fhw_subi_v2(Fhi,Fwi,
     &            Re,Pr,Sc)                  ! interfacial sublayer
                  
                  CALL Fhw_sub_update(Fhs,Fws,
     &            Fx_opt,d0,hd,z2,zeta2)     ! surface sublayer
      
                  Fh02       = (rhocp_rep / rhocp0) * Fhi + Fhs
                  Fw02       =   (rho_rep / rho0)   * Fwi + Fws
                  theta_star = Dtheta / Fh02                    ! <theta_star>_02
                  q_star     =     Dq / Fw02                    ! <q_star>_02
            
      
            else if (imode .eq. 2) then

                  CALL Fhw_sub_update(Fh12,Fw12,
     &            Fx_opt,d0,z1,z2,zeta2)     ! surface sublayer
      
            
                  theta_star = Dtheta / Fh12                    ! <theta_star>_12
                  q_star     =     Dq / Fw12                    ! <q_star>_12
      
      end if


c     Find TKE buoyancy production eb

            fac1 = (1d0 / theta_rep)       * ustar * theta_star
            fac2 = (b / (1d0 + b * q_rep)) * ustar * q_star
            eb   = -gacc * (fac1 + fac2)

c     Find temperature and humidity fluxes (Htheta,E)

            Htheta = -rhocp_rep * ustar * theta_star
            Eflux     = -rho_rep   * ustar * q_star
        
            

c     Don't adjust windspeed if atmosphere is basically neutral 
c     (Andreas et al., 2010 - SPEED subroutine)
            if (abs(L) .lt. 1000d0) then

                  if (L .lt. 0d0) then  
                        wstar = ustar*((-H_PBL/(vk*L))**0.333333d0)
                        u3 = sqrt(U3_og**2d0+ (beta*wstar)**2d0)

                  else if (L .gt. 0d0) then
c     Modified Windless Coefficient described in Eq. 2.6 of 
c     Andreas et al. (2010), based on Jordan et al. (1999)
                        u3 =  U3_og + beta2/cosh(U3_og)
                  end if
            
            end if

            CALL Fm_sub_update(Fm,
     &      Fx_opt,z0m,d0,z3,zeta3)

            ustar = u3 / Fm

c     Find TKE buoyancy production eb

            fac1 = (1d0 / theta_rep)       * ustar * theta_star
            fac2 = (b / (1d0 + b * q_rep)) * ustar * q_star
            eb   = -gacc * (fac1 + fac2)

c     Estimate Obukhov length L

            L = -ustar**3d0 / (vk * eb)
            
            !write(*,*) 'viscL: ', viscL
            !write(*,*) 'Re: ', Re
            !write(*,*) 'Pr: ', Pr
            !write(*,*) 'Sc: ', Sc
            !write(*,*) 'eb: ', eb
            
            !write(*,*) 'ustar: ', ustar
            !write(*,*) 'q_star: ', q_star
            !write(*,*) 'q_rep: ', q_rep
            !write(*,*) 'theta_rep: ', theta_rep
            !write(*,*) 'theta_star: ', theta_star
            
c     Find the bulk Richardson number Rib

            fac1 = gacc / (theta_rep * (z2 - z1)) * (z3 / u3)**2d0
            fac2 = Dtheta  + b * theta_rep / (1d0 + b * q_rep) * Dq
            Rib  = fac1 * fac2

c     Find surface values (theta0,q0) and improved viscosity nu if imode = 2


            if (imode .eq. 2) then

c     Set values at z = hd

                  thetad = theta(2)
                  rhd    = rh(2)
                  qd     = q(2)
                  rhod   = rho(2)
                  rhocpd = rhocp(2)

c     Estimate (theta0,q0) using (Pr,Sc,nu,Re,rho,rhocp) values at z = hd

                  CALL rho_sub_new(rhoa,rhow,dummy1,dummy2,boil,
     &                 dryair,volatile,P,thetad,rhd)
                
                  rho(1) = rhoa
                  rho(2) = rhow
                  
                  CALL mol_diffus_v4_sntherm(Pr,Sc,nu,
     & planet,thetad,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)


                  viscL           = nu / ustar
                  Re              = z0 / viscL

                  CALL Fhw_subi_v2(Fhi,Fwi,
     &            Re,Pr,Sc)                  ! interfacial sublayer

                  dtheta2          = -Fhi * Htheta / (ustar * rhocpd)
                  dq2              = -Fwi * Eflux      / (ustar * rhod)
                  theta0          = thetad - dtheta2
                  q0              =     qd - dq2

c     Iterate for improved estimates of (Pr,Sc,nu,Re,rho0,rhocp0,theta0,q0) at z = 0

                  theta0l = theta0
                  resid0   = 99999d0

                  do 270 while (abs(resid0) .gt. 1.0d-6)


                        CALL Rh_sub(dummy1,rh0,
     &                  dryair,volatile,P,theta0,q0)

                        CALL rho_sub_new(rhoa,rhow,rho0,dummy3,boil,
     &                  dryair,volatile,P,theta0,rh0)
            
                        CALL Cp_sub(dummy1,dummy2,Cp0,
     &                  dryair,volatile,theta0,q0)

                    
                        rhocp0             = rho0 * Cp0

                        rho(1) = rhoa
                        rho(2) = rhow
                  
                        CALL mol_diffus_v4_sntherm(Pr,Sc,nu,
     & planet,theta0,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

                    
                        viscL      = nu / ustar
                        Re         = z0 / viscL
                    
                        CALL Fhw_subi_v2(Fhi,Fwi,
     &                  Re,Pr,Sc)                  ! interfacial sublayer
                    
                        dtheta2     = -Fhi * Htheta / (ustar * rhocp0)
                        dq2         = -Fwi * Eflux      / (ustar * rho0)
                        theta0     = thetad - dtheta2
                        q0         =     qd - dq2

                        resid0     = theta0 - theta0l
                        theta0l    = theta0

270               continue
            end if


c     Update array values

            theta(1) = theta0
            rh(1)    = rh0
            q(1)     = q0
            rho(1)   = rho0
            rhocp(1) = rhocp0

c     Test for convergence
            TESTU = abs((ustar - ustar_old)/ustar)
            TESTT = abs((theta_star - theta_star_old)/theta_star)
            TESTQ = abs((q_star - q_star_old)/q_star)
            ic = ic + 1

c     For conditions implying that TESTT/U/Q is NaN,
c     set TESTT/U/Q to below CONV

            if (theta_star .eq. 0d0 .and. theta_star_old .eq. 0d0) then
                  TESTT = Conv-1d0
            elseif (ustar .eq. 0d0 .and. ustar_old .eq. 0d0) then
                  TESTU = Conv-1d0
            elseif (q_star .eq. 0d0 .and. q_star_old .eq. 0d0) then
                  TESTQ = Conv-1d0
            end if

223   continue      
c ---------- end of main iteration loop ---------------
  
c     Adjust Z(2) to track interfacial height hd
      Z(2) = hd                        


c     Find surface values

      CALL Cp_sub(Cpd0,Cpw0,dummy2,
     &                  dryair,volatile,theta0,q0)

    
      lambda0      = Cpw0 / Cpd0 - 1d0
      Hq           = (lambda0 * Cpd0 * theta0) * Eflux
      H0sen        = Htheta  + Hq
      E0           = Eflux
    
      CALL Ls_sub_sntherm(Ls,
     &     theta0,dliqvol)
    
      LE           = Ls*E0
    
      
      CALL es_sub_new(es,
     &       volatile,theta2)
      
      qlat_a1 = rh2*es/100d0
      
      CALL es_sub_new(es,
     &       volatile,theta1)
      
      qlat_a2 = rh1*es/100d0
      
      qlat = LE/(qlat_a1 - qlat_a2)
      qsen = H0sen/(theta2-theta1)

      return
      end
