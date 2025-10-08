c***********************************************************************
      subroutine condgas(K,
     &   beta,cv,m,T,Tc,viscosity)
c***********************************************************************

      implicit none

      double precision R,m,T,Tc,Tr
      double precision Z,cv,psi,alpha,beta
      double precision viscosity,K

c INPUTS:
c m         = molar mass in g/mol
c viscosity = dynamic viscosity in kg/(m s) or N s/m^2
c cv        = heat capacity at constant volume in J/ (mol K)

c OUTPUTS:
c K         = gas thermal conductivity  (W/(m K))

c CALLED FROM: mol_diffuse_v4_sntherm
c CALLS:       None

     
      R = 8.314462618d0 ! Gas constant, J/K/mol

      Tr = T/Tc       ! Reduced temperature

      Z = 2d0 + 10.5d0* Tr**2d0

      alpha = cv/R -1.5d0

      psi = 1d0 + alpha*((0.215d0 + 0.28288d0*alpha 
     &       - 1.061d0*beta + 0.26665d0*Z)
     &     /(0.6366d0 + beta*Z + 1.061d0*alpha*beta )   );


      K = viscosity*3.75d0*psi*R/(m*0.001d0) ! Thermal conductivity, W/(m K), Eq. 10-3.14, pg. 576


      return
      end
c***********************************************************************
      subroutine Cp_sub(Cpd,Cpw,Cp,dryair,volatile,T,q)                       
c***********************************************************************

      implicit none
      double precision Cpd,Cpw,Cp,T,q
      double precision T0,Ma,r,Mw
      integer dryair,volatile


c     Finds the specific heat of dry air, volatile species, and moist air.


c INPUTS:
c     dryair   = predominant content of dry air ('air','CO2','N2')
c     volatile = predominant volatilve('H2O','CH4')
c     T        = temperature (K)
c     q        = specific humidity (unitless)

c OUTPUTS:
c   Cpd      = specific heat of dry air (J/kg/K)
c   Cpw      = specific heat of water vapor or other volatile (J/kg/K)
c   Cp       = total specific heat (J/kg/K)

c Notes:

c     Cpd, Cpw values may be suspect for Titan (NIST equations don't go as low.


      T0  = 273.15d0


c     Specific heat of dry air

      if (dryair .eq. 1) then
     
        Cpd = 1005d0
     
      else if (dryair .eq. 2) then
     
        Cpd = 826.66d0 +(0.91315d0)*(T - T0)
     
      else if (dryair .eq. 3) then            
      ! https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas
        
        Ma = 0.0280d0
        
        Cpd = 28.98641d0 + 1.853978d0*(T/1000d0) 
     &   - 9.647459d0*(T/1000d0)**2d0
     &   + 16.63537d0*(T/1000d0)**3d0 + 0.000117d0/(T/1000d0)**2d0
        
        Cpd = Cpd/Ma ! to convert from J/mol/K to J/kg/K
      
      end if


c     Specific heat of water vapor (assumed to be ideal gas) or other volatile

      if (volatile .eq. 1) then
     
        r   = 2283d0/ T
     
        Cpw = 462.8d0 * (4d0 + r**2d0 *exp(-r))
     
      else if (volatile .eq. 2) then          
      ! % https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Units=SI&Mask=1#Thermo-Gas 
       
        Mw = 0.016043d0
       
        Cpw = -0.703029d0 + 108.4773d0*
     &   (T/1000d0) -42.52157d0*(T/1000d0)**2d0 
     &   + 5.862788d0*(T/1000d0)**3d0 + 0.678565d0/(T/1000d0)**2d0
       
        Cpw = Cpw/Mw ! to convert from J/mol/K to J/kg/K
      
      end if

c     Specific heat of moist air

      Cp = (1d0 - q) * Cpd + q * Cpw

      return
      end

c***********************************************************************
      subroutine es_sub_new(es,volatile,T)
c***********************************************************************
      implicit none
      double precision es,T
      double precision b_0,b_1,b_2,b_3,a_0,a_1,a_2,a_3,a_4
      integer volatile

c Saturation vapor pressure (es) of a volatile over the liquid and solid forms
c of the volatile.

c Notation:
c INPUTS:
c     volatile  = 'H2O' (1) or 'CH4' (2)
c     T              = temperature (K)

c OUTPUTS:
c     es             = saturation vapor pressure (Pa)



      if (volatile .eq. 1) then ! H2O
 

c Murphy & Koop (2005) Water Ice Coefficients
          b_0 = 9.550426d0
          b_1 = -5723.265d0
          b_2 = 3.53068d0
          b_3 = -0.00728332d0
c Sonntag (1990) Liquid Water Coefficients
          a_0 = 16.635764d0
          a_1 = -6096.9385d0
          a_2 = -2.711193d-2
          a_3 = 1.673952d-5
          a_4 = 2.433502d0


          
          if (T .le. 273.16d0) then

c Murphy & Koop (2005) formula over water ice (in Pa)

          es = exp(b_0+b_1/T+b_2*log(T)+b_3*T)

          
          else

c Sonntag (1990) formula over liquid water (in Pa)
   
          es = 100d0*exp(a_0 + a_1/T +a_2*T + a_3*T*T + a_4*log(T))

          end if

      else if (volatile .eq. 2) then ! CH4


          if (T .le. 90.66d0) then

c formula over methane ice

          es = 10d0**(7.69540d0 - 532.20d0/(T+1.842d0))*133.322368d0

          else

c formula over liquid methane

          es = 10d0**(6.61184d0 - 389.93d0/(T-7.16d0))*133.322368d0

          end if

      end if

      return
      end
c***********************************************************************
      subroutine Fhw_sub_update(Fh,Fw,
     &      Fx_opt,d0,zref,z,zeta)
c***********************************************************************

      double precision Fh,Fw,d0,z,zref,zeta,Zr,zetar
      double precision vk,beta_h,gamma_h,gamma_c,Pr0,a_h,b_h
      double precision Fh0,Fh_kansas,y,Psi_zeta,Psi_zetar
      double precision PI,lamb,lambr,y_2,Fhc
      integer Fx_opt
c 9/23/22: Updated code to include blend of Kansas-type function and free convective function
c so that zeta limit is extended to -infinity. alpha is assumed to be 1.
c 2/6/23 Changed gamma_h value to 13.3 from Dyer & Bradley (1982)
c 2/7/23 Adjusted convectice constant gamma_c to account for change in gamma_h

c Calculate the primitive functions for heat & water-vapor within the 
c surface sublayer, Fh(zref,z) & Fw(zref,z).

c Notation:

c INPUTS:
c     Fx_opt  Fx physics option [1 or 2]                        (scalar)
c     d0      zero-plane displacement                           (scalar)
c     zref    reference height                                  (scalar)
c     z       heights                                           (vector)
c     zeta    Monin-Obukov stability values for heights z       (vector)

c OUTPUTS:
c     Fh      primitive function Fh(zref,z)                     (vector)
c     Fw      primitive function Fw(zref,z)                     (vector)

      PI = 4.D0*DATAN(1.D0)

c critical parameters

      vk      = 0.40d0         ! Hogstrom (1996)
      beta_h  = 5.4d0          ! Cheng & Brutsaert (2005)
      gamma_h = 16d0           ! Brutsaert (1982)
      gamma_c = 34.15d0        ! Grachev et al. (2000)
      Pr0 = 0.98d0             ! Turbulent Prandtl number
      a_h = 5d0                ! Gryanik et al. (2020) 
      b_h = 0.4d0              ! Gryanik et al. (2020) 

c convert to the (Z = z-d0) coordinate system
      
      z = z - d0
      Zr = zref - d0

c find zeta at the reference height

      zetar = (Zr / z) * zeta

c find Fh(Zr,Z)
 
      Fh0 = log(z / Zr)

      if (zeta .gt. 0d0) then             ! stable conditions

            if (Fx_opt .eq. 1) then
      
                  Fh = Pr0*(Fh0 + beta_h * (zeta - zetar))
   
            else if (Fx_opt .eq. 2) then
                  
                  Fh = Pr0*Fh0 + (Pr0*a_h/b_h)*
     &             ( log(1d0+b_h*zeta) - log(1d0+b_h*zetar))


            end if

      else if (zeta .eq. 0d0) then          ! neutral conditions

            Fh = Pr0*Fh0

      else if (zeta .lt. 0d0) then          ! unstable conditions

         lambr = sqrt(1d0 - gamma_h * zetar)
         lamb  = sqrt(1d0 - gamma_h * zeta)
        
         Fh_kansas   = Fh0 + 2d0 * log((lambr + 1d0) / (lamb + 1d0))
         
         y = (1d0-gamma_c*zetar)**(1d0/3d0)
         
         y_2 =(1d0-gamma_c*zeta)**(1d0/3d0)
         
         Psi_zetar = 1.5d0*log((y**2d0 + y + 1d0)/3d0) 
     &    - sqrt(3d0)*atan((2d0*y + 1d0)/sqrt(3d0)) + pi/sqrt(3d0)
        
         Psi_zeta = 1.5d0*log((y_2**2d0 + y_2 + 1d0)/3d0) 
     &    - sqrt(3d0)*atan((2d0*y_2 + 1d0)/sqrt(3d0)) + pi/sqrt(3d0)  
        
         Fhc = Fh0 - Psi_zeta + Psi_zetar
         
         Fh = Pr0*(Fh_kansas + Fhc*zeta**2d0)/(1d0 + zeta**2d0)

       end if

      Fh = Fh/vk

c set Fh = 0 when z <= zref
 

      if (z .le. zref) then
            Fh = 0d0
      end if
      
c set Fw = Fh
      
      Fw = Fh

c Convert back z so it doesn't change the input variable
       z = z + d0
      return
      end
c***********************************************************************
      subroutine Fhw_subi_v2(Fh,Fw,
     &      Re,Pr,Sc)
c***********************************************************************

      implicit none
      double precision Fh,Fw,Re,Pr,Sc
      double precision Cs,Cr,Re1,Re2,Fh1,Fh2,Fw1,Fw2
      double precision x(2),y(2),y2(2),p(2),p2(2)
      double precision log_Re_transition,Re_transition(1866)
      integer i,idx
c Calculate the primitive functions for heat & water-vapor for the
c interfacial sublayer, Fh(0,hd) & Fw(0,hd).

      Cs = 1d0 / 13.6d0                ! Brutsaert (1975)
      Cr = 1d0 / 7.3d0                 ! Brutsaert (1975)


      if (Re .le. 0.135d0) then    ! aero smooth

            Fh = 1d0/Cs * Pr**(2d0/3d0)
            Fw = 1d0/Cs * Sc**(2d0/3d0)

      else if (Re .gt. 0.135d0 .and. Re .lt. 2d0) then ! aero transition

c     Perform log-log interpolation between smooth and rough regimes,
c     as suggested by Andreas (1987)

            Re1 = 0.135d0 
            Re2 = 2d0
            Fh1 = 1d0/Cs * Pr**(2d0/3d0)  
            Fh2 = 1d0/Cr * Pr**(1d0/2d0) * Re2**(1d0/4d0)
            Fw1 = 1d0/Cs * Sc**(2d0/3d0)
            Fw2 = 1d0/Cr * Sc**(1d0/2d0) * Re2**(1d0/4d0)
            
            x(1)  = Re1 
            x(2)  = Re2 
            y(1)  = Fh1
            y(2)  = Fh2
            y2(1) = Fw1
            y2(2) = Fw2


            do 101 i = 1,1866

                  Re_transition(i) = Re1 + (i-1)*0.001d0

101         continue


            do 123 i = 1,1866

                  if (Re_transition(i) .ge. Re) then
                  idx = i
                  exit
                  end if
123         continue

            
      CALL POLYFIT(LOG(x(1)), LOG(x(2)), LOG(y(1)), LOG(y(2)), p)
      CALL POLYFIT(LOG(x(1)), LOG(x(2)), LOG(y2(1)), LOG(y2(2)), p2)

            log_Re_transition = LOG(Re_transition(idx))

            Fh = EXP(p(1) * log_Re_transition + p(2))
            Fw = EXP(p2(1) * log_Re_transition + p2(2))


      else if (Re .gt. 2d0) then         ! aero rough

            Fh = 1d0/Cr * Pr**(1d0/2d0) * Re**(1d0/4d0)
            Fw = 1d0/Cr * Sc**(1d0/2d0) * Re**(1d0/4d0)

      end if

      return
      end


!_______________________________________________________________________
      SUBROUTINE POLYFIT(x1, x2, y1, y2, p)
      IMPLICIT NONE
      double precision x1, x2, y1, y2, p(2)
      double precision a, b

      a = (y2 - y1) / (x2 - x1)
      b = y1 - a * x1

      p(1) = a
      p(2) = b

      RETURN
      END
c***********************************************************************
      subroutine Fm_sub_update(Fm,
     &      Fx_opt,z0m,d0,z,zeta)
c***********************************************************************

      implicit none
      double precision Fm,z0m,d0,z,zeta,vk,zeta0m,yi,yi0m
      double precision Psi_zetam,Psi_zeta0m,xi0m,xi,t1,t2,Fmc
      double precision beta_m,gamma_m,gamma_c,a_m,b_m,Fm0,Fm_kansas
      double precision PI
      integer Fx_opt

c 9/23/22: Updated code to include blend of Kansas-type function and free convective function
c so that zeta limit is extended to -infinity.

c 2/6/23 Updated gamma_m to 28, based on Dyer & Bradley (1982) after Gary suggested it 
c - seems to reduce iterations by ~400!! Thanks Gary!

c 2/7/23 Adjusted convective constant gamma_c to account for change in gamma_m

c Calculate the primitive function for momentum, Fm(z0m,z).

c Notation:

c INPUTS:
c	Fx_opt  Fx physics option [1 or 2]			  (scalar)
c	z0m     momentum roughness length             (scalar)
c	d0      zero-plane displacement               (scalar)
c	z       heights                               (vector)
c	zeta    Monin-Obukhov stability values        (vector)

c OUTPUTS:
c	Fm      primitive function Fm(z0m,z)          (vector)

c 	Notes:

c (1) In it's current form, this function assumes that zeta corresponds to
c Monin-Obukhov stability parameter at the heights z.  Thus, the value
c zeta(i) corresponds to the zeta value at height z(i).

c (2) Fx_opt = 1	Cheng & Brutsaert (2005) for stable conditions
c            = 2	Gryanik et al (2020)           "
c ________________________________________________

      
      
      PI = 4.D0*DATAN(1.D0)

c critical parameters


      vk      = 0.40d0         ! Hogstrom (1996)
      beta_m  = 5.8d0          ! Cheng & Brutsaert (2005)
      gamma_m = 16d0           ! Brutsaert (1982)
      gamma_c = 10.15d0        ! Grachev et al. (2000)
      a_m = 5d0                ! Gryanik et al. (2020)
      b_m = 0.3d0              ! Gryanik et al. (2020)

c convert to the (Z = z-d0) coordinate system
            
      z = z - d0

c define zeta0m = z0m/L (z0m is already in the z-d0 coordinate system)

      zeta0m = (z0m / z) * zeta

c find Fm

      Fm0 = log(z / z0m)

      if (zeta .gt. 0d0) then             ! stable conditions

        if (Fx_opt .eq. 1) then
        
            Fm = Fm0 + beta_m * (zeta - zeta0m)
        
        else if (Fx_opt .eq. 2) then
            Fm = Fm0 + (3d0*a_m/b_m)*( (1d0 + b_m*zeta)**(1d0/3d0)
     &        - (1d0 + b_m*zeta0m)**(1d0/3d0) )
        end if

      else if (zeta .eq. 0d0) then       ! neutral conditions


        Fm = Fm0

      else if (zeta .lt. 0d0) then       ! unstable conditions

        xi0m = (1d0 - gamma_m * zeta0m)**0.25d0
      
        xi   = (1d0 - gamma_m * zeta)**0.25d0
      
        t1   = log((xi0m**2d0 + 1d0)*(xi0m + 1d0)**2d0
     &   / ((xi**2d0 + 1d0) * (xi + 1d0)**2d0))
      
        t2   = 2d0 * (atan(xi)- atan(xi0m))
      
        Fm_kansas   = Fm0 + t1 + t2
      
        yi0m = (1d0 - gamma_c * zeta0m)**(1d0/3d0)
      
        yi   = (1d0 - gamma_c * zeta)**(1d0/3d0)
      
        Psi_zeta0m = 1.5d0*log((yi0m**2d0 + yi0m + 1d0)/3d0) 
     &   - sqrt(3d0)*atan((2d0*yi0m + 1d0)/sqrt(3d0)) + PI/sqrt(3d0)
      
        Psi_zetam  = 1.5d0*log((yi**2d0 + yi + 1d0)/3d0)
     &    - sqrt(3d0)*atan((2d0*yi + 1d0)/sqrt(3d0)) + PI/sqrt(3d0)
      
        Fmc = Fm0 - Psi_zetam + Psi_zeta0m
    
        Fm = Fm_kansas/(1d0+zeta**2d0) + Fmc*zeta**2d0/(1d0+zeta**2d0)

      end if
      
      Fm = Fm / vk


c set Fm = 0 when z <= z0m

       if (z .le. z0m) then
       Fm = 0d0
       end if

c Convert back z so it doesn't change the input variable
       z = z + d0
       return
       end
c***********************************************************************
      subroutine Ls_sub_sntherm(Ls,T,dliqvol)         
c***********************************************************************

      implicit none
      double precision Ls,T,dliqvol

c     Find the latent heat of water.
c     (see: Fleagle & Businger, Intro to Atmos Physics)

c Notation:
c     INPUTS:
c     T       = temperature (K)
c     dliqvol = volume fraction of liquid water
c     OUTPUTS:
c     Ls      = latent heat of sublimation/evaporation

c     For wet snow (dliqvol > 0.02), use latent heat of evaporation
      
      if (T .lt. 273.16d0 .or. dliqvol .lt. 0.020d0) then

c     Sublimation

            Ls = (28.34d0 - 0.00149d0*(T - 273.15d0))*1.d05    
    
      else

c     Evaporation            

            Ls = (25.00d0 - 0.02274d0*(T - 273.15d0))*1.d05

      end if

      return
      end
c***********************************************************************
      subroutine mol_diffus_v4_sntherm(Pr,Sc,nu,
     & planet,T,rho,P,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)
c***********************************************************************
c CALLED FROM: 
c CALLS:       congas
      implicit none

      double precision md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b
      double precision md_a0,md_a1,md_a2,md_a3,md_a4,md_b0,md_b1,md_b2
      double precision md_b3,md_b4,md_beta_a,md_beta_b,md_gamma_a
      double precision md_gamma_b,md_kappa_a,md_kappa_b,md_m_ab
      double precision md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b
      double precision md_sum_v_a,md_sum_v_b
      double precision md_m_a

      double precision Pr,Sc,nu,T,P,rho(2),R_gas,mu_a,mu_b
      double precision moles,D,mumix,cp_a,cp_b,cv_a,cv_b,cpmix
      double precision y_a,y_b
      double precision Tr_a,Tr_b,lambda_tr_a,lambda_tr_b,lambda_tr_ab
      double precision lambda_tr_ba,epsilon,A_ab,A_ba,K_a,K_b,Kmix
      integer          planet


c INPUTS:
c    T      = temperature (K)                 (scalar)
c    rho    = gas density of each species     (vector)
c    mu     = dynamic viscosity
c    K      = thermal conductivity
c    Cp     = specific heat                   (J/kg/K)
c    D12    = binary diffusion coefficient    (m^2/s)

c OUTPUTS:
c    Pr     = Prandtl number                  (unitless)
c    Sc     = Schmidt number                  (unitless)
c    nu     = kinematic viscosity             (m^2/s)



      R_gas = 8.314462618d0     ! Gas constant, J/K/mol
      P     = P*1.0d-5        ! Convert from Pa to bar


      
      if (planet .eq. 1) then ! Earth
            
            moles = rho(1) / (md_m_a*0.001d0) + rho(2) 
     &       / (md_m_b*0.001d0)
            y_a  = (rho(1) / (md_m_a*0.001d0)) /moles ! mole fraction of A
            y_b  = (rho(2) / (md_m_b*0.001d0)) /moles ! mole fraction of B


            D = (0.00143d0*T**(1.75d0)) / ( P*sqrt(md_m_ab)* 
     &      (md_sum_v_a**(1d0/3d0) + md_sum_v_b**(1d0/3d0))**2d0 )! Eq. 11-4.4 (P must be in bars)
        
            D = D*0.0001d0 ! convert from cm^/s to m^2/s

c     Individual dynamic viscosities in kg/(m s) or N s/m^2 by setting mole fractions to 1 and 0, respectively
      

        !write(*,*) 'moles: ', moles
        !write(*,*) 'md_kappa_a: ', md_kappa_a

            CALL  viscosity_mix(mu_a,
     &  md_kappa_a,md_kappa_b,md_m_a,md_m_b,md_mu_a,md_mu_b,md_omega_a,
     &  md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a,md_Vc_b,1d0,0d0)

        
            CALL  viscosity_mix(mu_b,
     &  md_kappa_a,md_kappa_b,md_m_a,md_m_b,md_mu_a,md_mu_b,md_omega_a,
     &  md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a,md_Vc_b,0d0,1d0)

        
c     Eq. 9-5.24: Mixture dynamic viscosity in kg/(m s) or N s/m^2
        
              CALL  viscosity_mix(mumix,
     &  md_kappa_a,md_kappa_b,md_m_a,md_m_b,mu_a,mu_b,md_omega_a,
     &  md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a,md_Vc_b,y_a,y_b)
      
        
        
            cp_a = 0.001d0*md_m_a*(md_a0 + md_a1*T + md_a2*T**2d0 + 
     &      md_a3*T**3d0 + md_a4*T**4d0) ! Specific heat at constant pressure, J/ (mol K)
        
            cp_b = R_gas*(md_b0 + md_b1*T + md_b2*T**2d0 + 
     &       md_b3*T**3d0 + md_b4*T**4d0)


      else if (planet .eq. 2) then ! Mars

            moles = rho(1) / (md_m_a*0.001d0) + rho(2) /(md_m_b*0.001d0)
            y_a  = (rho(1) / (md_m_a*0.001d0)) /moles !  mole fraction of A
            y_b  = (rho(2) / (md_m_b*0.001d0)) /moles !  mole fraction of B

            D = (0.00143d0*T**(1.75d0)) / ( P*sqrt(md_m_ab)* 
     &            (md_sum_v_a**(1d0/3d0) 
     &            + md_sum_v_b**(1d0/3d0))**2d0 ) ! Eq. 11-4.4 (P must be in bars)
        
            D = D*0.0001d0 ! convert from cm^/s to m^2/s

c Individual dynamic viscosities in kg/(m s) or N s/m^2 by setting mole fractions to 1 and 0, respectively
            CALL viscosity_mix(mu_a,md_kappa_a,md_kappa_b,md_m_a,md_m_b,
     &  md_mu_a,md_mu_b,md_omega_a,md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a
     &  ,md_Vc_b,1d0,0d0)
      
            CALL viscosity_mix(mu_b,md_kappa_a,md_kappa_b,md_m_a,md_m_b,
     &  md_mu_a,md_mu_b,md_omega_a,md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a
     &  ,md_Vc_b,0d0,1d0)

c Eq. 9-5.24: Mixture dynamic viscosity in kg/(m s) or N s/m^2
            CALL viscosity_mix(mumix,
     &  md_kappa_a,md_kappa_b,md_m_a,md_m_b,mu_a,mu_b,
     &  md_omega_a,md_omega_b,md_Tc_a,md_Tc_b,T,md_Vc_a,md_Vc_b,y_a,y_b)


            cp_a = R_gas*(md_a0 + md_a1*T + md_a2*T**2d0 + md_a3*T**3d0 
     &      + md_a4*T**4d0) ! Specific heat at constant pressure, J/ (mol K)
        
            cp_b = R_gas*(md_b0 + md_b1*T + md_b2*T**2d0 + md_b3*T**3d0 
     &      + md_b4*T**4d0)

      else if (planet .eq. 3) then ! Titan

        !load('N2_CH4.mat');
      write(*,*)'mol_diffus_v4_sntherm: Titan has not been set up yet'

      end if


      cpmix = cp_a*y_a + cp_b*y_b ! Mixture specific heat, J/(mol K)
      cpmix = 1000d0*cpmix/(y_a*md_m_a + y_b*md_m_b)

      cv_a = cp_a - R_gas ! Specific heat at constant volume, J/(mol K)
      cv_b = cp_b - R_gas

      Tr_a = T/md_Tc_a !Reduced temperature
      Tr_b = T/md_Tc_b

      lambda_tr_a = (md_gamma_b* ( exp(0.0464d0*Tr_a) 
     &       - exp(-0.2412d0*Tr_a) )) ! Translational thermal conductivity

      lambda_tr_b = (md_gamma_a* 
     &   ( exp(0.0464d0*Tr_b) - exp(-0.2412d0*Tr_b)))

      lambda_tr_ba = lambda_tr_b/lambda_tr_a ! Ratio, Eq. 10-6.5, pg. 595
      
      lambda_tr_ab = lambda_tr_a/lambda_tr_b ! Ratio, Eq. 10-6.5, pg. 595


      epsilon = 1d0 ! Always taken as one, pg. 595


      A_ab = epsilon * (1d0 + sqrt(lambda_tr_ab) 
     & * (md_m_a/md_m_b)**(1d0/4d0) )**2d0
     & / sqrt(8d0*(1d0+md_m_a/md_m_b)) ! Eq. 10-6.2, pg. 595
      
      A_ba = epsilon * (1d0 + sqrt(lambda_tr_ba)
     &  * (md_m_b/md_m_a)**(1d0/4d0) )**2d0
     &  / sqrt(8d0*(1d0 + md_m_b/md_m_a))

      CALL condgas(K_a,
     &       md_beta_a,cv_a,md_m_a,T,md_Tc_a,mu_a) ! Thermal conductivity, W/(m K)

      CALL condgas(K_b,
     &       md_beta_b,cv_b,md_m_b,T,md_Tc_b,mu_b)

      Kmix = (y_a* K_a/(y_a + y_b*A_ab) + y_b*K_b/(y_a*A_ba + y_b) )! Mixture thermal conductivity, W/(m K), Eq. 10-6.1, pg. 594


c find Prandtl, Schmidt numbers, and kinematic viscosity

      
      Pr = (mumix / Kmix) * cpmix
      Sc = mumix / (sum(rho)*D) ! sum may not work
      nu = mumix / sum(rho)

      P     = P*100000d0        ! Convert back from bar to Pa

      
      return
      end
c***********************************************************************
      subroutine mol_setup(planet,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b,dryair,volatile,Ma,Mw)
c***********************************************************************

      implicit none
      
      double precision md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b
      double precision md_a0,md_a1,md_a2,md_a3,md_a4,md_b0,md_b1,md_b2
      double precision md_b3,md_b4,md_beta_a,md_beta_b,md_gamma_a
      double precision md_gamma_b,md_kappa_a,md_kappa_b, md_m_a,md_m_ab
      double precision md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b
      double precision md_sum_v_a,md_sum_v_b
      double precision Ma,Mw
      integer planet,dryair,volatile
      
      
c The following parameters are used to calculate the molecular properties
c of the gases involved
      
      if (planet .eq. 1) then ! Earth
      
c  %%% Air - H2O Properties for Molecular Diffusion %%%
c  %%% Data from Properties of Gases and Liquids, Poling et al., 2001 %%%
      
      ! a = Air
      ! b = H2O
      
          dryair   = 1       ! predominant dry air component (air)
          volatile = 1       ! predominant volatile species  (H2O)

          ! Appendix A, Critical Pressure, bar
          md_Pc_a   =   37.859999999999999d0
          md_Pc_b   =   2.206400000000000d2
          
          ! Critical temperature, K, from Lemmon et al. (2000), https://srd.nist.gov/jpcrdreprint/1.1285884.pdf
          
          md_Tc_a   =   1.325306000000000d2
          md_Tc_b   =   6.471400000000000d2
          
          ! Critical volume, cm^3/mol, https://www.ohio.edu/mechanical/thermo/property_tables/gas/Critical.html
          
          md_Vc_a   =   88.299999999999997d0
          md_Vc_b   =   55.950000000000003d0
          
          ! Air Specific Heat Constants
          md_a0     =   1.057500000000000d3
          md_a1     =   -0.448900000000000d0
          md_a2     =    0.001140700000000d0
          md_a3     =   -7.999900000000000d-7
          md_a4     =  1.932700000000000d-10
          
          ! H2O Specific Heat Constants, Appendix A, Section C
          
          md_b0     =     4.395000000000000d0
          md_b1     =    -0.004186000000000d0
          md_b2     =   1.405000000000000d-5
          md_b3     =  -1.564000000000000d-8
          md_b4     =  6.320000000000000d-12
          
          md_beta_a =     0.786200000000000d0 ! for non-polar, pg. 576
          md_beta_b =     0.780000000000000d0 ! Table III of Chung et al. (1984)
          
          ! Reduced inverse thermal conductivity [W/(m K)]^-1, Eq. 10-3.9, pg. 570
          
          md_gamma_a = 2.263400787457794d2
          md_gamma_b =  71.793205483616532d0
          
          ! Association factor, Table 9-1, pg. 474
          
          md_kappa_a =   0.0000d0
          md_kappa_b =    0.076000000000000d0
          
          ! Molar mass, g/mol
          
          md_m_a     =  28.965859999999999d0
          md_m_b     =  18.015280000000001d0
          md_m_ab    =  22.214364246623219d0 ! m_ab = 2 / [1/m_a + 1/m_b]
          
          ! Appendix A, dipole moment, debyes
          
          md_mu_a    =   0.0000d0
          md_mu_b    =   1.8000d0
          
          ! Appendix A, Pitzer acentric factor, acentric factors from perry's handbook for chemical engineers
          
          md_omega_a =   0.0000d0
          md_omega_b =   0.344000000000000d0
          
          ! Diffusion volume from Table 11.1, pg. 645
          
          md_sum_v_a =  19.699999999999999d0
          md_sum_v_b =  13.100000000000000d0

      else if (planet .eq. 2) then ! Mars
      
c %%% CO2 - H2O Properties for Molecular Diffusion %%%
c %%% Data from Properties of Gases and Liquids, Poling et al., 2001 %%%
      
      ! a = CO2
      ! b = H2O
           
          dryair   = 2       ! predominant dry air component (CO2)
          volatile = 1       ! predominant volatile species  (H2O) 
          
          ! Appendix A, Critical Pressure, bar
          
          md_Pc_a   =   73.739999999999995d0
          md_Pc_b   =   2.206400000000000d2
          
          ! Appendix A, Critical temperature, K
          
          md_Tc_a   =   3.041200000000000d2 
          md_Tc_b   =   6.471400000000000d2
          
          ! Critical volume, cm^3/mol, https://www.ohio.edu/mechanical/thermo/property_tables/gas/Critical.html! 
         
          md_Vc_a   =   94.069999999999993d0
          md_Vc_b   =   55.950000000000003d0
          
          ! CO2 Specific Heat Constants, coefficients from Table 12 of Lemmon et al. (2000)
          
          md_a0     =   3.259000000000000d0
          md_a1     =   0.001356000000000d0
          md_a2     =   1.502000000000000d-5
          md_a3     =   -2.374000000000000d-8
          md_a4     =   1.056000000000000d-11
          
          ! H2O Specific Heat Constants, Appendix A, Section C
          
          md_b0     =     4.395000000000000d0
          md_b1     =    -0.004186000000000d0
          md_b2     =   1.405000000000000d-5
          md_b3     =  -1.564000000000000d-8
          md_b4     =  6.320000000000000d-12
          
          md_beta_a =     0.692910500000000d0 ! for non-polar, pg. 576
          md_beta_b =     0.780000000000000d0 ! Table III of Chung et al. (1984)
          
          ! Reduced inverse thermal conductivity [W/(m K)]^-1, Eq. 10-3.9, pg. 570
          
          md_gamma_a =  2.054472758857511d2   
          md_gamma_b =  71.793205483616532d0
          
          ! Association factors, Table 9-1, pg. 474
         
          md_kappa_a =   0.0000d0             
          md_kappa_b =    0.076000000000000d0
          
          ! Molar mass, g/mol
       
          md_m_a     =  44.009999999999998d0
          md_m_b     =  18.015280000000001d0
          md_m_ab    =   25.565462108353238d0 ! 2 / [1/m_a + 1/m_b]
          
          ! Appendix A, dipole moment, debyes
          
          md_mu_a    =   0.0000d0             
          md_mu_b    =   1.8000d0             
          
          ! Appendix A, Pitzer acentric factor, acentric factors from perry's handbook for chemical engineers
          
          md_omega_a =   0.225000000000000d0 
          md_omega_b =   0.344000000000000d0  
          
          ! Diffusion volume from Table 11.1, pg. 645
          
          md_sum_v_a =  26.899999999999999d0  
          md_sum_v_b =  13.100000000000000d0  
      
      end if

c     Set molecular weight of dry air (kg/mol)

      if (dryair .eq. 1) then
          Ma = 0.028966d0  ! modern earth
      else if (dryair .eq. 2) then
          Ma = 0.0440098d0 ! mars
      else if (dryair .eq. 3) then
          Ma = 0.0280134d0 ! titan
      end if

c     Set molecular weight of volatile (kg/mol)

      if (volatile .eq. 1) then
          Mw = 0.01801528d0
      else if (volatile .eq. 2) then
          Mw = 0.016043d0
      end if
      
      
      return
      end

c***********************************************************************
      subroutine q_sub_sntherm(qo,Ma,Mw,P,T,Rh,es)
c***********************************************************************


      implicit none
      double precision rhoa,rhow,rho,qo
      double precision P,T,Rh
      double precision R,Ma,Mw,e,es
      integer boil

c     Find specific humidity

c INPUTS:
c     Ma       = molecular weight of dry air ('air','CO2','N2'; kg/mol)
c     MW       = molecular weight of volatile ('H2O','CH4'; kg/mol)
c     P        = pressure (Pa)
c     T        = temperature (K)
c     Rh       = relative humidity (0 - 1)

c OUTPUTS:
c     qo        = specific humidity (unitless)



! boolean variable indicating whether boiling is occurring
      boil = 0 

      R  = 8.31434d0  ! gas constant

c find water vapor pressure, e

       e  = Rh * es*100d0  ! to convert to Pa
      
      if (e > P) then

        write(*,*) 'rho_sub -- Boiling, for e: ',e, ', P: ',P,', P: ',T
    
        rhoa = 0d0
        rhow = 0d0
        boil = 1 ! boiling is occurring
        return
    
      end if

c find densities

      rhow = Mw * e / (R*T)
      rhoa = Ma * (P - e) / (R*T)
      rho  = rhoa + rhow
      qo    = rhow / rho
      
      return
      end
c***********************************************************************
      subroutine qsat_derivative(dqsat,Ma,Mw,P,T)
c***********************************************************************


      implicit none
      double precision dqsat,P,T,Ma,Mw,num_1,num_2,den_1,den_2
      double precision b_0,b_1,b_2,b_3,a_0,a_1,a_2,a_3,a_4
      double precision eso, es_derivative
      

c INPUTS:
c     Ma       = molecular weight of dry air ('air','CO2','N2'; kg/mol)
c     MW       = molecular weight of volatile ('H2O','CH4'; kg/mol)
c     P        = pressure (Pa)
c     T        = temperature (K)

c OUTPUTS:
c     dTqsat   = d(q_sat)/dT 

            
c Murphy & Koop (2005) Water Ice Coefficients
      b_0 = 9.550426d0
      b_1 = -5723.265d0
      b_2 = 3.53068d0
      b_3 = -0.00728332d0
c Sonntag (1990) Liquid Water Coefficients
      a_0 = 16.635764d0
      a_1 = -6096.9385d0
      a_2 = -2.711193d-2
      a_3 = 1.673952d-5
      a_4 = 2.433502d0
          
c Equilibrium (saturation) vapor pressure over water ice or liquid
           
      if (T .le. 273.16d0) then

c Murphy & Koop (2005) formula over water ice (in Pa)

          eso = exp(b_0 + b_1/T + b_2*log(T) + b_3*T)

          es_derivative = eso* (-b_1/T**2d0 + b_2/T + b_3)
          
      else

c Sonntag (1990) formula over liquid water (in Pa)
   
          eso = 100d0*exp(a_0 + a_1/T +a_2*T + a_3*T*T + a_4*log(T))
          
          es_derivative = eso* (- a_1/T**2d0 + a_2 + 2d0*a_3*T + a_4/T)

      end if
          
          

      num_1 = Mw*es_derivative

      den_1  = Ma*(P - eso)

      num_2 = Mw*Ma*eso*es_derivative

      den_2 =  (Ma*(P - eso))**2d0

      dqsat = num_1/den_1 + num_2/den_2

     
      return
      end
c***********************************************************************
      subroutine Rh_sub(e,Rh,
     & dryair,volatile,P,T,q)
c***********************************************************************

      implicit none
      double precision e,es,Rh,Ma,Mw,P,T,q,fac

      integer dryair,volatile

c     Find vapor pressure and relative humidity from the specific humidity.

c     Notation:

c     INPUTS:
c     dryair   = predominant content of dry air ('air','CO2','N2')
c     volatile = predominant volatilve('H2O','CH4')
c     P        = pressure (Pa)
c     T        = temperature (K)
c     q        = specific humidity

c     OUTPUTS:
c     e        = water vapor pressure (Pa)
c     Rh       = relative humidity (0 - 1)

c     Set molecular weight of dry air

      if (dryair .eq. 1) then
            Ma = 0.028966d0        ! Modern Earth
      else if (dryair .eq. 2) then
            Ma = 0.0440098d0       ! Mars
      else if (dryair .eq. 3) then
            Ma = 0.0280d0          ! Titan
      end if

c     Set molecular weight of volatile

      if (volatile .eq. 1) then
            Mw = 0.018016d0
      else if (volatile .eq. 2) then
            Mw = 0.016043d0
      end if

c     Saturation vapor pressure
     
      CALL es_sub_new(es,
     &       volatile,T)
 

c     Vapor pressure of volatile, e

      fac = 1d0 + (Mw / Ma) * (1d0/ q - 1d0)
      e   = P / fac

c     Relative humidity

      Rh = e / es

c     Reality check

      if (Rh .gt. 1d0) then
            Rh = 1d0
      end if
      
      return
      end
c***********************************************************************
      subroutine rho_sub_new(rhoa,rhow,rho,q,boil,
     &                       dryair,volatile,P,T,Rh)
c***********************************************************************


      implicit none
      double precision rhoa,rhow,rho,q
      double precision P,T,Rh
      double precision R,Ma,Mw,e,es
      integer boil,dryair,volatile

c     Find gas densities

c INPUTS:
c     dryair   = predominant content of dry air ('air','CO2','N2')
c     volatile = predominant volatile('H2O','CH4')
c     P        = pressure (Pa)
c     T        = temperature (K)
c     Rh       = relative humidity (0 - 1)

c OUTPUTS:
c     rhoa     = density of dry air (kg/m^3)
c     rhow     = water vapor density (kg/m^3)
c     rho      = total density of moist air (kg/m^3)
c     q        = specific humidity (unitless)



! boolean variable indicating whether boiling is occurring
      boil = 0 

! gas constant
      R  =  8.31446261815324d0

c     Set molecular weight of dry air

      if (dryair .eq. 1) then
        Ma = 0.028966d0  ! modern earth
      else if (dryair .eq. 2) then
        Ma = 0.0440098d0 ! mars
      else if (dryair .eq. 3) then
        Ma = 0.0280134d0 ! titan
      end if

c     Set molecular weight of volatile

      if (volatile .eq. 1) then
        Mw = 0.01801528d0
      else if (volatile .eq. 2) then
        Mw = 0.016043d0
      end if

c find water vapor pressure, e

      CALL es_sub_new(es,
     &       volatile,T)
      e  = Rh * es
      
      
      if (e > P) then

        write(*,*) 'rho_sub -- Boiling, for e: ',e, ', P: ',P,', P: ',T
    
        rhoa = 0d0
        rhow = 0d0
        boil = 1 ! boiling is occurring
        return
    
      end if

c find densities

      rhow = Mw * e / (R*T)
      rhoa = Ma * (P - e) / (R*T)
      rho  = rhoa + rhow
      q    = rhow / rho

      
      return
      end
c*********************************************************************** 
      SUBROUTINE SETUP_GRID(hd,Zmax,Z,N)
c***********************************************************************
      IMPLICIT NONE
      double precision hd, Zmax
      double precision Z(200)
      INTEGER N, i, nZ, nZ1, nZ2, nZ3, nZ4, nZ5
      double precision Z1grid(2), Z2grid(20), Z3grid(90)
      double precision Z4grid(10), Z5grid(10)
      double precision dZ

      ! Initialize Z array
      nZ = 0
      do i =1,200
            Z(i) = 0d0
      end do

      ! First segment
      Z1grid(1) = 0d0
      Z1grid(2) = hd
      nZ1 = 2

      if (hd .le. 1d0) then
            dZ = (1d0 - hd) / 4d0
            nZ2 = INT((1d0 - (hd + dZ)) / dZ + 1d0)
            do i = 1, nZ2
                  Z2grid(i) = hd + dZ + (i-1) * dZ
            end do
            nZ3 = 90
            do i = 1, nZ3
                  Z3grid(i) = 1d0 + 0.1d0 * i
            end do
      else if (hd .gt. 1d0 .AND. hd .le. 3d0) then
            dZ = (3d0 - hd) / 4d0
            nZ2 = INT((3d0 - (hd + dZ)) / dZ + 1d0)
            do i = 1, nZ2
                  Z2grid(i) = hd + dZ + (i-1) * dZ
            end do
            nZ3 = 70
            do i = 1, nZ3
                  Z3grid(i) = 3d0 + 0.1d0 * i
            end do
      else
            dZ = (7d0 - hd) / 4d0
            nZ2 = INT((7d0 - (hd + dZ)) / dZ + 1d0)
            do i = 1, nZ2
                  Z2grid(i) = hd + dZ + (i-1) * dZ
            end do
            nZ3 = 30
            do i = 1, nZ3
                  Z3grid(i) = 7d0 + 0.1d0 * i
            end do
      end if

      ! Combine Z1grid, Z2grid, and Z3grid into Z
      do i = 1, nZ1
            nZ = nZ + 1
            Z(nZ) = Z1grid(i)
      end do

      do i = 1, nZ2
            nZ = nZ + 1
            Z(nZ) = Z2grid(i)
      end do

      do i = 1, nZ3
            nZ = nZ + 1
            Z(nZ) = Z3grid(i)
      end do

      if (Zmax .gt. 10d0) then
            nZ4 = 10
            do i = 1, nZ4
                  Z4grid(i) = 10d0 + i
            end do

            do i = 1, nZ4
                  nZ = nZ + 1
                  Z(nZ) = Z4grid(i)
            end do

            if (Zmax .gt. 20d0) then
                  nZ5 = INT((Zmax - 22d0) / 2d0) + 1d0
                  do i = 1, nZ5
                        Z5grid(i) = 20d0 + 2d0 * i
                  end do

                  do i = 1, nZ5
                        nZ = nZ + 1
                        Z(nZ) = Z5grid(i)
                  end do
            end if
      end if

      N = nZ

      RETURN
      END !SUBROUTINE SETUP_GRID  REJ added "!" on 2025/02/24

c***********************************************************************
      subroutine Tqsat_derivative(dTqsat,Ma,Mw,P,T)
c***********************************************************************

      implicit none
      double precision dTqsat,P,T,Ma,Mw
      double precision b_0,b_1,b_2,b_3
      

c INPUTS:
c     Ma       = molecular weight of dry air ('air','CO2','N2'; kg/mol)
c     MW       = molecular weight of volatile ('H2O','CH4'; kg/mol)
c     P        = pressure (Pa)
c     T        = temperature (K)

c OUTPUTS:
c     dTqsat   = d(Tq_sat)/dT

            
c Murphy & Koop (2005) Water Ice Coefficients
          b_0 = 9.550426d0
          b_1 = -5723.265d0
          b_2 = 3.53068d0
          b_3 = -0.00728332d0


      dTqsat = (Mw * (b_2/T - b_1/T**2d0 + b_3) * T * 
     & exp(2d0 * b_2 * log(T) + 2d0 * b_3 * T
     & + (2d0 * b_1) /T + 2d0 * b_0))/
     & (Ma * (P - exp(b_2 * log(T) + b_3*T + b_1/T + b_0))**2d0)
     & + (Mw * (b_2/T - b_1/T**2d0 + b_3) * T * 
     & exp(b_2 * log(T) + b_3 * T + b_1/T + b_0))
     & /(Ma * (P - exp(b_2 * log(T) + b_3 * T + b_1/T + b_0)))
     & +(Mw * exp(b_2 * log(T) + b_3 * T + b_1/T + b_0))
     & /(Ma * (P - exp(b_2 * log(T) + b_3 * T + b_1/T + b_0)))

     
      return
      end
c***********************************************************************
      subroutine viscosity_mix(mumix,
     &   kappa_a,kappa_b,m_a,m_b,mu_a,mu_b,omega_a,
     &   omega_b,Tc_a,Tc_b,T,Vc_a,Vc_b,y_a,y_b)
c***********************************************************************


      implicit none

      double precision mumix,kappa_a,kappa_b,m_a,m_b,mu_a,mu_b,omega_a
      double precision omega_b,Tc_a,Tc_b,T,Vc_a,Vc_b,M_ab
      double precision y_a,y_b
      double precision A,B,C,D,E,F,eta_k_a,eta_k_b,eta_k_ab,eta_k_m
      double precision omega_ab,omega_m,omega_v,sigma_a,sigma_b,sigma_m
      double precision sigma_ab,kappa_m,kappa_ab,M_m,mu_m_4,mu_m
      double precision T_cm,V_cm,mu_rm,F_cm,T_star_m

      sigma_a = 0.809d0*Vc_a**(1d0/3d0) !
      sigma_b = 0.809d0*Vc_b**(1d0/3d0) ! Eq. 9-5.32

      sigma_ab = sqrt(sigma_a*sigma_b) ! % Eq. 9-5.33

      eta_k_a = Tc_a/1.2593d0 ! Eq. 9-5.34
      eta_k_b = Tc_b/1.2593d0 !

      eta_k_ab = sqrt(eta_k_a*eta_k_b) !  Eq. 9-5.35

      omega_ab = (omega_a + omega_b)/2d0 ! % Eq. 9-5.36

      kappa_ab = sqrt(kappa_a*kappa_b) ! % Eq. 9-5.38

      M_ab = 2d0*m_a*m_b/(m_a + m_b)     ! % Eq. 9-5.39

      !write(*,*) 'mu_a: ', mu_a
      
      sigma_m = (y_a)**2d0*(sigma_a)**3d0 + (y_b)**2d0*(sigma_b)**3d0 
     & + 2d0*y_a*y_b*(sigma_ab)**3d0 ! Eq. 9-5.25

      kappa_m = (y_a)**2d0*(kappa_a) + (y_b)**2d0*(kappa_b) 
     & + 2d0*y_a*y_b*(kappa_ab) ! % Eq. 9-5.27

      eta_k_m = ((y_a)**2d0*eta_k_a*(sigma_a)**3d0 + 
     & (y_b)**2d0*eta_k_b*(sigma_b)**3d0 
     & + 2d0*y_a*y_b*eta_k_ab*(sigma_ab)**3d0)/sigma_m ! Eq. 9-5.27


      

      M_m = ( (y_a**2d0 * eta_k_a * sigma_a**2d0 * sqrt(m_a) + 
     &     y_b**2d0 * eta_k_b * sigma_b**2d0 * sqrt(m_b) 
     & + 2d0 * y_a * y_b * eta_k_ab * sigma_ab**2d0 * 
     & sqrt(M_ab))/((eta_k_m)*sigma_m**(2d0/3d0)) )**2d0 ! % Eq. 9-5.40


      omega_m = ( y_a**2d0*omega_a * sigma_a**3d0 
     & + y_b**2d0*omega_b*sigma_b**3d0
     & + 2d0*y_a * y_b * omega_ab * sigma_ab**3d0  )/sigma_m ! Eq. 9-5.29

      mu_m_4 = ( y_a**2d0*mu_a**4d0/sigma_a**3d0  
     &     + y_b**2d0*mu_b**4d0/sigma_b**3d0 
     & +  2d0*y_a*y_b*mu_a**2d0*mu_b**2d0/sigma_ab**3d0  )*sigma_m

      mu_m = mu_m_4**(1d0/4d0) ! Eq. 9-5.30


      V_cm = sigma_m/(0.809d0)**3d0 !  Eq. 9-5.42

      T_cm = 1.2593d0* eta_k_m      ! Eq. 9-5.44

      mu_rm = 131.3d0 * mu_m / sqrt(V_cm*T_cm) ! Eq. 9-5.43

      F_cm = 1d0 - (0.275d0*omega_m) + 0.059035d0*mu_rm**4d0+ kappa_m

      T_star_m = T/eta_k_m


      A = 1.16145d0
      B = 0.14874d0
      C = 0.52487d0
      D = 0.77320d0
      E = 2.16178d0
      F = 2.43787d0

      omega_v = A*T_star_m**(-B) + C*exp(-D*T_star_m) + 
     & E*exp(-F*T_star_m)
      
      mumix = (26.69d0* F_cm *sqrt( M_m*T))/(sigma_m**(2d0/3d0)*omega_v)! micropoises

      mumix = mumix*1.0d-7 ! convert from micropoise to kg/(m s) or N s/m^2

      return
      end


          
