c***********************************************************************
c  SNTHERM24 TRUNCATED VERSION: Descending nodal order truncated
c  after 'END Radiative Solver Adding Doubling Method'
c  Second part of netsolarfluxes.f handled by TOTAL_SOLAR.f 
c  NETSOLARFLUXES models the spectrally integrated solar fluxes
c  based on the Delta-Eddington method modified to account for varying
c  refractive index boundaries and vertically inhomogeneous layers of
c  snow, firn, and glacier ice with impurities such as martian dust.
c  REFERENCES:
c  Joseph, J. H., Wiscombe, W. & Weinman, J. The delta-Eddington 
c  approximation for radiative flux transfer. (1976)
c  Wiscombe, W. J. & Warren, S. G. A model for the spectral albedo of 
c  snow.1. Pure snow. (1980)
c  Briegleb, P. & Light, B. A Delta-Eddington mutiple scattering 
c  parameterization for solar radiation in the sea ice component of the 
c  community climate system model. (2007)
c  Khuller, A. R., Christensen, P. R. & Warren, S. G. Spectral albedo of
c  dusty martian H2O snow and ice. (2021)
c  Whicker, C. A. et al. SNICAR-ADv4: a physically based radiative 
c  transfer model to represent the spectral albedo of glacier ice.(2022)
c  Liou, K. N. An introduction to atmospheric radiation. (2002)
c***********************************************************************
c   IMPORTANT note: routine returns Fabs_int (= dsol) instead of F_abs
cREJ      subroutine netsolarfluxes(Fabs_int,albedo,tau,omega,g,R_sfc,
cREJ_2025/04/26 subroutine netsolarfluxes(dsol,albedo,tau,omega,g,R_sfc,lyr_typ,
      subroutine netsolarfluxes_truncat(dsol,albedo,tau,omega,g,R_sfc,
     & lyr_typ,mu_not,rad,n,albedo_int,iy,jday,ihour,iter,solartest,
     & fdirup,fdifup,fdirdn,fdifdn,dlambda)  
     
c
c Called from Main (sntherm.f) 
c Calls total_solar 
c Note: Some items needed by total_solar.f are merely passed through 
c this routine--i.e., they are input and output, but not used.
c REJ is still working on the documentation
c
c INPUTS:
c albedo        Spectral albedo
c tau           Spectral optical depth of layers
c omega         Spectral single-scattering albedo of layers
c g             Spectral asymmetry parameter of layers
c R_sfc         Spectral albedo of surface underlying snow/ice
c lyr_typ       1 = not specular/Fresnel, 2 = specular/Fresnel
c mu_not        Cosine of zenith angle
c Fs0           Normalized spectral, direct-beam downward flux [W/m^2/band]
c rad           Broadband, downward solar flux [W/m^2] rad(1) from sntherm
c n             number of CVs/nodes
c iy,jday,hour  year(2 digits),julian day, and hour of simulation
c iter          iteration within the hourly subtime division
c solartest     true = solartest comparison, false = reg SNTHERM run
c porosity_rad  Porosity of ice  NOT USED??

c PASSED: through 'arrays_spectral' All dimensioned (nbr_wvl)

c lambda        Wavelength array [microns]
c mImag         Complex part of ice refractive index (spectral)
c n_ice         Real part of ice refractive index (spectral)
c Fd0           Normalized spectral, diffuse downward flux [W/m^2/band]
c FL_r_dif_a    Pre-calculated, diffuse Fresnel reflectivity above layer
c FL_r_dif_b    Pre-calculated, diffuse Fresnel reflectivity below layer

c PARAMETERS in 'const'
c nd = 100, nbr_wvl = 281, ncv = 50, nif = 51   Alotted array sizes
c These must eq/exceed the actual dimension. i.e., ncv must eq/exceed n

c OUTPUTS: (incomplete)
c albedo_int
c fdirup        
c fdifup
c fdirdn 
c fdifdn
c rad           Broadband, downward solar flux [W/m^2] rad(1) from sntherm    
c Fs0           Normalized spectral, direct-beam downward flux [W/m^2/band]
c Fs            Spectral, direct-beam downward flux [W/m^2/band]
c Fd            Spectral, diffuse-beam downward flux [W/m^2/band]


c***********************************************************************
      implicit none
      include 'const' !REJ 4/30/24
      include 'arrays_spectral' !REJ 4/30/24

c  Parameters
c    Sntherm
c      nd      = allotted storage for control vol. variables in 'const'
c    Netsolarfluxes 
c      nbr_wvl = allotted storage for wavelength intervals in 'const'
c      ncv     = allotted storage for control vol. variables in 'const'
c      nif     = allotted storage for interface variables in 'const'

c***********************************************************************      
c     Variable Declaration. Will complete descriptions later
c     ICVP = Inherent Control Volume Property
c     ITFP = Interface Property

      logical found_frsnl,found,check,do_write,solartest
      integer n, kfrsnl,i,j,iw,k,counter,ng,iy,jday,ihour,iter 
      double precision epsilon,exp_min,trmin,refk
      double precision rad,albedo_int,mu_not,PI,critical_angle
      double precision Fsum,Fsum2,dlambda

cREJ_2025/04/26      double precision Fs0    (nbr_wvl)          ! Wavelength
cREJ_2025/04/26      double precision Fs     (nbr_wvl)  
cREJ_2025/04/26      double precision Fd     (nbr_wvl)
      double precision Nreal  (nbr_wvl)
      double precision albedo (nbr_wvl)
      double precision R_sfc  (nbr_wvl)
      double precision temp1  (nbr_wvl) 
      double precision temp2  (nbr_wvl)

      double precision rdir   (nif)              ! Broadband Flux ITFP
      double precision rdif_a (nif)
      double precision rdif_b (nif)
      double precision tdif_a (nif)   
      double precision tdif_b (nif)
      double precision Fnet   (nif)
      double precision Fint   (nif)
      double precision trnlay (nif)

      double precision tdir(ncv)                 ! Broadband ICVP
      double precision Fabs_int(nd)              ! passes to dsol(nd) in sntherm
      double precision dsol(nd)
      integer lyr_typ(ncv),km1

      
      double precision tau   (nbr_wvl,ncv)       ! Spectral ICVP 
      double precision omega (nbr_wvl,ncv)
      double precision g     (nbr_wvl,ncv) 
      double precision tau0  (nbr_wvl,ncv)
      double precision omega0(nbr_wvl,ncv)
      double precision g0    (nbr_wvl,ncv)   

       
     
      double precision trndir (nbr_wvl,nif)      ! Spectral Fluxes ITFP
      double precision trntdr (nbr_wvl,nif)
      double precision trndif (nbr_wvl,nif)
      double precision rupdir (nbr_wvl,nif)
      double precision rupdif (nbr_wvl,nif)
      double precision rdndif (nbr_wvl,nif)
      double precision fdirup (nbr_wvl,nif)
      double precision fdirdn (nbr_wvl,nif)
      double precision fdifup (nbr_wvl,nif)
      double precision fdifdn (nbr_wvl,nif)
      double precision dfdir  (nbr_wvl,nif)
      Double precision dfdif  (nbr_wvl,nif)
      double precision F_up   (nbr_wvl,nif)
      double precision F_dwn  (nbr_wvl,nif)
      double precision F_net  (nbr_wvl,nif)
      !double precision F_abs  (nbr_wvl,ncv)

      complex refindx

c    P.24 Briegleb and Light, Coakley et al., 1983 
      double precision refkm1,refkp1, rintfc, tdndif, tdrrdir, puny
      double precision tautot, wtot, gtot, ftot, gwt, mu0, mu
      double precision ts, ws, gs, lm, ue, extins, ne
      double precision alp, gam, apg, amg, R1, T1, swt, smr, smt, trn
      double precision rdr, tdr, nr, mu0n, R2, T2, Rf_dif_b, Tf_dif_b
      double precision Rf_dir_a, Tf_dir_a, Rf_dif_a, Tf_dif_a

      double precision gauspt(8), gauswt(8) ! Guassian angle loop 

c***********************************************************************
c     Data Statements

      data epsilon/1d-5/,exp_min/1d-5/,trmin/1d-4/, puny/1d-10/
      data gauspt /0.9894009d0, 0.9445750d0, 0.8656312d0, 
     & 0.7554044d0,0.6178762d0, 0.4580168d0, 0.2816036d0, 
     & 0.0950125d0/
      data gauswt /0.0271525d0, 0.0622535d0, 0.0951585d0,  
     & 0.1246290d0,0.1495960d0, 0.1691565d0, 0.1826034d0, 
     & 0.1894506d0/
      data found/.false./,check/.false./,do_write/.false./
c***********************************************************************

      PI = 4.D0*DATAN(1.D0)
      open(unit =34, file ='test_netsolarfluxes_SN',status = "UNKNOWN") !REJ testing

      dlambda = lambda(2) - lambda(1) ! assuming lambda is equally-spaced
       
c      apg = 0d0 ! ARK_2025/2/14
c      amg = 0d0 ! ARK_2025/2/14
c      alp = 0d0 ! ARK_2025/2/14
c      gam = 0d0 ! ARK_2025/2/14
c      rdr = 0d0 ! ARK_2025/2/14
c      tdr = 0d0 ! ARK_2025/2/14
c      tautot = 0d0 

c  check input
c      check = .true.
      if(check)then
      write(34,*)'fl_r_dif_a'; write(34,10)fl_r_dif_a  !REJ test
      write(34,*)'fl_r_dif_b'; write(34,10)fl_r_dif_b  !REJ test
      write(34,*)'g'; write(34,10)g  !REJ test
      write(34,*)'lambda'; write(34,10)lambda  !REJ test
      write(34,*)'lyr_typ'; write(34,*)lyr_typ  !REJ test
      write(34,*)'omega'; write(34,10)omega  !REJ testF
      write(34,*)'R_sfc'; write(34,10)R_sfc  !REJ test
      write(34,*)'tau'; write(34,10)tau  !REJ test
      write(34,*)'rad',rad
      stop 'check input'
      endif
c      write(34,10)lyr_typ  !REJ test
10    format(6f13.8) !REJ test
c
c Initialize values      
c      do 110 k=1,n  !REJ Make top_n_down  5/5/24
c      write(*,*)'n,nbr_wvl,ncv,nif', n,nbr_wvl,ncv,nif
      do 110 k=n,1,-1
            rdir(k)   = 0d0 ! layer reflectivity to direct radiation
            rdif_a(k) = 0d0 ! layer reflectivity to diffuse radiation from above
            rdif_b(k) = 0d0 ! layer reflectivity to diffuse radiation from below
            tdir(k)   = 0d0 ! layer transmission to direct radiation (solar beam + diffuse)
            tdif_a(k) = 0d0 ! layer transmission to diffuse radiation from above
            tdif_b(k) = 0d0 ! layer transmission to diffuse radiation from below
            trnlay(k) = 0d0 ! solar beam transmission for layer (direct beam only)

            do 129 iw=1,nbr_wvl
                  Fd(iw) = 0d0
                  Fs(iw) = 0d0
                  Temp1(iw) = 0d0 ! ARK_2025/2/14
                  Temp2(iw) = 0d0 ! ARK_2025/2/14
                  nreal(iw) = 0d0 ! ARK_2025/2/14 

                  trndir(iw,k) = 0d0
                  trntdr(iw,k) = 0d0
                  trndif(iw,k) = 0d0

                  rupdir(iw,k) = 0d0
                  rupdif(iw,k) = 0d0
                  rdndif(iw,k) = 0d0
                  fdirup(iw,k) = 0d0
                  fdirdn(iw,k) = 0d0
                  fdifup(iw,k) = 0d0
                  fdifdn(iw,k) = 0d0

                  F_up(iw,k)   = 0d0
                  F_dwn(iw,k)  = 0d0
                  F_net(iw,k)  = 0d0
                  
                  dfdif(iw,k)  = 0d0 ! ARK_2025/2/14
                  dfdir(iw,k)  = 0d0 ! ARK_2025/2/14

                 ! ARK_2025/2/14 Needed to reset arrays
                  trndir(iw,n+1) = 1d0 ! solar beam down transmission from top
                  trntdr(iw,n+1) = 1d0 ! total transmission from layers above
                  trndif(iw,n+1) = 1d0 ! diffuse transmission for layers above
                  F_up(iw,n+1)   = 0d0 ! 
                  F_dwn(iw,n+1)  = 0d0 ! 
                  F_net(iw,n+1)  = 0d0 
                  rupdir(iw,n+1) = 0d0 ! reflectivity to direct radiation for layers below
                  rupdif(iw,n+1) = 0d0 ! reflectivity to diffuse radiation for layers below
                  rdndif(iw,n+1) = 0d0 ! reflectivity to diffuse radiation for layers above
                  fdirup(iw,n+1) = 0d0 ! ARK_2025/2/14 Needed to reset arrays
                  fdirdn(iw,n+1) = 0d0 ! ARK_2025/2/14 Needed to reset arrays
                  fdifup(iw,n+1) = 0d0 ! ARK_2025/2/14 Needed to reset arrays
                  fdifdn(iw,n+1) = 0d0 ! ARK_2025/2/14 Needed to reset arrays


129         continue
110   continue
   
      
      do 123 iw=1,nbr_wvl

c MULTIPLY BROADBAND SOLAR FLUX WITH SPECTRAL SOLAR FLUX

      if (Fs0(iw) .gt. 0d0) then  !REJ
      ! Use direct-beam
         Fs(iw) = rad * Fs0(iw)  !REJ
      endif

      ! Use diffuse beam                  
      if (Fd0(iw) .gt. 0d0) then
         Fd(iw) = rad * Fd0(iw)  !REJ
      endif

      temp1(iw) = n_ice(iw)**2d0 - mImag(iw)**2d0 
     &            + sin(acos(mu_not))**2d0
      temp2(iw) = n_ice(iw)**2d0 - mImag(iw)**2d0 
     &            - sin(acos(mu_not))**2d0

c     Ice adjusted refractive index  (Liou 2004 Eq. 5.4.18)   
      Nreal(iw) = (sqrt(2d0)/2d0) * ( temp1(iw) + (temp2(iw)**2d0 +
     &  4d0*n_ice(iw)**2d0 * mImag(iw)**2d0)**(0.5d0) )**0.5d0
      
123       continue
 
c     Find the first ice layer that occurs between the last snow layer 
c     and the first ice layer.
c     If the top layer is ice, total reflection will occur at high solar
c     zenith angles.
c     It is usually better to include a surface-scattering layer
c     to avoid aphysical results

      kfrsnl = 0
      found = .false.
c      do 1234 k = n,1,-1  ! REJ block.  Next 4
c        lyr_typ(k) = 1     
c1234  continue
c        lyr_typ(3) = 2  !REJ LTYPE  line 214
      do 234 k = 1, size(lyr_typ)    
         if (lyr_typ(k) .eq. 2 ) then
            kfrsnl = k
            found = .true.
            exit
         end if
234   continue

      if (.not. found) then
            kfrsnl = 0
      end if

c ----- BEGIN Radiative Solver Adding Doubling Method -----

      tau0    = tau
      g0      = g
      omega0  = omega
c  Proceed down one layer at a time; if the total transmission to
c  the interface just above a given layer is less than trmin, then no
c  Delta-Eddington computation for that layer is done.

      do 321 iw = 1,nbr_wvl      ! wavelengths
            do 323 counter = n,1,-1
c initialize all layer apparent optical properties to 0
            rdir(counter)   = 0d0
            rdif_a(counter) = 0d0
            rdif_b(counter) = 0d0 
            tdir(counter)   = 0d0 
            tdif_a(counter) = 0d0  
            tdif_b(counter) = 0d0  
            trnlay(counter) = 0d0  
323         continue
      
      ! begin main level loop
      
c      do 345 k = 1, n   ! number of layers  Reverse order to top down. REJ 5/5/24
      do 345 k = n, 1,-1   ! number of layers

c compute next layer Delta-eddington solution only if total transmission
c of radiation to the interface just above the layer exceeds trmin.
cREJ5/16/24  Important fix for descending index.  if (trntdr(iw,k) .gt. trmin ) then
            if (trntdr(iw,k+1) .gt. trmin )then
c initialize current layer properties to zero; only if total
c transmission to the top interface of the current layer exceeds the
c minimum, will these values be computed below:
               mu0 = mu_not

      ! ice adjusted refractive index (Liou 2002)
               nr = Nreal(iw)
cREJ               if (k .lt. kfrsnl .or. kfrsnl .eq. 0) then
               if (k .gt. kfrsnl .or. kfrsnl .eq. 0) then
      ! above Fresnel Layer, mu0 is unchanged
                  mu0n = mu0
cREJ               else if (k .ge. kfrsnl) then
               else if (k .le. kfrsnl) then
      ! mu0 under the Fresnel Layer
      ! Eq. (5.4.13) Liou 2002
                  mu0n = cos(asin(sin(acos(mu0))/nr))
               end if

! calculation over layers with penetrating radiation
               tautot = tau0(iw,k)
               wtot   = omega0(iw,k)
               gtot   = g0(iw,k)

               ftot   = g0(iw,k) * g0(iw,k)
c       write(*,*)'3 items',iw,k,tautot,wtot,gtot

c coefficient for delta eddington solution for all layers
c Eq. 50: Briegleb and Light 2007

      ! layer delta-scaled extinction optical depth
               ts   = (1d0-wtot*ftot) * tautot        
      ! layer delta-scaled single scattering albedo
               ws   = ((1d0-ftot)*wtot)/(1d0-wtot*ftot)
      ! layer delta-scaled asymmetry parameter 
               gs   = gtot/(1d0+gtot) 

      ! lambda
               lm   = sqrt(3d0*(1d0-ws)* (1d0-ws*gs))  
      ! u equation, term in diffuse reflectivity and transmissivity   
               ue   = 1.5d0 * (1d0-ws*gs)/lm 

      ! extinction, MAX function prevents error if exp(-lm*ts) is < 1e-5       
               extins = max(exp_min, exp(-lm*ts))

      ! N equation, term in diffuse reflectivity and transmissivity  
               ne     = (ue+1d0)**2d0/extins - (ue-1d0)**2d0*extins
c     first calculation of rdif, tdif using Delta-Eddington formulas
c     Eq.50 : Briegleb 1992; alpha and gamma for diffuse radiation
      
      !  R BAR = layer reflectivity to DIFFUSE radiation
               rdif_a(k) = (ue**2d0-1d0)*(1d0/extins - extins)/ne 
               
      ! T BAR layer transmissivity to DIFFUSE radiation
               tdif_a(k) = 4d0*ue/ne
            
c     evaluate rdir,tdir for direct beam
      
      ! transmission TOA was incorrect/Noted by Adi
                         
               trnlay(k) = max(exp_min, exp(-ts/mu0n)) 
c  Eq. 50: Briegleb and Light 2007; alpha and gamma for direct radiation
      ! %lp = alpha(ws,mu0n,gs,lm)
               alp = (0.75d0*ws*mu0n) * (1d0 + gs*(1d0-ws))/
     &               (1d0 - lm**2d0* mu0n**2d0 + epsilon) 

      !  gam = gamma(ws,mu0n,gs,lm)
               gam = (0.5d0 * ws) * (1d0 + 3d0*gs*mu0n**2d0*(1d0-ws))
     &               / (1d0 - lm**2d0* mu0n**2d0 + epsilon)

               apg = alp + gam
               amg = alp - gam
      ! layer reflectivity to DIRECT radiation
               rdir(k) = apg*rdif_a(k) + amg*(tdif_a(k)*trnlay(k) - 1d0)   
      ! layer transmissivity to DIRECT radiation
               tdir(k) = apg*tdif_a(k) + (amg* rdif_a(k)-apg+1d0)
     &                   *trnlay(k)
               
c     recalculate rdif,tdif using direct angular integration over rdir,tdir,
c     since Delta-Eddington rdif formula is not well-behaved (it is usually
c     biased low and can even be negative); use ngmax angles and gaussian
c     integration for most accuracy:
      
      ! use R1 as temporary
               R1 = rdif_a(k) 
      ! use T1 as temporary
               T1 = tdif_a(k)
               swt = 0d0
               smr = 0d0
               smt = 0d0

c     loop through the gaussian angles for the AD integral

               do 455 ng=1, size(gauspt)  ! gaussian angles (radians)

      ! solar zenith angles
                  mu  = gauspt(ng)
      ! gaussian weight       
                  gwt = gauswt(ng)
      ! sum of weights
                  swt = swt + mu*gwt     
      ! transmission          
                  trn = max(exp_min, exp(-ts/mu))

      ! alp = alpha(ws,mu0n,gs,lm)
                  alp = (0.75d0*ws*mu) * (1d0 + gs*(1d0-ws))
     &                  /(1d0 - lm**2d0 * mu**2d0 + epsilon)
      ! gam = gamma(ws,mu0n,gs,lm)
                  gam = (0.5d0* ws) * (1d0 + 3d0*gs*mu**2d0*(1d0-ws))
     &                  /(1d0-lm**2d0* mu**2d0 + epsilon)

                  apg = alp + gam
                  amg = alp - gam
                  rdr = apg*R1 + amg*T1*trn - amg
                  tdr = apg*T1 + amg*R1*trn - apg*trn + trn
      ! accumulator for rdif gaussian integration
                  smr = smr + mu*rdr*gwt 
      ! accumulator for tdif gaussian integration
                  smt = smt + mu*tdr*gwt 

455               continue  ! ng; gaussian angles for the AD integral

               rdif_a(k) = smr/swt
               tdif_a(k) = smt/swt

      ! homogeneous layer
               rdif_b(k) = rdif_a(k)
               tdif_b(k) = tdif_a(k)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Fresnel layer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (k .eq. kfrsnl) then
                  !write(*,*) 'Found a fresnel layer! woohoo!'
      ! ice complex index of refraction
                  refindx = complex(n_ice(iw),mImag(iw))

      ! critical angle where total internal reflection occurs
                  critical_angle = asin(refindx)

      ! compare incoming angle to critical angle
                  if (acos(mu_not) .lt. critical_angle) then

c     compute fresnel reflection and transmission amplitudes
c     for two polarizations: 1=perpendicular and 2=parallel to
c     the plane containing incident, reflected and refracted rays.

      ! Eq. (5.4.18a-b) of Liou (2002)

      ! reflection amplitude factor for perpendicular polarization
                     R1 = (mu0-nr*mu0n) / (mu0 + nr*mu0n)  
      ! reflection amplitude factor for parallel polarization
                     R2 = (nr*mu0 - mu0n) / (nr*mu0 + mu0n)
      ! transmission amplitude factor for perpendicular polarization
                     T1 = 2d0*mu0 / (mu0 + nr*mu0n)
      ! transmission amplitude factor for parallel polarization              
                     T2 = 2d0*mu0 / (nr*mu0 + mu0n)       

      ! unpolarized light for direct beam
      ! Eq. 21; Brigleb and light 2007
                     Rf_dir_a = 0.5d0 * ((R1**2d0) + (R2**2d0))
                     Tf_dir_a = 0.5d0 * (T1*T1 + T2*T2)*nr*mu0n/mu0
                     
                  else ! total internal reflection

                     Tf_dir_a = 0d0
                     Rf_dir_a = 1d0

                  endif ! critical angle check


c     precalculated diffuse reflectivities and transmissivities
c     for incident radiation above and below fresnel layer, using
c     the direct albedos and accounting for complete internal
c     reflection from below; precalculated because high order
c     number of gaussian points is required for convergence:

                  Rf_dif_a = FL_r_dif_a(iw)
                  Tf_dif_a = 1d0 - Rf_dif_a
                  Rf_dif_b = FL_r_dif_b(iw)
                  Tf_dif_b = 1d0 - Rf_dif_b
                  

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! the k = kfrsnl layer properties are updated to combined
      !! the fresnel (refractive) layer, always taken to be above
      !! the present layer k (i.e. be the top interface):

      ! denom interface scattering
                  rintfc   = 1d0 / (1d0-Rf_dif_b*rdif_a(k)) 
      
      ! layer transmissivity to DIRECT radiation
      ! Eq. B7; Briegleb & Light 2007
                  tdir(k)   = Tf_dir_a*tdir(k) +     
     &            Tf_dir_a*rdir(k) * Rf_dif_b*rintfc*tdif_a(k)
      
      ! layer reflectivity to DIRECT radiation
      ! Eq. B7; Briegleb & Light 2007
                  rdir(k)   = Rf_dir_a + Tf_dir_a*rdir(k) * 
     &                  rintfc*Tf_dif_b

      ! R BAR = layer reflectivity to DIFFUSE radiation (above)
      ! Eq. B9; Briegleb & Light 2007
                  rdif_a(k) = Rf_dif_a + Tf_dif_a*rdif_a(k) * 
     &                  rintfc*Tf_dif_b
      ! R BAR = layer reflectivity to DIFFUSE radiation (below)
      ! Eq. B10; Briegleb & Light 2007
                  rdif_b(k) = rdif_b(k) + tdif_b(k)*Rf_dif_b * 
     &                  rintfc*tdif_a(k)
      ! T BAR layer transmissivity to DIFFUSE radiation (above) 
      ! Eq. B9; Briegleb & Light 2007
                  tdif_a(k) = tdif_a(k)*rintfc*Tf_dif_a
      ! Eq. B10; Briegleb & Light 2007 
                  tdif_b(k) = tdif_b(k)*rintfc*Tf_dif_b

      ! update trnlay to include fresnel transmission
                trnlay(k) = Tf_dir_a*trnlay(k) 
               endif ! k = kfrsnl  End Fresnel Block

            endif ! trntdr(k,iw) > trmin
c     For the remainder of the routine, REJ chnged k-1 to k+1 and k+1 to k-1

c     Calculate the solar beam transmission, total transmission, and
c     reflectivity for diffuse radiation from below at interface k,
c     the top of the current layer k:
c
c              layers       interface
c
c               ---------------------  k+1
c                        k
c               ---------------------  k
c                        k-1
c               ---------------------  K-1

c     we ignore refraction between layers and underlying material:

c        
c                      layers       interface
c        
c              ---------------------  k+1
c                       k
c              ---------------------  k
c              \\\\\\\ underlying material \\\\\\\

      ! Eq. 51; Briegleb and Light 2007
c REJ 5/15/24.  For descending procedure, K and K+1 are interchanged 
c for interfaces in following block
c solar bean transmission from top
            trndir(iw,k) = trndir(iw,k+1)*trnlay(k)        
      ! interface multiple scattering for k
            refkm1         = 1d0/(1d0 - rdndif(iw,k+1)*rdif_a(k)) 
      ! direct tran times layer direct ref
            tdrrdir        = trndir(iw,k+1)*rdir(k)             
      ! total down diffuse = tot tran - direct tran
            tdndif         = trntdr(iw,k+1) - trndir(iw,k+1)   
      ! total transmission for layers above (Note:trndir is nodal) 
            trntdr(iw,k) = trndir(iw,k+1)*tdir(k) + 
     &      (tdndif + tdrrdir*rdndif(iw,k+1))*refkm1*tdif_a(k)

      ! Eq. 51; Briegleb and Light 2007

      ! reflectivity to diffuse radiation for layers above
            rdndif(iw,k) = rdif_b(k) + 
     &       (tdif_b(k)*rdndif(iw,k+1)*refkm1*tdif_a(k)) 
      ! diffuse transmission to diffuse beam for layers above  
            trndif(iw,k) = trndif(iw,k+1)*refkm1*tdif_a(k) ! Line 524      
345         continue  ! k   end main level loop; number of layers

c     compute reflectivity to direct and diffuse radiation for layers
c     below by adding succesive layers starting from the underlying
c     layers and working upwards:
c
c              layers       interface
c
c       ---------------------  k+1
c                 k
c       ---------------------  k
c                k-1
c       ---------------------  k-1

c     set the underlying ground albedo  
c     To avoid a 0 nodal suscript, use R_sfc(iw) for both direct and
c     diffuse reflectivity from layer/node = 0.
      ! reflectivity to direct radiation for layers below
                rupdir(iw,1) = R_sfc(iw) !REJ Ascend rupdir(iw,n+1)  
      ! reflectivity to diffuse radiation for layers below
                rupdif(iw,1) = R_sfc(iw) !REJ Ascend rupdir(iw,n+1)

      ! starts at the bottom and works its way up to the top layer

      
cREJ            do 567 k = n,1,-1            
C           do 567 k = 1,n    
C
C     ! Eq. B2; Briegleb and Light 2007
C     ! interface scattering
C         
CcREJ               refkp1        = 1d0/( 1d0 - rdif_b(k)*rupdif(iw,k+1))
C              refkp1        = 1d0/( 1d0 - rdif_b(k)*rupdif(iw,k))
C           
C    
C     ! dir from top layer plus exp tran ref from lower layer, interface
C     ! scattered and tran thru top layer from below, plus diff tran ref
C     ! from lower layer with interface scattering tran thru top from below
C
C               rupdir(iw,k+1) = rdir(k)
C    &       +  (  trnlay(k)  + rupdir(iw,k) 
C    &            + (tdir(k)-trnlay(k))* rupdif(iw,k))*refkp1*tdif_b(k)
C
C    
C     ! dif from top layer from above, plus dif tran upwards reflected and
C     ! interface scattered which tran top from below
C               rupdif(iw,k+1) = rdif_a(k) 
C    &          + tdif_a(k)*rupdif(iw,k)*refkp1*tdif_b(k) 
C       
      do 567 k = 2, n+1
        
          refkp1 = 1d0 / (1d0 - rdif_b(k - 1 ) * rupdif(iw,k - 1 ))

          rupdir(iw,k) = rdir(k - 1 )
     & + (trnlay(k - 1 ) * rupdir(iw,k - 1 )
     & + (tdir(k - 1 ) - trnlay(k - 1 )) 
     & * rupdif(iw,k - 1 )) * refkp1 * tdif_b(k - 1 )

          rupdif(iw,k) = rdif_a(k - 1 )
     & + tdif_a(k - 1 ) * rupdif(iw,k - 1 ) * refkp1 * tdif_b(k - 1 )

     
      !write(*,*) 'k = ', k
      !write(*,*) 'rupdir(iw,k) = ', rupdir(iw,k)
      !write(*,*) 'rupdif(iw,k) = ', rupdif(iw,k)

567         continue ! k

321   continue ! iw; number of wavelengths

       write(34,*)' END Radiative Solver Adding Doubling Method'
       write(34,*)'  iw   k     fdirup       fdifup      fdirdn  fdifdn'

      ! Fluxes at interfaces

      do 601 iw = 1,nbr_wvl
         
c5/24/2024         do 701 k =1,n+1
         do 701 k =n+1,1,-1
cREJ Changed all indices on FLUXES in following block from k to k+1
      !Eq. 52; Briegleb and Light 2007

      ! interface scattering

            refk          = 1d0/(1d0 - rdndif(iw,k)*rupdif(iw,k))
      ! dir tran ref from below times interface scattering, plus diff
      ! tran and ref from below times interface scattering
            fdirup(iw,k) = (trndir(iw,k)*rupdir(iw,k) + 
     &       (trntdr(iw,k)-trndir(iw,k)) 
     &       *rupdif(iw,k))*refk

      ! dir tran plus total diff trans times interface scattering plus
      ! dir tran with up dir ref and down dif ref times interface scattering
            fdirdn(iw,k) = trndir(iw,k) + (trntdr(iw,k) 
     &       - trndir(iw,k) + trndir(iw,k) 
     &      *rupdir(iw,k)*rdndif(iw,k))*refk

      ! diffuse tran ref from below times interface scattering
            fdifup(iw,k) = trndif(iw,k)*rupdif(iw,k)*refk
      
      ! diffuse tran times interface scattering
            fdifdn(iw,k) = trndif(iw,k)*refk
   
      ! dfdir = fdirdn - fdirup
            dfdir(iw,k) = trndir(iw,k) 
     &      + (trntdr(iw,k)-trndir(iw,k)) * (1d0 - rupdif(iw,k)) * refk 
     &      -  trndir(iw,k)*rupdir(iw,k)*(1d0 - rdndif(iw,k))*refk

            if (dfdir(iw,k) .lt. puny) then
               dfdir(iw,k) = 0d0
            endif

            dfdif(iw,k) = trndif(iw,k)*(1d0-rupdif(iw,k))*refk
            
            if (dfdif(iw,k) .lt. puny) then
               dfdif(iw,k) = 0d0
            endif

      If (k .ge. 25 .and. iw .ge. 8 .and. iw .lt. 50)then
11    format(2i4,6f13.8) !REJ test 
c      if( (float(iw)/10.0)-int(float(iw)/10.0).lt.dabs(0.009d0))  !REJ_2025/04/16
c     &  write(34,11)iw,k,fdirup(iw,k),fdifup(iw,k), 
c     &     fdirdn(iw,k),fdifdn(iw,k) !REJ_2025/04/16
      endif
701      continue ! k = number of layers

601   continue ! iw = number of wavelengths

c ----- END Radiative Solver Adding Doubling Method -----
cc       call total_solar(n,iy,jday,ihour,fdirup,fdifup,fdirdn,  !Restore
cc     & fdifdn,rad,Fs0,Fs,Fd0,Fd,mu_not,PI,Fabs_int,albedo_int,dlambda, !Restore
cc     & lambda,solartest)  !Restore

      return
      end