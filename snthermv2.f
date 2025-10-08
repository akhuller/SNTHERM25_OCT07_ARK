c***********************************************************************
      program snthermv2ARK
c      version sntherm24/sntherm24_June15  Updated to sntherm24_check3
c     Sept 09, 2025 Version with new precip
c***********************************************************************
c One-Dimensional Mass and Energy Balance Model for a Snowcover
c SNTHERM.89  Model developed by R. Jordan,  US Army Cold Regions
c                 Research and Engineering Lab, Hanover, NH
c             Model coded by R. Jordan, USA CRREL and
c                 J. Jones, Sparta Systems, Inc., Lexington, MA
c             SNTHERM24 solar absorption and turbulent transfer routines 
c                 developed and coded by A. Khuller, Jet Propulsion 
c                 Laboratory,  Pasadena, CA
c             Adaptation of these codes to SNTHERM by R. Jordan in 2024  
c***********************************************************************
c Complete TECHNICAL DOCUMENTATION can be found in Jordan, R.,"A
c One-Dimensional Temperature Model for a Snow Cover," USA-CRREL
c Special Report 91-16. A partial nomenclature list and user's
c instructions are contained in the file DOCUMENTATION. Instructions
c for running the code are contained in the report, "User's quide
c for CRREL one-dimensional snow temperature model," available upon
c request from R. Jordan
c DISCLAIMER: Internal documentation of the code is incomplete and 
c may contain errors.  No attempt has been made to reconcile the
c nomenclature in the code with that in the Technical Report.
c In some instances they agree, but in many there are minor
c inconsistacies in definition, so that care must be taken when
c comparing the equations in the Report with those in the code.
c
c RESTRICTIONS and LIMITATIONS on the code:  SNTHERM24 has much
c in common with the earlier release, sntherm89. The anticipated 
c SNTHERM2, that couples capillary tension with temperature depression 
c and adds capillary effects to water flow prediction, never reached 
c the level of public release. Instead, a temporary drain is  
c constructed which removes water from the system when it reaches the  
c snow-soil interface.F
c In the interim years since the its release, SNTHERM89 has been
c used extensively in polar regions. Along with updated solar absorption
c and turbulent transfer routines, the knowledge from this experience is 
c being added to sntherm24. 

C***********************************************************************
c  1. Declare arrays and variable types
c     Principal arrays and variables are declared in ARRAYS
c     ALL REAL VARIABLES ARE DOUBLE PRECISION.  An ending of 'o'
c     in a variable name denotes the previous (old) time step.
c***********************************************************************
      implicit none
      include 'const'
      include 'arrays'
      include 'arrays_spectral' ! NEW block of read-in spectral data
      include 'l_arrays.txt'    ! NEW common for logical H2Otype
c      
 !REJ 6/19/24
c       Note: The parameters ncv = 50, nif = 51, nwl = 281, nd = 100,
c       ld = 5 are sized in 'const'
c       Most SNTHERM CV/nodal arrays in common are sized with nd = 100 
c       and layertype (snow, frn, soil,etc.) ld = 5.  
c       Within MatLabSolar, CV arrays are dimensioned as n and interfaces
c       as n+1. Note that in SNTHERM, n and n+1 are the startup number of 
c       CVs and interfaces, which can change in the course of a run. n is
c       distinguished from the allotted storage dimension, nd, for CV arrays.
c       Memory for CV arrays within the SNTHERM MLS subroutines is currently 
c       dimensioned as ncv = 50 and nif = 51.  Utilizing ncv and nif retains
c       the distinction between CVs and interfaces.
c       nwl is simply a shorthand for the constant nbr_wvl,
c       the number of spectral bands between 0.2 and 30 microns.
 
c Local
      integer ido,iy,jday,ihour,iy2,jday2,ihour2,ipday,iphour,ithour  
      integer i,j,m,ii,ibasestep,icalcstep,iw !added iw REJ_2025/05/23
      integer isolarcalc,ircalc,islope,itimezone,itracks,itm
      integer newsno,iptype,ioutfiltrate, botnode,ifluxout
      integer ibounds(15,2),ititle,nzone,nz
      integer igood,ngoodmin,imetcalc,imin,istboff,iqturb
      integer writefmt,use_SSA,spyear,nc,istep_conv_out  !REJ_2025/09/04 
      integer iLAPs,iLAPs_type(nd) !REJ_25/01/09 impurities added
      double precision dtmin,dtmax,dssallowed,errtallowd,dtssmax,dtsum
      double precision dzn,dznm,pinv,tm,de0,bp,tmsg,rtmsq ! RJ 1/30/23 
      double precision snowdepth,difftemp,dtsmin,eta0
      double precision thsnow,tlsnow,ssisnow,bifall,bifallin,dmlimit
      double precision dlsdrw,overburden,dzinc,ci0,dzmin,dlvdrw
      double precision a1snow,rw,emsnow,albsnow,bwfall,blfall
      double precision rmsqt1,height(3)    
      double precision cdryair,dlatt,azslope,rhoair
      double precision bvi0,wmass(nd),elev,bvw0  !REJ_2025/05/09
      double precision ssminsnow,thkair,thkice,dsnowfall,unbarmax
      double precision dlongt,clearness,effceiling
      double precision frh(ld),g1,g2,r1,r2,t1,t2,wt,totaltime,dum
      double precision floo(nd),dzoo(nd),TsurfEst,Esurface,depth(nd) !REJ 1/7/2024
      double precision iLAPs_concen(nd) !REJ_2025/1/09 impurities
      double precision exp1  ! fmt and exp extra writes by ceretha 1996
      double precision meltflux,uosave,topflux,frvis,botdsol,botflux 
      double precision sum,timesum,Tkairsave,Total_Latent  !REJ_2025/01/22+06/05
      double precision bextnir,Totalmass,sum_vapor_dt,sum_delta_dz !REJ_2025/07/22 +08/27
      double precision snow_threshold,snow_dz_min,mass1,mass2,SSA1,SSA2 ! ARK_2025/09/11
      integer istart_eff,prev_iy ! ARK_2025/09/11

c    Debug Option
      logical debug
      integer idebug_step   
      character*12 debug_variable

c*m*m*m*m*m*m*m*m*m* MatLabSolar_Declarations *m*m*m*m*m*m*m*m*m*m*m*m*m
c     Rename "MatLabSolar" with something more descriptive !REJ_2025/1/09
 
      ! REJ  changed n to ncv, which is set in const
      double precision tau(nbr_wvl,ncv), omega (nbr_wvl,ncv) 
      double precision g (nbr_wvl,ncv)      
      double precision t_net_lyr(nbr_wvl), w_net_lyr (nbr_wvl)
      double precision g_net_lyr (nbr_wvl)
      double precision R_sfc(nbr_wvl), albedo(nbr_wvl), Fnet(ncv)
      double precision F_abs(nbr_wvl, ncv) !Fs added by REJ.Moved to common 8/1/24
      double precision r_d, rho_d, porosity_rad, ssArea, rbbl 
      double precision lyr_typ(ncv), mu_not  !Fs0 out RBJ Moved to common 8/1/24
      double precision fdirup(nbr_wvl,ncv), fdifup(nbr_wvl,ncv) !REJZ
      double precision fdirdn(nbr_wvl,ncv), fdifdn(nbr_wvl,ncv) !REJZ
      double precision albedo_int,dlambda,fabs_int(nd),pi,time0,time 
     &    !REJ_2025/07/24
      logical solartest,do_Mie,MLS_run_once
      logical print,istop,recross,repeat,fullnode,
     &    thinnode,ICVPfixd,Neumann,prnt_hour,spinup !REJ2025/05/24
      logical writematlabflag,writefluxmatlabflag ! ARK 7/29/24
      logical writematlablastyearflag ! ARK 12/26/24  
      logical usemeasuredbottomtemp   ! ARK_2025/10/06
      character*160 dust_mie_dir, mie_dir, miefile
      character*160 filename,x1,fnm(nfiles),fname 
      character*160 spfilename

      common /insolc/  r1(4,4),r2(4,4),t1(4,4),t2(4,4),wt(4,6)

c Function declarations
      integer nmelt
      double precision fvapri,thrk,fvaprw,fliquid,fgrain,fd_init
         !REJ_2025/03/19
      integer century,planet ! ARK 12/26/24
c*tt*tt*tt*tt*tt*tt*tt*tt*KCTURB_Declarations*tt*tt*tt*tt*tt*tt*tt*tt*tt
      logical USEKCTURB ! ARK_2025/1/13
      double precision H_PBL,ZMAX,BETA,BETA2 ! ARK_2025/1/11
      integer FX_OPT,dryair,volatile ! ARK 12/26/24 [century & planet declared above]
      double precision md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b
      double precision md_a0,md_a1,md_a2,md_a3,md_a4,md_b0,md_b1,md_b2
      double precision md_b3,md_b4,md_beta_a,md_beta_b,md_gamma_a
      double precision md_gamma_b,md_kappa_a,md_kappa_b, md_m_a,md_m_ab
      double precision md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b
      double precision md_sum_v_a,md_sum_v_b
      double precision q2,Fm,rho_rep,rhocp_rep,qo,qairo
      double precision Cpd0,Fh02,Fw02,H0sen,lambda0,Ls,LE
      double precision Ma,Mw,dqsat ! ARK_2025/09/19 Removed ustar
c END KCTURB_Declarations
c*tt*tt*tt*tt*tt*tt*tt*tt*KCTURB_Declarations*tt*tt*tt*tt*tt*tt*tt*tt*tt
C***********************************************************************
c    1b. Data statements
C**********************************************************************
c Data statements
      data ibasestep,ii,igood,ititle/4*0/ 
      data rmsqt1/0d0/,tm,tmsg/2*0d0/,snowdepth/0d0/
      data dzinc/.04d0/,de0/0.9d-4/,albsnow/0.78d0/,emsnow/0.99d0/ ! ARK_2025/09/11 changed emsnow to 0.99 from 0.97
cREJ 8/5/24      data rw/461.296d0/,a1snow/100D0/,dzmin/2d-3/,ci0/2117d0/
      data rw/461.296d0/,a1snow/400d0/,dzmin/2d-3/,ci0/2117d0/ !REJ_2025/08/31
      data ssisnow/0.04d0/,thkair/2.30d-2/
      data dzn/1.66667d-2/,dznm/3.33333d-2/,cdryair/1.d3/,rhoair/1.276/
      data ido/0/,ithour/0/,newsno/0/
      data icalcstep/0/,recross/.false./,repeat/.false./
      data g1/2.0d-8/,g2/4.0d-12/ ! ARK_2025/09/10 changed g1
      data botflux/0.0d0/,botdsol/1d-7/,frvis/0.666667d0/ !REJ 1/07/24
      PI = 4.D0*DATAN(1.D0) ! ARK 7/26/24
      data istart_eff/0/
c_______________________________________________________________________
c Debug input feature.  Comment-out if not used
      debug = .false.; idebug_step = 0 
         

c NEW INPUTS THAT NEED TO BE INTEGRATED INTO INPUT FILE ARK_2025/11/1     
      ! Read in Cosine of Zenith Angle or set to fixed value ARK 7/26/24
      mu_not  = 1d0    ! ARK 12/27/24
      planet  = 1      ! ARK 12/27/24 (1 = Earth, 2= Mars)
      century = 19     ! ARK 12/27/24 (for Earth solar zenith angle calc - NOT USED ANYMORE!)
      ZMAX    = 120d0  ! Max grid size for turbulence calculations (m) 
      BETA    = 1.25d0 ! Gustiness parameter (can be measured)
      BETA2   = 0d0    ! Windless coefficient (set to zero)
      FX_OPT  = 2      ! Stability function type (1 = old, 2 = new)
      
      ! To be read from GETMET since this changes hourly
      H_PBL   = 600d0  ! ARK_2025/1/11 Height of boundary layer (m)     
c END NEW INPUTS 
     
c     Next are logical declarations. Tailor to individual needs.
c     Note: ICVPfixd & Neumann must only be declared in this one location
c     REJ Note: ICVPfixd replaces Khuller
      data ICVPfixd/.false./,Neumann/.false./,prnt_hour/.false./
     &    !REJ_2025/08/05 ICVPfixd = .F. !2025/04/12 Chng fixd back
      data fullnode/.false./ !11/07/23 Was triggering dt = dtmin ??
      data usemeasuredbottomtemp /.true./ !ARK_2025/10/06
      snow_dz_min = 1d-4 ! ARK_2025/09/11
      snow_threshold = 9.25d-6  ! ARK_2025/09/11 

c         Recheck!
      Total_Latent = 0d0 !REJ_2025/06/05
c*m*m*m*m*m*m*m*m*m* MatLabSolar_Data *m*m*m*m*m*m*m*m*m*m*m*m*mm*m*m*      
      ! Read in Cosine of Zenith Angle or set to fixed value ARK 7/26/24
      
      mu_not = 1d0

c*m*m*m*m*m*m*m*m*m* MatLabSolar_Data *m*m*m*m*m*m*m*m*m*m*m*m*mm*m*m*
cREJ_2025/01/10      data do_Mie/.true./,spinup/.false./ cf
      data do_Mie/.true./,spinup/.false./  !REJ_2025/02/13 Change
cREJ_2025/01/10      data do_Mie/.false./,spinup/.false./ !
      data solartest/.false./,MLS_run_once/.false./  !REJ_2025/04/29
c      data writematlabflag/.true./,writefluxmatlabflag/.true./ ! ARK 12/19/24
c      data writematlablastyearflag/.true./ ! ARK 12/19/24
      data writematlabflag/.true./,writefluxmatlabflag/.true./ ! ARK 7/29/24 ./,REJ_2025/02/26
      data writematlablastyearflag/.false./ ! ARK 12/19/24,REJ_2025/02/26
c*tt*tt*tt*tt*tt*tt*tt*tt*TURBULANCE SWITCH*tt*tt*tt*tt*tt*tt*tt*tt*tt*
       data USEKCTURB/.true./ ! true = Use KCTURB, false = Use QTURB
c     &    !ARK ~2025/11/1
c      data USEKCTURB/.false./ ! true = Use KCTURB, false = Use QTURB 
cS*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*
c OPTIONAL SPIN-UP. Data Section: Set spinup/.true./ 
c                                 Set ICVP.fixd/.true./ 
c                                 Set Neumann/.true./
c                                 Specify spyear below 
c                                 set 2-digit version ID line 
c                             ICVP = Inherent Control Volume Properties
      nsoil = 0  !REJ_2025/01/16  !TEMPORARY(perhaps not needed)  FIX LATER 
      if(.not. spinup)then
        spyear=1  !Standard run. NO SPIN-UP
      else
        if(.not.ICVPfixd)stop 'Use Fixed ICVPs for spinup'
        if(.not. Neumann)stop 'Use flux bottom BC for spinup'
        spyear=5  !For large spyear, limit optional printouts  !REJ_2025/01/25 Chgd from 20
c        if(do_Mie .and. spyear .gt. 12) stop 'increase spyear with care'
      endif

      if(ICVPfixd)then ! ARK_2025/2/10
         write(*,*) 'Using Fixed Inherent Control Volume Properties'
      end if

c----START SPIN_UP DO-LOOP----------------------------------------------

      DO ii=1,spyear ! End of this loop is at the bottom of MAIN
      if(spinup) write(*,*)'Begin spin-up year', ii
     
cS*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S*
c
c***********************************************************************
c  2. Open files and initialize optional spin-up runs
c***********************************************************************
      If(ii .eq. 1)Then
c  Read in file names from FILENAME. (see DOCUMENTATION)
      open(81,file='FILENAME',status='old')
      read(81,*) (fnm(i),i=1,nfiles)
      close(81)
c
c  Check for repeat of file names
      do 1 i=1,nfiles
        fname=fnm(i)
        do 2 j=1,nfiles
           istop = i.ne. j .and. fname .eq. fnm(j)
           if(istop) then
             write(*,*)'  **WRITE:_I/O File name repeated.**'
             stop '**Error in I/O file names. Execution Terminated **'
           endif
 2      continue
 1    continue
      Endif
      open(90,file=fnm(1),status='old')       !Input file (getinput.f)
      open(88,file=fnm(2),status='old')       !Met file (getmet.f)
c     open(89,file=fnm(3),status='old')       !Opt Measured Temperature
      open(80,file=fnm(4),status='unknown')   !Main output file
      open(7, file=fnm(5),status='unknown')   !Opt flux file (flux.f)
      open(99,file=fnm(6),status='unknown')   !Opt filtrate file (filt.out)
      open(180, file = 'Tsurfest_error',status='unknown') !REJ_2025/02/11. temp write
      open(20, file='output.txt', status='unknown') ! matlab output file
      open(49, file='fluxoutput.txt',status='unknown') ! matlab flux output file
c     Notes: Met file must contain 1 year of data for spinups 
c            Caution. Filtrate optional output can be very large!
c     Most of the following files are used for testing and debugging
      open(95, file= 'test_out',    status='unknown')!REJ_2025/03/07/ Dedicated test_out file
      open(115,file= 'massbal_out', status='unknown') !REJ2025_05/28 Dedicated MB testing
      open(120,file= 'converge_out',status='unknown') !REJ2025_05/28 ditto converge testing
      open(125,file= 'vapor _bal',  status='unknown') !REJ2025_07/22 ditto vapor_mass balance
      open(140,file= 'combo_out',   status='unknown') !REJ_2025/08/04
      open(45, file= 'summit_2007_precip', status='unknown')!REJ_2025/04/30
      open(46,file = 'summit_2007_compact',status='unknown')!REJ_2025/04/3

c  Check that the correct layerin file is used for solartest comparison.
cREJ ~Aug/24     if(solartest .and. fnm(1) .ne. 'MatLabSolar_solartest.in')
cREJ ~Aug/24     & stop  'Not the solartest layer.in file. Check FILENAME'
      if(.not. solartest .and. fnm(1) .eq. 'MatLabSolar_solartest.in')
     & stop  'Solartest layer.in file found basic run. Check FILENAME'

c*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m

      If(do_mie)Then
          dust_mie_dir = "sntherm_Mars_Dust_Mie_bin/"
          mie_dir = "sntherm_Ice_Mie_bin_v2/"
c     MATLABSOLAR INCLUSION. REJ 4/26/24   
         include 'csv.txt'  !Content is function of wavelength only 
c        csv.txt opens and closes units 10 thru 16
c        csv.txt reads:  lambda(nbr_wvl),   g_data(nbr_wvl),
c                        mImag(nbr_wvl),    n_ice(nbr_wvl),
c                        FL_r_dif_a(nbr_wvl),FL_r_dif_b(nbr_wvl),
c                        Fd0(nbr_wvl) or Fsd0(nbrwvl) !REJ ~Aug/24

      ! ARK 7/29/24 Moved this here since these do not ever change
      ! Initialize both Fs0 and Fd0. REJ. Ask Adi about next lines. Fs0 now
      !                                   read in csv.txt. REJ_2025/03/10
              do i = 1,nbr_wvl
              R_sfc(i) = 0d0 ! spectral albedo of underlying surface
              !Fs0(i)   = 0d0 ! direct-beam downward flux [W/m^2/band] ARK 7/26/24
              !Fs0(i) = Fs0(i)/PI ! ARK 7/26/24- removed on 12/26/24
              Fd0(i) = 0d0  !REJ 7/31/24
              end do

      EndIf
      
c     KCTURB - LOAD GAS MOLECULAR INFORMATION. ARK_2025/1/13
      If(USEKCTURB)Then
      call mol_setup(planet,
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b,dryair,volatile,Ma,Mw)

      end if
c***********************************************************************
c  3. Read in snowpack/soil parameters and initial data from layer.in
c     Important notes: layers and elements number from the bottom up.
c***********************************************************************
       call getinput(ifluxout,isolarcalc,ircalc,islope,itracks,
     &    itm,ioutfiltrate,imetcalc,ngoodmin,itimezone,iy2,jday2,ihour2,
     &    pinv,bp,frvis,albsnow,height,dtmin,dtsmin,dtmax,dtssmax, !RJ
     &    dssallowed,errtallowd,emsnow,dlatt,elev,dlongt,azslope,dzmin,
cREJ_2025/05/24 dzn,dznm,fnm,ssisnow,frh,bifallin,dmlimit,istboff,iqturb,eta0, 
     &    dzn,dznm,fnm,ssisnow,frh,bwfall,dmlimit,istboff,iqturb,eta0,
     &    exp1) 
        bifallin = bwfall  !REj_2025/05/25 Temporary.read as this. bwfall in main

      planet = 1   ! ARK 12/27/24
      century = 19 ! ARK 12/27/24    

       if(frvis.gt. 1d0)stop 'Inpt Bext is now frvis, frac of vis solar'
 
       do 12 i=nsoil+1,n
         binew(i)=dmlimit/1.15d0
c REJ 2025/1/9  Next is temporary specifications for impurities.  
c               Move to getinput later
         iLAPs_type(i) = 1
         iLAPs_concen(i) = 1d-7 !Warren(1984) for Antarctic soot
12     continue

       writefmt = ifluxout/10  !two lines added by Ceretha 1996
       ifluxout = ifluxout - (writefmt * 10)  !REJ needs to check this

cS*S*S*S*S*S*S*S*S*S*SPINUP FOR OUT-YEARS*S*S*S*S*S*S*S*S*S*S*S*S*S*S*
      if(ii .ge. 2)then !New annual spinup input !REJ+4 12/29/23
         read(92,*) ! Read over heading REJ 1/7/2024
         do 15 i=1,n        
           read(92,*)j,depth(i),to(i),dzo(i),bwo(i),do(i)
15       continue 
         close(92)
       endif   
c***********************************************************************
c  4. Calculate constant parameters
c***********************************************************************
      call calconstant(a1snow,bp,bvi0,cdryair,ci0,dlsdrw,
     & height,rw,snowdepth,ssisnow,ssminsnow,thsnow,
     & tlsnow,bvw0,dlvdrw)

c***********************************************************************
c  5.  Initialize density; other items.
c***********************************************************************
c      initial=1 !REJ 9/19/24
      dt=dtmin
      dt=dmax1(dt,1d0)
      dto=dtmin
      totaltime=dtmin
      istep_conv_out = 550 !REJ_2025/09/04
      do 550 i=1,n
         m=k(i)
c        **NEW**  H2OLayers include snow/firn/ice  !REJ
         H2Olayer(i) = .false. 
         if(  ltype(k(i)) .eq. 1 .or. 
     &      ( ltype(k(i)) .ge. 90 .and. ltype(k(i)).le. 99 ) )
     &   H2Olayer(i) = .true.    
         if(do(i) .le. dtol1 .and. H2Olayer(i) )then
            do(i)=fd_init(bi(i))  ! estimate grain size  !RBJ 8/1/24
c           for the initial step, bmelt is used to flag unavailable
c           grain diameter data.
            bmelt(i)=999d0
         endif
         d(i)=do(i)
         r(i)=(d(i))/2d0 !Use for Mie computations REJ 4/27/24       
         call density(to(i),bwo(i),td(i),td13(i),flo(i),blo(i),bi(i),
     1        bt(i),dmass(i),bdjp(m),a243(m),bd(m),a1(m),
     2        0,dzo(i),flgo(i),sso(i),dice,ssi(i),porosity(i),
     &        dmvol(m),ltype(m),impermeable(i),idelete(i),
     &        solidporosity(i),ipond(i),dicevol(i),dliqvol(i), 
     &        rhowater,H2Olayer(i)) !REJ 6/24/24
        
           dmasso(i)=bwo(i)*dzo(i) !REJ temp
     
         if(sso(i) .gt. .95d0)then
           print *,'Initial value of effective saturation illegally',
     &     ' exceeded 0.95 for node',i,'.  Lower water content slightly'
     &     ,' and restart program.'
           stop 'sso exceeds 0.95 in Section 5. MAIN'
         endif
         bw(i)=bwo(i)  !REJ_2025/05/30  Next 4 needed?
         bl(i)=blo(i)
         ss(i)=sso(i)
         dz(i) = dzo(i)!RJ added 12/13/23 Review exactly what density does
         melt(i)=nmelt(to(i),th(m),tl(m))
         floo(i)=flo(i)
         if(ICVPfixd)then  !11/3/23 RJ Initialize to zero. 
     &     ! should remain 	
           !unbar(i) = 0d0 ARK 7/25/24 experiment
           !uvapor(i) = 0d0 ARK 7/25/24 experiment
           pdzdtc(i) = 0d0  
         endif
 550  continue
      sbt3o=sb*to(n)*to(n)*to(n)
      
      if (.not. USEKCTURB) then ! ARK_2025/1/15
      !  water vapor pressure of air (mbar) 
          eso=fvapri(to(n),1.0d2,e0)
      else
          CALL es_sub_new(eso,volatile,to(n)) 
          eso = eso/100d0 ! convert Pa to mbar
      end if
      
      es=eso
      ct(1)=((cl*flo(1)+ci(1)*(1.-flo(1)))*bwo(1)+bdcd(k(1)))/bt(1)

c***********************************************************************
ct*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t
c     Begin loop over time steps. This is the main loop!  New met
c     parameters are read at the basic time step rate (usually once per
c     hour). Iteration time steps range between dtmin and dtmax seconds
c     and are automatically determined by the code in order that
c     convergence criteria are met.  Met parameters are linearly
c     interpolated between past and current read-in values.
ct*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t
c***********************************************************************
c      write(45,*)'      Input precip for Summit_2007 boosted by x20'
      write(45,*)
     & '      Date-Time  ihour  imin   SWE(m/hr)     Depth(mm/hr)   '
     & ,      '    D(m)       CV' !REJ
      write(45,*)'T',T(n)

1001  icalcstep=icalcstep+1
      IF(iter .eq. 0)THEN
         ibasestep=ibasestep+1
c***********************************************************************
c  6a. Read in met data
c      At start of precip event, istart (passed in arrays) is set to 1
c      Reset to 0 in old.f
c***********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'6a VARIABLE = ',TOPFLUXV 
c      REJ: see above idebug/debug to use this debug feature.
c      Replace TOPFLUXV with variable you want to trace, also upper case.
c      Do a global replace, making sure that it is case sensitive.   
c***********************************************************************
         call getmet(*999,imetcalc,isolarcalc,islope,iy,jday,ihour,imin,
c     & ibasestep,ido,iptype,ircalc,itimezone,ititle,itm,new_sno, !REJ_2025/06/15
     & ibasestep,ido,iptype,ircalc,itimezone,ititle,itm,snow_age, !REJ_2025/06/15
     & a1snow,azslope,bifall,bifallin,bp,clearness,dlatt,dlongt,
     & dsnowfall,effceiling,elev,tlsnow,tm,tmsg,albsnow,bwfall,blfall,
     & exp1,ICVPfixd,time0) !REJ_2025_July added time0
      ENDIF
c    REJ_2025/07/31 Tprecip set to wetbulb T in GETMET. Passed thru common
c    REJ_2025/08/19 Boosting precip for Summit moved to GETMET

      ! Calculate cosine of solar zenith angle ! ARK 12/27/24
c         call zenith(mu_not,dlatt,azslope,century,ihour,iy,jday,planet,
c     & elev)
         mu_not =1D0
c       REJ_2025  Somewhere before this I declared dsnowfall as 1d-3.
c       I need to find this and locate it. But until then declare:
         dsnowfall = 120d-6 ! Fresh snow diameter ARK_2025/09/11
         if(usemeasuredbottomtemp) To(1) =tmsg ! ARK_2025/10/06

c         istart_eff = 0 ! ARK_2025/09/11
c***********************************************************************
c
c  6b. Subdivide basic time interval and interpolate read-in met values
c     Initial Flag =1 for first calc step and then set to 0
c     If new node initialized in this iteration, use minimum time step 
c     dtmin
c***********************************************************************
c123   if(istart.eq.1 .or. fullnode) dt=dtmin !REJ_2025 currently = 5 sec
c123   if(istart.eq.1 .or. fullnode) dt=5d0 !REJ_2025/09/04 currently

123   if((istart.eq.1 .and. snowrate .gt. snow_threshold) .or. fullnode)
     &            dt=5d0 ! ARK_2025/09/16 


cREJ      if(imin .le. 0)call subtime(repeat)  !12/26/24. This caused a problem in Tkair??
      if(initial .le. 0)call subtime(repeat,timesum,tkairsave,icalcstep)!REJ_2025/06/10
      time = dfloat((jday-1)*24) + dfloat(ihour) + (dt/3600d0) 
     &  + dfloat(imin)/60d0- time0

c        write(*,*)'MAIN',iy,jday,ihour,imin,dt,icalcstep,prcp ! ARK_2025/09/16

      if (icalcstep .eq. 1) then 
            prev_iy = iy ! store first year ! ARK_2025/09/16
            write(*,*)'Year = ',iy
      endif
      if (iy .gt. prev_iy) then ! ARK_2025/09/16
             write(*,*)'Previous year = ',prev_iy
             write(*,*)'Next year = ',iy
             prev_iy=iy
      endif       
c
c  Next can be re-instated to see if arrays are being reset properly
c  You must also instate the DEBUG section within reset.f	
c  Set the flux ouput option to 0, since this file is used as a
c  scratch file.
c     m=i+1
c     do 241 iar=1,narray1
c      write(7,*)(array1(i,iar),i=1,n)
c241    continue
c     do 257 iar=1,narray2
c      write(7,*)(array2(i,iar),i=1,n)
c257    continue
c     do 261 iar=1,niarray
c     write(7,*)(iarray(i,iar),i=1,n)
c261   continue
c      rewind 7
c***********************************************************************
      if(istart .ne. 1)snow_age(n) = snow_age(n) + dt/3600d0
c      write(115,*)'MAIN 6b snow_age(n)',snow_age(n)
c***********************************************************************
c  7.  Snowfall, sleet or rain ponding on impermeable top element
c***********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'7. VARIABLE = ',TOPFLUXV
c***********************************************************************
c  If the top node is impermeable and rainfall occurs, initiate a new
c  top node.
c      ddzdtp(n) = 0d0 !REJ_2025/05/09 Initialized 2 lines now in newsnow?  RECHECK
c      us(n) = 0d0  RECHECK
      if(istart .eq. 1)
     &   snow_age(n) = dtmin/3600d0 !Age of CV in hours
     &       !REJ_2025/06/15 Replaces 'newso'
c  Ponding water branch
      if((impermeable(n).eq.1.or.sso(n).gt.1d0).and.rainrate.gt.0d0)then
c     Add ponding water to an impermeable or water saturated node
c     Ponding is not implemented in SNTHERM.89
         call newsnow(rainrate,rainrateo,1d-2,1d3,ssminsnow,1d-3,
     1     repeat,0d0,1d3,fullnode,a1snow,thinnode)!REJ_2025/08/31
     
      end if
c Snowfall branch
c      ICVPfixd = .false. !REJ_2025_05/08. Redundant,  Declared above
      if(snowrate.gt.0.0.and. ICVPfixd)then!REJ_2025/04/29
      
         stop 'No snowfall in the ICVPfixd version'  
      
      elseif(snowrate .gt. 0.0)then ! Add new CV with index n+1 ! REJ 11/10/23
         
            !prec_buf =  snowrate * dt   ! [m snow depth]  ARK_2025/09/11        
         

         
            if(snowrate .lt. snow_threshold)then ! ARK_2025/09/16       
      ! Don't add a new layer for very low snowfall rates

      ! Add thickness of snowfall to top layer
                  mass1  = bwo(n)*dzo(n) ! original mass [kg]
                  mass2  = snowrate*dt*rhowater ! mass of snowfall [kg]

                  dz(n)  =  dzo(n) + mass2/bwfall ! updated thickness [m]

                  bwo(n) = (mass1 + mass2)/dz(n) ! updated density [m]
                  dzo(n) = dz(n)
                  bw(n)  = bwo(n)
      
      ! Mass-weighted average of new grain diameter                 
                  SSA1   = 6d0/(dice*do(n)) ! Specific surface area of original [m^2/kg]
                  SSA2   = 6d0/(dice*dsnowfall) ! Specific surface area of snowfall [m^2/kg]
                  do(n)  = 6d0/(917d0*((mass1*SSA1 + mass2*SSA2)
     &                                 /(mass1 + mass2)))  
                  d(n)   = do(n)

                  istart = 0 ! treat the rest of the code as if no snowfall occurred
            else

      ! Add a new layer for moderate-to-high snowfall rates                 

            call newsnow(snowrate,snowrateo,dzinc,bwfall,ssminsnow 
c     1   ,dsnowfall,repeat,ssisnow,blfall,fullnode,a1snow)REJ2025/05/24
     1           ,dsnowfall,repeat,ssisnow,blfall,fullnode,a1snow !REJ2025/05/24
     2           ,thinnode)!REJ2025/06/03
     
          
            endif
                 
            
          !dz(n) = 0d0  !REJ-2025/05/25 This no longer needed? ! ARK_2025/09/11
         
      endif

c  Nodal array size check
      if(n .gt. nd-1)stop 'Increase nodal array size nd in "const"' 
c**********************************************************************
c     8.  Compaction rate for snow
c     IMPORTANT NOTE: COMPACT CURRENTLY ONLY FOR SNOW. NOT FIRN/ICE
c     ALSO, Compact only runs for the first iteration in a time step.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'8. VARIABLE = ',TOPFLUXV
c**********************************************************************
c         ICVPfixd = .false.  !REJ REMOVE sec 8. Uncomment to activate
c  REJ_2025/08/05 Above now redundant. Declared above

c  Natural compaction and metamorphosis.  The compaction rate
c  is recalculated for every new basestep or continuously when the snow
c  is wet or during new snowfall and for 72 hours afterwards.
c  **Note that the rate is negative and fractional. Thus, the change 
c  in nodalthickness is pdzdtc[s-1]*dz[m]*dt[s] !REJ_2025/06/02
c  Wet snow copaction needs revision.  Meteo France may have an 
c  algorithm for this.

      IF(nsoil .lt. n)THEN !REJ 6/23/24
      iwet=0
      do 370 i=nsoil+1,n
        if(flo(i) .gt. 1d-3)iwet=1
370   continue
      ENDIF
      IF(.not. ICVPfixd .and. n .gt. nsoil)THEN !REJ 11/3/2023 6/22/24      
       if(iter .eq. 1 .or. iwet .eq. 1 .or. newsno .ge. 1 .or.
     &    snow_age(n) .lt. 48d0)then 
c         pdzdtc(i) = 0d0  !Initialize REJ_2025/05/17 undid since
c           compact is not called for each iteration. 
         overburden=dmass(n)/2d0 !REJ_2025/05/09. Surface node
         do 377 i=n,nsoil+1,-1
           if(i .ge. n-1)overburden = overburden+dmass(i) !REJ_2025/05/09
c           if(ltype(k(i)) .eq.1 .or. ltype(k(i)) .eq. 90) !7/14/24 REJ
           if(H2Olayer(i))  !REJ_2025/01/20
     &     call compact(bi(i),to(i),blo(i),overburden,pdzdtc(i),ss(i),
     &     dice,bwo(i),flo(i),floo(i),dto,unbar(i),dzo(i),dzoo(i),
cREJ_2025/06/09     &     melt(i),dmlimit,eta0,binew(i),i,n,ICVPfixd,icalcstep)
     &     melt(i),dmlimit,eta0,bifallin,i,n,ICVPfixd,icalcstep)
 377     continue
       endif
      ENDIF    
 
c  Phase change related density changes are treated in Sec. 14.
c
c REJ_2025/08/05  ICVPfixd = .true. !REJ uncomment if above set .false.
c**********************************************************************
c  9.  Optical parameters and solar extinction
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'9. VARIABLE = ',TOPFLUXV
c**********************************************************************
c  Calculate dsol if there is solar radiation, otherwise dsol = 0d0.
      !sdown= 0d0 ! ARK NO SUN TEST 12/14/24

45    IF(sdown .gt. 0d0 .and. H2Olayer(n))THEN 
c       Temporary fix
c        if(istart .eq. 1)stop 'main sec. 9'

!     ARK_2025/09/16 Commented out next if-statement
c        if(istart .eq. 1)then
c         write(*,*)'hit precip(istart .eq. 1. tempFix)in sec. 9: n=',n
c         dz(n) = -us(n)*dt/bwfall; bw(n) = bwfall
c         write(*,*)'MAIN:dz(n),bw(n),d(n)',dz(n),bw(n),d(n),us(n),
c     &    bwfall
c        elseif(dabs(us(n)) .gt. 0d0)then  ! WORK IN VAPOR
c         dz(n) = dzo(n)-us(n)*dt/bwfall ; bw(n) = bwo(n)
c        else
c         continue
c        endif

        if(istart .eq. 1) then ! ARK_2025/09/17
            write(*,*)'MAIN:dz(n),bw(n),d(n)',dz(n),bw(n),d(n),us(n)
        endif

        if(dz(n) .le. 0d0 .or. bw(n) .le. 0d0 .or. d(n) .le. 0d0) then ! ARK_2025/09/17
            write(*,*)'MAIN:dz(n),bw(n),d(n)',dz(n),bw(n),d(n),us(n)
        endif
c       Next c'd out for now      
c        if(istart .eq. 1)do_Mie = .false.  !REJ_2025/06/27 Use simpler form for new CV
        If(do_Mie)Then
c*m*m*m*m*m*m*m*m*m*m*m*m*m*m* START *m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m  

         If(.not. MLS_run_once .or. .not. ICVPfixd)Then  
            MLS_run_once = .true.  
            MLS_run_once = .false. !REJ_2025/05/07 Forces execution of
c                                   Mie routine for all timesteps

c          Calculate optical properties for each layer from Mie Theory.
c          Wavelenth loop is within the subroutine.

           j = 1 !set firstwavelength. 
            do 432 i = n,1,-1 ! CV loop. 
c           do 432 i = 1,n  !  CV loop.  REJ_2025/05/22
c          Call opt_properties(ltype(k(i)),(d(i))/2d0, rbbl, bw(i), r_d, ! ARK_2025/2/14
           Call opt_properties(ltype(k(i)),(do(i))/2d0, rbbl, bw(i), ! REJ_2025/03/09
     &      r_d,rho_d,dz(i),concen(i),dust_mie_dir,mie_dir,use_SSA, ! REJ_2025/03/09
     &      porosity_rad,ssArea,tau(j,i),omega(j,i),g(j,i),i)!REJ_2025/03/04 Added i for testing
432       continue
c           write(*,*)(iw,omega(iw,n),iw=1,281)
c           stop 'opt-P'
c       Computes spectral fluxes for all interfaces

cREJ_2025_04_12 Call netsolarfluxes(dsol,albedo,tau,omega,g,R_sfc,lyr_typ,
           Call netsolarfluxes_truncat(dsol,albedo,tau,omega,g,R_sfc,
c     &     mu_not,Fs0,rad(1),n,albsnow,iy,jday,ihour,iter,solartest, !REJ Aug 2024
     &     lyr_typ,mu_not,rad(1),n,albsnow,iy,jday,ihour,iter,
     &     solartest,fdirup,fdifup,fdirdn,fdifdn,dlambda)
         Endif

c      Computes broadband arrays 
         call total_solar(n,iy,jday,ihour,fdirup,fdifup,fdirdn,
c     &   fdifdn,rad(1),Fs0,Fs,mu_not,PI,dsol,albsnow,dlambda,        !REJ Aug 2024
     &   fdifdn,rad(1),mu_not,PI,dsol,albsnow,dlambda,
c     &   solartest,spinup) !REJ_2025/03/13
     &   solartest,spinup,icalcstep)
      
c        Note: In sntherm rad(3) = sdown,sup,dirdown
         if(dabs(rad(2)-9999d0) .lt. 1d0 .or. rad(2) .lt. 1d-4)then
           sup   = albsnow*sdown
           solar = sdown-sup   !Note: solar = rad(1)*(1d0-albsnow)
         else 
           write(*,*)'Sdown and Sup are both measured.' 
           write(*,*)'Need to adjust dsol for the Mie case'   
           stop 'MAIN Section 9' 
         endif

c*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m*m END m*m*m*m*m*m*m*m*m*m*m*m*m*m*m

        Else 
c        sdsol needs to be updated following implementation of
c        MLS routine
         solar = sdown*(1d0-albsnow) !REJ_2025/06/27
         bextnir = 60d0
         nsoil = 0
         call sdsol(dsol,d,dmass,fext,botnode,nsoil,n,solar,bextnir)
cREJ         Call sdsol(dsol,d,bw,dz,botdsol,nsoil,n,solar,frvis
cREJ     &    ,jday,ihour,depth,ltype,k) 
        Endif
c      Reset do_Mie to .true.
        do_Mie = .true.
      ELSEIF(solar .le. 0d0)THEN
        do 112 i=1,n+1
           dsol(i)=0d0
112     continue
      ELSE   !Soil
         dsol(n)=solar  ! All solar absorbed in top soil node
         dsol(n+1)=0d0  !??
         do 113 i=1,n-1
           dsol(i)=0d0
113      continue
      ENDIF
C***********************************************************************
c  10.  Top fluxes and stability functions. Option QTURB or KCTURB.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'10. VARIABLE = ',TOPFLUXV
c      REJ note: thinnode is now actually a minimum mass (80d0*dzmin)
c      required for a (surface) CV before energy balance is vald.
c      Now updated in Section 16. dzmin currently set to 2mm
c***********************************************************************
c   REJ comment:  Recheck integration of this section for initial timestep
C     
c      USEKCTURB = .false. ! Use 2001 QTURB method ! ARK_2025/09/20 Set in the beginning already
      IF (.not. USEKCTURB) THEN !ARK_2025/1/13
       !Update vapor pressure in air (mbar)
       !Use fvapri if rh is defined relative to ice saturation, fvaprw
       !if relative to water saturation. Depends on Instrumentation
c       ea=fvapri(tkair,rh,e0)
        ea=fvaprw(tkair,rh,e0)

        m=k(n) ! m = the layer type: snow, soil, etc.

c    First estimate new surface temperature, for a better computation
c    of qsen and qlat coefficients. Added in Nov, 1996.
c    Note: ARK_2025/1/15 qsen and qlat are not used by KCTURB 
       topfluxk = 0d0; topfluxv = 0d0; topflux=0d0 !REJ_2025/09/06      
       If(icalcstep.le.1 .or. melt(n).eq.1 .or. dabs(u(n+1)).gt.0d0
     &  .or. thinnode)Then !Skip Tsurf estimate
         if(thinnode .and. prcp .gt. 0d0)then           
           if(istart .eq. 1)then
            TsurfEst=(0.2*Tprecip+0.8d0*To(n-1))
            To(n)   = TsurfEst
           else
c            TsurfEst=(-us(n)*dt*Tprecip + dmasso(n)*To(n))/
c     &               (-us(n)*dt + dmasso(n)) !REJ Tried and abandoned this
            TsurfEst = To(n)
           endif
         else
           TsurfEst=To(n)     
         endif
       Else ! Proceed with TsurfEst estimate.
         ! First update qsen and qlat for wind speed and add windless exchange  
         ! coefficients, csk(m) for sensible and ck(m) latent heat. 
         ! Bulk transfer coefficients Ch and Ce are not updated yet.
         topfluxk = 0d0; topfluxv=0d0  !REJ_2025/09/06
         if(wspo .ne. 0d0)then  
            qsen=(qsen-csk(m))*wsp/wspo + csk(m)
            qlat=(qlat-ck(m))*wsp/wspo + ck(m)
         endif

         !Compute linear estimate for incoming surface heat flux as a f[T(n)]

         ! ARK 9/6/24 correction. The surface emissivity, em(m), 
         ! should not multiply dirdown (downward longwave)
      
c        Next if-then block added by !REJ_2025/08/03
         if(u(n+1) .lt. 0d0)then      !rainfall
           convect = -cl*u(n+1) !Recheck cl and the n+1 index
         elseif(us(n) .lt. 0d0)then   !snowfall
           convect = -ci(n)*us(n)
         else                         !No precip
           convect = 0d0
         endif 

c         topfluxk=(em(m)*dirdown+qsen*tkair+qlat*(ea+frh(m)*eso*
c     &    21.452d0)+3d0*em(k(n)) *sbt3o*to(n))/2d0 ! ARK 9/6/24
c     &    REJ_2025/09/08 added -To(n) to convect term 3 lines below
c         topfluxk=(dirdown+qsen*tkair+qlat*(ea+frh(m)*eso* ! ARK 9/6/24
c     &    21.452d0)+3d0*em(k(n))*sbt3o*To(n) +
c     &   convect*(Tprecip-To(n)))/2d0 !REJ_2025/08/03 + convect*Tprecip
         topfluxk=(dirdown+qsen*tkair+qlat*(ea+frh(m)*eso* !REJ_2025/09/08
     &    21.452d0)+3d0*em(k(n))*sbt3o*To(n) +
     &   convect*(Tprecip-To(n)) )/2d0

         topfluxv=-(qsen+frh(k(n))*22.452d0*eso*qlat/to(n)+
     &   4d0*em(k(n))*sbt3o)/2d0      !REJ_2025/09/08 removed -convect 
c     &    4d0*em(k(n))*sbt3 - convect )/2d0 !REJ_2025/08/03 -convect !REJ_2025/08/29
c    &   ! was sbt3

         ! Solve surface CV energy eq for Tsurfest
         dum=qs(n)*dto/dt !REJ qs(n)=cxmass/dt. Should be (qso+qs)/2 ?     
         TsurfEst=(dum*To(n)+.5d0*dsol(n)+topfluxk+bbo(n)
     &   -heatfluxbtop)/(dum-topfluxv) 
c-----------------------------------------------------------------------
c optional testing 
         if(icalcstep .ge.999999)then ! precip starts at 560
           write(120,*)'step topo,top,compute',topfluxo,topflux,
     &     topfluxv*Tsurfest,topfluxk,topfluxv,Tsurfest,
     &    'tprecip',tprecip,'convect',convect,'us',us(n),
     &    'convect*Tprecip',convect*(Tprecip-to(n))
c           if(dabs(topfluxk+topfluxv*Tsurfest-topflux).gt. 8d0)
c     &      stop 'topflux error'  !REJ Reinstate this later
         endif
c-----------------------------------------------------------------------
c        Limit temperature change to 5 degree/hour.

         if(Tsurfest .gt. To(n))TsurfEst=dmin1(TsurfEst,To(n)+dt/720.)
         if(Tsurfest .lt. To(n))TsurfEst=dmax1(TsurfEst,To(n)-dt/720.)

         Endif !End TsurfEst
     
         Esurface=fvapri(TsurfEst,1.0d2,e0)
         m=k(n) !REJ_2025/02/21  Redundant? Check

         call qturb(height,Tkair,TsurfEst,Ea,Esurface,wsp,wspo,qlat,
     1   qsen,dliqvol(n),cdryair,Rw,Ck(m),Csk(m),snowdepth,ltype(m),
     2   znaught(m),Cd(m),rhoair,bp,dlogTt(m),dlogQq(m),dlogWw(m),
     3   n,icalcstep) !REJ_2025/05/21 Added n and icalcstep

c       Update topfluxk and topfluxv with updated qsen and qlat.
      
        topfluxk=(dirdown+qsen*tkair+qlat*(ea+frh(k(n))*eso*  !REJ
     &      21.452d0)+3d0*em(k(n))  ! REJ removed em from dirdown
     &     *sbt3o*to(n)+convect*(Tprecip-To(n)))/2d0 !REJ_2025/09/08 
                                            ! REJ added -To(n) above 
        if(n .le.nsoil)
     &   topfluxk=topfluxk+cl*u(n+1)*(Tprecip-to(n))/2d0 !Recheck this

         topfluxv=-(qsen+frh(k(n))*22.452d0*eso*qlat/to(n)+
     &    4d0*em(k(n))*sbt3o )/2d0  !REJ_2025/09/08 removed "- convect" 
 

c         if(Tsurfest .lt. To(n))TsurfEst=dmax1(TsurfEst,To(n)-dt/720.)  !REJ Check this!!
c         REJ_2025 REMOVE above c (=comment)
      
c       End of QTURB Section

      ELSE ! Use KCTURB method.  ARK added the following section ca. 9/6/24

        m=k(n) ! m = the layer type: snow, soil, etc.

c      Compute vapor pressure in air with respect to water(mbar)

       CALL es_sub_new(ea,volatile,tkair) 
       ea = ea*(rh/100d0)/100d0 ! convert Pa to mbar and account for relative humidity
       
      topfluxk = 0d0; topfluxv = 0d0; topflux=0d0 ! ARK_2025/09/22  

c    First compute a new surface temperature estimate TsurfEst

      If(icalcstep.le.1 .or. melt(n).eq.1 .or. dabs(u(n+1)).gt.0d0
     &  .or. thinnode)Then !Skip estimate
        if(thinnode .and. prcp .gt. 0d0)then ! ARK_2025/09/20 Copied from QTURB branch 
           if(istart .eq. 1)then
            TsurfEst=(0.2*Tprecip+0.8d0*To(n-1))
            To(n)   = TsurfEst
           else
c            TsurfEst=(-us(n)*dt*Tprecip + dmasso(n)*To(n))/
c     &               (-us(n)*dt + dmasso(n)) !REJ Tried and abandoned this
            TsurfEst = To(n)
           endif
         else
           TsurfEst=To(n)     
         endif
      Else !Proceed with TsurfEst estimate.'
      
c        Next if-then block added by ! ARK_2025/09/22
         if(u(n+1) .lt. 0d0)then      !rainfall
           convect = -cl*u(n+1) !Recheck cl and the n+1 index
         elseif(us(n) .lt. 0d0)then   !snowfall
           convect = -ci(n)*us(n)
         else                         !No precip
           convect = 0d0
         endif 

c     Calculate some old values (REJ Question. Could the next 3 be moved to Sec. 20?)

      ! calculate old surface specific humidity for use in topfluxk and topfluxv

      CALL q_sub_sntherm(qo,Ma,Mw,100d0*bp,to(n),1d0,eso)
      
      ! calculate old air equilibrium (saturated) vapor pressure
      CALL es_sub_new(eao,volatile,tkairo) 
          eao = eao/100d0 ! convert mbar to Pa
      
      ! calculate old air specific humidity for use in topfluxk and topfluxv    
      CALL q_sub_sntherm(qairo,Ma,Mw,100d0*bp,tkairo,(rh/100d0),eao)

      
      ! calculates derivative of d(q_sat)/dT for use in topfluxv
      CALL qsat_derivative(dqsat,Ma,Mw,100d0*bp,to(n))

C     Compute linear estimate for incoming surface heat flux as a f[T(n)]
      !Incoming surface heat flux estimate = topfluxk + topfluxv*T(n)

C     Here is the KCTURB version of topfluxk ARK_2025/2/3 
         topfluxk =(dirdown + 

      ! replacing the qsen term
     &  (rhocp_rep/(Fm*Fh02) )* ( (tkairo-to(n))*wsp + wspo*to(n) )
     &  + (rho_rep*Cpd0*lambda0*tkair/(Fm*Fw02))  *
     & (wsp *(qairo - frh(k(n))*qo) + wspo * frh(k(n))* dqsat*to(n))   
      ! qsen term ends 

      ! replacing the qlat term
     &    + (rho_rep*Ls/(Fm*Fw02))*( wsp * (qairo - frh(k(n))*qo)
     &    + wspo * frh(k(n)) * dqsat * to(n)) 
      ! qlat term ends

cREJ_2025/09/08     &    + 3d0*em(k(n))*sbt3o*to(n) + convect*Tprecip)/2d0 !REJ_2025/08/05
     &    + 3d0*em(k(n))*sbt3o*To(n) + convect*(Tprecip-To(n)))/2d0
     &   !REJ_2025/09/08 
     
c	   Here is the KCTURB version of topfluxv ARK_2025/2/3                       
          topfluxv = -(
      
      ! replacing the qsen term
     &    rho_rep*Cpd0*lambda0*tkair*frh(k(n))*wspo*dqsat/(Fm*Fw02) 
     &     + rhocp_rep*wspo/(Fm*Fh02)  
      ! qsen term ends
      
      ! replacing the qlat term
     &    + rho_rep*Ls*frh(k(n))*wspo*dqsat/(Fm*Fw02)
      ! qlat term ends

c    &     +4d0*em(k(n))*sbt3o - convect)/2d0 !REJ_2025/08/05 !REJ_2025/09/08
     &     +4d0*em(k(n))*sbt3o)/2d0 !REJ_2025/09/08

      ! Solve surface CV energy eq for Tsurfest
         dum=qs(n)*dto/dt
         
         TsurfEst=(dum*To(n)+.5d0*dsol(n)+topfluxk+bbo(n)
     &   -heatfluxbtop)/(dum-topfluxv)
     
c        Limit temperature change to 5 degree/hour - why?

         if(Tsurfest .gt. To(n))TsurfEst=dmin1(TsurfEst,To(n)+dt/720.)
         if(Tsurfest .lt. To(n))TsurfEst=dmax1(TsurfEst,To(n)-dt/720.)

       Endif  !End TsurfEst

          CALL es_sub_new(esurface,volatile,TsurfEst)  
          esurface = esurface/100d0 

          m=k(n) !ARK_2025/09/22  Redundant? Check
c      Call Khuller & Clow (2024) Turbulence Model 
c REJ_2025/02/20 question.  where are qsen and qlat defined? 
c ARK_2025/09/20 kcturb does not compute or use qlat or qsen, so SNTHERM should always check for .not. usekcturb to use qlat
          CALL kcturb(Cpd0,Fh02,Fw02,H0sen,lambda0,LE,Ls,qlat,qsen,q2,
c REJ_2025/05/12 Not sure when Fm was replaced with ustar. Chgd back     
c     & rho_rep,rhocp_rep,ustar,
     & rho_rep,rhocp_rep,Fm,
     & planet,znaught(m),Zmax,H_PBL,bp,beta,beta2, ! ARK_2025/09/20
     & Fx_opt,0d0,height(1),height(2),tsurfest,tkair,1d0,rh/100d0,
     & wsp,frh(m),ltype(m),snowdepth,dliqvol(n),
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b) 

c      recompute topflukv and topfluxv after updates from call to kcturb

         ! ARK_2025/2/3  
         if(initial .eq. 1)then ! we don't have tkairo,wspo yet
c         REJ Comment.  Move this to start of routine?  First check if TT section skipped
c         for initial?    
              tkairo = tkair
              wspo   = wsp
         end if
      
       ! calculate old surface specific humidity for use in topfluxk and topfluxv
      CALL q_sub_sntherm(qo,Ma,Mw,100d0*bp,to(n),1d0,eso)
      
      ! calculate old air equilibrium (saturated) vapor pressure
      CALL es_sub_new(eao,volatile,tkairo) 
          eao = eao/100d0
      
      ! calculate old air specific humidity for use in topfluxk and topfluxv    
      CALL q_sub_sntherm(qairo,Ma,Mw,100d0*bp,tkairo,(rh/100d0),eao)
      
      ! calculates derivative of d(q_sat)/dT for use in topfluxv
      CALL qsat_derivative(dqsat,Ma,Mw,100d0*bp,to(n))
            
c     Here is the KCTURB version of topfluxk ARK_2025/2/3  
       topfluxk =(dirdown +
     
      ! replacing the qsen term
     &    (rhocp_rep/(Fm*Fh02) )* ( (tkairo-to(n))*wsp + wspo*to(n) )
     &    + (rho_rep*Cpd0*lambda0*tkair/(Fm*Fw02))  *
     &    (wsp *(qairo - frh(k(n))*qo) + wspo * frh(k(n))* dqsat *to(n))   
      ! qsen term ends

      ! replacing the qlat term
     &    + (rho_rep*Ls/(Fm*Fw02))*(wsp*(qairo - frh(k(n))*qo)
     &    + wspo*frh(k(n))*dqsat*to(n)) 
      ! qlat term ends

c    & + 3d0*em(k(n))*sbt3o*to(n) + convect*Tprecip)/2d0
c     &!REJ_2025/08/05 !REJ_2025/09/08
     & + 3d0*em(k(n))*sbt3o*to(n) + convect*(Tprecip-To(n)))/2d0
c    &!REJ_2025/08/05 !REJ_2025/09/08     
      
c     Here is the KCTURB version of topfluxv ARK_2025/2/3                       
          topfluxv = -(
      
      ! replacing the qsen term
     &    rho_rep*Cpd0*lambda0*tkair*frh(k(n))*wspo*dqsat/(Fm*Fw02) 
     &     + rhocp_rep*wspo/(Fm*Fh02)  
      ! qsen term ends
      
      ! replacing the qlat term
     &    + rho_rep*Ls*frh(k(n))*wspo*dqsat/(Fm*Fw02)
      ! qlat term ends
      
c     & +4d0*em(k(n))*sbt3o - convect)/2d0 !REJ_2025/08/05 !REJ_2025/09/08
     & +4d0*em(k(n))*sbt3o)/2d0 !REJ_2025/08/05 !REJ_2025/09/08

cREJ_2025/02/24      end if
      
      if(n .le. nsoil)topfluxk=topfluxk+cl*u(n+1)*to(n)/2d0  !REJ recheck this!

      ENDIF !End of KCTURB Section

c Next added on March 28, 1995. No evaporation for dry soil.
c     if(n .le. nsoil .and. bwo(n) .le. 1d-10)qlat=0d0 ! Bare soil case
c Jan 29,2001.  Shuts down latent heat flux for low of high water content
      if(n .le. nsoil)then ! Bare soil case
         frh(k(n))=1d0
         if(bwo(n).ge.0.95 *1d3*porosity(n)  .and. ea .gt.eso)qlat=0d0
     &    ; LE =0D0    ! ARK_2025/09/20
         if(bwo(n) .le. bdjp(m) + 1d-1 .and. ea .lt. eso) qlat=0d0
     &    ; LE =0D0    ! ARK_2025/09/20
      endif
      if(initial .ge. 1)then
c         dlong= em(m)*(dirdown-sbt3o*to(n)) 
         dlong= dirdown-em(m)*(sbt3o*to(n)) ! ARK 9/6/24
          
          if (.not. USEKCTURB) then ! ARK_2025/1/13 because KCTURB doesn't use qlat or qsen in the same manner
              sensheat=qsen*(tkair-to(n))
              dlatheat=qlat*(ea-frh(m)*eso)
          else
              sensheat = -H0sen  ! sensible heat flux (W/m^2)
              dlatheat = -LE     ! latent heat flux (W/m^2)
         end if
         
         convect=0d0
      endif
      u(n+1)=-1.0d3*rainrate
      
c***********************************************************************
c  11.  Water infiltration from rain or melt
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'11. VARIABLE = ',TOPFLUXV
c***********************************************************************
c  Do this section if there is water flow.
c  Skip if a new node has been initialized during this iteration.
c  Snowcover/soil elements are grouped into nzone blocks where there
c  is contiguous water flow.  Block limits are designated in routine
c  NLIMIT.
c      goto 1111  !skipped infiltration for now ! ARK_2025/09/18
       IF(.not. ICVPfixd)THEN !10/30/23  RJ
       call nlimit(ibounds,nzone)
       if(nzone .gt. 0)then
        do 65 i=1,nzone
         nz=i
         if(ibounds(nz,2).gt.ibounds(nz,1).or.ibounds(nz,1).lt.1)stop
     &   '**Execution halted in MAIN sec 11.  Incorrect node bounds**'
          call filtrate(ibounds(nz,2),ibounds(nz,1),*1001,dzmin,
     &         Neumann)  !added Neumann on 12/10/23 
65      continue
        if(ioutfiltrate.eq.1)call outfiltrate(ititle,ibounds,ihour,
     1  dtsum,nzone,jday)
       endif  
       unbarmax=0d0
       umax=0d0
       do 451 i=1,n
         umax=dmin1(u(i),umax)
         if(sso(i) .lt. 1d0 .or. i.eq. n)then
           unbar(i)=(uo(i+1)-uo(i)+u(i+1)-u(i))/2d0
cTemp.     Next used in drain. Remove later
           if (i .le. nsoil) unbar(i) = 0d0 !10/22/23 REJ
         else
           unbar(i)=u(i+1)-u(i)
cTemp.     Next used in drain. Remove later
           if (i .le. nsoil) unbar(i) = 0d0
cTempend      
           if(dabs(unbar(i)).gt.1d-10)stop '**Net flow for saturated 
     &     element is not zero. Something must be wrong in FILTRATE  
     &     routine**'
         endif
         unbarmax=dmin1(unbar(i),unbarmax)
         if(dabs(unbar(i)) .gt. 0d0)then !Next 4. REJ_25/01/16
c            write(*,*)'jday',jday,'i',i, 'unbar',unbar(i)
c           stop 'unbar(i)'  !Did not hit this.  Thus, no water flow.
         endif
 451  continue
       ENDIF  ! End of skipped filtrate computation !10/30/23  REJ
c1111   continue ! ARK_2025/09/18
c***********************************************************************
c  12.  Sublimation and diffusion of water vapor in snow
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'12. VARIABLE = ',TOPFLUXV
c***********************************************************************
c  Based on old temperatures, calculate mass diffused vapor flux.
c  For surface element, change thickness dzo to correspond to mass
c  evaporative losses/gains.
c  REJ 1/5/2024 Next 3. Changed and changed back on 1/6/2024
c  REJ_2025/08/05        ICVPfixd = .true. !2025 May 18? 
c  REJ_2025/08/05      IF(ICVPfixd)then 
        call diffusion(de0,bp,bvi0,dlsdrw,frh,dlvdrw,bvw0,ICVPfixd,
     &     H2Olayer,icalcstep)!REJ_2025/06/04
c REJ_2025/08/05      ENDIF  
c***********************************************************************
c  12b  Compute new effective grain size for snow or firn ice?
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'12b. VARIABLE = ',TOPFLUXV
c***********************************************************************
c
c      ICVPfixd = .false.  
c REJ_2025/08/05      write(*,*)'ICVPfixd',ICVPfixd
      IF(.not. ICVPfixd .or. n .eq. nsoil)THEN !REJ_2025/05/05 omitted comment  
c                                               Why "n .eq. nsoil"?     
       do 42 i=nsoil+1,n
c         if(ltype(k(i)) .eq. 1 .or. ltype(k(i)) .ge. 90)then !REJ 6/20/24
         if(H2Olayer(i) )then !REJ 6/20/24
           d(i)=fgrain(g1,g2,ufvapor(i),do(i),dliqvol(i),dt)
c REJ3 1/23/24  if(d(i) .gt. 5d-3)d(i)=5d-3  increased max to 5d-2 ??
c REJ_2025/05/05 Omit or change algorithm for ice.
           if(d(i) .gt. 5d-2)d(i)=5d-2 
         endif 
42     continue
      ELSE  !REJ added this branch 6/20/24
       do 43 i = nsoil+1,n
           d(i)=do(i)
43     continue      
      ENDIF 
c REJ_2025/08/05      ICVPfixd = .false.  !REJ REMOVE sec 12b
c***********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'13. VARIABLE = ',TOPFLUXV
c***********************************************************************
c  13. Calculate new values for dmass, bw, bt, dz and expansion/
c  contraction of water saturated elements with phase change.  
c  WE STILL NEEED TO MAKE SOME ADJUSTMENTS FOR ICE LAYERS
c  REJ_2025/05/17 Entire Sec. 13 replaced with call to massbal.f
c***********************************************************************
c      write(115,*)'MAIN 13.jday,ihour,imin,icalstep,dt',
c     &  jday,ihour,imin,icalcstep,dt

      if(.not.ICVPfixd .or. prcp .gt. 0d0)
     & Call massbal(icalcstep,ICVPfixd,dtmin,dzmin,bwfall,sum_vapor_dt,
     &  Totaltime,sum_delta_dz,thinnode)!REJ_2025/07/22 & 24 & 08/27
     

          if(icvpfixd .and. dabs(unbar(n)) .gt.  1d-10)then
            write(*,*)'UNBAR:STEP',unbar(n),jday,ihour,icalcstep
            stop 'MAIN: Section 13'
          endif

          Totalmass = Totalmass -uvapor(42)*dt   !REJ testing: remove later

c      Set a minimum mass
          thinnode = .false.  !REJ_2025/08/30
c         if(dz(n) .lt. dzmin .or. dmass(n) .lt. 80d0*dzmin) !REJ_2025/08/30
          if(dmass(n) .lt. 80d0*dzmin) ! changed limit to just mass.
     &    thinnode = .true. !REJ_2025/08/30

c        if(prcp .gt. 0d0) !A comparison of delta(mass(n)) with snowfall
c     &  write(120,*)'MB  bw(n)*dz(n)',bw(n)*dz(n),'mass(n)',dmass(n),
c     &  dmass(n)-dmasso(n),-us(n)*dt
        
c**********************************************************************
c**********************************************************************
c*****    END MASS BALANCE   ***      START THERMAL BALANCE       *****
c______________________________________________________________________
c  IMPORTANT PROGRAMMING NOTE:  Only the 'old' designations should be
c  used for variables BWO and DZO prior to Sec. 13 and for variables TO
c  and BLO prior to the final thermal balance.
c***********************************************************************
c  14.  Calculation  of thermal parameters
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'14. VARIABLE = ',SNOWRATE
c***********************************************************************
      H2Olayer(n) = .true.
      do 37 i=1,n
          if(H2Olayer(i))then  !REJ 6/24/24
c  Effective thermal conductivity of snow, including vapor diffusion.
c  RJ note:  Vapor contribution is 0d0 for now.  Dec 7, 2023 Update on
c   June 29, 2025.  Need to check current state of this.
c        thkice=2.290d0  ! And REJ turned this off 1/5/2024
         thkice=780d0/to(i)-0.615d0 ! REJ Turned back on temp dependence 
         thk(i) = thkair+(7.75d-5 *bw(i)+ 1.105d-6*bw(i)*bw(i))*(thkice
     1   -thkair)
        if(dliqvol(i) .le.0.02)then
          thk(i)=thk(i)+dls*df(i)
        else
          thk(i)=thk(i)+dlv*df(i)
        endif
      else
c  thermal conductivity for soil.  recalculate for unfrozen soil only
c  when the water content changes.
         m=k(i)
         if(ibasestep .eq. 1)goto 557  !RECHECK THIS!! REJ 6/21/24
         if(dabs(bwo(i)- bw(i)).lt.dtol1 .and. to(i).gt.273.15)goto 37
 557     thk(i)=thrk(to(i),dicevol(i),dliqvol(i),dkmg(m),
     1        icoarse(m),dkdry(m),dksatu(m),porosity(i),flo(i))
      endif
 37   continue

c RESTORE 08/03/25     call thparam(Neumann,357,18,prnt_hour) 
      call thparam(Neumann,357,18,prnt_hour,120,icalcstep)!REJ_2025/08/10
c      write(120,*)'Aft THP',icalcstep,bt(n),ct(n),dz(n),
c     &     bb(n),ci(n)      
c***********************************************************************
c  15.  Do next if this is the initial iteration
c        Skip energy balance and go to print-out.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'15. VARIABLE = ',TOPFLUXV
c***********************************************************************
      if(initial .eq. 1)then
c         hg is the basal upward thermal from soil to snow 
         if(nsoil .eq. 0)then
           hg = 0d0 !In keeping with the 0d0 bottom flux for Khuller project
c          However, will not be 0d0 for shallow stack
         elseif(nsoil.lt.n)then 
          hg=-qk(nsoil+1)*(to(nsoil+1)-to(nsoil)) 
         endif              
         call old(ido,topflux)

c         write(120,*)'bef initial fbb:topfluxk,topfluxv,to(n),EST', !T_REMOVE
c     &      topfluxk,topfluxv,to(n), topfluxk+topfluxv*to(n) !T_REMOVE

         call fbb(qk,to,dsolo,qf,topfluxk+topfluxv*to(n),botflux,bbo,n, ! RJ 
     &      nsoil,jday,ihour,Neumann,prnt_hour,H2Olayer,heatfluxbtop, ! REJ_2025/02/09 
     &      icalcstep,bbo,120) !! REJ_2025/08/17+09/04
c         call fbb(qk,to,dsolo,qf,0.5d0*(topfluxk+topfluxv*to(n)),botflux,
c     &      ,bbo,n, ! REJ_2025/08/17 Changed back 08/29
cREJ_2025/02/09     &      nsoil,jday,ihour,Neumann,prnt_hour,H2Olayer) ! RJ12/17+6/24/24

c         write(120,*)'MAIN Aft fbb call bb,bbo',bb(n),bbo(n)

         dtsum=0d0
         print= .true.
         iter=0
         topfluxo=topfluxk+to(n)*topfluxv

c July 10 1996 - added Bert's flux format (ceretha)
           
         if(ifluxout .ne. 0)call flux(jday,ihour,ibasestep,iptype,
     &    height,iy,clearness,imin,isolarcalc,effceiling,islope,
     &    ifluxout)
          
         initial=0
         goto 1000  ! Print-out for this basestep
      end if

c***********************************************************************
c  16.  Thermal balance
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'16. VARIABLE = ',TOPFLUXV
c***********************************************************************
      call thermal(dtmin,dzmin,recross,*1001,errtallowd,Neumann,
     &     357,18,prnt_hour,icalcstep,120,thinnode,TsurfEst)!10/22/2023 RJ12/17!REJ_2025/08/11
c      write(*,*)'bw(n) aft thermal',icalcstep,bw(n),ci(n),dz(n),bb(n)
     
c  Calculate top fluxes using newly calculated top node temperature
      m=k(n)
      
      if (.not. USEKCTURB) then ! ARK_2025/1/13
          es=fvapri(t(n),1.0d2,e0)
          
          call qturb(height,Tkair,T(n),Ea,es,wsp,wspo,qlat,
     1 qsen,dliqvol(n),cdryair,Rw,Ck(m),Csk(m),snowdepth,ltype(m),
     2 znaught(m),Cd(m),rhoair,bp,dlogTt(m),dlogQq(m),dlogWw(m),n,
     3 icalcstep) ! REJ_2025/05/21 Added n and icalcstep
      
      sbt3=sb*t(n)*t(n)*t(n)
      
      else
      
          CALL es_sub_new(es,volatile,t(n))
          es = es/100d0

c      Call Khuller & Clow (2024) Turbulence Model  

 
      CALL kcturb(Cpd0,Fh02,Fw02,H0sen,lambda0,LE,Ls,qlat,qsen,q2,
     & rho_rep,rhocp_rep,Fm, ! ARK_2025/09/19
     & planet,znaught(m),Zmax,H_PBL,bp,beta,beta2, ! ARK_2025/09/20
     & Fx_opt,0d0,height(1),height(2),T(n),tkair,1d0,rh/100d0, ! ARK_2025/09/20
     & wsp,frh(m),ltype(m),snowdepth,dliqvol(n),
     & md_Pc_a,md_Pc_b,md_Tc_a,md_Tc_b,md_Vc_a,md_Vc_b,md_a0,md_a1,
     & md_a2,md_a3,md_a4,md_b0,md_b1,md_b2,md_b3,md_b4,md_beta_a,
     & md_beta_b,md_gamma_a,md_gamma_b,md_kappa_a,md_kappa_b,
     & md_m_a,md_m_ab,md_m_b,md_mu_a,md_mu_b,md_omega_a,md_omega_b,
     & md_sum_v_a,md_sum_v_b)

      endif
      sbt3=sb*t(n)*t(n)*t(n) !REJ_2025/05/16 added
      dlong= dirdown-em(m)*(sbt3*t(n)) ! ARK 9/6/24
      
      if(n .le. nsoil)then
         if(bwo(n).ge.0.95 *1d3*porosity(n)  .and. ea .gt.eso)then
         qlat=0d0; LE = 0d0 ! ARK_2025/09/20 
           endif
         if(bwo(n) .le. bdjp(m)+ 1d-1.and. ea .lt. eso) qlat=0d0
     &   ; LE = 0d0 ! ARK_2025/09/20 
      endif
      
      if(nsoil .eq. n)then
        hg = 0d0  ! See comment above
      elseif(nsoil .lt.n)then
        hg=-qk(nsoil+1)*(t(nsoil+1)-t(nsoil))
      endif
      
      if(.not. USEKCTURB)then ! ARK_2025/1/13     
       sensheat=qsen*(tkair-t(n))
       dlatheat=qlat*(ea-frh(m)*es) 
      else
       sensheat = -H0sen
       dlatheat = -LE     
      endif 

      dlong= dirdown-em(m)*sbt3*t(n)! REJ_2025/08/17   
  
      convect = -ci(n)*us(n)!REJ_2025/08/02 REJ_2025/08/02.  Redo this when rain is back in
c  Out for now  if(n .le. nsoil)convect=convect+cl*u(n+1)*t(n)!REJ6/23 Rethink this
cREJ_2025/09/08 topflux=(dlong+sensheat+dlatheat+convect*(Tprecip-T(n)))/2d0 !REJ_2025/08/17
      topflux=(dlong+sensheat+dlatheat+convect*(Tprecip-To(n)))/2d0 !REJ_2025/09/08
     &    !REJ_2025/09/09
      if(icalcstep .ge. istep_conv_out)
     &   write(120,*)'3 bef FBB: topflux & Est',
     &                       topflux,topfluxk +topfluxv*T(n)!REJ_2025/09/04                   
c      write(*,*)' bef FBB',icalcstep,bw(n),ci(n),dz(n),bb(n)
c-----------------------------------------------------------------------------
c OK to delete this block 
c REJ_2025/08/17      call fbb(qk,t,dsol,qf,topfluxo,botflux,bb,n,nsoil, ! RJ 11/03/23 
cREJ_2025/08/28      call fbb(qk,t,dsol,qf,0.5d0*topflux,botflux,bb,n,nsoil, ! RJ 11/03/23+REJ/08/17
c      call fbb(qk,t,dsol,qf,0.5d0*topflux,botflux,bb,n,nsoil, ! RJ 11/03/23+REJ/08/17
c     &      jday,ihour,Neumann,prnt_hour,H2Olayer)  ! RJ 11/03/23 
c-----------------------------------------------------------------------------

      call fbb(qk,t,dsol,qf,topflux,botflux,bb,n,nsoil, ! RJ 11/03/23+REJ/08/17
     &      Jday,ihour,Neumann,prnt_hour,H2Olayer,heatfluxbtop, !REJ_2025/02/09
     &      icalcstep,bbo,120) !REJ_2025/08/17 Changed back on 08/29/2025;+ 09/04
      
      if(istart .eq. 1 .and. dabs(bbo(n)) .le. 0d0)bbo(n)=bb(n)
      
c-----------------------------------------------------------------------------
c Ditto this block
c      write(*,*)'AFTER FBB',icalcstep,
c     & 'topfluxo,topflux,EQ',topfluxo,topflux,topfluxk+topfluxv*T(n),
c     &   topfluxk,topfluxv,'no precip',dlong+sensheat+dlatheat,'precip',
c     &   convect*us(n)*(tprecip-T(n))
c       if(prcp .gt. 0d0)write(*,*)'bb(n),bbo(n)',
c     &    bb(n),bbo(n)
c       if(prcp .gt. 0d0)write(120,*)'bb(n),bbo(n)',bb(n),bbo(n)
c      write(180,142)jday,ihour,Tsurfest,T(n),To(n),Tsurfest-T(n) !REJ_2025  Testing
c      if(melt(n) .gt. 0)stop 'melt:  Main' !REJ_2025/02/13  melt(n) at Tsurfest
c142   format(2i6,4f10.4)
c-----------------------------------------------------------------------------
       Total_Latent = Total_Latent + (dlatheat*dt/dls) !REJ_2025 testing
c       write(*,*)'total_latent',Total_Latent  !REJ_2025 ditto
C***********************************************************************
c  17.  Final adjustments.  Adjust liquid water fraction after thermal
c       balance, determine melt state and element midpoint position.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'17. VARIABLE = ',TOPFLUXV
c***********************************************************************
      qf(n+1)=.5*cl*u(n+1)
      unbarmax=0d0
      umax=0d0
      do 40 i=1,n
         m=k(i)
      if(i .gt. nsoil)then !Don't see why flo for soil is excluded here??
        floo(i)=flo(i)
        dzoo(i)=dzo(i)
      endif
c      if(i .eq. n)write(*,*)'bw(n) bef density',icalcstep,bw(n),ci(n),
c     &   dz(n),bb(n)
         call density(t(i),bw(i),td(i),td13(i),flo(i),bl(i),bi(i),bt(i),
     1        dmass(i),bdjp(m),a243(m),bd(m),a1(m),melt(i),dz(i),
     2        flgo(i),ss(i),dice,ssi(i),porosity(i),dmvol(m),
     &        ltype(m),impermeable(i),idelete(i),solidporosity(i), !REJ 6/24/24
     &        ipond(i),dicevol(i),dliqvol(i),rhowater,H2Olayer(i)) !REJ 6/24/24
c      if(i .eq. n)write(*,*)'bw(n) aft density',icalcstep,bw(n),ci(n),
c     &   dz(n),bb(n)
         if(idelete(i) .eq. 1)then
             if(i .gt. 1)iskip(i-1)=1
             iskip(i)=1
         endif
         ci(i)=-13.3+7.8*t(i)
         melt(i)=nmelt(t(i),th(m),tl(m))
         umax=dmin1(u(i),umax)
         unbar(i)=(uo(i+1)-uo(i)+u(i+1)-u(i))/2d0
c Next added on Jan 7, 2002
         if(i .le. nsoil)unbar(i)=0d0 !no waterflow in soil
         unbarmax=dmin1(unbar(i),unbarmax)
 40   continue
      if(thinnode)melt(n)=0
c      write(*,*)'dz(n),icalcstep',dz(n),icalcstep
      call fz(z,dz,snowdepth,nsoil,n)  
C***********************************************************************
c  18.  Check if convergence criteria are met. Estimate next time step.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'18. VARIABLE = ',TOPFLUXV
C***********************************************************************
c      write(*,*)'bw(n) bef converge',bw(n),ci(n),dz(n),bb(n)
      call converge(repeat,igood,ngoodmin,dtsum,dssallowed,
c     1     errtallowd,dtmin,dtmax,dtsmin,dtssmax,print) 11/05/23/  RJ
     1     errtallowd,dtmin,dtmax,dtsmin,dtssmax,print,Neumann !REJ_2025/06/23
     2     ,icalcstep,120) !REJ_2025/06/23
c      if(icalcstep .gt. 800)write(*,*)icalcstep,'dt = ',dt
c      if(dt .lt. 0.04d0)stop 'MAIN after converge'
      if(dt .lt. 0.01d0)stop 'MAIN after converge'  !REJ_2025/08/07
c***********************************************************************
c  19.  Optionally print-out summary of met conditions and top fluxes.
c***********************************************************************
c July 10 1996 - added Bert's flux format (ceretha)
      if((iter.eq.0 .or. print) .and. ifluxout .ne. 0)call flux
     &(jday,ihour,ibasestep,iptype,height,iy,clearness,imin,isolarcalc,
     & effceiling,islope,ifluxout)
c***********************************************************************
c  20.  Establish old values for mass, thermal and met parameters.
c       Or if convergence criteria not met, repeat past iteration.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'20. VARIABLE = ',TOPFLUXV
c***********************************************************************
      if(.not. repeat)then
c Sept 9 1996 - added water output variables (ceretha)
       uosave=uo(nsoil+1) !REJ. Right now, don't see how this fits in
c                         6/23/24. See line labeled 1000
cREJ_2025/08/17         call old(ido)
         call old(ido,topflux) !REJ_2025/08/17 
         totaltime=totaltime+dt
      else
         goto 1001 ! Return to start of MAIN LOOP
      endif
c***********************************************************************
c  21.  Combine or divide thin or thick snow elements
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'21. VARIABLE = ',ICVPfixd
C***********************************************************************
c
c  After precip has stopped, divide thick snow or ponded water elemnts
      if((prcp .le.0d0.and.(dzo(n).gt.dzn .or. dzo(n-1) .gt. dznm)))then
        nold=n
c April 4, 1995      call subdivide(dzn,dznm,a1snow,dzinc)
c    NOTE: Several debug lines removed here  RJ12/17
        call subdivide(dzn,dznm,a1snow,dzinc,floo,45,icalcstep)
     &     !REJ_2025/08/20  dt=1d0
CR           if(icalcstep .gt. 8000)stop 'MAIN after combo'
      end if
c Unless it is snowing or water is ponding on top, combine thin elements 
      call combinenodes(dzmin,a1snow,thsnow,tlsnow,ssisnow,icalcstep,
cREJ_2025/08/20     &140) !REJ_2025/08/04 Added iunit 140
     &45) !REJ_2025/08/04 Added iunit 140
c      if(icalcstep .gt. 100000)stop 'MAIN after combo' ! ARK_2025/09/10

c  Reset idelete and iskip flags
      do 3000 i=1,n+1
        iskipo(i)=iskip(i)
        iskip(i)=0
        idelete(i)=0
3000  continue
      nold=n 
C***********************************************************************
c   22.  Print-out for this basic time step..
c        If measured surface temps are provided, compute rms error.
c**********************************************************************
       if(debug .and. icalcstep .ge. idebug_step)
     & write(*,*)icalcstep,'22. VARIABLE = ',ICVPfixd
c***********************************************************************
c Sept 9 1996 - added variables for water output (ceretha)

1000  meltflux=meltflux + 0.5 * (uosave + uo(nsoil+1)) * dto  !REJ check this

cREJ      IF(ibasestep .eq.0 .or. print)THEN
      IF(print)THEN !print set /.true./ in initial and in converge
         if(itm .ge. 1)then
cREJ  6/26/24  if(ltype(k(n)).eq.1 .and.tm .gt. 273.1.and.tm.lt.373.15)
            if(H2Olayer(n) .and.tm .gt. 273.1.and.tm.lt.373.15)then
               tm=273.15d0
            endif
            if(tm.lt.373.15 .and. ibasestep .gt. 1)then
               rmsqt1=rmsqt1+(tm-t(n))**2
               difftemp=t(n)-tm
               rtmsq=dsqrt(rmsqt1/dble(ibasestep))
            endif
         endif

      if(writefluxmatlabflag) then ! ARK 7/29/24
              call writefluxmatlab(iy,jday,ihour,ibasestep,iptype,
     &        height,clearness,imin,isolarcalc,effceiling,islope,
     &        ifluxout)
          endif
      
      if (writematlabflag) then
        call writematlab(pinv,ibasestep,ihour,iy,jday,tm,
     &        de0,bp,difftemp,rtmsq,icalcstep,fnm,tmsg,itm,
     &        height,frvis,dtmin,dtmax,dtsmin,
     &        dtssmax,dssallowed,errtallowd,ngoodmin,imin,albsnow,
     &        islope,isolarcalc,ircalc,elev,azslope,
     &        bifallin,dmlimit,istboff,ssisnow,iqturb,eta0,ICVPfixd,
     &        Neumann,depth,H2Olayer) !ARK 7/29/24
       endif      
      
c
c July 8 1996 added option to write shortened output format (ceretha)
c
      writefmt = 0  ! standard write output. 
c     writefmt = 2  ! >= 2 ceretha water output file

       if (writefmt .eq. 0)then
c        if(jday .eq. 185 .and. ihour .eq. 0)then  !Just for Summit testing.  REMOVE LATER
c           do 661 i=1,n
c             write(90,'(3F11.3,F11.6)')t(i),dz(i),bw(i),d(i)
c661        continue
c           stop 'after write'
c        endif
cS          write(*,*)'BW42',bw(42)
          call write(pinv,ibasestep,ihour,iy,jday,tm,
     &        de0,bp,difftemp,rtmsq,icalcstep,fnm,tmsg,itm,
     &        height,frvis,dtmin,dtmax,dtsmin,
     &        dtssmax,dssallowed,errtallowd,ngoodmin,imin,albsnow,
     &        islope,isolarcalc,ircalc,elev,azslope,
     &        bifallin,dmlimit,istboff,ssisnow,iqturb,eta0,ICVPfixd,
     &        Neumann,depth,H2Olayer,snowdepth) !REJ 6/24/24
       else if (writefmt .ge. 2)then
        call write2(pinv,ibasestep,imin,ihour,iy,jday,icalcstep,
     &	      meltflux,writefmt)
       endif
         meltflux=0d0
      endif

c***********************************************************************
ct*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t
c END OF MAIN LOOP.
      goto 1001  ! Return to start of MAIN LOOP

 999  continue   ! Code jumps to 999 when the end of Met File is found
  
ct*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t
c***********************************************************************
      if(itm .ge.1)then  
         write(80,2000)rtmsq
2000    format(/,5x,'root mean square error of surface temp = ',f10.4,/)    
      end if
cS*S*S*S*S*S*S*S*S*WRITE SPINUP OUTPUT FILES*S*S*S*S*S*S*S*S*S*S*S*S*S*S
c REJ First generate the file name 
      If(spinup)Then
       write(*,*)'End spin-up year', ii
c     Generate names for spin-up output files
       spfilename = fnm(4)
       nc = len_trim(spfilename)-3
       spfilename(nc:nc+1)='v4'  !REJ 8_2_24
       nc = nc + 2  !RBJ 8_2_24
       if(ii .lt. 10)then
         spfilename(nc:nc+8) = '_spinyr0'
         write(spfilename(nc+8:nc+9),'(i1)')ii 
       else
         spfilename(nc:nc+7) = '_spinyr'
         write(spfilename(nc+7:nc+9),'(i2)')ii
       endif
       write(*,*)'spfilename  ',spfilename

c REJ Block to write short spin-up output files
        if(ICVPfixd)then ! 10/22/23 RJ Write final data to short output 
          open(92, file = spfilename ,status='unknown')
          write(92,*)
     &   'Node   Depth(m)   Temp(K)   dz(m) Density(kg/m3) Diameter(m)'
          do 700 i = 1,n         
             write(92,55)i,depth(i),t(i),dz(i),bw(i),d(i)
55           format(i4,4F11.3,F11.6) 
700       continue
          rewind(92)
        EndIf

         close (90); close(80); close(88)

c    Next for initializing between years during spinup
         n=0;ibasestep=0;igood=0;ititle=0        !REJ 12/29/23

cS*S*S*S*S*S*S*S*S*S*S*END SPIN-UP WRITE*S*S*S*S*S*S*S*S*S*S*S*S*S*S*S
      ENDIF

      close(180)  !REJ_2025/02/11. temp write close  !REMOVE REJ_2025/02/11
c      stop  'temp stop at end of main spinyr = 1'     !REMOVE REJ_2025/02/11
      END DO  ! END OF YEARLY SPIN_UP (ii) LOOP  REJ 12/29/23

c***********************************************************************
c   23.  Closure
c***********************************************************************
      stop '** Execution Completed ** '
      end
