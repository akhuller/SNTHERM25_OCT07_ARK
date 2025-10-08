c  Printed on 03/10/2025
c  sntherm24/sntherm_June15 2024 version with noted updates
c  Several c_d out testing prints removed on 03/10/2025.  Version WITH
c  testing saved to opt_properties_03_10_2025.f  This is somewhat cleaned
c  up version
c  ARK fixed dust-snow and dust-ice mixing formulae on 9/22/2025

c***********************************************************************
c  OPT_PROPERTIES calculates optical properties of a single layer
c  for a mixture of snow/firn/glacier ice and impurities.
c  
c  REFERENCES:
c  Wiscombe, W. J. & Warren, S. G. A model for the spectral albedo of 
c  snow, I: Pure snow. (1980)
c  Warren, S. G., & Wiscombe, W. J. A model for the spectral 
c  albedo of snow. II: Snow containing atmospheric aerosols. (1980)
c  Mullen, P. C. & Warren, S. G. Theory of the optical properties of 
c  lake ice. (1988)
c  Khuller, A. R., Christensen, P. R. & Warren, S. G. Spectral albedo of
c  dusty martian H2O snow and ice. (2021)
cc***********************************************************************

      subroutine opt_properties(ltype,radius,rbbl,bw,r_dust, rho_d,
c     &   dz,concen,lambda,mImag,g_data,dust_mie_dir,mie_dir,
     &   dz,concen,dust_mie_dir,mie_dir,
     &   use_SSA,porosity_rad,ssArea,t_net_lyr,w_net_lyr,g_net_lyr,ii) !REJ_2025/03/04
c
c Called from sntherm (MAIN) Section 9
c
c REJ note on 2025/03/07. *Several inputs in microns converted to meters [MKS]
c                          All units now in MKS, except wavelength remains as microns.
c
c INPUTS:
c ltype         sntherm layer type. i.e.,snow=1 soil 2 to 4,ice/firn 90 to 94
c                 (see specifics in get_input.f) 
c radius        Radius of snow/firn/glacier ice grain [microns]. *Now [meters]
c rbbl          Radius of bubbles within firn/glacier ice [microns) *Now [meters]
c bw            Bulk water density of snow/firn/glacier ice [kg/m^3]
c r_dust        Radius of dust grain [microns] *Now [meters]
c rho_d         Density of dust [kg/m^3]
c dz            thicknessness of CV (20m => infinite optical depth)
c concen        Dust concentration [ppmw] *Now [ppw]
c lambda        Wavelength [microns] micron wavelengths for snow (radius<500e-6 m)
c mImag         Complex part of ice refractive index (spectral)
c g_data        Air bubble asymmetry parameter in non-absorbing ice 
c dust_mie_dir  Directory with pre-calculated Mie outputs for 0.2-3.0 
c               micron wavelengths for impurities
c mie_dir       Directory with pre-calculated Mie outputs for 0.2-3.0 
c use_SSA       0 = Use porosity_rad to calculate SSA, 1 = use input SSA
c porosity_rad  Porosity of ice
c ssArea        Specific surface area (SSA): the area of air-ice interfaces                                 per unit massof ice [m2/kg]

c OUTPUTS:
c t_net_lyr     Net optical depth of layer
c w_net_lyr     Net single-scattering albedo of layer
c g_net_lyr     Net asymmetry parameter of layer

      implicit none
      include 'const'
      include 'arrays_spectral'  !spectral csv input
      integer nx,ny
      double precision rho_Ice
      parameter (nx=281, ny = 5, rho_ice = 917d0)
      integer fid,fid_d,iter,ltype,use_SSA,i,io_stat  
c     Note: ltype returns the layertype (snow,ice,firn,soiltype,etc) for node i
      ! include next parameter line in const file
      double precision PI,ssArea         
c      double precision lambda(nx), g_data(nx), mImag(nx) REJ In arrays_spectral
      double precision bw, dz, concen, rho_d, r_dust
      double precision beta_abs_glacier(nx), beta_abs_ice(nx)
      double precision beta_sca_glacier, kSca_glacier
      double precision kAbs_glacier(nx), rbbl
      double precision sigma_s_abs(nx), sigma_s_sca(nx),g_ice (nx)
      double precision sigma_d_abs(nx), sigma_d_sca(nx),g_d (nx)
      double precision Q_abs(nx), Q_sca(nx), kAbs_snow(nx)
      double precision Q_d_abs(nx), Q_d_sca(nx)
      double precision kSca_snow(nx), kAbs_dust(nx), kSca_dust(nx)
      double precision kExt_snow(nx), kExt_dust(nx), kExt_glacier(nx) ! ARK_2025/09/22
      double precision t_net_lyr(nx), w_net_lyr (nx), g_net_lyr (nx)
      double precision V_s, A_s, vAir ,vIce, porosity_rad,radius
      double precision Cdust, Csnow, Cglacier,V_d,A_d
      character*160 mie_dir,dust_mie_dir, filename, filename_d, x1, x2
      integer i1,ii !REJ_2025/03/10 
      double precision mie_data(nx*5), mie_d_data(nx*5)

      PI = 4.D0*DATAN(1.D0)
      ! test code follows
c      test = 2.1
c      write (x1,'(F3.1)') test
c      filename= trim(mie_dir) // trim(x1) // 'f.bin'
      
c      open(NEWUNIT = fid, FILE = filename, STATUS = "OLD", 
c     &   ACCESS = "STREAM", FORM = "UNFORMATTED", IOSTAT = io_stat)
      
c        do i =1,nx*5   
c          read(fid, IOSTAT = io_stat) mie_data(i)
c          if (i .ge. nx*2+1 .and. i .le. nx*3) then
            !g_ice(i-nx*2)     = mie_data(i)
            !write(*,*) g_ice(i-nx*2)
c          else if (i .ge. nx*3+1 .and. i .le. nx*4) then
            !Q_abs(i-nx*3)     = mie_data(i)
            !write(*,*)  Q_abs(i-nx*3)
c          else if (i .ge. nx*4+1 .and. i .le. nx*5) then
            !Q_sca(i-nx*4)     = mie_data(i)
            !write(*,*)  Q_sca(i-nx*4)
c          end if
c         end F
c        close(fid)

c      radius = radius*1.0d-6  ! convert from microns to m   !REJ Now input in meters
      concen = concen/1.0d6   ! convert from ppmw to ppw REJ: MUST CHNG LATER
      Cdust = concen          ! Mass fraction of dust
      Csnow = 1.0 - Cdust       ! Mass fraction of snow
      Cglacier = 1.0 - Cdust    ! Mass fraction of firn/glacier ice
      fid   = 5
      fid_d = 6

c      bw = bw  ! ARK_2025/2/14  REJ Note: rhos replaced by bw
      
      if (concen .eq. 0.0) then ! no need to calculate dust properties
         

         do i=1,nx ! loop over wavelengths !REJ_2025/03/04 
c            write(*,*)i,nx
c dust mass absorption coeff. [m^2/kg]
            kAbs_dust(i) = 0.0
c dust mass scattering coeff. [m^2/kg]            
            kSca_dust(i) = 0.0
c dust asymmetry parameter            
            g_d(i)       = 0.0
c dust mass extinction coeff. [m^2/kg]           ! ARK_2025/09/22 
            kExt_dust(i) = kAbs_dust(i) + kSca_dust(i)
         end do
      else ! need to calculate dust properties
        
c  Load pre-calculated impurity Mie outputs for 0.205 to 3.0 microns

        write (x2,'(F3.1)') r_dust ! convert double to string
c        filename_d = trim(dust_mie_dir)//'Mars_dust_'//trim(x2)//'0.bin' !ARK_2025/2/10S

        filename_d =
     &  trim(dust_mie_dir)// '/' //'Mars_dust_'//trim(x2)//'0.bin' !ARK_2025/2/10
        ! note that I had to manually add a zero at the end above
        
        
        r_dust     = r_dust*1.0d-6 ! convert from microns to m

c Extract Mie data

      open(NEWUNIT = fid_d, FILE = filename_d, STATUS = "OLD", 
     &  ACCESS = "STREAM", FORM = "UNFORMATTED", IOSTAT = io_stat)
      
        do i =1,nx*3 ! since there are 3*nx rows of data divided by nx
          read(fid_d, IOSTAT = io_stat) mie_d_data(i)
          
          if (i  .le. nx) then
c dust asymmetry parameter                   
            Q_d_abs(i)     = mie_d_data(i)
            !write(*,*) Q_d_abs(i) 
          else if (i .ge. nx+1 .and. i .le. nx*2) then
c dust absorption efficiency
            Q_d_sca(i-nx)     = mie_d_data(i)
            !write(*,*)  Q_d_sca(i-nx)
          else if (i .ge. nx*2+1 .and. i .le. nx*3) then
c dust scattering efficiency            
            g_d(i-nx*2)     = mie_d_data(i)
            !write(*,*)  g_d(i-nx*2)
          end if
         end do
        close(fid_d)
        
c dust sphere volume [m^3]            
            V_d         = (4./3.)*PI*r_dust**3.0

c dust sphere cross-sectional area [m^2]            
            A_d         = PI*r_dust**2.0

          do 414 iter=1,nx ! loop over wavelengths

c dust absorption cross-section [m^2]
            sigma_d_abs(iter) = Q_d_abs(iter)*A_d
            
c dust scattering cross-section [m^2]
            sigma_d_sca(iter) = Q_d_sca(iter)*A_d

c dust mass absorption coeff. [m^2/kg]
            kAbs_dust(iter)   = sigma_d_abs(iter)/(rho_d*V_d)

c dust mass scattering coeff. [m^2/kg]
            kSca_dust(iter)   = sigma_d_sca(iter)/(rho_d*V_d)

c dust mass extinction coeff. [m^2/kg]           ! ARK_2025/09/22 
            kExt_dust(iter) = kAbs_dust(iter) + kSca_dust(iter)   
414       continue
      
      end if  !REJ End Dust Block 
      
c  If using bubble radii, trigger firn/glacier case

      if (rbbl .gt. 0) then 
        radius    = 2000d0  ! some arbitrary radius which won't be used ARK/2025/2/14
        rbbl = rbbl*1.0d-6
      end if
       i = 1  ! RECHECK THIS REJ_2025/03/07
c Check if we are modeling snow vs. firn/glacier ice
       
      IF (radius .lt. 500*1.0d-6) THEN ! Start Snow REJ_2025/03/04
c        if(ltype .eq. 1)then !REJ See getinput.f for layer type codes
c snow sphere volume [m^3]            
        V_s         = (4./3.)*PI*radius**3
        
c snow sphere cross-sectional area [m^2] 
        A_s         = PI*radius**2.0

c Load pre-calculated Mie outputs for 0.2 to 3.0 micron wavelengths
c       In directory sntherm_ice_bin_v2  REJc_2025/03/03    
        if (radius*1.0d6 .lt. 5.0) then  !radius in meters, converted to microns
          write (x1,'(F3.1)') radius*1.0d6 ! convert double to string
c         0.1fbin to 4.9fbin  REJc_2025/03/03
          filename= trim(mie_dir) // '/' // trim(x1) // 'f.bin' ! ARK_2025/2/10
        else
          write (x1,'(I3.0)') nint(radius*1.0d6)
c         5.bin to 500.bin: REJc_2025/03/03
          filename= trim(mie_dir) // '/' // trim(adjustl(x1)) // '.bin' ! ARK_2025/09/12 to ensure leading spaces don't cause issues for <100
        end if
 
          
c Extract Mie data for snow

        open(NEWUNIT = fid, FILE = filename, STATUS = "OLD", 
     &  ACCESS = "STREAM", FORM = "UNFORMATTED", IOSTAT = io_stat)
      
        do i =1,nx*5 ! since there are 5*nx rows of data divided by nx
          read(fid, IOSTAT = io_stat) mie_data(i)
c          if (mie_data(i) .gt. 0d0) stop 'not zero' !REJc_2025/03/10 Ask Adi 
          if (i .ge. nx*2+1 .and. i .le. nx*3) then
c snow asymmetry parameter            
            g_ice(i-nx*2)     = mie_data(i)
            !write(*,*) g_ice(i-nx*2)
          else if (i .ge. nx*3+1 .and. i .le. nx*4) then
c snow absorption efficiency            
            Q_abs(i-nx*3)     = mie_data(i)
            if(i-nx*3 .le. 20) then
c              write(*,*)i,nx,i-nx*3, Q_abs(i-nx*3)
            endif
          else if (i .ge. nx*4+1 .and. i .le. nx*5) then
c snow scattering efficiency            
            Q_sca(i-nx*4)     = mie_data(i)
            !write(*,*)  Q_sca(i-nx*4)
          end if
         end do
        close(fid)

        do 789 iter=1,nx ! loop over wavelengths

c snow absorption cross-section [m^2]
          sigma_s_abs(iter) = Q_abs(iter)*A_s
          
c snow scattering cross-section [m^2]
          sigma_s_sca(iter) = Q_sca(iter)*A_s

c snow mass absorption coeff. [m^2/kg]
          kAbs_snow(iter) = sigma_s_abs(iter)/(rho_ice*V_s) 

c snow mass scattering coeff. [m^2/kg]
          kSca_snow(iter) = sigma_s_sca(iter)/(rho_ice*V_s)
          
c snow mass extinction coeff. [m^2/kg]           ! ARK_2025/09/22 
          kExt_snow(iter) = kAbs_snow(iter) + kSca_snow(iter)
            
c REJ_2025/03/06 Next 5 for Testing
c          write(*,*)'ksca_snow',  Test shortened version
c     &     ksca_snow(iter),Q_sca(iter)*0.75d0/(917d0*radius)
c          write(*,*)'kAbs_snow',
c     &     kabs_snow(iter),Q_abs(iter)*0.75d0/(917d0*radius) 
c        write(95,*)iter,lambda(iter),kAbs_snow(iter),kSca_snow(iter)
  
c------ snow-dust mixing---------c
 
c net single-scattering albedo [unitless]
          w_net_lyr(iter)        = (Csnow*kSca_snow(iter) +   ! ARK_2025/09/22
     &    Cdust*kSca_dust(iter))/(Csnow*kExt_snow(iter) 
     &    + Cdust*kExt_dust(iter)) 

c net asymmetry parameter [unitless]
          g_net_lyr(iter)        = (Csnow*g_ice(iter)*kSca_snow(iter) + 
     &    Cdust * g_d(iter) * kSca_dust(iter))/
     &    (Csnow * kSca_snow(iter) + Cdust * kSca_dust(iter))

c net optical depth [unitless]     
          t_net_lyr(iter)         = dz*bw*(Csnow*kExt_snow(iter) ! ARK_2025/09/22
     &     + Cdust*kExt_dust(iter))
          
          
789       continue

      ELSE !  Firn/glacier ice  
         ! stop 'firn'

c fraction of air (porosity)
      if (use_SSA .eq. 0) porosity_rad = 1.0 - bw/917.
      
          if (rbbl .gt. 0)  then ! if input as bubble radius (microns)
            ssArea = (3.0*porosity_rad)/(bw*rbbl) ! [m^2/kg]
          else        ! if input as ice grain radius 
            ssArea = 3.0/(radius*917.) ! [m^2/kg]
          end if

c fraction of air (porosity)
          vAir     = porosity_rad  ! [unitless]

c fraction of ice 
          vIce     = 1.0 - vAir ! [unitless]

c glacier ice scattering coefficient [m^-1]
          beta_sca_glacier = (1 - porosity_rad) * 917.0 * ssArea * 0.5

c mass scattering coeff. of the glacier ice [m^2/kg]
          kSca_glacier     = beta_sca_glacier/bw
          
c      write(*,*)'beta,bw,kSca_glac',beta_sca_glacier,bw,kSca_glacier
c     & ,ssArea/2d0
      
c next REJ 2025/04/10.. write opt_properties output
        open(145, file = 'o_prop_out', status = 'unknown')
        do 790 iter = 1 , nx ! loop over wavelengths

c pure ice scattering coefficient [m^-1]
          beta_abs_ice(iter)     = 4.0 * PI * mImag(iter) 
     &                             /(lambda(iter)*1.0d-6) 
     
c glacier ice absorption coefficient [m^-1]
          beta_abs_glacier(iter) = vIce*beta_abs_ice(iter)    

c mass absorption coeff. of the glacier ice [m^2/kg]          
          kAbs_glacier(iter)     = beta_abs_glacier(iter)/bw

c mass absorption coeff. of the glacier ice [m^2/kg]          
          kExt_glacier(iter)     = kAbs_glacier(iter) 
     &                            + kSca_glacier ! ARK_2025/09/22
          
c          write(95,*)lambda(iter),kabs_glacier(iter),
c     &     kSca_glacier
 
          
c------ firn/glacier ice-dust mixing---------c

c net single-scattering albedo of layer [unitless]
          w_net_lyr(iter)    = (Cglacier*kSca_glacier 
     &     + Cdust*kSca_dust(iter))/(Cglacier*kExt_glacier(iter)
     &     + Cdust*kExt_dust(iter)) ! ARK_2025/09/22

c net optical depth of layer [unitless]

          t_net_lyr(iter)    = dz*bw*(Cglacier*kExt_glacier(iter)
     &                        + Cdust*kExt_dust(iter)) ! ARK_2025/09/22
           
c net asymmetry parameter of layer
          g_net_lyr(iter)    = (Cglacier*g_data(iter)*kSca_glacier +
     &                         Cdust*g_d(iter)*kSca_dust(iter))/
     &                         (Cglacier*kSca_glacier + 
     &                         Cdust*kSca_dust(iter))
     
 110     format(i3) !next writes REJ_2025_4/10
c         if(iter .ge. 10 .and. iter .le. 20) 
c     &    write(145,*)ii,iter,w_net_lyr(iter),t_net_lyr(iter),
c     &  g_net_lyr(iter)

790      continue

      ENDIF ! End Firn/Ice
c      stop 'end opt'


      return
      end
      




    




