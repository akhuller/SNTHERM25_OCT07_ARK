c***********************************************************************
c SUBROUTINE TOTAL_SOLAR.f
c     written by REJ on 6/25/2024.
c     ARK 12/26/24 Removed mu_not*pi multiplication 
c     Excerpted from Khuller MatLab routine, netsolarfluxes.f
c     First called after initial step, when both the F_net and Fabs_int 
c     files are printed, for comparison with solartest output. Fint
c     and albedo_int are now included in the Fabs_int file.
c     Thereafter, the Fabs_int file only prints out once per hour.
c***********************************************************************
      subroutine total_solar(n,iy,jday,ihour,fdirup,fdifup,fdirdn,
cREJ_2025/04/26&fdifdn,rad,Fs0,Fs,mu_not,PI,Fabs_int,albedo_int,dlambda,
     &fdifdn,rad,mu_not,PI,Fabs_int,albedo_int,dlambda,
     &solartest,spinup,icalcstep)!REJ_2025/04/26 added icalcstep
c***********************************************************************
      implicit none
      include 'const' !REJ 4/30/24
      include 'arrays_spectral' !REJ 4/30/24
c INPUTS:
c n		number of CVs
c iy,jday,ihour  year(2 digits),julian day, and hour of simulation
c fdirup        
c fdifup
c fdirdn 
c fdifdn
c rad           Broadband, downward solar flux [W/m^2] rad(1) from sntherm    
c Fs0           Normalized spectral, direct-beam downward flux [W/m^2/band]
c Fs            Spectral, direct-beam downward flux [W/m^2/band]
c Fd            Spectral, diffuse-beam downward flux [W/m^2/band]

c mu_not        Cosine of zenith angle; PI
c PASSED FROM ARRAYS SPECTRAL
c Fd0           Normalized spectral, diffuse-beam downward flux [W/m^2/band]

c OUTPUTS:
c albedo_int    broadband albedo
c Fnet_int      broadband net radiative flux
c Fabs_int      broadband absorbed radiation

c iter          iteration within the hourly subtime division

      integer n,iy,jday,ihour,iter,iw,i,k,io_stat,icalcstep  !REJ_2025/04/26
      integer idum  !REJ_2025 'REMOVE LATER!'
      logical solartest,solarwrite,do_write,spinup
      double precision fdirup(nbr_wvl,ncv), fdifup(nbr_wvl,ncv)
      double precision fdirdn(nbr_wvl,ncv), fdifdn(nbr_wvl,ncv)
      double precision F_dn  (nbr_wvl,ncv), F_up  (nbr_wvl,ncv)
      double precision F_net (nbr_wvl,ncv), F_abs  (nbr_wvl,ncv)
      double precision Fabs_int(nd), dsol(nd), Fnet_int(ncv)
cREJ_2025_04_26      double precision Fs0(nbr_wvl), Fs(nbr_wvl),Fd(nbr_wvl)
      double precision albedo(nbr_wvl),dlambda,mu_not
      double precision albedo_int,rad,Fsum,Fsum2,pi

      data solarwrite/.false./

      
          
      
      ! Write out net spectral fluxes for each interface = layer+1 
      ! (format is (n+1)*nbr_wvl in one row)
      ! I use matlab to transpose the row into an array of size
      ! [n+1 nbr_wvl]

c  Testing :remove block later!!
      idum = 630

      do 55 iw = 6,25
        do 56 k = n,n-4,-1
cF        if(icalcstep .gt. idum)write(*,*)iw,k,fdirup(iw,k)
56      continue
cF        if(icalcstep .gt. idum)write(*,*)'rad*Fs0(iw)',rad*Fs0(iw)
55     continue
cF       if(icalcstep .gt. idum+10)stop 'total_solar'
c  End testing block

      if(solartest)then  !REJ
        open(unit=10, file='F_net_solartest.txt', status='unknown')
      else
        open(unit=10, file='F_net_SN.txt', status='unknown')
      endif

      ! Write out spectrally-integrated absorbed fluxes

      if(solartest)then  !REJ
        open(unit=11, file='Fabs_int_solartest.txt',status='unknown')     
      else
        open(unit=11, file='Fabs_int_SN.txt', status='unknown')
      endif

c MULTIPLY BROADBAND SOLAR FLUX WITH SPECTRAL SOLAR FLUX

      do 123 iw=1,nbr_wvl
      Fs(iw) = 0d0
      if (Fs0(iw) .gt. 0d0) then  !REJ
      ! Use direct-beam
         Fs(iw) = rad * Fs0(iw)  !REJ
      endif

      ! Use diffuse beam 
      Fd(iw) = 0d0                 
      if (Fd0(iw) .gt. 0d0) then
         Fd(iw) = rad * Fd0(iw)  !REJ
      endif

123   continue

cccccc Added next section from netsolarfluxes (since this subroutine was broken into two by REJ - ARK_2025/2/14 ccc     c  Since 
        do 110 k=n+1,1,-1
          dsol(k) = 0d0
          Fabs_int(k) = 0d0 
          do 129 iw=1,nbr_wvl
              F_up(iw,k)   = 0d0
              F_dn(iw,k)  = 0d0
              F_net(iw,k)  = 0d0
              F_abs(iw,k)  = 0d0 
129	      continue
110   continue
cccccc ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc ccc
      
      
      
      do 702 iw = 1,nbr_wvl
cREJ         do 703 k = 1,n+1
         do 703 k = n+1,1,-1
c            F_up(iw,k) = (fdirup(iw,k)*(Fs(iw)*mu_not*PI) ARK 12/26/24
c     &                   + fdifup(iw,k)*Fd(iw))               
c            F_dn(iw,k) = (fdirdn(iw,k)*(Fs(iw)*mu_not*PI) ARK 12/26/24
c     &                   + fdifdn(iw,k)*Fd(iw))
     
            F_up(iw,k) =  (fdirup(iw,k)*(Fs(iw)) 
     &                   + fdifup(iw,k)*Fd(iw))

            F_dn(iw,k) = (fdirdn(iw,k)*(Fs(iw)) 
     &                   + fdifdn(iw,k)*Fd(iw))
cF          if((k .eq. n+1 .or. k .eq. n) .and. iw  .lt. 20)write(*,*) !REJ2025  Remove after testing'
cF     &    'tot_sol'.   k,iw,fdirdn(iw,k),Fs(iw) ,fdifdn(iw,k),Fd(iw)     
      ! Net spectral flux at each interface
            
c           F_net(iw,k) = F_up(iw,k) - F_dn(iw,k)
            
            F_net(iw,k) = F_dn(iw,k) - F_up(iw,k) !REJ SNTHERM convention
            
            
            if(.not. spinup .and. solartest)then !Override with care. Extreme output!!
             write(10, '(F6.2,i4,F23.20)')lambda(iw),k,F_net(iw, k) !REJ chged format .
            endif
            
            
703      continue 
      ! Spectral albedo
        albedo(iw) = F_up(iw,n+1)/F_dn(iw,n+1)
702   continue
      close(unit=10) 
      
c REJ Next loop calculates the integrated surface albedo

      Fsum = 0d0 
      Fsum2 = 0d0
      k= n+1
      do 712 iw = 1,nbr_wvl ! Sum to get the spectrally-integrated albedo
           Fsum  = Fsum  +  F_dn(iw,k)*dlambda*1000d0
           Fsum2 = Fsum2 +  F_up(iw,k) *dlambda*1000d0
712   continue

      albedo_int = Fsum2/Fsum 
cF      if(icalcstep.ge.idum)write(*,*)'albedo_int2',icalcstep,Fsum2/Fsum,
cF     &      Fsum  !REJ_2025 delete this after testing
cF      if(icalcstep .eq. idum +10)stop 'total_solar'  !REJ_2025 delete this after testing

c      'check:Fsum and rad should be equal',Fsum,rad,ihour,jday
      if(dabs(Fsum-rad).gt.1d0)then !Note: Initial Fd0 must sum to 1d0
c        write(*,*)'Fsum,sdown',Fsum,rad
c        stop
      endif

      do 704 iw = 1,nbr_wvl
cREJ         do 705 k = 1,n
         do 705 k = n,1,-1  
      ! Absorbed spectral flux in each layer:  Note definition in Matlab of 
c       flux positive upwards?  I am changing it to positive downwards.

         F_abs(iw,k) = F_net(iw,k+1) - F_net(iw,k)
         
705      continue 
704   continue 

cREJ      do k = 1,n+1
      do k = n+1,1,-1
         Fsum = 0d0
         do iw =1,nbr_wvl
            Fsum = Fsum + F_net(iw,k)*dlambda*1000d0
         enddo
cREJ         Fint(k) = -Fsum
             Fnet_int(k) = Fsum
      enddo
      
      do k = n,1,-1
         Fsum = 0d0
         do iw =1,nbr_wvl
            Fsum = Fsum + F_abs(iw,k)*dlambda*1000d0
         enddo
         Fabs_int(k) = dmax1(0d0,Fsum) !REJ_2025/07/24
         dsol(k)     = dmax1(0d0,Fsum) !REJ_2025/07/24
      enddo   
      
c     WRITE Broadband Fabs_int, Fint, and albedo_int
      if(spinup)do_write = .false.   
      IF(do_write)THEN  
c      Note: solartest is specified in MAIN
       If(.not. solartest)Then ! regular Sntherm run
          solarwrite = .false.  !REJ added on 2025/05/05
c         Prints out Fabs and Fint once per day at Noon.
          if(iter .le. 1 .and. ihour .eq. 12) 
     &    solarwrite=.true.
       Else !Compare sntherm embedded MatLabSolar with solartest
c       Prints Fabs_int and Fint for jday=1,ihour=0 and stops.
c       Use input files for solartest in FILENAME
c       In solartest.f,chng thickness of top 3 CVs to .01,.03, and .06m
c                      to get accurate comparison,
          if(jday .eq. 1 .and. ihour .eq. 0 .and. iter .le. 1) 
     &    solarwrite=.true.
       EndIf
      ENDIF
      solarwrite = .false.   ! REMOVE by end of Sat.
      IF(solarwrite)THEN 
      write(11,*)'Year=',iy, 'Jday=',jday, 'Hour=', ihour
      write(11,*)
      write(11,*)' Node/               Abs                    NetFlux'  
      write(11,*)' BotIntFace','                        TopIntFace = '
     &  ,n+1
      write(11,'(F66.20)')Fnet_int(n+1)
      do k = n,1,-1  
         write(11,'(i6,2F30.20)')k,Fabs_int(k),Fnet_int(k)
      enddo
      write(11,*)
      write(11,*)'albedo_int = ', albedo_int, 'sdown (=rad)', rad 
      write(11,*)
CREJ      close(unit =11) !Note: file is continuosly overwritten
      if(solartest)
     &   stop 'Fabs_int_solartest.txt written. **Execution complete**'
      ENDIF  
c***********************************************************************
      return
      end