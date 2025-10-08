c***********************************************************************
      subroutine write_solar(n,iy,jday,ihour,F_abs,Fabs_int,Fint,
     %    albedo_int,rad,dlambda)
c     written by REJ on 6/25/2024.
c     Excerpted from Khuller MatLab routine, netsolarfluxes.f
c     First called after initial step, when both the F_net and Fabs_int 
c     files are printed, for comparison with solartest output. Fint
c     and albedo_int are now included in the Fabs_int file.
c     Thereafter, the Fabs_int file only prints out once per hour.
c***********************************************************************
      implicit none
      include 'const'
      integer n,k,iy,jday,ihour,iw
      double precision F_abs(nbr_wvl,ncv),Fabs_int(ncv),Fint(ncv)
      double precision albedo_int,rad,dlambda,Fsum
      write(17,*)'writing to unit 17'
      write(18,*)'Year=',iy, 'Jday=',jday, 'Hour=', ihour
      write(18,*)
      write(18,*)' Node/               Abs                    NetFlux'  
      write(18,*)' BotIntFace','                        TopIntFace = '
     &  ,n+1
      write(18,'(F66.20)')Fint(n+1)
      do k = n,1,-1
         Fsum = 0d0
         do iw =1,nbr_wvl
            Fsum = Fsum + F_abs(iw,k)*dlambda*1000d0
         enddo
         Fabs_int(k) = Fsum
         write(18, '(i6,2F30.20)')k,Fabs_int(k),Fint(k)
      enddo
      write(18,*)
      write(18,*)'albedo_int = ', albedo_int, 'sdown (=rad)', rad 
      write(18,*)
      close(18)

      return
      end