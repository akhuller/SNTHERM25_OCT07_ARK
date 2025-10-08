***********************************************************************
c  WRITEMATLAB generates principal output file. Parameters for the top
c  node are printed following the completion of a basic time step (in 
c  sync with time periods when met values are provided).  Full snow/soil
c  profiles are printed at an optional time interval specified in the
c  input parameter IPRNT
c***********************************************************************
      subroutine writematlab(pinv,ibasestep,ihour,iy,jday,tmn,
     &     de0,bp,difftemp,rtmsq,icalcstep,fnm,tmsg,itm,
     &     height,frvis,dtmin,dtmax,dtsmin,
     &     dtssmax,dssallowed,errtallowd,ngoodmin,imin,albsnow,islope,
     &     isolarcalc,ircalc,elev,azslope,bifallin,dmlimit,istboff,
     &     ssisnow,iqturb,eta0,Khuller,Neumann,snowdepth,H2Olayer) ! ARK_2025/2/12



c Called from MAIN
c
      implicit none
      include 'const'
      include 'arrays'
c arguments
c
c pinv : major print-out interval [hr]
c ibasestep: Current number of basestep from start of problem
c ihour : Hour of simulation, (read from met time-hack) 
c imin: Minute of simulation, (read from met time-hack)
c iy : Year of simulation: last 2 digits, (read from met time-hack)
c jday : Julian day of simulation, (read from met time-hack)
c tmn : measured surface temperature (in readtm)
c de0 : effective diffusion coefficient for water vapor in snow  
c bp : approximate average barometric pressure [mb]
c difftemp : difference between measured and calculated surface temperature [k]
c rtmsq : rms error of (calculated - measured surface temperature)
c icalcstep : current number of calcultions from start of problem
c fnm(nfiles) : list of i/o file names
c tmsg : measured snow/ground interface temperature [k]
c itm : read measured temperature flag 1=yes 0=no
c height : height above ground of measured met values
cRJ 12/16/23 bext : snow extinction coefficient for near-IR 
c Fraction of solar spectrum  0.4 to 1.12 microns  cRJ 12/16/23
c dtmin : minimum allowed time step [s]
c dtmax : maximum allowed time step [s]
c dtsmin:  Minimum allowed time step when water flow is present [s]
c dtssmax : maximum time step for saturation criteria [s]
c dssallowed : maximum allowed nodal change in saturation per
c errtallowd: maximum allowed linearization error per time-step in heat 
c             balance eq, expressed in units of temperature [K]
c ngoodmin : minimum number of consecutive good calculations required
c
      integer ibasestep,ihour,iy,jday,istboff,iqturb
      integer itm,icalcstep,ngoodmin,imin,islope,isolarcalc,ircalc
      double precision pinv,bp,de0,tmn,difftemp,rtmsq,tmsg,eta0
cRJ 12/16/23      double precision height(3),bext,dtmin,dtmax,dtsmin
      double precision height(3),frvis,dtmin,dtmax,dtsmin !RJ 12/16/23
      double precision dtssmax,dssallowed,errtallowd,albsnow
      double precision elev,azslope,bifallin,dmlimit,ssisnow,snowdepth ! ARK_2025/2/12
      logical Khuller,Neumann,H2Olayer(nd)  ! 10/22/23 + 6/24/24 REJ
c 
c Local
c
c dzsnow: Thickness of snow in temperature profile.
c fnm(nfiles): list of i/o file names.
c i: looping index variable.
c irem: irem = mod(jday,100).
c m: m=k(i).
c profile: name of output file for temperature profile of snow.
c tsg: t(sg)  = Predicted Snow/Ground Interface Temperature (K).
c precip : Precip type. S=snow, R=rain (character string)
c phase : Phase state, F, T or M (character string)
c
c Passed through common (incomplete)
c
c snowrateo : Old snow rate [m/s]

c rainrateo : Old Rain Rate [m/s]
c n : Number of nodes in model
c k : Layer type assoicated with node
c t : Nodal temperature [K]
c th(ld) : Upper meltzone temperature limit [K]
c tl(ld) : Lower meltzone temperature limit [K]
c
      integer i,m,irem
      double precision tsg,dzsnow,swe
      logical do_profiles  !RJ 12/15/23
      character*45 fnm(nfiles) !RJ 12/16/23
      character*12 profile !RJ 12/16/23
      character*1 precip,phase(nd)

      
      if(icalcstep .eq. 1)then ! ARK_2025/2/12
          write(20,60) (iy,jday,ihour,imin,i,z(n)+(dz(n)/2d0)-z(i),
     &    to(i),bt(i),bw(i),bl(i),ct(i),thk(i),do(i)/2d0,i=n,1,-1)
60    format(i2,1x,i3,1x,i3,1x,i3,1x,i3,1x,f8.5,f8.2,f8.2,f8.2,f8.2, ! ARK_2025/09/17
     & f8.2,f8.2,f8.5)  
      else
          write(20,70) (iy,jday,ihour,imin,i,z(n)+(dz(n)/2)-z(i),
     &    t(i),bt(i),bw(i),bl(i),ct(i),thk(i),d(i)/2d0,i=n,1,-1)
70    format(i2,1x,i3,1x,i3,1x,i3,1x,i3,1x,f8.5,f8.2,f8.2,f8.2,f8.2, ! ARK_2025/09/17
     & f8.2,f8.2,f8.5)  
      end if
      
 
      return
      end
