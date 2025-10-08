c***********************************************************************
c  TRIDIAGN. Thomas Tri-Diagonal-Matrix algorithm for a Neumann lower BC 
c  References: Wikipedia. Tridiagonal Matrix Algorithm. 
c  H.I-I Conte, S.D., and deBoor, C., Elementary Numerical Analysis, 
c  McGraw-hill, New York. FORTRAN code W. H. Mason
c  https://archive.aoe.vt.edu/mason/Mason_f/CAtxtAppH1.pdf
c  Note typo in above document T(n) = bp(n)/ap(n,2), Not bp(n)/a(n,2)
c  Look for a more standard reference
c***********************************************************************
c     An entirely new routine.  11/6/2023 RJ
      subroutine tridiagN(n,a,T,B,Neumann,melt)

c arguments
c
c a(nd,3) : Coefficient array in linear equation matrix
c b(nd) : Constant array in linear equation matrix.  See subroutne fbb.f
c n : Number of nodes in model
c t(nd) : Nodal temperature [K]
c
c Called from filtrate and thermal
c
      implicit none
      include 'const'
      integer n,melt(nd)
      double precision a(nd,3),T(nd),B(nd),To(nd)
      logical Neumann
c
c local
c
c i: looping index.
c m: array index (m=n-i).
c n1: n1=n-1. index used in downsweep loop.
c ap(nd,3),bp(nd)  Altered "prime" array
c xm is a parameter in the solution
c do balchk: if true, prints out balance chk for all nodes, pauses

      integer i,n1,m
      double precision ap(nd,3),Bp(nd),xm,x
      logical dobalchk

      if(.not. Neumann)stop 'Use tridiag for fixed T bottom BC'

      dobalchk = .false. ! Set as true for balchk print

c     Upsweep
c     Solves simultaneous equations for nodes i and i+1 to
c     eliminate T(i-1). Replacs matrices a(nd,3) and b(nd)
c     with the 'prime' matrices ap(nd,3) and bp(nd) as  
c     i sweeps upward from 2 to n-1. Each nodal equation is
c     left with 2 unknowns, T(i) and T(i+1).

      ap(1,3)=a(1,3); ap(1,2)=a(1,2); bp(1)=b(1) 
      do 15 i=2,n
           xm = a(i,1)/ap(i-1,2)
           ap(i,3)= a(i,3)
           ap(i,2)= a(i,2)-xm*ap(i-1,3)
           ap(i,1)=0d0
           bp(i) = b(i)-xm*bp(i-1)
30    format(i4,5f13.9)
 15   continue



c     Top node
c     When the top node is reached, only two unknowns,
c     T(n) and T(n-1), remain. Simultaneous solution of the  
c     two top equations for T(n) follows.

      T(n) = bp(n)/ap(n,2)

c     Downsweep.
c     Since T(i+1) is now known, the equation for T(i) can 
c     be solved for directly as i sweeps downward from node
c     n-1 to node 1.
      n1=n-1
      do 10 i=1,n1
         m=n-i
         t(m)= (bp(m) - ap(m,3)*T(m+1))/ap(m,2)
10    continue

      if(dobalchk)then  !dobalchk = .true. for print of 
     &    !thermal balance after Neumann tridiag'
        write(*,*)'BalChk: i,A*T,B',      !Bottom node
     &     1,a(1,3)*T(2)+a(1,2)*T(1),b(1)
        do i = 2,n-1                        
           write(*,*)'BalChk: i,A*T,B',   !Interior nodes
     &     i,a(i,3)*T(i+1)+a(i,2)*T(i)+a(i,1)*T(i-1),b(i)
        enddo
       write(*,*)'BalChk: i,A*T,B',       !Top node
     &     n,a(n,2)*T(n)+a(n,1)*T(n-1),b(n)
       stop 'TRIDIAGN. After BALCHK'
      endif

      return
      end
