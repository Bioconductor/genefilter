c By R Gray, March 19, 2000, DFCI
c Copyright (C) 2000 Robert Gray
c Distributed under the GNU public license
c
c t-statistics
c first ng1 columns of d assumed to be group 1, other ng-ng1 assumed to be
c group2.  Note: single precision stats
c
c  Modified by R. Gentleman, 2004, just extracted the ttest stats and 
c  computed a ratio on demand - or fold-change

      subroutine fastt(d,n,ng,ng1,z,dm,eqv,ratio)
      real d(n,ng),z(n),dm(n)
      integer n,ng,ng1,ng2,eqv,ratio
c  initialize
      ng2=ng-ng1
      do 61 i=1,n
         call tst2GM(d(i,1),ng1,ng2,n,z(i),dm(i), eqv, ratio)
 61   continue
      return
      end

      subroutine tst2GM(d,ng1,ng2,n,tst,dm,eqv, ratio)
c columns 1 to ng1 in group 1, ng1+1 to ng1+ng2 in group 2
      real d(n,ng1+ng2),tst,dm
      double precision dm1,dm2,dss1,dss2
      integer ng1,ng2,n,i,eqv, ratio
      dm1=0
      dm2=0
      dss1=0
      dss2=0
      do 10 i=1,ng1
         dm1=dm1+d(1,i)
 10   continue
      dm1=dm1/ng1
      do 11 i=1,ng1
         dss1=dss1+(d(1,i)-dm1)**2
 11   continue
      do 12 i=1,ng2
         dm2=dm2+d(1,ng1+i)
 12   continue
      dm2=dm2/ng2
      do 13 i=1,ng2
         dss2=dss2+(d(1,ng1+i)-dm2)**2
 13   continue
      if( ratio.eq.0) then
         dm=dm1-dm2
      endif
      if( ratio.eq.1) then
         dm=dm1/dm2
      endif 
      if (dss1.eq.0.and.dss2.eq.0) then
         tst=0
         return
      endif
c intermediate calculations in dp, so stats with many ties give same sp result
c regardless of order of calculations
      if( eqv .eq. 1 ) then
         tst=(dm1-dm2)/sqrt((1.d0/ng1+1.d0/ng2)*(dss1+dss2)/(ng1+ng2-2))
         return
      endif
      tst=(dm1-dm2)/sqrt(dss1/((ng1-1)*ng1)+dss2/((ng2-1)*ng2))
      end

