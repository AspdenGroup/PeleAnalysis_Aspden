#include "CONSTANTS.H"
#include "ArrayLim.H"
#include "BC_TYPES.H"
#include "SPACE.H"
#include "xtra_F.H"

#define DIM BL_SPACEDIM

c
c     Higher order extrap
c
      subroutine FORT_PUSHVTOGHO(dlo, dhi, U, DIMS(U), nc)
      implicit none
      integer nc
      integer dlo(DIM), dhi(DIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U),nc)
      integer i,j,k,n
      integer bc(DIM,2)
      REAL_T xlo(DIM), dx(DIM)
c
c     Make up something for these that gets what we want
c
      do n=1,DIM
         bc(n,1) = HOEXTRAP
         bc(n,2) = HOEXTRAP
         dx(n) = one
         xlo(n) = one
      end do

      do n=1,nc
         call filcc(U(D_DECL(ARG_L1(U),ARG_L2(U),ARG_L3(U)),n),DIMS(U),dlo,dhi,dx,xlo,bc)
      end do

      end
c
c     First order extrap
c
      subroutine FORT_PUSHVTOGFO(dlo, dhi, U, DIMS(U), nc)
      implicit none
      integer nc
      integer dlo(DIM), dhi(DIM)
      integer DIMDEC(U)
      REAL_T U(DIMV(U),nc)
      integer i,j,k,n
      integer bc(DIM,2)
      REAL_T xlo(DIM), dx(DIM)
c
c     Make up something for these that gets what we want
c

      do n=1,DIM
         bc(n,1) = FOEXTRAP
         bc(n,2) = FOEXTRAP
         dx(n) = one
         xlo(n) = one
      end do

      do n=1,nc
         call filcc(U(D_DECL(ARG_L1(U),ARG_L2(U),ARG_L3(U)),n),DIMS(U),dlo,dhi,dx,xlo,bc)
      end do

      end
