/* -------------------------------------------------------------------- */
/*  CoarsenPLT_F.F                                                      */
/* -------------------------------------------------------------------- */
#include <REAL.H>
#include <ArrayLim.H>
#include <CONSTANTS.H>
#include <AmrDeriveCoarsen_F.H>
#include <BC_TYPES.H>

#if (BL_SPACEDIM==3)
      
#define SDIM 3

      subroutine FORT_AVGDN(lo, hi, c, DIMS(c), f, DIMS(f), ratio, nc)
      implicit none
      integer lo(SDIM), hi(SDIM), ratio, nc
      integer DIMDEC(c)
      integer DIMDEC(f)
      REAL_T c(DIMV(c),nc)
      REAL_T f(DIMV(f),nc)
      integer i,j,k,n,i2,j2,k2

      if (ratio.ne.2) then
         write(6,*) 'Need to generalize avgdown to handle r != 2'
         stop
      end if

      do k=lo(3),hi(3)
         k2 = 2*k
         do j=lo(2),hi(2)
            j2 = 2*j
            do i=lo(1),hi(1)
               i2 = 2*i
               do n=1,nc
                  c(i,j,k,n) = 0.125d0*
     &                 (f(i2,j2,  k2,  n)+f(i2+1,j2,  k2,  n)
     &                 +f(i2,j2+1,k2,  n)+f(i2+1,j2+1,k2,  n)
     &                 +f(i2,j2,  k2+1,n)+f(i2+1,j2,  k2+1,n)
     &                 +f(i2,j2+1,k2+1,n)+f(i2+1,j2+1,k2+1,n))
               end do
            end do
         end do
      end do
      end
      
#else

#define SDIM 2

      subroutine FORT_AVGDN(lo, hi, c, DIMS(c), f, DIMS(f), ratio, nc)
      implicit none
      integer lo(SDIM), hi(SDIM), ratio, nc
      integer DIMDEC(c)
      integer DIMDEC(f)
      REAL_T c(DIMV(c),nc)
      REAL_T f(DIMV(f),nc)
      integer i,j,n,i2,j2

      if (ratio.ne.2) then
         write(6,*) 'Need to generalize avgdown to handle r != 2'
         stop
      end if

      do j=lo(2),hi(2)
         j2 = 2*j
         do i=lo(1),hi(1)
            i2 = 2*i
            do n=1,nc
               c(i,j,n) = 0.25d0*
     &              (f(i2,j2,  n)+f(i2+1,j2,  n)
     &              +f(i2,j2+1,n)+f(i2+1,j2+1,n))
            end do
         end do
      end do

      end

#endif
      
