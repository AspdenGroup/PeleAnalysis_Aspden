#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDeriveProgIntegration_F.H"
      
#if (BL_SPACEDIM==2)
#define SDIM 2
      subroutine FORT_INTEGRATE(dat,nv,DIMS(dat),DIMS(grid),delta,turb,
     &     isize,jsize,tv)
   
      implicit none
      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer isize, jsize, nv, tv
      REAL_T  delta(2)
      REAL_T  dat(DIMV(dat),0:nv)
      REAL_T  turb(0:jsize*tv-1)
      integer i, j, v, n
      integer ilo, jlo
      integer ihi, jhi
      REAL_T  area, isect
      REAL_T  dx, dy
      
      ilo = grid_l1
      ihi = grid_h1
      jlo = grid_l2
      jhi = grid_h2
 
      dx = delta(1)
      dy = delta(2)

      area = dx

      do j = jlo, jhi
         do i = ilo, ihi
c     Check intersect flag
            if (dat(i,j,nv).gt.1e-8) then
               do n = 0, tv-1
                  if (dat(i,j,0).ge.0.7) then
                     turb(j*tv+n)=turb(j*tv+n)+area*1 
                  endif
               end do
            endif
         end do
      end do
      end
#else
#define SDIM 3

      subroutine FORT_INTEGRATE(dat,nv,DIMS(dat),DIMS(grid),delta,turb,
     &                          isize,jsize,ksize,tv)

      implicit none
      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer isize, jsize, ksize, nv, tv
      REAL_T  delta(3)
      REAL_T  dat(DIMV(dat),0:nv)
      REAL_T  turb(0:ksize*tv-1)
      
      integer i, j, k, v, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      
      REAL_T  area, isect
      REAL_T  dx, dy, dz
      
      ilo = grid_l1
      ihi = grid_h1
      jlo = grid_l2
      jhi = grid_h2
      klo = grid_l3
      khi = grid_h3
      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      area = dx*dy

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
c     Check intersect flag
               if (dat(i,j,k,nv).gt.1e-8) then
                  do n = 0, tv-1
                     if (dat(i,j,k,0).ge.0.9) then
                        turb(k*tv+n)=turb(k*tv+n)+area*1
                     endif
                  end do
               endif
            end do
         end do
      end do
      end
#endif
