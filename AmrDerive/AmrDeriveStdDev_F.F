#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDeriveStdDev_F.H"

#define SDIM 3

      subroutine FORT_INTEGRATE(dat,nv,DIMS(dat),DIMS(grid),delta,turb,
     &                          isize,jsize,ksize,tv,axis,size)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer isize, jsize, ksize, nv, tv, axis, size
      REAL_T  delta(3)
      REAL_T  dat(DIMV(dat),0:nv)
      REAL_T  turb(0:size*tv-1)

      integer i, j, k, v, n, iteraxis
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
      if (axis == 0) then
        area = dy*dz
      else if (axis == 1) then
        area = dx*dz
      else if (axis == 2) then
        area = dx*dy
      endif

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi

      !c     Check intersect flag
               if (dat(i,j,k,nv).gt.1e-8) then
                 if (axis == 0) then
                   iteraxis = i
                 else if (axis == 1) then
                   iteraxis = j
                 else if (axis == 2) then
                   iteraxis = k
                 endif
                 do n = 0, nv-1
                    turb(iteraxis*tv+n)=turb(iteraxis*tv+n)+area*dat(i,j,k,n)
                 end do

               endif
            end do
         end do
      end do

      end

      subroutine FORT_STDDEV(dat,nv,DIMS(dat),DIMS(grid),delta,turb,
     &                          isize,jsize,ksize,tv,axis,size,mean)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer isize, jsize, ksize, nv, tv, axis, size
      REAL_T  delta(3)
      REAL_T  dat(DIMV(dat),0:nv)
      REAL_T  turb(0:size*tv-1)
      REAL_T  mean(0:size*tv-1)


      integer i, j, k, v, n, iteraxis
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
      if (axis == 0) then
        area = dy*dz
      else if (axis == 1) then
        area = dx*dz
      else if (axis == 2) then
        area = dx*dy
      endif

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
              if (dat(i,j,k,nv).gt.1e-8) then
               if (axis == 0) then
                 iteraxis = i
               else if (axis == 1) then
                 iteraxis = j
               else if (axis == 2) then
                 iteraxis = k
               endif

               do n = 0, nv-1
                  turb(iteraxis*tv+n)=turb(iteraxis*tv+n)+(dat(i,j,k,n)-mean(iteraxis*tv+n))*(dat(i,j,k,n)-mean(iteraxis*tv+n))
               end do
             endif
            end do
         end do
      end do

      end
