#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDerive_udot_grad2u_F.H"

#define SDIM 3

      subroutine FORT_INTEGRATE(dat,nv,DIMS(dat),DIMS(grid),delta,turb,
     &                          isize,jsize,ksize,tv)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(grid)
      integer isize, jsize, ksize, nv, tv
      REAL_T  delta(3)    !array of cell sizes (dx dy dz)
      REAL_T  dat(DIMV(dat),0:nv)
      REAL_T  turb(0:ksize*tv-1)

      integer i, j, k, v, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi

      REAL_T  area, isect
      REAL_T  dx, dy, dz
      REAL_T  gradux, gradvx, gradwx
      REAL_T  graduy, gradvy, gradwy
      REAL_T  graduz, gradvz, gradwz
      REAL_T  lapu, lapv, lapw


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

      do k = klo, khi !z
         do j = jlo, jhi !y
            do i = ilo, ihi !x
               gradux = (dat(i+1,j,k,0) - (2 * dat(i,j,k,0)) + dat(i-1,j,k,0)) / (dx * dx)
               graduy = (dat(i,j+1,k,0) - (2 * dat(i,j,k,0)) + dat(i,j-1,k,0)) / (dy * dy)
               graduz = (dat(i,j,k+1,0) - (2 * dat(i,j,k,0)) + dat(i,j,k-1,0)) / (dz * dz)
               gradvx = (dat(i+1,j,k,1) - (2 * dat(i,j,k,1)) + dat(i-1,j,k,1)) / (dx * dx)
               gradvy = (dat(i,j+1,k,1) - (2 * dat(i,j,k,1)) + dat(i,j-1,k,1)) / (dy * dy)
               gradvz = (dat(i,j,k+1,1) - (2 * dat(i,j,k,1)) + dat(i,j,k-1,1)) / (dz * dz)
               gradwx = (dat(i+1,j,k,2) - (2 * dat(i,j,k,2)) + dat(i-1,j,k,2)) / (dx * dx)
               gradwy = (dat(i,j+1,k,2) - (2 * dat(i,j,k,2)) + dat(i,j-1,k,2)) / (dy * dy)
               gradwz = (dat(i,j,k+1,2) - (2 * dat(i,j,k,2)) + dat(i,j,k-1,2)) / (dz * dz)
               lapu = gradux + graduy + graduz
               lapv = gradvx + gradvy + gradvz
               lapw = gradwx + gradwy + gradwz

               turb(k*tv) = turb(k*tv)+(area * ((dat(i,j,k,0) * lapu) + (dat(i,j,k,1) * lapv) + (dat(i,j,k,2) * lapw)))

            end do
         end do
         print *, area*(dat(ihi,jhi,k,0)*lapu + dat(ihi,jhi,k,1)*lapv + dat(ihi,jhi,k,2)*lapw)
      end do

      end
