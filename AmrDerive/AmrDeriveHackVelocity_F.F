#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDeriveHackVelocity_F.H"

#if (BL_SPACEDIM==3)
      
      subroutine FORT_HACKVELOCITY(dat,DIMS(dat),outDat,DIMS(outDat),xlo,xhi,
     &                             delta,nVars,fSize)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(outDat)
      integer nVars, fSize
      REAL_T  xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      REAL_T  delta(BL_SPACEDIM)
      REAL_T  dat(DIMV(dat),nVars)
      REAL_T  outDat(DIMV(outDat),nVars)
      REAL_T  r, x, y
      
      integer i, j, k, n
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer iproblo, jproblo, kproblo
      integer iprobhi, jprobhi, kprobhi
      integer isioproc, nRho
      REAL_T  vol, sum
      
      call bl_pd_is_ioproc(isioproc)
      
      iproblo = fSize
      jproblo = fSize
      kproblo = fSize

      iprobhi = (xhi(1)-xlo(1))/delta(1) - fSize - 1
      jprobhi = (xhi(2)-xlo(2))/delta(2) - fSize - 1
      kprobhi = (xhi(3)-xlo(3))/delta(3) - fSize - 1

      ilo = max(outDat_l1,iproblo)
      jlo = max(outDat_l2,jproblo)
      klo = max(outDat_l3,kproblo)
      
      ihi = min(outDat_h1,iprobhi)
      jhi = min(outDat_h2,jprobhi)
      khi = min(outDat_h3,kprobhi)

      do n=1, nVars
         do k = klo, khi
c            z = xlo(3) + (float(k)+0.5d0)*delta(3)
            do j = jlo, jhi
c               y = xlo(2) + (float(j)+0.5d0)*delta(2)
               do i = ilo, ihi
c                  x = xlo(1) + (float(i)+0.5d0)*delta(1)

                  outDat(i,j,k,n) = dat(i,j,k,n)

               end do
            end do
         end do
      end do
      
      end

#else

c     2D version
      
      subroutine FORT_HACKVELOCITY(dat,DIMS(dat),outDat,DIMS(outDat),xlo,xhi,
     &                             delta,nVars,fSize)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(outDat)
      integer nVars, fSize
      REAL_T  xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      REAL_T  delta(BL_SPACEDIM)
      REAL_T  dat(DIMV(dat),nVars)
      REAL_T  outDat(DIMV(outDat),nVars)
      REAL_T  r, x, y, r0
      REAL_T  rlo_i, rhi_i, rlo_j, rhi_j
      
      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      integer iproblo, jproblo
      integer iprobhi, jprobhi
      integer isioproc, nRho
      REAL_T  vol, sum
      
      call bl_pd_is_ioproc(isioproc)



      


      iproblo = fSize
      jproblo = fSize

      iprobhi = (xhi(1)-xlo(1))/delta(1) - fSize - 1
      jprobhi = (xhi(2)-xlo(2))/delta(2) - fSize - 1

      ilo = outDat_l1
      ihi = outDat_h1
      jlo = outDat_l2
      jhi = outDat_h2



      
      
      do n=1, nVars
         do j = jlo, jhi
            y = xlo(2) + (float(j)+0.5d0)*delta(2)
            do i = ilo, ihi
               x = xlo(1) + (float(i)+0.5d0)*delta(1)
               rlo_i = 0.25 * fSize
               rhi_i = (xhi(1)-xlo(1))/delta(1) - (0.25*fSize) - 1
               rlo_j = 0.25 * fSize
               rhi_j = (xhi(2)-xlo(2))/delta(2) - (0.25*fSize) - 1


               

               
               outDat(i,j,n) = dat(i,j,n) * ((1.0d0 * tanh((i - rlo_i) / 25) * tanh((j - rlo_j) / 25)) * (1.0d0 * tanh((i - rlo_i) / 25) * tanh((rhi_j - j) / 25)) *  (1.0d0 * tanh((rhi_i - i) / 25) * tanh((j -rlo_j) / 25)) * (1.0d0 * tanh((rhi_i - i) / 25) * tanh((rhi_j - r0) / 25)))
             
                 
               if (tanh(i - rlo_i) .LT. 0) then
                  outDat(i,j,n) = 0
               end if
               if (tanh(j - rlo_j) .LT. 0) then
                  outDat(i,j,n) = 0
               end if
               if (tanh(rhi_i - i) .LT. 0) then
                   outDat(i,j,n) = 0
               end if
               if (tanh(rhi_j - j) .LT. 0) then
                   outDat(i,j,n) = 0
               end if




c     Hack to only do top and bottom, no grad
               
c               rlo_j = fSize
c               rhi_j = (xhi(2)-xlo(2))/delta(2) - (fSize) - 1
c               outDat(i,j,n) = dat(i,j,n)

c               if (j .LT. rlo_j) then
c                  outDat(i,j,n) = 0
c               end if
c               if (j .GT. rhi_j) then
c                  outDat(i,j,n) = 0
c               end if
               
            end do
         end do
      end do
      

      
      end

#endif

      
