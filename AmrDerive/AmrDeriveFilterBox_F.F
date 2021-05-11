#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDeriveFilter_F.H"

      subroutine FORT_FILTER(dat,DIMS(dat),outDat,DIMS(outDat),xlo,xhi,
     &                       delta,nVars,fSize,favre,denVar)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(outDat)
      integer nVars, fSize, denVar
      integer favre(nVars)
      REAL_T  xlo(BL_SPACEDIM), xhi(BL_SPACEDIM)
      REAL_T  delta(BL_SPACEDIM)
      REAL_T  dat(DIMV(dat),nVars)
      REAL_T  outDat(DIMV(outDat),nVars)

      integer i, j, k, n
      integer ii, jj, kk
      integer ilo, jlo, klo
      integer ihi, jhi, khi
      integer isioproc, nRho
      integer alternate
      REAL_T  vol, numFiltersx, numFiltersy, numFiltersz

      call bl_pd_is_ioproc(isioproc)

      ilo = outDat_l1
      ihi = outDat_h1
      jlo = outDat_l2
      jhi = outDat_h2
      klo = outDat_l3
      khi = outDat_h3

      vol = 1

c     adjust for fortran indexing
      nRho = denVar + 1

      do n=1, nVars

         if (favre(n)>0) then
c     do Favre averaging for this variables
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                    numFiltersx = abs(i)/fSize
                    numFiltersy = abs(j)/fSize
                    numFiltersz = abs(k)/fSize
                    alternate = floor(numFiltersx) + floor(numFiltersy) + floor(numFiltersz)
                    outDat(i,j,k,n) = mod(alternate,2)
c     do i, j, k
                  end do
               end do
            end do
c     favre
         else
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                    numFiltersx = abs(i)/fSize
                    numFiltersy = abs(j)/fSize
                    numFiltersz = abs(k)/fSize
                    alternate = floor(numFiltersx) + floor(numFiltersy) + floor(numFiltersz)
                    outDat(i,j,k,n) = mod(alternate,2)
c     do i, j, k
                  end do
               end do
            end do
c     favre
         endif
c     do n
      end do

c     Evaluate mean
c     do Favre averaging first (i.e. before changing density)
      do n=1, nVars
         if (favre(n).gt.0) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     outDat(i,j,k,n) = outDat(i,j,k,n) / outDat(i,j,k,nRho)
                  end do
               end do
            end do
         endif
      end do
c     now do non-Favre quantities
      do n=1, nVars
         if (favre(n).eq.0) then
            do k = klo, khi
               do j = jlo, jhi
                  do i = ilo, ihi
                     outDat(i,j,k,n) = outDat(i,j,k,n) / vol
                  end do
               end do
            end do
         endif
      end do

      end





























#if 0
      do n=1, nVars
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi

c     time-saving looping strategy
                  if (i.eq.ilo) then
                     if (j.eq.jlo) then
                        if (k.eq.klo) then

c     (0,0,0)
c     Very first cell, compute full filter
                           sum = 0.0d0
                           do kk = k-fSize,k+fSize
                              do jj = j-fSize,j+fSize
                                 do ii = i-fSize,i+fSize
                                    sum = sum + dat(ii,jj,kk,n)
                                 end do
                              end do
                           end do
c     if k=0
                        else
c     (0,0,k)
c     loop over i-j plane
                           sum = outDat(i,j,k-1,n)
                           do jj = j-fSize,j+fSize
                              do ii = i-fSize,i+fSize
                                 sum = sum - dat(ii,jj,k-fSize-1,n) + dat(ii,jj,k+fSize,n)
                              end do
                           end do
c     if k=0
                        end if
c     if j=0
                     else
c     (0,j,k)
c     loop over i-k plane
                        sum = outDat(i,j-1,k,n)
                        do kk = k-fSize,k+fSize
                           do ii = i-fSize,i+fSize
                              sum = sum - dat(ii,j-fSize-1,kk,n) + dat(ii,j+fSize,kk,n)
                           end do
                        end do
c     if j=0
                     end if
c     if i=0
                  else
c     (i,j,k)
c     loop over j-k plane
                     sum = outDat(i-1,j,k,n)
                     do kk = k-fSize,k+fSize
                        do jj = j-fSize,j+fSize
                           sum = sum - dat(i-fSize-1,jj,kk,n) + dat(i+fSize,jj,kk,n)
                        end do
                     end do
c     if k=0
                  end if

c     post integration - copy sum to outDat
                  outDat(i,j,k,n) = sum

c     do i
               end do
c     do j
            end do
c     do k
         end do
c     do n
      end do

#endif
