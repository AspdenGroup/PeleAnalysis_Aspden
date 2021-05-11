#include "REAL.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"

#include "AmrDeriveHit_F.H"

#define SDIM 3

      subroutine FORT_HIT (dat,DIMS(dat),DIMS(grid),delta,Avgs,PdfData,PdfMin,PdfMax)

      implicit none

      integer DIMDEC(dat)
      integer DIMDEC(grid)
      REAL_T  delta(SDIM)
#ifdef DO_FORCING
      REAL_T  dat(DIMV(dat), SDIM*2)
#else
      REAL_T  dat(DIMV(dat), SDIM)
#endif
      REAL_T  Avgs(0:BARVARS-1)
      integer  PdfData(0:12*PDFBINS-1)
      REAL_T  PdfMin(0:11)
      REAL_T  PdfMax(0:11)

      integer ilo, ihi, jlo, jhi, klo, khi
      integer i, j, k, isioproc
      REAL_T  dx, dy, dz, vol
      REAL_T  dx2, dy2, dz2
      REAL_T  ux, uy, uz
      REAL_T  fx, fy, fz
      REAL_T  D, D2
      REAL_T  SXX, SYY, SZZ, SXY, SYZ, SZX
      REAL_T  UXX, UYY, UZZ
      REAL_T  VXX, VYY, VZZ
      REAL_T  WXX, WYY, WZZ
      REAL_T  q(0:11)
      integer bin, qi, bin2

c      write (*,*) "In FORT_HIT"

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

c      write (*,*) "delta: ",dx,dy,dz

      dx2= dx*dx
      dy2= dy*dy
      dz2= dz*dz
      
      vol = dx*dy*dz

      ilo = ARG_L1(grid)
      ihi = ARG_H1(grid)
      jlo = ARG_L2(grid)
      jhi = ARG_H2(grid)
      klo = ARG_L3(grid)
      khi = ARG_H3(grid)

c     Calculate integrals
      do i = 0, BARVARS-1
         Avgs(i) = zero
      end do

      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
               ux = dat(i,j,k,1)
               uy = dat(i,j,k,2)
               uz = dat(i,j,k,3)
#ifdef DO_FORCING
               fx = dat(i,j,k,4)
               fy = dat(i,j,k,5)
               fz = dat(i,j,k,6)
#endif

c     Strain rate tensor components (the long way - should cancel lots of 'two's)
               SXX =   (dat(i+1,j,k,1)-dat(i-1,j,k,1))/(two*dx)
               SYY =   (dat(i,j+1,k,2)-dat(i,j-1,k,2))/(two*dy)
               SZZ =   (dat(i,j,k+1,3)-dat(i,j,k-1,3))/(two*dz)
               SXY = ( (dat(i,j+1,k,1)-dat(i,j-1,k,1))/(two*dy) + (dat(i+1,j,k,2)-dat(i-1,j,k,2))/(two*dx) ) / two
               SYZ = ( (dat(i,j,k+1,2)-dat(i,j,k+1,2))/(two*dz) + (dat(i,j+1,k,3)-dat(i,j+1,k,3))/(two*dy) ) / two
               SZX = ( (dat(i+1,j,k,3)-dat(i-1,j,k,3))/(two*dx) + (dat(i,j,k+1,1)-dat(i,j,k+1,1))/(two*dz) ) / two
               D   = ( SXX*SXX + SYY*SYY + SZZ*SZZ ) + two * ( SXY*SXY + SYZ*SYZ + SZX*SZX ) 
c     Laplacians
               UXX = (dat(i+1,j,k,1)-two*dat(i,j,k,1)+dat(i-1,j,k,1)) / dx2
               UYY = (dat(i,j+1,k,1)-two*dat(i,j,k,1)+dat(i,j-1,k,1)) / dy2
               UZZ = (dat(i,j,k+1,1)-two*dat(i,j,k,1)+dat(i,j,k-1,1)) / dz2
               VXX = (dat(i+1,j,k,2)-two*dat(i,j,k,2)+dat(i-1,j,k,2)) / dx2
               VYY = (dat(i,j+1,k,2)-two*dat(i,j,k,2)+dat(i,j-1,k,2)) / dy2
               VZZ = (dat(i,j,k+1,2)-two*dat(i,j,k,2)+dat(i,j,k-1,2)) / dz2
               WXX = (dat(i+1,j,k,3)-two*dat(i,j,k,3)+dat(i-1,j,k,3)) / dx2
               WYY = (dat(i,j+1,k,3)-two*dat(i,j,k,3)+dat(i,j-1,k,3)) / dy2
               WZZ = (dat(i,j,k+1,3)-two*dat(i,j,k,3)+dat(i,j,k-1,3)) / dz2
               D2  = ux*(UXX+UYY+UZZ) + uy*(VXX+VYY+VZZ) + uz*(WXX+WYY+WZZ)
c     Averages
               Avgs(UXBAR)  = Avgs(UXBAR)  + ux
               Avgs(UYBAR)  = Avgs(UYBAR)  + uy
               Avgs(UZBAR)  = Avgs(UZBAR)  + uz
               Avgs(UX2BAR) = Avgs(UX2BAR) + ux*ux
               Avgs(UY2BAR) = Avgs(UY2BAR) + uy*uy
               Avgs(UZ2BAR) = Avgs(UZ2BAR) + uz*uz
               Avgs(DBAR)   = Avgs(DBAR)   + D
               Avgs(D2BAR)  = Avgs(D2BAR)  + D2
#ifdef DO_FORCING
               Avgs(FBAR)   = Avgs(FBAR)   + ux*fx + uy*fy + uz*fz
#endif

c     Pdfs
c     Velocities
               q(0)  = ux
               q(1)  = uy
               q(2)  = uz
c     Gradient of u_x
               q(3)  = (dat(i+1,j,k,1)-dat(i-1,j,k,1))/(two*dx)
               q(4)  = (dat(i,j+1,k,1)-dat(i,j-1,k,1))/(two*dy)
               q(5)  = (dat(i,j,k+1,1)-dat(i,j,k-1,1))/(two*dz)
c     Gradients of u_y
               q(6)  = (dat(i+1,j,k,2)-dat(i-1,j,k,2))/(two*dx)
               q(7)  = (dat(i,j+1,k,2)-dat(i,j-1,k,2))/(two*dy)
               q(8)  = (dat(i,j,k+1,2)-dat(i,j,k-1,2))/(two*dz)
c     Gradients of u_z
               q(9)  = (dat(i+1,j,k,3)-dat(i-1,j,k,3))/(two*dx)
               q(10) = (dat(i,j+1,k,3)-dat(i,j-1,k,3))/(two*dy)
               q(11) = (dat(i,j,k+1,3)-dat(i,j,k-1,3))/(two*dz)

               do qi = 0, 11
                  bin = max(0,min(PDFBINS-1,int( (q(qi)-PdfMin(qi))*(dfloat(PDFBINS)/(PdfMax(qi)-PdfMin(qi))) )))
                  PdfData(qi*PDFBINS+bin) = PdfData(qi*PDFBINS+bin) + 1
               enddo

            end do
         end do
      end do

      end

















#if 0
      call bl_pd_is_ioproc(isioproc)
     
      if (isioproc) then
         write (*,*) "lo: ", ilo, jlo, klo
         write (*,*) "hi: ", ihi, jhi, khi
         write (*,*) "Grid is valid region"
         write (*,*) "grid1: ", ARG_L1(grid), ARG_H1(grid)
         write (*,*) "grid2: ", ARG_L2(grid), ARG_H2(grid)
         write (*,*) "grid3: ", ARG_L3(grid), ARG_H3(grid)
         write (*,*) "Dat is array dimension"
         write (*,*) "dat1: ", ARG_L1(dat), ARG_H1(dat)
         write (*,*) "dat2: ", ARG_L2(dat), ARG_H2(dat)
         write (*,*) "dat3: ", ARG_L3(dat), ARG_H3(dat)
      endif
#endif




#if 0
      REAL_T  dudx, dudy, dudz
      REAL_T  dvdx, dvdy, dvdz
      REAL_T  dwdx, dwdy, dwdz
      REAL_T  SZY, SYX, SXZ
      
               dudx = ( dat(i+1,j,k,1) - dat(i-1,j,k,1) ) / (two*dx)
               dudy = ( dat(i,j+1,k,1) - dat(i,j-1,k,1) ) / (two*dy)
               dudz = ( dat(i,j,k+1,1) - dat(i,j,k-1,1) ) / (two*dz)
               dvdx = ( dat(i+1,j,k,2) - dat(i-1,j,k,2) ) / (two*dx)
               dvdy = ( dat(i,j+1,k,2) - dat(i,j-1,k,2) ) / (two*dy)
               dvdz = ( dat(i,j,k+1,2) - dat(i,j,k-1,2) ) / (two*dz)
               dwdx = ( dat(i+1,j,k,3) - dat(i-1,j,k,3) ) / (two*dx)
               dwdy = ( dat(i,j+1,k,3) - dat(i,j-1,k,3) ) / (two*dy)
               dwdz = ( dat(i,j,k+1,3) - dat(i,j,k-1,3) ) / (two*dz)
               SXX  = ( dudx + dudx ) / two
               SXY  = ( dvdx + dudy ) / two
               SXZ  = ( dwdx + dudz ) / two
               SYX  = ( dudy + dvdx ) / two
               SYY  = ( dvdy + dvdy ) / two
               SYZ  = ( dwdy + dvdz ) / two
               SZX  = ( dudz + dwdx ) / two
               SZY  = ( dvdz + dwdy ) / two
               SZZ  = ( dwdz + dwdz ) / two
               D    = SXX*SXX + SXY*SXY + SXZ*SXZ + SYX*SYX + SYY*SYY + SYZ*SYZ + SZX*SZX + SZY*SZY + SZZ*SZZ
#endif
