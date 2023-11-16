module theta_module
  
  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  use amrex_bc_types_module
  use amrex_filcc_module
  
  implicit none

  public

contains

  subroutine FORT_INTEGRATE(dat,dlo,dhi,lo,hi,xlo,delta,turb,rsize,levRsize,ksize,levKsize,nv,tv) bind(c,name='FORT_INTEGRATE')
    implicit none

    integer, intent(in) :: rsize, levRsize,ksize,levKsize, nv, tv, lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(in) :: xlo(3)
    real(amrex_real), intent(in) :: delta(3)
    real(amrex_real), intent(in) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nv+1)
    real(amrex_real), intent(inout) :: turb(0:ksize*rsize*tv-1)
    
    real(amrex_real) :: gridDx
    
    integer :: i, j, k, v
    integer :: ilo, jlo, klo
    integer :: ihi, jhi, khi
    
    real(amrex_real) ::  rho, ux, uy, uz, YH2, YO2, YN2
    real(amrex_real) ::  drhodx, duxdx, duydx, duzdx, dtempdx, dYH2dx, dYO2dx, dYN2dx
    real(amrex_real) ::  drhody, duxdy, duydy, duzdy, dtempdy, dYH2dy, dYO2dy, dYN2dy
    real(amrex_real) ::  drhodz, duxdz, duydz, duzdz, dtempdz, dYH2dz, dYO2dz, dYN2dz
    real(amrex_real) :: dx, dy, dz, tdx, tdy, tdz, dxt, dyt, dzt
    real(amrex_real) :: d2uxdx2, d2uydy2, d2uzdz2
    
    integer :: ii,jj,kk,kn,nn,numSplit,kklo,kkhi,ridx, idx
    real(amrex_real) :: dxnn,dynn,dznn
    real(amrex_real) :: xlow,ylow,zlow,xctr,yctr,zctr,xx,yy,zz,rr
    real(amrex_real) :: rhoctr, uxctr, uyctr, uzctr, tempctr
    real(amrex_real) ::YH2ctr, YO2ctr, YN2ctr
    
    real(amrex_real) :: deltax, deltay, deltaz
    real(amrex_real) :: ur, ut, temp
      
    real(amrex_real) :: mdx, mdy
    real(amrex_real) :: drx, dry, dtx, dty
    real(amrex_real) :: xxp, yyp, xxm, yym, rrp, rrm
    real(amrex_real) :: uxp, uyp, uxm, uym, urp, urm, utp, utm
    real(amrex_real) :: durdr, dutdr, durdt, dutdt
    real(amrex_real) :: vol
    integer :: idx_rho, idx_ux, idx_uy, idx_uz, idx_temp, idx_YH2, idx_YO2, idx_YN2

    integer :: isioproc

    real(amrex_real) :: rhoMin, rhoMax
      
    !call bl_pd_is_ioproc(isioproc)

    !isioproc=0

    !if (isioproc.eq.1) then 
    !   write (6,*) "---"
    !   write (6,*) grid_l1,grid_h1
    !   write (6,*) grid_l2,grid_h2
    !   write (6,*) grid_l3,grid_h3
    !   write (6,*) dat_l1,dat_h1
    !   write (6,*) dat_l2,dat_h2
    !   write (6,*) dat_l3,dat_h3
    !   write (6,*) rsize, ksize
    !   write (6,*) levRsize, levKsize
    !   write (6,*) nv,tv
    !endif

    ilo = lo(1)
    ihi = hi(1)
    jlo = lo(2)
    jhi = hi(2)
    klo = lo(3)
    khi = hi(3)
      
    dx = delta(1)
    dy = delta(2)
    dz = delta(3)
    
    tdx = 2.0*dx
    tdy = 2.0*dy
    tdz = 2.0*dz
    
    dxt = dx*dx
    dyt = dy*dy
    dzt = dz*dz
    
    numSplit = 2
    kn   = ksize/levKsize
    nn   = numSplit*kn
    dxnn = dx/float(nn)
    dynn = dy/float(nn)
    !     Note this is right, it really should be kn not nn
    dznn = dz/float(kn)
    !     Mod dx
    mdx = sqrt(dxnn*dxnn+dynn*dynn)
    vol = dxnn*dynn*dznn
    gridDx=dx/dfloat(kn)
    
    idx_ux   = 1
    idx_uy   = 2
    idx_uz   = 3
    idx_rho  = 4
    idx_temp = 5
    idx_YH2 = 6
    idx_YO2 = 7
    idx_YN2 = 8
    
    
    rhoMin = 1.d100
    rhoMax = -rhoMin
    do k = dlo(3), dhi(3)
       do j = dlo(2), dhi(2)
          do i = dlo(1), dhi(1)
             rho = dat(i,j,k,idx_rho)
             if (rho.lt.rhoMin) then
                rhoMin = rho
             endif
             if (rho.gt.rhoMax) then
                rhoMax = rho
             endif
          enddo
       enddo
    enddo
    do k = klo, khi
!     Bounds for vertical numerical integration - confused about this
       kklo=k*kn
       kkhi=kklo+kn-1
       !     Calculate left hand edge of cell
         zlow = xlo(3) + dz*dfloat(k-klo)
!     And x centre/er
         zctr = zlow + 0.5*dz
                     
         do j = jlo, jhi
!     Calculate front of cell
            ylow = xlo(2) + dy*dfloat(j-jlo)
!     And y centre/er
            yctr = ylow + 0.5*dy
            
            do i = ilo, ihi
!     Calculate left hand edge of cell
               xlow = xlo(1) + dx*dfloat(i-ilo)
!     And x centre/er
               xctr = xlow + 0.5*dx
               
!     Check intersect flag
               if (dat(i,j,k,nv+1).gt.1e-8) then
                  uxctr   = dat(i,j,k,idx_ux)
                  uyctr   = dat(i,j,k,idx_uy)
                  uzctr   = dat(i,j,k,idx_uz)
                  rhoctr  = dat(i,j,k,idx_rho)
                  tempctr = dat(i,j,k,idx_temp)
                  !YH2ctr = dat(i,j,k,idx_YH2)
                  !YO2ctr = dat(i,j,k,idx_YO2)
                  !YN2ctr = dat(i,j,k,idx_YN2)
                  duxdx   = (dat(i+1,j,k,idx_ux)   - dat(i-1,j,k,idx_ux)  )/tdx
                  duydx   = (dat(i+1,j,k,idx_uy)   - dat(i-1,j,k,idx_uy)  )/tdx
                  duzdx   = (dat(i+1,j,k,idx_uz)   - dat(i-1,j,k,idx_uz)  )/tdx
                  drhodx  = (dat(i+1,j,k,idx_rho)  - dat(i-1,j,k,idx_rho) )/tdx
                  dtempdx = (dat(i+1,j,k,idx_temp) - dat(i-1,j,k,idx_temp))/tdx
                  !dYH2dx  = (dat(i+1,j,k,idx_YH2)  - dat(i-1,j,k,idx_YH2) )/tdx
                  !dYO2dx  = (dat(i+1,j,k,idx_YO2)  - dat(i-1,j,k,idx_YO2) )/tdx
                  !dYN2dx  = (dat(i+1,j,k,idx_YN2)  - dat(i-1,j,k,idx_YN2) )/tdx

                  duxdy   = (dat(i,j+1,k,idx_ux)   - dat(i,j-1,k,idx_ux)  )/tdy
                  duydy   = (dat(i,j+1,k,idx_uy)   - dat(i,j-1,k,idx_uy)  )/tdy
                  duzdy   = (dat(i,j+1,k,idx_uz)   - dat(i,j-1,k,idx_uz)  )/tdy
                  drhody  = (dat(i,j+1,k,idx_rho)  - dat(i,j-1,k,idx_rho) )/tdy
                  dtempdy = (dat(i,j+1,k,idx_temp) - dat(i,j-1,k,idx_temp))/tdy
                  !dYH2dy  = (dat(i,j+1,k,idx_YH2)  - dat(i,j-1,k,idx_YH2) )/tdy
                  !dYO2dy  = (dat(i,j+1,k,idx_YO2)  - dat(i,j-1,k,idx_YO2) )/tdy
                  !dYN2dy  = (dat(i,j+1,k,idx_YN2)  - dat(i,j-1,k,idx_YN2) )/tdy

                  duxdz   = (dat(i,j,k+1,idx_ux)   - dat(i,j,k-1,idx_ux)  )/tdz
                  duydz   = (dat(i,j,k+1,idx_uy)   - dat(i,j,k-1,idx_uy)  )/tdz
                  duzdz   = (dat(i,j,k+1,idx_uz)   - dat(i,j,k-1,idx_uz)  )/tdz
                  drhodz  = (dat(i,j,k+1,idx_rho)  - dat(i,j,k-1,idx_rho) )/tdz
                  dtempdz = (dat(i,j,k+1,idx_temp) - dat(i,j,k-1,idx_temp))/tdz
                  !dYH2dz  = (dat(i,j,k+1,idx_YH2)  - dat(i,j,k-1,idx_YH2) )/tdz
                  !dYO2dz  = (dat(i,j,k+1,idx_YO2)  - dat(i,j,k-1,idx_YO2) )/tdz
                  !dYN2dz  = (dat(i,j,k+1,idx_YN2)  - dat(i,j,k-1,idx_YN2) )/tdz

                  d2uxdx2 = (dat(i+1,j,k,idx_ux)-2.0*dat(i,j,k,idx_ux)+dat(i-1,j,k,idx_ux))/dxt
                  d2uydy2 = (dat(i,j+1,k,idx_uy)-2.0*dat(i,j,k,idx_uy)+dat(i,j-1,k,idx_uy))/dyt
                  d2uzdz2 = (dat(i,j,k+1,idx_uz)-2.0*dat(i,j,k,idx_uz)+dat(i,j,k-1,idx_uz))/dzt
                  
!     Integrate numerically
                  do jj=1,nn
                     yy = ylow + dynn*(dfloat(jj)-0.5)
                     
                     do ii=1,nn
                        xx = xlow + dxnn*(dfloat(ii)-0.5)
                        rr = sqrt(xx*xx+yy*yy)
                        
                        ridx = int(rr/(sqrt(2.0)*gridDx))
                        
                        !write (6,*) ridx

                        if (ridx.lt.rsize) then
                           
                           do kk=kklo,kkhi
                              
                              zz = zlow + dznn*(dfloat(kk-kklo)+0.5)
!     How far are we from the centre of the cell?
                              deltax = xx-xctr
                              deltay = yy-yctr
                              deltaz = zz-zctr

!     Do slope reconstruction
                              ux   = uxctr   + deltax*duxdx   + deltay*duxdy   + deltaz*duxdz
                              uy   = uyctr   + deltax*duydx   + deltay*duydy   + deltaz*duydz
                              uz   = uzctr   + deltax*duzdx   + deltay*duzdy   + deltaz*duzdz
                              rho  = rhoctr  + deltax*drhodx  + deltay*drhody  + deltaz*drhodz
                              temp = tempctr + deltax*dtempdx + deltay*dtempdy + deltaz*dtempdz
                              !YH2  = YH2ctr  + deltax*dYH2dx  + deltay*dYH2dy  + deltaz*dYH2dz
                              !YO2  = YO2ctr  + deltax*dYO2dx  + deltay*dYO2dy  + deltaz*dYO2dz
                              !YN2  = YN2ctr  + deltax*dYN2dx  + deltay*dYN2dy  + deltaz*dYN2dz
                              
                              ur = (ux*xx+uy*yy)/rr
                              ut = (uy*xx-ux*yy)/rr

!     Get radial derivative of radial and azimuthal velocity (more complicated than it sounds)
                              drx = mdx*xx/rr
                              dry = mdx*yy/rr
                              xxp = xx + drx
                              yyp = yy + dry
                              rrp = sqrt(xxp*xxp+yyp*yyp)
                              xxm = xx - drx
                              yym = yy - dry
                              rrm = sqrt(xxm*xxm+yym*yym)
                              uxp = uxctr + (xxp-xctr)*duxdx + (yyp-yctr)*duxdy + deltaz*duxdz
                              uyp = uyctr + (xxp-xctr)*duydx + (yyp-yctr)*duydy + deltaz*duydz
                              uxm = uxctr + (xxm-xctr)*duxdx + (yym-yctr)*duxdy + deltaz*duxdz
                              uym = uyctr + (xxm-xctr)*duydx + (yym-yctr)*duydy + deltaz*duydz

                              urp = (uxp*xxp+uyp*yyp)/rrp
                              utp = (uyp*xxp-uxp*yyp)/rrp
                              urm = (uxm*xxm+uym*yym)/rrm
                              utm = (uym*xxm-uxm*yym)/rrm

                              durdr = (urp-urm)/(2.0*mdx)
                              dutdr = (utp-utm)/(2.0*mdx)
!     Get azimuthal derivative of radial and azimuthal velocity (more complicated than it sounds)
                              dtx =-mdx*yy/rr
                              dty = mdx*xx/rr
                              xxp = xx + dtx
                              yyp = yy + dty
                              rrp = sqrt(xxp*xxp+yyp*yyp)
                              xxm = xx - drx
                              yym = yy - dry
                              rrm = sqrt(xxm*xxm+yym*yym)
                              uxp = uxctr + (xxp-xctr)*duxdx + (yyp-yctr)*duxdy + deltaz*duxdz
                              uyp = uyctr + (xxp-xctr)*duydx + (yyp-yctr)*duydy + deltaz*duydz
                              uxm = uxctr + (xxm-xctr)*duxdx + (yym-yctr)*duxdy + deltaz*duxdz
                              uym = uyctr + (xxm-xctr)*duydx + (yym-yctr)*duydy + deltaz*duydz
                              urp = (uxp*xxp+uyp*yyp)/rrp
                              utp = (uyp*xxp-uxp*yyp)/rrp
                              urm = (uxm*xxm+uym*yym)/rrm
                              utm = (uym*xxm-uxm*yym)/rrm
                              durdt = (urp-urm)/(2.0*mdx)
                              dutdt = (utp-utm)/(2.0*mdx)
                              
                              !idx = ridx*tv
                              idx = (kk*rsize+ridx)*tv
!     write (6,*) idx
!     0 - Let's keep track of how many cells contribute to the integral
                              turb(idx) = turb(idx) + vol
!     1 - Density
                              idx = idx + 1
                              turb(idx) = turb(idx) + rho*vol
!     2 - Flow vars not density-weighted
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + ur*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + ut*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + uz*vol
                              idx = idx + 1
                              turb(idx) = turb(idx) + temp*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + YH2*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + YO2*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + YN2*vol
!     9 - Flow vars density-weighted
                              idx = idx + 1
                              turb(idx) = turb(idx) + rho*ur*vol
                              idx = idx + 1
                              turb(idx) = turb(idx) + rho*ut*vol
                              idx = idx + 1
                              turb(idx) = turb(idx) + rho*uz*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + rho*temp*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + rho*YH2*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + rho*YO2*vol
                              !idx = idx + 1
                              !turb(idx) = turb(idx) + rho*YN2*vol
!     15 - Second order (ur)
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*ur*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*ut*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*uz*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*temp*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*YH2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*YO2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ur*YN2*vol
!      22 - second order (ut)                       
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*ut*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*uz*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*temp*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*YH2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*YO2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*ut*YN2*vol
!      28 - second order (uz and scalar)
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*uz*uz*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*uz*temp*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*temp*temp*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*uz*YH2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YH2*YH2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*uz*YO2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YO2*YO2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*uz*YN2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YN2*YN2*vol
!     37 - scalar gradients
!                              turb(idx) = turb(idx) + rho*temp*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YH2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YO2*vol
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*YN2*vol
!     41 - scalar dissipation
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + (dYH2dx*dYH2dx + dYH2dy*dYH2dy + dYH2dz*dYH2dz)*vol                    
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + (dYO2dx*dYO2dx + dYO2dy*dYO2dy + dYO2dz*dYO2dz)*vol                              
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + (dYN2dx*dYN2dx + dYN2dy*dYN2dy + dYN2dz*dYN2dz)*vol

!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*(dYH2dx*dYH2dx + dYH2dy*dYH2dy + dYH2dz*dYH2dz)*vol                    
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*(dYO2dx*dYO2dx + dYO2dy*dYO2dy + dYO2dz*dYO2dz)*vol                              
!                              idx = idx + 1
!                              turb(idx) = turb(idx) + rho*(dYN2dx*dYN2dx + dYN2dy*dYN2dy + dYN2dz*dYN2dz)*vol
                                
                              
                              !if (idx-(kk*rsize+ridx)*tv).ge.tv) then
                              !   write (*,*) idx-(kk*rsize+ridx)*tv
                              !   call bl_error('TURB VARS TOO SMALL')
                              !endif

!     kk
                           end do
!     r .lt. rsize
                        endif
                        
!     ii
                     end do
!     jj
                  end do
                  
!     isect
               endif
               
!     i
            end do
!     j
         end do
!     k
      end do
      
    end subroutine FORT_INTEGRATE
!
!     Higher order extrap
!
    subroutine pushvtog(lo, hi, dlo, dhi, U, U_lo, U_hi, nc) bind(c,name='pushvtog')
      implicit none
      integer, intent(in) :: nc, lo(3),  hi(3), dlo(3), dhi(3)
      integer, intent(in) :: U_lo(3), U_hi(3)
      real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),U_lo(3):U_hi(3),nc)
      
      integer :: n
      real(amrex_real) :: xlo(3), dx(3)
      
      ! Make up something for these that gets what we want
      dx(1:3)  = 1._amrex_real
      xlo(1:3) = 0._amrex_real
      
      do n = 1,nc
         call hoextraptocc(U(:,:,:,n),U_lo(1),U_lo(2),U_lo(3),U_hi(1),U_hi(2),U_hi(3),lo,hi,dx,xlo)
      enddo
      
    end subroutine pushvtog

  end module theta_module
