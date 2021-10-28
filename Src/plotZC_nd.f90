module stream_module

  use amrex_fort_module, only : amrex_real, dim=>amrex_spacedim
  
  implicit none

  public

contains

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

  
  
  subroutine pushvtog2d(lo, hi, dlo, dhi, U, U_lo, U_hi, nc) bind(c,name='pushvtog2d')
    implicit none
    integer, intent(in) :: nc, lo(2),  hi(2), dlo(2), dhi(2)
    integer, intent(in) :: U_lo(2), U_hi(2)
    real(amrex_real), intent(inout) :: U(U_lo(1):U_hi(1),U_lo(2):U_hi(2),nc)

    integer :: n
    real(amrex_real) :: xlo(2), dx(2)

    ! Make up something for these that gets what we want
    dx(1:2)  = 1._amrex_real
    xlo(1:2) = 0._amrex_real

    do n = 1,nc
       call hoextraptocc(U(:,:,n),U_lo(1),U_lo(2),U_hi(1),U_hi(2),lo,hi,dx,xlo)
    enddo
      
  end subroutine pushvtog2d

 
end module stream_module
