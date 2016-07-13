module polar

    use utility

contains


! caclulate polarization vector using edge pixels of a cell
subroutine getMWPolar2( c, p, rSim, sigma, xcell)
    implicit none
    integer,  intent(in), dimension(:,:)  :: sigma, xcell
    integer,  intent(in), dimension(2)    :: rSim
    real(b8), intent(in), dimension(:,:)  :: c
    real(b8), intent(inout), dimension(2) :: p
    real(b8), dimension(2) :: q, qtot
    integer,  dimension(2) :: nn
    integer :: edge, i, nl
    real(b8) :: ms, s, rx, ry

    qtot(:) = 0.0_b8
    nl = 1
    ! iterate over all pixels in a cell
    do while( xcell(nl,1) /= 0 )
        ! check whether pixel xcell(nl,:) is on an edge
        q(:) = 0.0_b8
        edge = 0
        do i = 1, 4
            call nnGet( i, nn, rSim, xcell(nl,:))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(xcell(nl,1),xcell(nl,2)) /= sigma(nn(1),nn(2)) )then
                ! ls is an edge
                q(:) = real(nn(:)) - real(xcell(nl,:)) + q(:)
                edge = edge + 1
            endif
        enddo
        ! normalize q vector and sample signal
        if( dot_product(q,q) /= 0.0 )then
            ! write(*,*) ' q =',q, ' c =',c( xcell(nl,1), xcell(nl,2)), xcell(nl,1), xcell(nl,2)
            q = q / sqrt( dot_product(q,q) )
            q = c( xcell(nl,1), xcell(nl,2)) * q
            qtot = qtot + q
            ! write(*,*) 'qtot =',qtot
        endif

        nl = nl + 1
    enddo
    ! normalize qtot vector
    if( dot_product(qtot,qtot) /= 0.0 )then
        qtot = qtot / sqrt( dot_product(qtot,qtot) )
    endif
    ! write(*,*) 'qtot =', qtot
    ! update polarization vector
    p = qtot

end subroutine getMWPolar2


! ! update polarization vector for Many Wrongs (MW) mechanism
! subroutine getMWPolar( p, c, xcell)
!     implicit none
!     real(b8), intent(in) :: plrR
!     integer,  intent(in),  dimension(:,:) :: xcell
!     real(b8), intent(in),  dimension(:)   :: cellCOM
!     real(b8), intent(inout), dimension(2) :: p
!     real(b8), dimension(2) :: q, qpixel
!     integer :: i, j, nl
!     real(b8) :: ms, s, rx, ry, qmag
!
!     q(:) = 0.0
!     qmag = 1.0
!     call occupyCount( nl, xcell(:,:) )
!     ! iterate over all pixels within one cell
!     do i = 1, nl
!         qpixel(:) = 0.0
!         ! calculate pixel vector
!         qpixel(1) = rx - cellCOM(1)
!         qpixel(2) = ry - cellCOM(2)
!         if( dot_product(qpixel,qpixel) /= 0.0 )then
!             qpixel = qpixel / sqrt( dot_product(qpixel,qpixel) )
!             qpixel = c( xcell(nl,1), xcell(nl,2)) * qpixel
!         endif
!         q = q + s * qpixel
!     enddo
!     if( dot_product(q,q) /= 0.0 )then
!         q = q / sqrt( dot_product(q,q) )
!     endif
!     p = q
!
! end subroutine getMWPolar


! output total cluster polarization to fort.141
subroutine wrtPlrTotal( nRun, N, p, tstep)
    implicit none
    integer,  intent(in) :: N, nRun, tstep
    real(b8), intent(in), dimension(:,:) :: p
    real(b8) :: px, py, pT
    integer  :: i, j

    px = 0.0_b8
    py = 0.0_b8
    pT = 0.0_b8
    do i = 1, N
        px = px + p(i,1)
        py = py + p(i,2)
    enddo
    ! pT = sqrt( px**2 + py**2)
    write(100+nRun,*) px, py, tstep

end subroutine wrtPlrTotal


end module
