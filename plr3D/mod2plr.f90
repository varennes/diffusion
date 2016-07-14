module polar

    use utility

contains


! caclulate polarization vector using edge pixels of a cell
subroutine getMWPolar2( c, p, rSim, sigma, xcell)
    implicit none
    integer,  intent(in), dimension(:,:,:) :: sigma
    integer,  intent(in), dimension(:,:)   :: xcell
    integer,  intent(in), dimension(3)     :: rSim
    real(b8), intent(in), dimension(:,:,:) :: c
    real(b8), intent(inout), dimension(3)  :: p
    real(b8), dimension(3) :: q, qtot
    integer,  dimension(3) :: nn
    integer  :: edge, i, nl
    real(b8) :: rx, ry

    qtot = 0.0_b8
    nl = 1
    ! iterate over all pixels in a cell
    do while( xcell(nl,1) /= 0 )
        ! check whether pixel xcell(nl,:) is on an edge
        q(:) = 0.0_b8
        edge = 0
        do i = 1, 6
            call nnGet( i, nn, rSim, xcell(nl,:))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(xcell(nl,1),xcell(nl,2),xcell(nl,3)) /= sigma(nn(1),nn(2),nn(3)) )then
                ! ls is an edge
                q(:) = real(nn(:)) - real(xcell(nl,:)) + q(:)
                edge = edge + 1
            endif
        enddo
        ! normalize q vector and sample signal
        if( dot_product(q,q) /= 0.0 )then
            ! write(*,*) ' q =',q, ' c =',c( xcell(nl,1), xcell(nl,2)), xcell(nl,1), xcell(nl,2)
            q = q / sqrt( dot_product(q,q) )
            q = c( xcell(nl,1), xcell(nl,2), xcell(nl,3)) * q
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


! update polarization vector for Many Wrongs (MW) mechanism
subroutine getMWPolar( p, c, xcell)
    implicit none
    integer,  intent(in),  dimension(:,:) :: xcell
    real(b8), intent(in), dimension(:,:,:)  :: c
    real(b8), intent(inout), dimension(3) :: p
    real(b8), dimension(3) :: com, q, qtot
    real(b8) :: qmag
    integer  :: i, nl

    qtot = 0.0_b8
    nl = 1
    ! iterate over all pixels within one cell
    do while( xcell(nl,1) /= 0 )
        ! calculate pixel vector
        q    = 0.0_b8
        qmag = 0.0_b8
        com  = 0.0_b8
        call calcCellCOM( xcell, com)
        do i = 1, 3
            q(i) = real(xcell(nl,i)) - com(i)
            qmag = q(i)*q(i) + qmag
        enddo
        qmag = sqrt(qmag)
        if( qmag /= 0.0 )then
            q = q / qmag
            q = c( xcell(nl,1), xcell(nl,2), xcell(nl,3)) * q
        endif
        qtot = qtot + q
        nl = nl + 1
    enddo
    if( dot_product(qtot,qtot) /= 0.0 )then
        qtot = qtot / sqrt( dot_product(qtot,qtot) )
    endif
    p = qtot

end subroutine getMWPolar


! calculate center of mass of a single cell
subroutine calcCellCOM( xcell, com)
    ! x = array of all the cell lattice sites
    ! com = center of mass of a single cell
    implicit none
    integer, intent(in),  dimension(:,:) :: xcell
    real(b8),    intent(out), dimension(:) :: com
    integer :: i, j, nl
    com = 0.0_b8
    call occupyCount( nl, xcell(:,:))
    do i = 1, nl
        do j = 1, 3
            com(j) = com(j) + real(xcell(i,j))
        enddo
    enddo
    com = com / real(nl)
end subroutine calcCellCOM


! output total cluster polarization to fort.141
subroutine wrtPlrTotal( nRun, N, p, tstep)
    implicit none
    integer,  intent(in) :: N, nRun, tstep
    real(b8), intent(in), dimension(:,:) :: p
    real(b8) :: px, py, pz, pT
    integer  :: i, j

    px = 0.0_b8
    py = 0.0_b8
    pz = 0.0_b8
    pT = 0.0_b8
    do i = 1, N
        px = px + p(i,1)
        py = py + p(i,2)
        pz = pz + p(i,3)
    enddo
    ! pT = sqrt( px**2 + py**2)
    write(100+nRun,*) px, py, pz, tstep

end subroutine wrtPlrTotal


end module
