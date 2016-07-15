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


! update polarization vector for Many Wrongs (MW) mechanism
subroutine getMWPolar( p, c, xcell)
    implicit none
    integer,  intent(in),  dimension(:,:) :: xcell
    real(b8), intent(in), dimension(:,:)  :: c
    real(b8), intent(inout), dimension(2) :: p
    real(b8), dimension(2) :: com, q, qtot
    integer :: i, j, nl

    qtot(:) = 0.0_b8
    call occupyCount( nl, xcell(:,:) )
    ! iterate over all pixels within one cell
    do i = 1, nl
        q = 0.0_b8
        ! calculate pixel vector
        com = 0.0_b8
        call calcCellCOM( xcell, com)
        q(1) = real(xcell(i,1)) - com(1)
        q(2) = real(xcell(i,2)) - com(2)
        if( dot_product(q,q) /= 0.0 )then
            q = q / sqrt( dot_product(q,q) )
            q = c( xcell(nl,1), xcell(nl,2)) * q
        endif
        qtot = qtot + q
    enddo
    if( dot_product(qtot,qtot) /= 0.0 )then
        qtot = qtot / sqrt( dot_product(qtot,qtot) )
    endif
    p = qtot

end subroutine getMWPolar


! EC mechanism for assigning cell polarization vectors
! repulsion vectors are weighted by difference of local concentration to that at
! the center of the cluster. Polarization vectors are adaptive.
subroutine getECPolar( iCell, N, c, g, p, rSim, sigma, x)
    implicit none
    integer,  intent(in) :: iCell, N
    real(b8), intent(in) :: g
    real(b8), intent(in), dimension(:,:)   :: c
    real(b8), intent(inout), dimension(2)  :: p
    integer,  intent(in), dimension(2)     :: rSim
    integer,  intent(in), dimension(:,:)   :: sigma
    integer,  intent(in), dimension(:,:,:) :: x
    real(b8), dimension(N,2) :: com
    real(b8), dimension(2) :: clstCOM, q, qtmp
    integer,  dimension(N) :: nnL
    integer  :: i, j, k, nl, nlSum, P1, P2
    real(b8) :: msCOM, s

    ! calculate COM of cluster
    nlSum   = 0
    clstCOM = 0.0_b8
    do i = 1, N
        nl = 1
        do while ( x(i,nl,1) /= 0 )
            do k = 1, 2
                clstCOM(k) = clstCOM(k) + real(x(i,nl,k))
            enddo
            nl = nl + 1
        enddo
        nlSum = nlSum + nl
    enddo
    clstCOM = clstCOM / real(nlSum)
    ! get mean signal at cluster COM
    msCOM  = c(1,1) + g * (clstCOM(1) - 1.0_b8)
    write(*,*) 'clstCOM =',clstCOM,' msCOM =',msCOM

    ! calulcate com of each cell
    do i = 1, N
        call calcCellCOM( x(i,:,:), com(i,:))
    enddo

    ! calculate cell perimeter and enclusure
    call getContactL( iCell, N, nnL, rSim, sigma, x(iCell,:,:))
    P1 = int(perimCalc(rSim, sigma, x(iCell,:,:)))
    P2 = sum(nnL)
    if( P1 /= P2 )then
        ! calculate repulsion vector
        do j = 1, N
            if( nnL(j) /= 0 )then
                qtmp = com(i,:) - com(j,:)
                qtmp = qtmp / sqrt( dot_product( qtmp, qtmp))
                q(:) = q(:) + real(nnL(j)) * qtmp
            endif
        enddo
        if( dot_product( q, q) /= 0.0 )then
            q = q / sqrt( dot_product( q, q))
        endif
        ! get signal value at cell COM
        s = c( nint(com(iCell,1)), nint(com(iCell,2)))
        q = (s - msCOM) * q / msCOM
    endif
    ! update polarization vector
    p = q
end subroutine getECPolar


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
        do j = 1, 2
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
