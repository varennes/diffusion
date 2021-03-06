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


! EC mechanism for assigning cell polarization vectors
! repulsion vectors are weighted by difference of local concentration to that at
! the center of the cluster. Polarization vectors are adaptive.
! local concentration is average over all pixels within the cell
subroutine getECPolar2( iCell, N, c, g, p, rSim, sigma, x)
    implicit none
    integer,  intent(in) :: iCell, N
    real(b8), intent(in) :: g
    real(b8), intent(in), dimension(:,:,:)   :: c
    real(b8), intent(inout), dimension(3)  :: p
    integer,  intent(in), dimension(3)     :: rSim
    integer,  intent(in), dimension(:,:,:) :: sigma, x
    real(b8), dimension(N,3) :: com
    real(b8), dimension(3) :: clstCOM, q, qtmp
    integer,  dimension(N) :: nnL
    integer  :: i, j, k, nl, nlSum, P1, P2
    real(b8) :: msCOM, s

    ! calculate COM of cluster
    nlSum   = 0
    clstCOM = 0.0_b8
    do i = 1, N
        nl = 1
        do while ( x(i,nl,1) /= 0 )
            do k = 1, 3
                clstCOM(k) = clstCOM(k) + real(x(i,nl,k))
            enddo
            nl = nl + 1
        enddo
        nlSum = nlSum + (nl-1)
    enddo
    clstCOM = clstCOM / real(nlSum)
    ! get mean signal at cluster COM
    msCOM  = c(1,1,1) + g * (clstCOM(1) - 1.0_b8)
    ! write(*,*) 'clstCOM =',clstCOM,' msCOM =',msCOM

    ! calulcate com of each cell
    do i = 1, N
        call calcCellCOM( x(i,:,:), com(i,:))
    enddo

    ! calculate cell perimeter and enclusure
    call getContactL( iCell, N, nnL, rSim, sigma, x(iCell,:,:))
    q  = 0.0_b8
    P1 = int(perimCalc(rSim, sigma, x(iCell,:,:)))
    P2 = sum(nnL)
    if( P1 /= P2 )then
        ! calculate repulsion vector
        do j = 1, N
            if( nnL(j) /= 0 )then
                qtmp = com(j,:) - com(i,:)
                qtmp = qtmp / sqrt( dot_product( qtmp, qtmp))
                q(:) = q(:) + real(nnL(j)) * qtmp
            endif
        enddo
        if( dot_product( q, q) /= 0.0 )then
            q = q / sqrt( dot_product( q, q))
        endif
        ! get signal value by averaging over all of the cells' pixels
        s  = 0.0_b8
        nl = 1
        do while ( x(iCell,nl,1) /= 0 )
            s  =  s + c( x(iCell,i,1), x(iCell,i,2), x(iCell,i,3))
            nl = nl + 1
        enddo
        s = s / real(nl-1)
        q = (s - msCOM) * q / msCOM
        ! write(*,*) 'com', com(iCell,:)
        ! write(*,*) 's=',s,'msCOM',msCOM,'clstCOM',clstCOM,  'q =',q
    endif
    ! update polarization vector
    p = q
end subroutine getECPolar2


! EC mechanism for assigning cell polarization vectors
! repulsion vectors are weighted by difference of local concentration to that at
! the center of the cluster. Polarization vectors are adaptive.
subroutine getECPolar( iCell, N, c, g, p, rSim, sigma, x)
    implicit none
    integer,  intent(in) :: iCell, N
    real(b8), intent(in) :: g
    real(b8), intent(in), dimension(:,:,:)   :: c
    real(b8), intent(inout), dimension(3)  :: p
    integer,  intent(in), dimension(3)     :: rSim
    integer,  intent(in), dimension(:,:,:) :: sigma, x
    real(b8), dimension(N,3) :: com
    real(b8), dimension(3) :: clstCOM, q, qtmp
    integer,  dimension(N) :: nnL
    integer  :: i, j, k, nl, nlSum, P1, P2
    real(b8) :: msCOM, s

    ! calculate COM of cluster
    nlSum   = 0
    clstCOM = 0.0_b8
    do i = 1, N
        nl = 1
        do while ( x(i,nl,1) /= 0 )
            do k = 1, 3
                clstCOM(k) = clstCOM(k) + real(x(i,nl,k))
            enddo
            nl = nl + 1
        enddo
        nlSum = nlSum + (nl-1)
    enddo
    clstCOM = clstCOM / real(nlSum)
    ! get mean signal at cluster COM
    msCOM  = c(1,1,1) + g * (clstCOM(1) - 1.0_b8)
    ! write(*,*) 'clstCOM =',clstCOM,' msCOM =',msCOM

    ! calulcate com of each cell
    do i = 1, N
        call calcCellCOM( x(i,:,:), com(i,:))
    enddo

    ! calculate cell perimeter and enclusure
    call getContactL( iCell, N, nnL, rSim, sigma, x(iCell,:,:))
    q  = 0.0_b8
    P1 = int(perimCalc(rSim, sigma, x(iCell,:,:)))
    P2 = sum(nnL)
    if( P1 /= P2 )then
        ! calculate repulsion vector
        do j = 1, N
            if( nnL(j) /= 0 )then
                qtmp = com(j,:) - com(i,:)
                qtmp = qtmp / sqrt( dot_product( qtmp, qtmp))
                q(:) = q(:) + real(nnL(j)) * qtmp
            endif
        enddo
        if( dot_product( q, q) /= 0.0 )then
            q = q / sqrt( dot_product( q, q))
        endif
        ! get signal value at cell COM
        s = c( nint(com(iCell,1)), nint(com(iCell,2)), nint(com(iCell,2)))
        q = (s - msCOM) * q / msCOM
        ! write(*,*) 'com', com(iCell,:)
        ! write(*,*) 's=',s,'msCOM',msCOM,'clstCOM',clstCOM,  'q =',q
    endif
    ! update polarization vector
    p = q
end subroutine getECPolar


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
    write(100+nRun,"(E16.8)", advance="no") px
    write(100+nRun,"(E17.8)", advance="no") py
    write(100+nRun,"(E17.8)", advance="no") pz
    write(100+nRun,"(I10)", advance="no") tstep
    write(100+nRun,*) ''

end subroutine wrtPlrTotal


end module
