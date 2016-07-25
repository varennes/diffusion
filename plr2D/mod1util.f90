module utility

    ! b8 will be used to define reals with 14 digits
    integer, parameter:: b8 = selected_real_kind(14)

contains


! find the system size
subroutine getClusterLength( l, sysSize, sigma)
    implicit none
    integer, intent(out) :: l
    integer, intent(in),  dimension(2)     :: sysSize
    integer, intent(in),  dimension(:,:) :: sigma
    integer, dimension(2) :: lmin, lmax
    integer :: i

    lmin(:) = 0
    lmax(:) = 0
    ! find minimum poistion of cluster
    do i = 1, sysSize(1)
        if ( sum( sigma(i,:) ) > 0 ) then
            lmin(1) = i
            exit
        end if
    enddo
    do i = 1, sysSize(2)
        if ( sum( sigma(:,i) ) > 0 ) then
            lmin(2) = i
            exit
        end if
    enddo
    ! find maximum poistion of cluster
    do i = sysSize(1), 1, -1
        if ( sum( sigma(i,:) ) > 0 ) then
            lmax(1) = i
            exit
        end if
    enddo
    do i = sysSize(2), 1, -1
        if ( sum( sigma(:,i) ) > 0 ) then
            lmax(2) = i
            exit
        end if
    enddo
    l = 0
    do i = 1, 2
        if ( (lmax(i)-lmin(i)+1) > l ) then
            l = lmax(i) - lmin(i) + 1
        end if
    enddo
end subroutine getClusterLength


! create an array of all the cells lattice sites
subroutine makeX( N, rSim, sigma, x)
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    ! lsCount = array that keeps track the number of lattice sites
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(2)     :: rSim
    integer, intent(in),  dimension(:,:)   :: sigma
    integer, intent(out), dimension(:,:,:) :: x
    integer :: cellIndex, i, nx, ny
    integer, allocatable :: lsCount(:)

    allocate( lsCount(N) )
    lsCount(:) = 0
    x = 0

    do nx = 1, rSim(1)
        do ny = 1, rSim(2)
            if( sigma(nx,ny) /= 0 )then

                cellIndex = sigma(nx,ny)
                lsCount(cellIndex) = lsCount(cellIndex) + 1

                i = lsCount(cellIndex)
                x(cellIndex,i,:) = [nx,ny]
            endif
        enddo
    enddo

    deallocate( lsCount )
end subroutine makeX


! initialize sigma array such that cells are randomly clumped together
subroutine itlSigmaRandom( N, r0, rSim, sigma)
    ! r0 = initial dimnesions of individual cells
    ! rSim = dimensions of cell initialization space
    ! sigma = array of cell labels
    implicit none
    integer, intent(in) :: N
    integer, intent(in),    dimension(2)   :: r0, rSim
    integer, intent(inout), dimension(:,:) :: sigma
    integer, allocatable :: cellGrid(:,:), ocpyGrid(:,:)
    integer, dimension(4,2) :: dr
    integer, dimension(2) :: center, cell1, cellk
    integer :: count, i, j, k
    real(b8) :: r

    allocate( cellGrid(N,2) )
    allocate( ocpyGrid(N+1,N+1) )

    sigma    = 0
    cellGrid = 0
    ocpyGrid = 0
    ocpyGrid(N/2+1,N/2+1) = 1
    cellGrid(1,:) = [ N/2+1, N/2+1]

    dr(1,:) = [ 1, 0]
    dr(2,:) = [ 0, 1]
    dr(3,:) = [-1, 0]
    dr(4,:) = [ 0,-1]

    center(1) = rSim(1)/2 + 1
    center(2) = rSim(2)/2 + 1
    cell1(1)  = center(1) - r0(1)/2
    cell1(2)  = center(2) - r0(2)/2

    sigma( cell1(1):cell1(1)+r0(1)-1, cell1(2):cell1(2)+r0(2)-1) = 1

    i = 1
    k = 1
    do i = 1, N-1
        if( k >= N )then
            exit
        endif

        call random_number(r)
        j = 1 + floor(r*4)

        do count = 1, 4
            if( k >= N )then
                exit
            endif

            if( ocpyGrid( cellGrid(i,1)+dr(j,1), cellGrid(i,2)+dr(j,2)) == 0 )then
                k = k + 1
                cellGrid(k,1) = cellGrid(i,1)+dr(j,1)
                cellGrid(k,2) = cellGrid(i,2)+dr(j,2)

                ocpyGrid( cellGrid(k,1), cellGrid(k,2)) = k
            endif

            j = j + 1
            if( j > 4 )then
                j = 1
            endif
        enddo
    enddo

    do k = 2, N
        i = cellGrid(1,1) - cellGrid(k,1) ! grid space distance from cell 1
        j = cellGrid(1,1) - cellGrid(k,2)

        cellk(1) = cell1(1) + i*r0(1)
        cellk(2) = cell1(2) + j*r0(2)

        if( cellk(1) < 2 .OR. (cellk(1)+r0(1)-1) > rSim(1)+1 )then
            write(*,*) 'Cells out of bounds! You done goofed!'
            exit
        endif
        if( cellk(2) < 2 .OR. (cellk(2)+r0(1)-1) > rSim(2)+1 )then
            write(*,*) 'Cells out of bounds! You done goofed!'
            exit
        endif

        sigma( cellk(1):cellk(1)+r0(1)-1, cellk(2):cellk(2)+r0(2)-1) = k
    enddo

end subroutine itlSigmaRandom




subroutine itlEdge( edge, ne, rSim, sigma)
    ! edge = array of lattice sites with diff. sigma neighbors
    implicit none
    integer, intent(out) :: ne
    integer, intent(out), dimension(:,:) :: edge
    integer, intent(in),  dimension(2)   :: rSim
    integer, intent(in),  dimension(:,:) :: sigma
    integer, dimension(2) :: nn
    integer :: i1, i2, j, k

    ne = 1
    edge = 0

    do i1 = 1, rSim(1) + 2
    do i2 = 1, rSim(2) + 2
        do k = 1, 4
            call nnGet( k, nn, rSim, [i1,i2])

            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= sigma(i1,i2) )then
                edge(ne,:) = [i1,i2]
                ne = ne + 1
                exit
            endif
        enddo
    enddo
    enddo
    ne = ne - 1

end subroutine itlEdge


! calculate center of mass of the group of cells
subroutine calcXCOM( N, x, xCOMt)
    ! N = total number of cells
    ! x = array of all the cells lattice sites
    ! xCOMt = array of center of mass of all cells at some time
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(:,:,:) :: x
    real(b8),    intent(out), dimension(:)     :: xCOMt
    integer :: i, j, k, nl, nlSum

    nlSum = 0
    xCOMt = 0.0

    do i = 1, N
        call occupyCount( nl, x(i,:,:))

        do j = 1, nl
            do k = 1, 2
                xCOMt(k) = xCOMt(k) + real(x(i,j,k))
            enddo
        enddo

        nlSum = nlSum + nl
    enddo

    xCOMt = xCOMt / real(nlSum)

end subroutine calcXCOM


! Output coordinates of the nearest neighbor (nn)
subroutine nnGet( i, nn, rSim, x)
    ! i indicates the nn we are interested in
    ! i = 1 -- nn up
    ! i = 2 -- nn right
    ! i = 3 -- nn down
    ! i = 4 -- nn left
    ! nn = array of nearest neighbor coordinates
    ! rSim = dimensions of simulation space
    ! x = coordinates of point
    implicit none
    integer, intent(in) :: i
    integer, dimension(2), intent(in) :: rSim, x
    integer, dimension(2), intent(out) :: nn

    if( i == 1 )then
        nn(1) = x(1) + 1
        nn(2) = x(2)

    elseif( i == 2 )then
        nn(1) = x(1)
        nn(2) = x(2) + 1

    elseif( i == 3 )then
        nn(1) = x(1) - 1
        nn(2) = x(2)

    elseif( i == 4 )then
        nn(1) = x(1)
        nn(2) = x(2) - 1

    endif

    if( nn(1) > (rSim(1) + 2) .OR. nn(1) < 1 )then
        nn = [ 0, 0]
    elseif( nn(2) > (rSim(2) + 2) .OR. nn(2) < 1 )then
        nn = [ 0, 0]
    endif

end subroutine nnGet


! Count the number of occupied lattice points
subroutine occupyCount( nl, xcell )
    ! xcell =  array of lattice sites occupied by one cell x(i,:,:)
    ! nl = number of lattice sites contained in xcell
    implicit none
    integer, dimension(:,:), intent(in) :: xcell
    integer, intent(out) :: nl

    ! count how many lattice sites one cell contains
    nl = 1
    do while( xcell(nl,1) /= 0 )
        nl = nl + 1
    enddo
    nl = nl - 1
end subroutine occupyCount


! calculate cell k's contanct lengths with all its neighbors
subroutine getContactL( k, N, nnLk, rSim, sigma, xcell)
    implicit none
    integer, intent(in) :: N
    integer, intent(out), dimension(:)   :: nnLk
    integer, intent(in),  dimension(:,:) :: sigma
    integer, intent(in),  dimension(:,:) :: xcell
    integer, dimension(2) :: nn, rSim
    integer :: i, inn, j, k, nl

    nnLk = 0

    call occupyCount( nl, xcell)

    do i = 1, nl
        do inn = 1, 4
            call nnGet( inn, nn, rSim, xcell(i,1:2))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= k .AND. sigma(nn(1),nn(2)) /= 0 )then
                nnLk( sigma(nn(1),nn(2))) = nnLk( sigma(nn(1),nn(2))) + 1
            endif
        enddo
    enddo

end subroutine getContactL


! calculate perimeter of a cell
real(b8) function perimCalc(rSim, sigma, xcell)
    implicit none
    integer,  intent(in), dimension(:,:) :: sigma
    integer,  intent(in), dimension(:,:) :: xcell
    integer,  intent(in), dimension(2) :: rSim
    integer, dimension(2) :: nn
    integer :: i, inn, k, nl
    real(b8) :: P

    P = 0.0
    k = sigma( xcell(1,1), xcell(1,2)) ! cell label
    call occupyCount( nl, xcell )
    do i = 1, nl
        do inn = 1, 4
            ! write(*,*) xcell(i,1:2)
            call nnGet( inn, nn, rSim, xcell(i,1:2))
            if( nn(1) == 0 )then
                cycle
            endif
            if( sigma(nn(1),nn(2)) /= k )then
                P = P + 1.0
            endif
        enddo
    enddo
    perimCalc = P
end function perimCalc


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


end module
