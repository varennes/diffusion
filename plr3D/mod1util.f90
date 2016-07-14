module utility

    ! b8 will be used to define reals with 14 digits
    integer, parameter:: b8 = selected_real_kind(14)

contains


! create an array of all the cells lattice sites
subroutine makeX( N, rSim, sigma, x)
    ! L = number of lattice sites along one dimension
    ! N = total number of cells
    ! sigma = array of cell labels
    ! x = array of all the cells lattice sites
    ! lsCount = array that keeps track the number of lattice sites
    implicit none
    integer, intent(in) :: N
    integer, intent(in),  dimension(3)     :: rSim
    integer, intent(in),  dimension(:,:,:)   :: sigma
    integer, intent(out), dimension(:,:,:) :: x
    integer :: cellIndex, i, nx, ny, nz
    integer, allocatable :: lsCount(:)

    allocate( lsCount(N) )
    lsCount(:) = 0
    x = 0

    do nx = 1, rSim(1)
        do ny = 1, rSim(2)
            do nz = 1, rsim(3)
                if( sigma(nx,ny,nz) /= 0 )then
                    cellIndex = sigma(nx,ny,nz)
                    lsCount(cellIndex) = lsCount(cellIndex) + 1

                    i = lsCount(cellIndex)
                    x(cellIndex,i,:) = [nx,ny,nz]
                endif
            enddo
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
    integer, intent(in),    dimension(3)   :: r0, rSim
    integer, intent(inout), dimension(:,:,:) :: sigma
    integer, allocatable :: cellGrid(:,:), ocpyGrid(:,:,:)
    integer, dimension(6,3) :: dr
    integer, dimension(3) :: center, cell1, cellk
    integer :: count, i, j, k, l
    real(b8) :: r

    allocate( cellGrid(N,3) )
    allocate( ocpyGrid(N+1,N+1,N+1) )

    sigma    = 0
    cellGrid = 0
    ocpyGrid = 0
    ocpyGrid(N/2+1,N/2+1,N/2+1) = 1
    cellGrid(1,:) = [ N/2+1, N/2+1, N/2+1]

    dr(1,:) = [ 1, 0, 0]
    dr(2,:) = [ 0, 1, 0]
    dr(3,:) = [-1, 0, 0]
    dr(4,:) = [ 0,-1, 0]
    dr(5,:) = [ 0, 0, 1]
    dr(6,:) = [ 0, 0,-1]

    center(1) = (rSim(1)-2)/2 + 2
    center(2) = rSim(2)/2 + 1
    center(3) = rSim(3)/2 + 1
    cell1(1)  = center(1) - r0(1)/2
    cell1(2)  = center(2) - r0(2)/2
    cell1(3)  = center(3) - r0(3)/2

    sigma( cell1(1):cell1(1)+r0(1)-1, cell1(2):cell1(2)+r0(2)-1, cell1(3):cell1(3)+r0(3)-1) = 1
    ! write(*,*) 'cellGrid(1) =', cellGrid(1,:)
    ! write(*,*) 'x',cell1(1),':',cell1(1)+r0(1)-1
    ! write(*,*) 'y',cell1(2),':',cell1(2)+r0(2)-1
    ! write(*,*) 'z',cell1(3),':',cell1(3)+r0(3)-1

    i = 1
    k = 1
    do i = 1, N-1
        if( k >= N )then
            exit
        endif

        call random_number(r)
        j = 1 + floor(r*6.0)

        do count = 1, 6
            if( k >= N )then
                exit
            endif
            if( ocpyGrid( cellGrid(i,1)+dr(j,1), cellGrid(i,2)+dr(j,2), cellGrid(i,3)+dr(j,3)) == 0 )then
                k = k + 1
                cellGrid(k,1) = cellGrid(i,1)+dr(j,1)
                cellGrid(k,2) = cellGrid(i,2)+dr(j,2)
                cellGrid(k,3) = cellGrid(i,3)+dr(j,3)
                ! write(*,*) 'dr =',dr(j,:)
                ! write(*,*) 'k =', k,'cellGrid =',cellGrid(k,:)

                ocpyGrid( cellGrid(k,1), cellGrid(k,2), cellGrid(k,3)) = k
            endif

            j = j + 1
            if( j > 6 )then
                j = 1
            endif
        enddo
    enddo

    do k = 2, N
        i = cellGrid(k,1) - cellGrid(1,1) ! grid space distance from cell 1
        j = cellGrid(k,2) - cellGrid(1,2)
        l = cellGrid(k,3) - cellGrid(1,3)

        cellk(1) = cell1(1) + i*r0(1)
        cellk(2) = cell1(2) + j*r0(2)
        cellk(3) = cell1(3) + l*r0(3)

        ! write(*,*) 'grid dist',i,j,l,' cellk',cellk

        if ( cellk(1) < 1 .OR. (cellk(1)+r0(1)-1) > rSim(1) ) then
            write(*,*) 'Cells out of bounds! You done goofed!'
            write(*,*) '1', cellk(1),cellk(1)+r0(1)-1, 'rSim(1)=',rSim(1)
            exit
        elseif ( cellk(2) < 1 .OR. (cellk(2)+r0(2)-1) > rSim(2) ) then
            write(*,*) 'Cells out of bounds! You done goofed!'
            write(*,*) '2', cellk(2),cellk(2)+r0(2)-1, 'rSim(2)=',rSim(2)
            exit
        elseif ( cellk(3) < 1 .OR. (cellk(3)+r0(3)-1) > rSim(3) ) then
            write(*,*) 'Cells out of bounds! You done goofed!'
            write(*,*) '3', cellk(3),cellk(3)+r0(3)-1, 'rSim(3)=',rSim(3)
            exit
        endif

        sigma( cellk(1):cellk(1)+r0(1)-1, cellk(2):cellk(2)+r0(2)-1, cellk(3):cellk(3)+r0(3)-1) = k
    enddo

end subroutine itlSigmaRandom


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
    ! i = 5 -- nn front
    ! i = 6 -- nn back
    ! nn = array of nearest neighbor coordinates
    ! rSim = dimensions of simulation space
    ! x = coordinates of point
    implicit none
    integer, intent(in) :: i
    integer, dimension(3), intent(in) :: rSim, x
    integer, dimension(3), intent(out) :: nn

    if( i == 1 )then
        nn(1) = x(1) + 1
        nn(2) = x(2)
        nn(3) = x(3)
    elseif( i == 2 )then
        nn(1) = x(1)
        nn(2) = x(2) + 1
        nn(3) = x(3)
    elseif( i == 3 )then
        nn(1) = x(1) - 1
        nn(2) = x(2)
        nn(3) = x(3)
    elseif( i == 4 )then
        nn(1) = x(1)
        nn(2) = x(2) - 1
        nn(3) = x(3)
    elseif( i == 5 )then
        nn(1) = x(1)
        nn(2) = x(2)
        nn(3) = x(3) + 1
    elseif( i == 6 )then
        nn(1) = x(1)
        nn(2) = x(2)
        nn(3) = x(3) - 1
    endif

    if( nn(1) > rSim(1) .OR. nn(1) < 1 )then
        nn = [ 0, 0, 0]
    elseif( nn(2) > rSim(2) .OR. nn(2) < 1 )then
        nn = [ 0, 0, 0]
    elseif( nn(3) > rSim(3) .OR. nn(3) < 1 )then
        nn = [ 0, 0, 0]
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


end module
