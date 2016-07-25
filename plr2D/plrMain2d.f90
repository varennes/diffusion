program plrMW2D
! diffusion on a 2d lattice
! a constant gradient is imposed by special boundary conditions
! cells are polarized using MW
use utility
use polar

implicit none

real(b8), parameter :: d  = 10.00_b8  ! diffusion coefficient
real(b8), parameter :: g  =  1.00_b8  ! concentration gradient
real(b8), parameter :: dt = 0.01_b8 ! time-step size

integer :: nRunTotal ! total number of instances
integer :: ncell     ! total number of cells

integer :: i, j, k, l, n, nx, ny, nTmod, nTfinal, nRun
integer :: sysSize(2), r0(2)
integer,  allocatable :: sigma(:,:), xCell(:,:,:)
real(b8), allocatable :: c(:,:), cDelta(:,:)
real(b8), allocatable :: p(:,:)

call init_random_seed()
! read on a few parameters
open(unit=10, file="param.dat", action="read")
read(10,*) ncell
read(10,*) nRunTotal
write(*,*) 'ncell =', ncell, 'nRunTotal =', nRunTotal
! set system size
! additional lattice sites needed to create gradient
nx = 120
ny = 120
sysSize(1) = nx + 2 ! this is the gradient direction
sysSize(2) = ny
! set cell parameters
r0 = [ 10, 10] ! cell dimensions

! allocate arrays
allocate( c( sysSize(1), sysSize(2)))
allocate( cDelta( sysSize(1), sysSize(2)))
allocate( sigma( sysSize(1), sysSize(2)))
allocate( xCell( ncell, r0(1)*r0(2)*2, 2))
allocate( p( ncell, 2))

! iterate over number of instances
do nRun = 1, nRunTotal

    ! initialize polarization
    p = 0.0_b8
    ! initialize concentration
    c(:,:)       =  0.0_b8
    cDelta(:,:)  =  0.0_b8
    ! initalize gradient
    do i = 2, sysSize(1)
        do j = 1, sysSize(2)
            c(i,j) = c(i,j) + g * real(i-1)
        enddo
    enddo

    ! do i = 1, sysSize(1)
    !     write(*,*) c(i,:), i
    ! enddo
    ! write(*,*)

    ! initialize cluster of cells
    xcell = 0
    call itlSigmaRandom( ncell, r0, sysSize, sigma)
    call makeX( ncell, sysSize, sigma, xCell)

    ! calculate cluster length and nTfinal
    call getClusterLength( l, sysSize, sigma)
    nTfinal = 10 * int( real(l)**2 / (d*dt) )
    nTmod   = floor( real(nTfinal) / 1000.0)
    if ( nTmod == 0 ) then
        nTmod = floor( real(nTfinal) / 100.0)
        if ( nTmod == 0 ) then
            nTmod = 10
        endif
    endif
    nTfinal = 50
    nTmod   = 5
    write(*,*) 'nTfinal =', nTfinal, 'l =', l
    write(*,*) 'nTmod =', nTmod
    write(*,*)


    ! time evolution of chemical concentration
    do n = 1, nTfinal
        ! calculate cDelta for each lattice site
        do i = 2, sysSize(1)-1
            do j = 1, sysSize(2)
                cDelta(i,j) = dt * getcDelta( i, j, c, sysSize)
            enddo
        enddo
        ! update c for each lattice site
        do i = 2, sysSize(1)-1
            do j = 1, sysSize(2)
                if ( (cDelta(i,j) + c(i,j)) > 0.0 ) then
                    c(i,j) = cDelta(i,j) + c(i,j)
                end if
                ! write(*,*) i,j, 'c =', c(i,j)
                ! if ( c(i,j) /= c(i,j) ) then
                !     write(*,*) i,j, 'c =', c(i,j)
                ! end if
            enddo
        enddo
        ! write(*,*)
        ! do i = 1, sysSize(1)
        !     write(*,*) c(i,:), i
        ! enddo
        ! write(*,*)
        ! do j = 1, sysSize(2)
        !     write(110,*) c(:,j) , n
        ! enddo
        ! update polarization of each cell
        if ( mod( n, nTmod) == 0 ) then
            do i = 1, ncell
                call getECPolar2( i, ncell, c, g, p(i,:), sysSize, sigma, xCell)
                ! call getECPolar( i, ncell, c, g, p(i,:), sysSize, sigma, xCell)
                ! call getMWPolar2( c, p(i,:), sysSize, sigma, xCell(i,:,:))
                ! call getMWPolar( p(i,:), c, xCell(i,:,:))
            enddo
            ! output total polarization
            call wrtPlrTotal( nRun, ncell, p, n)
        endif
    enddo

    write(*,*) 'instance', nRun, 'complete'
enddo ! ends instances loop

! do i = 1, sysSize(1)
!     do j = 1, sysSize(2)
!         write(*,*) sum(cTime(:,i,j)) / real(nfinal)
!     enddo
! enddo
close(10)

contains
    ! calculate cDelta
    real(b8) function getcDelta( i, j, c, sysSize)
        implicit none
        integer,  intent(in) :: i, j, sysSize(2)
        real(b8), intent(in) :: c(:,:)
        real(b8) :: cd, cv, ec
        cd = 0.0_b8
        cv = 0.0_b8
        ! check for periodic boundaries
        if ( i == 1 ) then
            cd = c( i+1, j) + c( sysSize(1), j) + cd
        elseif ( i == sysSize(1) ) then
            cd = c( 1, j) + c( i-1, j) + cd
        else
            cd = c( i+1, j) + c( i-1, j) + cd
        end if
        if ( j == 1 ) then
            cd = c( i, j+1) + c( i, sysSize(2)) + cd
        elseif ( j == sysSize(2) ) then
            cd = c( i, 1) + c( i, j-1) + cd
        else
            cd = c( i, j+1) + c( i, j-1) + cd
        end if
        ! calculate variance and sample noise
        cv = d * (cd + (4.0_b8) * c( i, j))
        ec = normal( 0.0_b8, sqrt(cv))
        ! caluclate cDelta
        cd = d * (cd - (4.0_b8) * c( i, j))

        getcDelta = cd + ec
    end function getcDelta


    ! returns random number between 0 - 1
    function ran1()
        implicit none
        real(b8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    ! returns a random number from a normal distribution
    function normal(mean,sigma)
        implicit none
        real(b8) normal,tmp
        real(b8) mean,sigma
        integer flag
        real(b8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0_b8
            do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
                r1=2.0_b8*ran1()-1.0_b8
                r2=2.0_b8*ran1()-1.0_b8
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0_b8*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal

end program


! initialize random_seed
subroutine init_random_seed()
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
end subroutine init_random_seed
