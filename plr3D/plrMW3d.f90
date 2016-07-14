program plrMW2D
! diffusion on a 3d lattice
! a constant gradient is imposed by special boundary conditions
! cells are polarized using MW
use utility
use polar

implicit none

real(b8), parameter :: d  = 1.00_b8  ! diffusion coefficient
real(b8), parameter :: g  = 1.00_b8  ! concentration gradient
real(b8), parameter :: dt = 0.01_b8 ! time-step size

integer,  parameter :: nRunTotal = 1 ! total number of instances
integer,  parameter :: ncell = 5    ! total number of cells

integer :: i, j, k, n, n1, n2, nTfinal, nRun
integer :: sysSize(3), r0(3)
integer,  allocatable :: sigma(:,:,:), xCell(:,:,:)
real(b8), allocatable :: c(:,:,:), cDelta(:,:,:)
real(b8), allocatable :: p(:,:)

call init_random_seed()

! set system size
n1 = 6
n2 = 6
! additional lattice sites needed to create gradient
! system is symmetric perpendicular to gradient
sysSize(1) = n1 + 2 ! this is the gradient direction
sysSize(2) = n2
sysSize(3) = n2
! set cell parameters
r0 = [ 2, 2, 2] ! cell dimensions

! number of time-steps to iterate over
nTfinal = 10 * int( (real(r0(1))*(real(ncell)**(0.333)))**2 / (d*dt) )
nTfinal = 10
write(*,*) 'nTfinal =', nTfinal
write(*,*)


! allocate arrays
allocate( c( sysSize(1), sysSize(2), sysSize(3)))
allocate( cDelta( sysSize(1), sysSize(2), sysSize(3)))
allocate( sigma( sysSize(1), sysSize(2), sysSize(3)))
allocate( xCell( ncell, r0(1)*r0(2)*r0(3)*2, 3))
allocate( p( ncell, 3))

! iterate over number of instances
do nRun = 1, nRunTotal

    ! initialize polarization
    p = 0.0_b8
    ! initialize concentration
    c(:,:,:)      = 10.0_b8 - g
    cDelta(:,:,:) =  0.0_b8
    ! initalize gradient
    do i = 2, sysSize(1)
        do j = 1, sysSize(2)
        do k = 1, sysSize(3)
            c(i,j,k) = c(i,j,k) + g * real(i-1)
        enddo
        enddo
    enddo

    ! do i = 1, sysSize(1)
    !     write(*,*) c(i,:,1), i
    ! enddo
    ! write(*,*)

    ! initialize cluster of cells
    xcell = 0
    call itlSigmaRandom( ncell, r0, sysSize, sigma)
    call makeX( ncell, sysSize, sigma, xCell)
    ! write(*,*) 'sigma'
    ! do i = 1, sysSize(1)
    !     write(*,*) ' x =', i,'plane'
    !     do j = 1, sysSize(2)
    !         write(*,*) sigma(i,j,:)
    !     enddo
    !     write(*,*)
    ! enddo
    ! write(*,*)
    ! write(*,*) ' coordinates of cell pixels'
    ! do i = 1, ncell
    !     write(*,*) 'cell',i
    !     do j = 1, r0(1)*r0(2)*r0(3)
    !         write(*,*) xcell(i,j,:)
    !     enddo
    !     write(*,*)
    ! enddo

    ! time evolution of chemical concentration
    do n = 1, nTfinal
        ! calculate cDelta for each lattice site
        do i = 2, sysSize(1)-1
            do j = 1, sysSize(2)
            do k = 1, sysSize(3)
                cDelta(i,j,k) = getcDelta( i, j, k, c, sysSize)
            enddo
            enddo
        enddo
        ! update c for each lattice site
        do i = 2, sysSize(1)-1
            do j = 1, sysSize(2)
            do k = 1, sysSize(3)
                c(i,j,k) = dt * cDelta(i,j,k) + c(i,j,k)
                ! write(*,*) i,j, 'c =', c(i,j)
                ! if ( c(i,j) /= c(i,j) ) then
                !     write(*,*) i,j, 'c =', c(i,j)
                ! end if
            enddo
            enddo
        enddo
        ! write(*,*)
        ! do i = 1, sysSize(1)
        !     write(*,*) c(i,:), i
        ! enddo
        ! write(*,*)

        ! ! concentration averaged over 3rd dimension
        ! do j = 1, sysSize(2)
        !     do i = 1, sysSize(1)
        !         write(110,"(F8.3)", advance="no") sum(c(i,j,1:sysSize(3)))/real(sysSize(3))
        !     enddo
        !     write(110,"(I3)", advance="no") n
        !     write(110,*) ''
        ! enddo
        ! ! concentration averaged over 2nd dimension
        ! do k = 1, sysSize(3)
        !     do i = 1, sysSize(1)
        !         write(120,"(F8.3)", advance="no") sum(c(i,1:sysSize(2),k))/real(sysSize(2))
        !     enddo
        !     write(120,"(I3)", advance="no") n
        !     write(120,*) ''
        ! enddo

        ! update polarization of each cell
        if ( mod( n, 10) == 0 ) then
            do i = 1, ncell
                call getMWPolar2( c, p(i,:), sysSize, sigma, xCell(i,:,:))
                ! call getMWPolar( p(i,:), c, xCell(i,:,:))
            enddo
            ! output total polarization
            call wrtPlrTotal( nRun, ncell, p, n)
        endif
    enddo

    write(*,*) 'instance', nRun, 'complete'
enddo ! ends instances loop

open(unit=10, file="param.dat", action="write")
write(10,*) ncell
write(10,*) nRunTotal
close(10)

contains
    ! calculate cDelta
    real(b8) function getcDelta( i, j, k, c, sysSize)
        implicit none
        integer,  intent(in) :: i, j, k, sysSize(3)
        real(b8), intent(in) :: c(:,:,:)
        real(b8) :: cd, cv, ec
        cd = 0.0_b8
        cv = 0.0_b8
        ! check for periodic boundaries
        if ( i == 1 ) then
            cd = c( i+1, j, k) + c( sysSize(1), j, k) + cd
        elseif ( i == sysSize(1) ) then
            cd = c( 1, j, k) + c( i-1, j, k) + cd
        else
            cd = c( i+1, j, k) + c( i-1, j, k) + cd
        end if
        if ( j == 1 ) then
            cd = c( i, j+1, k) + c( i, sysSize(2), k) + cd
        elseif ( j == sysSize(2) ) then
            cd = c( i, 1, k) + c( i, j-1, k) + cd
        else
            cd = c( i, j+1, k) + c( i, j-1, k) + cd
        endif
        if ( k == 1 ) then
            cd = c( i, j, k+1) + c( i, j, sysSize(3)) + cd
        elseif ( k == sysSize(3) ) then
            cd = c( i, j, 1) + c( i, j, k-1) + cd
        else
            cd = c( i, j, k+1) + c( i, j, k-1) + cd
        end if
        ! calculate variance and sample noise
        cv = d * (cd + (6.0_b8) * c( i, j, k))
        ec = normal( 0.0_b8, sqrt(cv))
        ! write(100,*) ec, cv
        ! caluclate cDelta
        cd = d * (cd - (6.0_b8) * c( i, j, k))

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
