program diff
! diffusion on a 2d lattice
! a constant gradient is imposed by special boundary conditions
! tracks time evolution of chemical concentration

implicit none

integer,  parameter :: b8 = selected_real_kind(14)

real(b8), parameter :: d  = 0.01_b8  ! diffusion coefficient
real(b8), parameter :: g  = 1.00_b8  ! concentration gradient
real(b8), parameter :: dt = 0.001_b8 ! time-step size

integer :: i, j, k, n, nx, ny, nfinal
integer :: sysSize(2)
real(b8), allocatable :: c(:,:), cDelta(:,:), cTime(:,:,:)

! set system size
! additional lattice sites needed to create gradient
nx = 3
ny = 3
sysSize(1) = nx + 2
sysSize(2) = ny + 2

! number of time-steps to iterate over
nfinal = 1000

! allocate arrays
allocate( c( sysSize(1), sysSize(2)))
allocate( cDelta( sysSize(1), sysSize(2)))
allocate( cTime( nfinal, sysSize(1), sysSize(2)))

! initialize concentration
c(:,:)       = 10.0_b8
cDelta(:,:)  =  0.0_b8
cTime(:,:,:) =  0.0_b8
! initalize boundary conditions
do j = 1, sysSize(2)
    c(sysSize(1),j) = c(1,j) + g * real(sysSize(1)-1)
    write(*,*) c(sysSize(1),j)
enddo

call init_random_seed()

! time evolution of chemical concentration
do n = 1, nfinal
    ! calculate cDelta for each lattice site
    do i = 2, sysSize(1)-1
        do j = 2, sysSize(2)-1
            cDelta(i,j) = getcDelta( i, j, c, sysSize)
        enddo
    enddo
    ! update c for each lattice site
    do i = 2, sysSize(1)-1
        do j = 2, sysSize(2)-1
            c(i,j) = dt * cDelta(i,j) + c(i,j)
            cTime(n,i,j) = c(i,j)
            ! write(*,*) c(i,j)
        enddo
    enddo
enddo

do i = 1, sysSize(1)
    do j = 1, sysSize(2)
        write(*,*) sum(cTime(:,i,j)) / real(nfinal)
    enddo
enddo

write(*,*)
! write(*,*) ' time series'
! do i = 1, nfinal
!     write(*,*) cTime(i,1,1)
! enddo


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
        write(100,*) ec, cv
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
