module utils

    public::swap                 ! Swap values of an array
    public::iswap                ! Swap values of an integer array
    public::quicksort            ! Quick sorting recursively
    public::iquicksort           ! Quick sorting recursively (integers)
    public::combination          ! Calculate combination
    public::factorial            ! Calculate factorial
    public::triangular_rec       ! Calculate triangular number recursively
    public::triangular_n         ! Calculate triangular number
    public::position_triang      ! Calculate index of array form indexes of a matrix
    public::subdsyevv            ! Diagonalize matrix
    public::remove_ext           ! Remove extension from a file name

contains

    !***********************
    !*  Look for extension in file name and remove it
    !*  -in
    !*  input     : Input name
    !*  -out
    !*  trimmed   : Input without extension
    !***********************
    function remove_ext(input) result(trimmed)
        implicit None
        character(*),intent(in) :: input
        character(100)          :: trimmed
        integer                 :: index

        index=0
        trimmed = trim(input)
        index = scan(trimmed,'.') ! Look for the index of a dot
        if (index>1) then ! If it has an extension
            trimmed = trimmed(1:index-1)
        end if
    end function remove_ext

    !***********************
    !*  This subroutine is to do the diagonalization of a symmetric matrix
    !*  From Stefano Evangelisti
    !*  mat : Input matrix
    !*  d   : Dimension of the matrix
    !*  W   : Store the eigenvalues
    !*  VR  : Store the eigenvectors
    !***********************
    subroutine subdsyevv(mat,W,d,VR)
        implicit none

        integer ::   k, l,LW,INFO,d
        real(kind=8), dimension(d,d) :: mat
        real(kind=8), dimension(d) :: W,VR(d,d)
        real(kind=8), allocatable :: ap(:,:)
        real(kind=8), allocatable :: WR(:)

        allocate (ap(d,d))
        allocate (WR(1))
        ap(:,:)=mat(:,:)
        LW=-1
    
        ! Give the value in the upper part of the mat to AP matrix
        do k=1,d
            ap(k,k)=mat(k,k)
            do l=k,d
                ap(k,l)=mat(k,l)
            enddo
        enddo
    
        ! Diagonalization of the mat matrix
        call DSYEV('V','U',d,ap,d,W,WR,lw,INFO)
        lw=int(WR(1))
        deallocate (WR)
        allocate (WR(lw))
        call DSYEV('V','U',d,ap,d,W,WR,lw,INFO)
        vr(:,:)=ap(:,:)
        return
    end subroutine subdsyevv

    !***********************
    !*  Subroutine to swap values of an array
    !*  -inout
    !*  array     : Input array with assumed dimensions
    !*  -in
    !*  i,j       : Indexes of the values to swap
    !***********************
    subroutine swap(array,a,b)
        implicit None
        integer,intent(in)                     :: a,b
        real(kind=8),intent(inout)             :: array(*)
        real(kind=8)                           :: temp

        temp = array(a)     ! Store value in temp
        array(a) = array(b) ! Change values
        array(b) = temp     ! Restore value in temp

    end subroutine swap

    !***********************
    !*  Subroutine to swap values of an array for integers 
    !*  -inout
    !*  array     : Input array with assumed dimensions
    !*  -in
    !*  i,j       : Indexes of the values to swap
    !***********************
    subroutine iswap(array,a,b)
        implicit None
        integer,intent(in)                :: a,b
        integer,intent(inout)             :: array(*)
        integer                           :: temp

        temp = array(a)     ! Store value in temp
        array(a) = array(b) ! Change values
        array(b) = temp     ! Restore value in temp

    end subroutine iswap

    !***********************
    !*  Recursive subroutine for quick sorting
    !*  -inout
    !*  arr        : Array to sort
    !*  -in
    !*  first,last : Indexes of where to start and finish the sorting
    !*  order      : Order in which sort the array
    !*      1  : Increasing order
    !*     -1  : Decreasing order
    !***********************
    recursive subroutine quicksort(arr,first,last,order)
        implicit none
        integer,intent(in)                      :: first,last,order
        real(kind=8),intent(inout)              :: arr(*) 
        real(kind=8)                            :: x
        integer                                 :: i,j

        ! Initialize value and indexes
        x = arr((first+last)/2)
        if (order==1) then
            i = first; j = last

            do while (.true.) ! Main loop
                do while (arr(i) < x) ! Look for smaller values
                    i=i+1
                end do
                do while (x < arr(j)) ! Look for bigger values
                    j=j-1
                end do
                if (i >= j) exit ! If ordered, exit
                call swap(arr,i,j) ! Swap values
                i=i+1; j=j-1 ! Change indexes
            end do
            ! Call quicksort again
            if (first < i-1) call quicksort(arr, first, i-1, 1)
            if (j+1 < last)  call quicksort(arr, j+1, last, 1)
        else if (order==-1) then
            i = last; j = first
            do while (.true.) ! Main loop
                do while (arr(i) > x) ! Look for smaller values
                    i=i+1
                end do
                do while (x > arr(j)) ! Look for bigger values
                    j=j-1
                end do
                if (i >= j) exit ! If ordered, exit
                call swap(arr,i,j) ! Swap values
                i=i+1; j=j-1 ! Change indexes
            end do
            ! Call quicksort again
            if (last < i-1) call quicksort(arr, i-1, last, -1)
            if (j+1 < first)  call quicksort(arr, first, j+1, -1)
        end if
    end subroutine quicksort

    !***********************
    !*  Recursive subroutine for quick sorting for integers
    !*  -inout
    !*  arr        : Array to sort
    !*  -in
    !*  first,last : Indexes of where to start and finish the sorting
    !*  order      : Order in which sort the array
    !*      1  : Increasing order
    !*     -1  : Decreasing order
    !***********************
    recursive subroutine iquicksort(arr,first,last, order)
        implicit none
        integer,intent(in)                      :: first,last,order
        integer,intent(inout)                   :: arr(*) 
        real(kind=8)                            :: x
        integer                                 :: i,j

        ! Initialize value and indexes
        x = arr((first+last)/2)
        
        if (order==1) then
            i = first; j = last
            do while (.true.) ! Main loop
                do while (arr(i) < x) ! Look for smaller values
                    i=i+1
                end do
                do while (x < arr(j)) ! Look for bigger values
                    j=j-1
                end do
                if (i >= j) exit ! If ordered, exit
                call iswap(arr,i,j) ! Swap values
                i=i+1; j=j-1 ! Change indexes
            end do
            ! Call quicksort again
            if (first < i-1) call iquicksort(arr, first, i-1, 1)
            if (j+1 < last)  call iquicksort(arr, j+1, last, 1)
        else if (order==-1) then
            i = last; j = first
            do while (.true.) ! Main loop
                do while (arr(i) > x) ! Look for smaller values
                    i=i+1
                end do
                do while (x > arr(j)) ! Look for bigger values
                    j=j-1
                end do
                if (i >= j) exit ! If ordered, exit
                call iswap(arr,i,j) ! Swap values
                i=i+1; j=j-1 ! Change indexes
            end do
            ! Call quicksort again
            if (last < i-1) call iquicksort(arr, i-1, last, -1)
            if (j+1 < first)  call iquicksort(arr, first, j+1, -1)
        end if
    end subroutine iquicksort

    !***********************
    !*  Calculate the combination of N over m
    !*  Using integer(kind=8) for factorials over 32 bits
    !*  -in
    !*  N,m : Integers
    !*  -out
    !*  f   : Result
    !***********************
    function combination(N,m) result(f)
        implicit None

        integer,intent(in) :: N,m
        integer(kind=8)    :: N_64,m_64
        integer            :: f
        
        N_64 = int8(N); m_64 = int8(m)
        f=int((factorial(N_64))/(factorial(m_64)*factorial(N_64-m_64)))
        
    end function combination

    !***********************
    !*  Calculate the factorial of n (n!), 64 bits version
    !*  -in
    !*  n : Integers
    !*  -out
    !*  f : Result
    !***********************
    recursive function factorial(n)  result(f)
        integer(kind=8), intent(in) :: n
        integer(kind=8)             :: f
        if(n == 0) then
            f=1
        else
            f=n*factorial(n-1)
        end if
    end function factorial

    !***********************
    !*  Calculate the sum of the n-1 previous numbers
    !*  Example: 
    !*      triangular_rec(5) = 4+3+2+1 = 10
    !*      triangular_rec(9) = 9+8+7+6+5+4+3+2+1 = 45
    !*  Deprecated, not efficient enough, see triangular_n
    !*  Difference in efficiency
    !*       n = 125000
    !*              Normal     Recursive
    !*  Result    7812437500   7812437500
    !*  Time (s)   0.000002     0.008050
    !*
    !*  -in
    !*  n : Integer
    !*  -out
    !*  f : Result
    !***********************
    recursive function triangular_rec(n)  result(f)
        integer, intent(in) :: n
        integer             :: f
        if(n <= 1) then
            f=0
        else
            f=(n-1)+triangular_rec(n-1)
        end if
    end function triangular_rec

    !***********************
    !*  Calculate the sum of the n-1 previous numbers
    !*  Example: 
    !*      triangular_m(5) = 4+3+2+1 = 5*4/2 = 10
    !*      triangular_m(10) = 9+8+7+6+5+4+3+2+1 = 10*9/2 = 45
    !*  -in
    !*  n : Integer
    !*  -out
    !*  f : Result
    !***********************
    function triangular_n(n)  result(f)
        integer, intent(in) :: n
        integer             :: f

        f = n*(n-1)/2

    end function triangular_n

    !***********************
    !*  Calculate the index of an array corresponding
    !*  to the index of a lower triangular matrix
    !*  -in
    !*  n,m : Indexes of the matrix
    !*  -out
    !*  f   : Index of the array
    !***********************
    function position_triang(n,m) result(f)
        integer,intent(in) :: n,m
        integer            :: f
        if (n>m)then
            f=triangular_n(n-1)+m
        else if (m>n)then
            f=triangular_n(m-1)+n
        end if
    end function position_triang   

    !***********************
    !*  Print coupling array as a triangular matrix
    !*  -in
    !*  unit      : Unit to write
    !*  couplings : Coupling array
    !*  n_sites   : Number of sites
    !***********************
    subroutine print_couplings(unit,couplings,n_sites)
        implicit none

        integer,intent(in)                    :: unit,n_sites
        real(kind=8), allocatable, intent(in) :: couplings(:)
        
        integer                               :: i,pos,k
        character(50)                         :: frmt
        
        ! Print header
        write(unit,'(a6)',advance="no")''
        do i=1,n_sites-1
            write(unit,'(2x,i7)',advance="no")i
        end do
        write(unit,*)
        
        ! Print array as a triangular matrix
        pos=1
        do i=2,n_sites
            write(frmt,'(a,i4,a)')'(x,i4,x,',i-1,'(2x,f7.3))'
            write(unit,frmt)i,(couplings(k),k=pos,pos+i-2)
            pos=pos+i-1
        end do

    end subroutine print_couplings  

    !***********************
    !*  Read input
    !*  -in
    !*  input     : Input name
    !*  -out
    !*  verb      : Verbose (debug) output
    !*  n_sites   : Number of sites
    !*  n_coupl   : Number of couplings
    !*  couplings : Coupling array
    !***********************
    subroutine read_input(input,verb,n_sites,n_coupl,couplings)
        implicit none

        character(:), allocatable, intent(in)  :: input
        
        integer,intent(out)                    :: n_sites, n_coupl
        logical,intent(out)                    :: verb
        real(kind=8), allocatable, intent(out) :: couplings(:)
            
        integer                                :: iost,i,j,size,pos
        real(kind=8)                           :: coupl
        character(1)                           :: dummy_bool
        
        open(77, file=input)

        ! Read debug mode boolean
        read(77,*) dummy_bool
        if (dummy_bool=='T')then
            verb=.true.
        else
            verb=.false.
        end if

        !Read number of sites and couplings
        read(77,*) n_sites
        read(77,*) n_coupl

        ! Get size of coupling array
        size=triangular_n(n_sites)
        
        ! Allocate and initialize couplings array
        allocate(couplings(size))
        couplings = 0.d0

        do while (.true.)
            ! Read coupling
            read(77,*,iostat=iost) i,j,coupl
            ! If end of file, exit while
            if (iost==-1)then
                exit
            else
                ! Get position in the array
                pos=position_triang(i,j)
                ! Store coupling value in array
                couplings(pos)=coupl
            end if
        end do
    end subroutine read_input

    !***********************
    !*  Parse arguments
    !*  -out
    !*  input :  Input name
    !*  output:  Output name
    !***********************
    subroutine parse_args(input,output)
        implicit none

        character(:), allocatable, intent(out) :: input, output
        character(100)                         :: trimmed = 'output'
        integer                                :: inp,i

        ! If no arguments, ask for input
        if (command_argument_count()<1) then
            write(*,*) 'No input specified'
            write(*,*) 'Please, specify an input'
            allocate(character(100) :: input)
            read(*,*) input
        ! Else, get input name
        else
            call get_command_argument(1, length=inp)
            allocate(character(inp) :: input)
            call get_command_argument(1, input)
        end if

        ! If more than one argument, get output name from second argument
        if (command_argument_count()>1) then
            call get_command_argument(2,length=i)
            allocate(character(i) :: output)
            call get_command_argument(2,output)
            if (output=='-o') then ! If flag
                deallocate(output)
                if (command_argument_count()>2) then ! If output name provided, read
                    call get_command_argument(3,length=i)
                    allocate(character(i) :: output)
                    call get_command_argument(3,output)
                else ! If not output name provided, generate one from input
                    allocate(character(inp+8) :: output) 
                    trimmed = remove_ext(input)
                    output = trim(trimmed)//'.out'
                end if
            end if
        ! Else, no output
        else
            allocate(character(4) :: output)
            output='None'
        end if
    end subroutine parse_args

end module utils