module main_sub
    use utils
    use print_f

    public::iterate            ! Main iterate loop
    public::assign_ms          ! Assign spin
    public::two_minimum        ! Find minimum eigenvalues
    public::create_ham         ! Create Hamiltonian
    public::diagonal           ! Calculate diagonal eigenvalues
    public::compare            ! Compare determinants
    public::deter              ! Generate determinants
    public::Integ_to_Bool      ! Transform integer determinants to boolean arrays
    
contains

    !***********************
    !*  Subroutine for iterating over all ms spaces and get eigenvalues
    !*  -in
    !*  unit      : Unit to write output
    !*  verb      : Boolean flag to write debug info
    !*  n_sites   : Number of sites
    !*  couplings : Couplings array
    !*  -out
    !*  f_eigenv  : Array containing the ms and the eigenvalues for each subspace
    !*  ms_spaces : Number of subspaces
    !*  arr_dets  : Array containing the number of determinants for each subspace
    !***********************
    subroutine iterate(unit,verb,n_sites,couplings,f_eigenv,ms_spaces,arr_dets)
        implicit None
        integer,intent(in)                   :: unit,n_sites
        real(kind=8), allocatable,intent(in) :: couplings(:)
        logical,intent(in)                   :: verb
        integer,intent(out)                  :: ms_spaces
        integer,allocatable,intent(out)      :: arr_dets(:)
        real(kind=8),allocatable,intent(out) :: f_eigenv(:,:)

        integer                              :: i,j,k,accu,compt,alpha,beta
        integer                              :: min_number_alpha,ndeter,max_deter_total
        integer,allocatable                  :: tab(:)
        real(kind=8)                         :: ms
        real(kind=8),allocatable             :: ham(:,:), evect(:,:), evals(:)
        logical*1,allocatable                :: det(:,:)
        character(50)                        :: frmt

       
        ms_spaces = n_sites/2+1 ! Number of subspaces
        min_number_alpha=n_sites-ms_spaces+1 ! Min number of alpha spins
        max_deter_total=combination(n_sites,min_number_alpha) ! Max number of determinants

        call writesep(unit,0,0) ! Separator

        ! Print basic info
        write(unit,'(x,a,/)')'Starting calculation'
        write(unit,104)'Number of sites =',n_sites
        write(unit,104)'Number of ms spaces =',ms_spaces
        write(unit,104)'Max number of determinants =',max_deter_total
        write(unit,104)'Min number of alpha e =',min_number_alpha

        call writesep(unit,0,0) ! Separator

        alpha=n_sites ! Number of alpha spins at start (max)
        allocate(f_eigenv(ms_spaces,max_deter_total+1)) ! Allocate results matrix
        f_eigenv=0.d0 ! Initialize
        allocate(arr_dets(ms_spaces)) ! Allocate array of number of determinants

        ! Main iteration
        do i=1,ms_spaces
            if (allocated(tab))then
                deallocate(tab)
                deallocate(det)
                deallocate(ham,evals,evect)
            end if

            beta=n_sites-alpha ! Calculate beta spins
            ms=dble(alpha-beta)/2.d0 ! Calculate ms value
            ndeter=combination(n_sites,beta) ! Calculate number of determinants

            ! Write down info
            write(unit,105)'Number of alpha/beta =',alpha,'/',beta
            write(unit,106)'Ms value =',ms
            write(unit,104)'Number of determinants =',ndeter

            call writesep(unit,0,1) ! Separator
            
            allocate(tab(ndeter)) ! Allocate tab array
            accu=0;compt=0 ! Initialize values for determinant generation
            call deter(unit,n_sites,alpha,tab,ndeter,accu,compt) ! Generate determinants
            call iquicksort(tab,ndeter,1,-1) ! Reverse the order of the determinants

            ! Print determinants as integers
            write(frmt,102)'(x,',ndeter,'(2x,i4))'
            write(unit,103)'Integer value of determinants: '
            call iarrout(unit,ndeter,tab) ! Print 1D array of integers
            
            call writesep(unit,0,1) ! Separator

            allocate(det(ndeter,n_sites)) ! Allocate boolean determinant array
            call Integ_to_Bool(n_sites,ndeter,tab,det) ! Transform integer to boolean arrays

            ! Print determinants with its logical value if they are not too many
            if ((ndeter<=100).or.(verb)) then
                write(frmt,102)'(x,',n_sites,'(x,L1))'
                write(unit,103)'Logical value of the determinants: '
                do j=1,ndeter
                    write(unit,frmt)(det(j,k),k=1,n_sites)
                end do
            else
                write(unit,'(x,a)')'The determinats are too many to print'
            end if

            allocate(ham(ndeter,ndeter)) ! Allocate Hamiltonian matrix

            ! Create Hamiltonian
            call create_ham(unit,verb,n_sites,ham,ndeter,det,couplings)

            call writesep(unit,0,1) ! Separator

            ! Print Hamiltonian matrix if size is not big
            if ((ndeter<=100).or.(verb)) then
                write(frmt,102)'(x,',ndeter,'(2x,f7.3))'
                write(unit,103)'Hamiltonian:'
                call matout(unit,ndeter,ndeter,ham,11,4) ! Print Hamiltonian matrix
            else
                write(unit,'(x,a)')'The Hamiltonian matrix is too big to print'
            end if

            allocate(evect(ndeter,ndeter),evals(ndeter)) ! Allocate arrays for eigenvectors and eigenvalues
            call subdsyevv(ham,evals,ndeter,evect) ! Diagonalize

            call writesep(unit,0,1) ! Separator

            ! Print eigenvalues and eigenvectors if they are not many
            if ((ndeter<=100).or.(verb)) then
                write(unit,103)'Eigenvalues and eigenvectors ='
                call eigenout(unit,n_sites,ndeter,evect,evals,det,11,4)
            else
                write(unit,'(x,a)')'The eigenvalues and eigenvectors are too many to print'
            end if

            arr_dets(i) = ndeter ! Store number of determinants
            f_eigenv(i,1) = ms   ! Store ms value
            f_eigenv(i,2:ndeter+1) = evals(1:ndeter) ! Store eigenvalues

            alpha=alpha-1 ! One alpha spin less

            call writesep(unit,0,0) ! Separator
        end do
        
        102 format (a,i5,a)
        103 format (/,x,a)
        104 format (x,a30,2x,i4)
        105 format (x,a30,2x,i3,x,a1,x,i3)
        106 format (x,a30,2x,f6.1)
    end subroutine iterate

    !***********************
    !*  Assign spin state multiplicity
    !*  -in
    !*  unit      : Unit to write output
    !*  verb      : Boolean flag to write debug info
    !*  f_eigenv  : Matrix with eigenvalues for each subspace
    !*  ms_spaces : Number of subspaces
    !*  arr_dets  : Array with number of determinant
    !*  -out
    !*  final_m   : Final ordered matrix with ms and S
    !***********************
    subroutine assign_ms(unit,verb,f_eigenv,ms_spaces,arr_dets,final_m)
        implicit none
        integer,intent(in)                   :: unit,ms_spaces
        integer,allocatable,intent(in)       :: arr_dets(:)
        real(kind=8),allocatable,intent(in)  :: f_eigenv(:,:)
        logical,intent(in)                   :: verb
        real(kind=8),allocatable,intent(out) :: final_m(:,:)

        integer                              :: i,j,k
        real(kind=8),allocatable             :: first(:),second(:)
        real(kind=8),parameter               :: thresh= 1e-6
        logical                              :: notfound
        
        ! Allocate and initialize final matrix
        allocate(final_m(ms_spaces+1,arr_dets(ms_spaces)+1))
        final_m(:,:)=0.d0

        final_m(1:ms_spaces,1)=f_eigenv(:,1) ! Assign ms values on first column of final array
        
        ! Assign S values to last row of final array
        do i=ms_spaces,1,-1
            final_m(ms_spaces+1,2:arr_dets(i)+1)=f_eigenv(i,1)
        end do
        
        ! Assign eigenvalues of first subspace on first row of final matrix
        do i=2,arr_dets(1)+1
            final_m(1,i)=f_eigenv(1,i)
        end do
        
        do i=1,ms_spaces-1  ! Iterate over the subspaces

            ! Allocate arrays
            if (allocated(first)) then
                deallocate(first)
            end if
            if (allocated(second)) then
                deallocate(second)
            end if

            ! Allocate first and second subspace arrays
            allocate(first(arr_dets(i)),second(arr_dets(i+1)))

            ! Assign subspaces
            ! First subspace is assigned from the final ordered matrix
            ! Second subspace is assigned from the unordered result
            first(:)=final_m(i,2:arr_dets(i)+1)
            second(:)=f_eigenv(i+1,2:arr_dets(i+1)+1)

            call quicksort(second,1,arr_dets(i+1),1) ! Sort second array, just in case

            if (verb) then ! Debug output
                write(unit,*)
                write(unit,'(x,a16,i4,a4,f4.1,a6,i4)')'Comparing space=',i,' ms=',f_eigenv(i,1),' size=',arr_dets(i)
                write(unit,'(x,a16,i4,a4,f4.1,a6,i4)')'      and space=',i+1,' ms=',f_eigenv(i+1,1),' size=',arr_dets(i+1)

                write(unit,'(/,x,a)')'First subspace:'
                call arrout(unit,size(first),first,11,4)
                write(unit,'(/,x,a)')'Second subspace:'
                call arrout(unit,size(second),second,11,4)

                call writesep(unit,0,1)
            end if
            
            notfound=.false. ! Initialize values for the algorithm
            j=1;k=1
            do while (j<arr_dets(i)+1)
                ! If values are equal
                if (abs(first(j)-second(k))<thresh) then
                    if (notfound) then ! If previously the value was not equal
                        call swap(second,j,k) ! Swap to maintain increasing order
                        notfound=.false. ! Set boolean to false
                        ! Increase the index of the first subspace,
                        ! and set index of the second equal to the first
                        j=j+1;k=j
                    else
                        ! If the last pair of values were equal, increase both indexes
                        j=j+1;k=k+1
                    end if
                else ! If values are not equal
                    ! Swap values if the indexes are not equal
                    if (j/=k) call swap(second,j,k)
                    notfound=.true. ! Set boolean value to true
                    k=k+1 ! Increase the index of the second subspace
                end if
            end do
            ! Set the i+1 row of the final matrix equal to the ordered array of the second subspace
            final_m(i+1,2:arr_dets(i+1)+1)=second(:)
        end do

        deallocate(first,second)

        ! Print the final ordered matrix
        write(unit,'(/,x,a,/)')'Final energy assignation:'
        call s_print(unit,ms_spaces,arr_dets,final_m,9,4)

    end subroutine assign_ms

    !***********************
    !*  Subroutine to find the ground state, the first excited state and their S
    !*  -in
    !*  unit      : Unit to write output
    !*  array     : Array to find states
    !*  -out
    !*  minv      : Ground state energy
    !*  smin      : S value of the ground state
    !*  secminv   : First excited state energy
    !*  ssecmin   : S value of the first excited state
    !*  gap       : Gap between states
    !***********************
    subroutine two_minimum(unit,array,minv,smin,secminv,ssecmin,gap)
        implicit None
        integer,intent(in)                   :: unit
        real(kind=8),allocatable,intent(in)  :: array(:,:)
        real(kind=8),intent(out)             :: minv,smin,secminv,ssecmin,gap
        real(kind=8),allocatable             :: energies(:),spins(:)
        integer                              :: dim1,dim2,minl,secminl
        logical,allocatable                  :: mask(:)
        real(kind=8),parameter               :: thresh = 1d-12
        
        ! Get dimensions of input array
        dim1 = size(array,1)
        dim2 = size(array,2)

        ! Allocate several arrays
        allocate(mask(dim2-1),energies(dim2-1),spins(dim2-1))

        ! Assign energy and spin arrays
        energies(:) = array(dim1-1,2:)
        spins(:) = array(dim1,2:)
        
        mask=.true. ! Set mask to true

        
        minv = minval(energies,1,mask)          ! Look for min val
        minl = minloc(energies,1,mask)          ! Get location of min value
        smin = spins(minl)                      ! Set S value using location of min val
        
        mask(minl) = .false. ! Set location of the min val as false in mask array
        
        do while (.true.)
            secminv = minval(energies,1,mask)       ! Look for second min val
            secminl = minloc(energies,1,mask)       ! Get location of second min val
            ssecmin = spins(secminl)                ! Set S value using location of min val
            if (abs(minv-secminv)<thresh) then      ! Avoid getting a degenerated state
                mask(secminl) = .false.
            else
                exit
            end if
        end do

        gap = secminv - minv ! Calculate gap

        if(abs(minv) < thresh) minv = 0d0       ! If minv is smaller than threshold, set as 0.d0
        if(abs(secminv) < thresh) secminv = 0d0 ! If secminv is smaller than threshold, set as 0.d0

        ! Print result
        write(unit,102)'Ground state energy =',minv,'(S =',smin,')'
        write(unit,102)'First excited state energy =',secminv,'(S =',ssecmin,')'
        write(unit,103)'Gap =',gap

        102 format (x,a30,2x,f9.3,2x,a4,f4.1,a1)
        103 format (x,a30,2x,f9.3)

    end subroutine two_minimum

    !***********************
    !*  Construct the Hamiltonian
    !*  -in
    !*  unit      : Unit to write output
    !*  verb      : Boolean flag to write debug info
    !*  n_sites   : Number of sites
    !*  ndeter    : Number of determinants
    !*  det       : Determinant array
    !*  couplings : Couplings array
    !*  -out
    !*  ham       : Hamiltonian matrix
    !***********************
    subroutine create_ham(unit,verb,n_sites,ham,ndeter,det,couplings)
        implicit none

        integer,intent(in)                   :: unit,ndeter,n_sites
        real(kind=8), allocatable,intent(in) :: couplings(:)
        logical,intent(in)                   :: verb
        logical*1,allocatable,intent(in)     :: det(:,:)
        real(kind=8),intent(out)             :: ham(ndeter,ndeter)

        integer                              :: i,j,trues
        integer,allocatable                  :: locs(:)
        logical*1                            :: diff(n_sites)
        logical*1,allocatable                :: det1(:),det2(:)

        character(100)                       :: frmt

        ! Allocate two arrays for not creating a ghost copy
        allocate(det1(size(det,2)),det2(size(det,2)))
        
        ham=0.d0 ! Initialize the Hamiltonian matrix
        
        if (verb) then ! Write separator if debug mode
            call writesep(unit,0,1)
        end if

        ! Loop over pairs of determinants
        do i=1,ndeter
            do j=i,ndeter
                ! Assign the determinants to variables (avoid ghost copy)
                det1(:)=det(i,:)
                det2(:)=det(j,:)

                if (i==j)then ! If diagonal term

                    if (verb) then ! Debug output
                        write(unit,'(x,a)')'Diagonal term'
                        write(frmt,'(a,i4,a,i4,a)')'(x,a5,',n_sites,'(L1),2x,a5,',n_sites,'(L1))'
                        write(unit,frmt)'Det1=',det1,'Det2=',det2
                    end if

                    ham(i,j)=diagonal(n_sites,det1,couplings) ! Set diagonal term

                else ! If not a diagonal term
                    ! Compare both determinants
                    call compare(n_sites,det1,det2,diff,locs,trues)

                    if (verb) then ! Debug output
                        write(unit,'(x,a)')'Extradiagonal term'
                        write(frmt,'(a,i4,a,i4,a,i4,a)')'(x,a5,',n_sites,'(L1),2x,a5,',n_sites,'(L1),2x,a4,',n_sites,'(L1))'
                        write(unit,frmt)'Det1=',det1,'Det2=',det2,'XOR=',diff
                        write(frmt,'(a,i4,a)')'(x,a11,x,',size(locs),'(i4))'
                        write(unit,frmt)'Locations =',locs
                    end if

                    if (trues==2)then ! If 2 and only 2 differences
                        ! Set extradiagonal terms i,j and j,i
                        ham(i,j)=-couplings(position_triang(locs(1),locs(2)))
                        ham(j,i)=-couplings(position_triang(locs(1),locs(2)))
                    end if
                end if
                ! Write separator if debug
                if ((i/=ndeter.or.j/=ndeter).and.verb)then
                    call writesep(unit,0,2)
                end if            
            end do
        end do
        
        deallocate(det1,det2) ! Deallocate arrays
        
    end subroutine create_ham

    !***********************
    !*  Calculate the diagonal terms of the Hamiltonian
    !*  -in
    !*  n_sites   : Number of sites
    !*  couplings : Couplings array
    !*  det       : Determinant of the diagonal
    !*  -out
    !*  f         : Diagonal term of the hamiltonian
    !***********************
    function diagonal(n_sites,det,couplings) result(f)
        implicit none
        integer,intent(in)                   :: n_sites
        real(kind=8), allocatable,intent(in) :: couplings(:)
        logical*1,intent(in)                 :: det(n_sites)

        integer                              :: i,j
        real(kind=8)                         :: f

        f=0.d0
        ! Check pair of sites and add the coupling if not equal
        do i=1,n_sites-1
            do j=i+1,n_sites
                if (det(i).neqv.det(j))then
                    f=f+couplings(position_triang(i,j))
                end if
            end do
        end do
    end function diagonal

    !***********************
    !*  Compare two determinants
    !*  -in
    !*  n_sites   : Number of sites
    !*  couplings : Couplings array
    !*  det1,det2 : Determinants to compare
    !*  -out
    !*  diff      : Array with booleans for different values
    !*  locations : Locations of the different values
    !*  n_trues   : Number of different values
    !***********************
    subroutine compare(n_sites,det1,det2,diff,locations,n_trues)
        implicit None
        
        integer,intent(in)                :: n_sites
        logical*1,intent(in)              :: det1(n_sites),det2(n_sites)
        integer,intent(out)               :: n_trues
        logical*1,intent(out)             :: diff(n_sites)
        integer,allocatable,intent(inout) :: locations(:)

        integer                           :: i,j
        
        ! Check if allocated
        if (allocated(locations))then
            deallocate(locations)
        end if

        ! Check intersection of two determinants
        do i=1,n_sites
            diff(i) = xor(det1(i),det2(i))
        end do
        
        n_trues =count(diff) ! Count number of trues in array

        allocate(locations(n_trues)) ! Allocate locations array
        
        ! Loop to store indexes of the values that are different
        j=1
        do i=1,n_sites
            if (diff(i))then
                locations(j)=i
                j=j+1
            end if
        end do
    end subroutine compare

    !***********************
    !*  This subroutine generate, recursively, 
    !*  all the determinants of a given number of sites and value of ms
    !*  -in
    !*  unit   : Unit to write
    !*  n      : Number of sites
    !*  m      : Number of alpha spins
    !*  ndeter : Number of determinants
    !*  -out
    !*  tab    : Determinants of dimension n (m "1" values and n-m "0" values)
    !*  accu   : Accumulator variables
    !*  compt  : Counter of the number of determinants generated 
    !***********************
    recursive subroutine deter(unit,n,m,tab,ndeter,accu,compt)
        implicit none
        integer, intent(in)  :: unit,n,m,ndeter
        integer, intent(out) :: tab(1:ndeter),accu,compt
        integer              :: accu1,i
        
        if (m>n) then ! If number of alpha spins > number of sites, stop
            write(unit,*) 'ERROR'
            write(unit,*) 'Number of alpha spins cannot be bigger than the number of sites.'
            stop
        endif

        ! If there are no more alpha spins then the next ones are beta
        ! Then you have to fill, in binary, with 0's
        ! That means that "you add 0" to the corresponding integer 
        if (m == 0) then       
            compt=compt+1      ! 1 determinant generated
            tab(compt)=accu    ! Determinant is accu
        ! If the number of remaining sites is equal to the number of alpha
        ! Then there are only alpha spins, in binary you have to fill with 1's
        else if (n == m) then 
            do i=0,n-1  
                accu=accu+2**i ! Max integer
            enddo
            compt=compt+1      ! 1 determinant generated  
            tab(compt)=accu    ! Integer value equal acc
        ! Else        
        else
            accu1=accu         ! Store accu value in accu1
            accu=accu+2**(n-1) ! Accu equal to 2^(n-1)
            ! Add a beta spin and move to the next site (i.e. one site less, same number of alpha)
            ! Call deter again for n-1 sites and m alpha spins
            call deter(unit,n-1,m,tab,ndeter,accu1,compt)       
            ! Add an alpha, move to the next site and reduce by 1 the number of remaining alpha spins
            ! Call deter again for n-1 sites and m-1 alpha spins
            call deter(unit,n-1,m-1,tab,ndeter,accu,compt)
        endif     
    end subroutine deter

    !***********************
    !*  Convert integers to boolean arrays
    !*  -in
    !*  nsites    : Number of sites
    !*  ndeter    : Number of determinants
    !*  tab       : Contains all the determinants in integer format
    !*  det       : Contains all the determinants in boolean format
    !***********************
    subroutine Integ_to_Bool(nsites,ndeter,tab,det)
        implicit none
        
        integer,intent(in)    :: nsites, ndeter,tab(ndeter)
        logical*1,intent(out) :: det(ndeter,nsites)

        logical*1             :: tmpdeter(nsites)
        integer               :: i, j, k

        ! Initialize tmpdeter array
        do k=1,nsites  
            tmpdeter(k)=.False.
        enddo
        
        do i=1,ndeter ! Iterate over all integers
            j=tab(i) ! Save integer as j
            do k=1,nsites ! Iterate over nsites
                if (mod(j,2)==0) then ! If module of j/2 is 0
                    tmpdeter(nsites-k+1)=.false. ! Set as false
                else ! If module is not 0
                    tmpdeter(nsites-k+1)=.true. ! Set as true
                end if
                j=j/2
                det(i,nsites-k+1)=tmpdeter(nsites-k+1) ! Set 
            enddo
        enddo
    end subroutine Integ_to_Bool

end module main_sub