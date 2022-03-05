module print_f

    use utils

    public::writeresume          ! Print a resume
    public::s_print              ! Print final matrix with assigned ms and S
    public::eigenout             ! Print eigenvalues and eigenvectors
    public::matout               ! Print matrix of dimension NxM
    public::arrout               ! Print array of dimension N
    public::iarrout              ! Print integer array of dimension N
    public::print_couplings      ! Print coupling matrix
    public::writehead            ! Write header
    public::writesep             ! Write separator
    public::writeend             ! Write ending of program

contains

    !***********************
    !*  Print resume
    !*  -in
    !*  input     : Input name
    !*  n_sites   : Number of sites
    !*  n_coupl   : Number of couplings
    !*  couplings : Coupling array
    !*  ms_spaces : Number of sites
    !*  arr_dets  : Number of determinants for each subspace
    !*  final_m   : Final ordered matrix with ms and S
    !*  minv      : Ground state energy
    !*  smin      : S value of the ground state
    !*  secminv   : First excited state energy
    !*  ssecmin   : S value of the first excited state
    !*  gap       : Gap between states
    !*  elapsed_t : Elapsed time
    !***********************
    subroutine writeresume(input, output,n_sites,n_coupl,couplings,ms_spaces,arr_dets, &
        & final_m,minv,smin,secminv,ssecmin,gap,elapsed_t)
        integer,intent(in)                        :: n_sites,n_coupl,ms_spaces
        integer,allocatable,intent(in)            :: arr_dets(:)
        real(kind=8),allocatable,intent(in)       :: couplings(:),final_m(:,:)
        real(kind=8),intent(in)                   :: minv,smin,secminv,ssecmin,gap,elapsed_t
        !integer                                   :: index
        character(*),intent(in)                   :: input,output
        character(100)                            :: resume

        ! Generate resume file name
        if (output=='None') then
            resume = remove_ext(input)
            resume = trim(resume)//'.resume.out' ! Add an extension
        else
            resume = remove_ext(output)
            resume = trim(resume)//'.resume.out' ! Add an extension
        end if

        ! Generate resume
        open(50,file=resume)

        write(50,*)'--------------------------  RESUME  --------------------------'
        call writesep(50,0,0)
        write(50,100)'Input =',input
        call writesep(50,0,0)
        write(50,101)'Number of sites =',n_sites
        write(50,101)'Number of couplings =',n_coupl
        write(50,101)'Number of Ms subspaces =',ms_spaces
        write(50,101)'Max number of determinants =',arr_dets(ms_spaces)
        call writesep(50,0,0)
        write(50,104)'Coupling matrix ='
        call print_couplings(50,couplings,n_sites)
        call writesep(50,0,0)
        write(50,104)'Final energy assignation ='
        call s_print(50,ms_spaces,arr_dets,final_m,9,4)
        call writesep(50,0,0)
        write(50,102)'Ground state energy =',minv,'(S =',smin,')'
        write(50,102)'First excited state energy =',secminv,'(S =',ssecmin,')'
        write(50,103)'Gap =',gap
        call writesep(50,0,0)
        write(50,105)'Time elapsed =',elapsed_t,' s'
        call writesep(50,0,0)

        close(50)

        write(6,'(a,a,a)')'Resume file saved in ',trim(resume),'.'

        100 format (x,a28,2x,a)
        101 format (x,a28,2x,i4)
        102 format (x,a28,2x,f9.4,2x,a4,f4.1,a1)
        103 format (x,a28,2x,f9.4)
        104 format (x,a28,/)
        105 format (x,a28,2x,f10.5,a)
    end subroutine

    !***********************
    !*  Print the energies classified by ms and S with precision Fent.dec
    !*  -in
    !*  unit      : Unit to write
    !*  ms_spaces : Number of sites
    !*  arr_dets  : Number of determinants for each subspace
    !*  final_m   : Final matrix with
    !*  ent       : Total size of float
    !*  dec       : Number of decimal places
    !***********************
    subroutine s_print(unit,ms_spaces,arr_dets,final_m,ent,dec)
        implicit none

        integer,parameter             :: ncol = 6
        real(kind=8),parameter        :: thresh = 1d-12

        integer,intent(in)            :: unit,ms_spaces,arr_dets(ms_spaces),ent,dec
        real(kind=8),intent(in)       :: final_m(ms_spaces+1,arr_dets(ms_spaces)+1)

        integer                       :: ilower,iupper,num,i,j,spaces,n,k,l,t
        real(kind=8)                  :: B(ncol),ms
        character(100)                :: frmt1,frmt2,frmt3

        ! Set the format
        spaces = 2+ent-6
        write(frmt1,'(a,i2,a,i2,a)')'(3x,10(',spaces,'x,i6))'
        write(frmt2,'(a,i2,a,i2,a)')'(x,f3.1,10(2x,f',ent,'.',dec,'))'
        write(frmt3,'(a,i2,a)')'(4x,10(2x,f',ent,'.1),a2)'

        n = arr_dets(ms_spaces)
        do ilower=1,n,ncol
            if (ilower/=1)then ! Print space between loops
                write(unit,*)''
            end if
            iupper = min(ilower + ncol - 1,n) ! Last column to print in this cycle
            num = iupper - ilower + 1
            write(unit,frmt1) (j,j=ilower,iupper) ! Write column numbers
            write(unit,'(x,a2)')'Ms'
            do i=1,ms_spaces ! Iterate over rows
                ms = final_m(i,1)
                k=arr_dets(i)+1 ! Set number of columns for this subspace
                if (k<ilower)then ! If all columns have been printed, print only the ms
                    write(unit,'(x,f3.1)')ms
                else
                    l = min(iupper+1,k) ! Set maximum columns to print for each subspace
                    do j=ilower+1,l
                        B(j-ilower) = final_m(i,j)
                    enddo
                    t = l-ilower
                    do j=1,t ! Change small values with 0.d0
                        if(abs(B(j)) < thresh) B(j) = 0d0
                    enddo
                    write(unit,frmt2) ms,(B(j),j=1,t) ! Print columns for this row
                end if
            enddo
            ! Print S values
            do j=ilower,iupper
                B(j-ilower+1) = final_m(ms_spaces+1,j+1)
            enddo
            write(unit,frmt3,advance='no') (B(j),j=1,num)
            write(unit,*)' S'
        enddo
    end subroutine s_print

    !***********************
    !*  Print the eigenvalues and eigenvectors with precision Fent.dec
    !*  -in
    !*  unit    : Unit to write
    !*  n_sites : Number of sites
    !*  ndets   : Number of determinants
    !*  evec    : Eigenvectors
    !*  eval    : Eigenvalues
    !*  dets    : Determinant array
    !*  ent     : Total size of float
    !*  dec     : Number of decimal places
    !***********************
    subroutine eigenout(unit,n_sites,ndets,evec,eval,dets,ent,dec)
        implicit none

        integer,parameter             :: ncol = 5
        real(kind=8),parameter        :: thresh = 1d-12

        integer,intent(in)            :: unit,n_sites,ndets,ent,dec
        real(kind=8),intent(in)       :: evec(ndets,ndets),eval(ndets)
        logical*1,intent(in)          :: dets(ndets,n_sites)

        integer                       :: ilower,iupper,num,i,j,k
        integer                       :: spaces,spaces_2,spaces_3
        real(kind=8)                  :: B(ncol)
        character(100)                :: frmt1,frmt2,frmt3,frmt4

        ! Set the format
        spaces = ent-4
        spaces_2 = n_sites-3
        spaces_3 = n_sites+1
        
        write(frmt2,'(a,i2,a,i2,a,i2,a)')'(x,a6,',spaces_2,'x,10(2x,f',ent,'.',dec,'))'
        write(frmt3,'(a,i2,a,i2,a,i2,a)')'(x,',n_sites,'(l1),3x,10(2x,f',ent,'.',dec,'))'

        do ilower=1,ndets,ncol
            if (ilower/=1)then ! Print space between loops
                write(unit,*)''
            end if
            iupper = min(ilower + ncol - 1,ndets) ! Last column to print in this cycle
            num = iupper - ilower + 1
            write(frmt1,'(a,i2,a,i2,a,i2,a)')'(',spaces_3+3,'x,',num,'(',spaces-1,'x,"E("i4")"))'
            write(unit,frmt1) (j,j=ilower,iupper) ! Write column numbers
            do j=ilower,iupper
                B(j-ilower+1) = eval(j)
            enddo
            do j=1,num ! Change small values with 0.d0
                if(abs(B(j)) < thresh) B(j) = 0d0
            enddo
            write(unit,frmt2)'Energy',(B(j),j=1,num)
            write(frmt4,'(a,i2,a,i2,a,i2,a)')'(',spaces_3+3,'x,',num,'(',spaces-1,'x,"Ψ("i4")"))'
            ! Print wave number
            write(unit,frmt4)(j,j=ilower,iupper)
            do i=1,ndets ! Iterate over rows
                do j=ilower,iupper
                    B(j-ilower+1) = evec(i,j)
                enddo
                do j=1,num ! Change small values with 0.d0
                    if(abs(B(j)) < thresh) B(j) = 0d0
                enddo
                ! Print determinant and columns for this row
                write(unit,frmt3) (dets(i,k),k=1,n_sites),(B(j),j=1,num) 
            enddo
        enddo
    end subroutine eigenout

    !***********************
    !*  Print the MxN array A with precision Fent.dec
    !*  -in
    !*  unit : Unit to write
    !*  m,n  : Dimensions of the array
    !*  A    : Array to print
    !*  ent  : Total size of float
    !*  dec  : Number of decimal places
    !***********************
    subroutine matout(unit,m,n,A,ent,dec)
        implicit none

        integer,parameter             :: ncol = 5
        real(kind=8),parameter        :: thresh = 1d-12

        integer,intent(in)            :: unit,m,n,ent,dec
        real(kind=8),intent(in)       :: A(m,n)

        integer                       :: ilower,iupper,num,i,j,spaces
        real(kind=8)                  :: B(ncol)
        character(100)                :: frmt1,frmt2

        ! Set the format
        spaces = 2+ent-6
        write(frmt1,'(a,i2,a,i2,a)')'(3X,10(',spaces,'x,i6))'
        write(frmt2,'(a,i2,a,i2,a)')'(i5,10(2x,f',ent,'.',dec,'))'

        do ilower=1,n,ncol
            iupper = min(ilower + ncol - 1,n) ! Last column to print in this cycle
            num = iupper - ilower + 1
            write(unit,frmt1) (j,j=ilower,iupper) ! Write column numbers
            do i=1,m ! Iterate over rows
                do j=ilower,iupper
                    B(j-ilower+1) = A(i,j)
                enddo
                do j=1,num ! Change small values with 0.d0
                    if(abs(B(j)) < thresh) B(j) = 0d0
                enddo
                write(unit,frmt2) i,(B(j),j=1,num) ! Print columns for this row
            enddo
        enddo
    end subroutine matout

    !***********************
    !*  Print the N array A with precision Fent.dec
    !*  -in
    !*  unit : Unit to write
    !*  n    : Dimension of the array
    !*  A    : Array to print
    !*  ent  : Total size of float
    !*  dec  : Number of decimal places
    !***********************
    subroutine arrout(unit,n,A,ent,dec)
        implicit none

        integer,parameter             :: ncol = 5
        real(kind=8),parameter        :: thresh = 1d-12

        integer,intent(in)            :: unit,n,ent,dec
        real(kind=8),intent(in)       :: A(n)

        integer                       :: ilower,iupper,num,j,spaces
        real(kind=8)                  :: B(ncol)
        character(100)                :: frmt1,frmt2

        ! Set the format
        spaces = 2+ent-6
        write(frmt1,'(a,i2,a,i2,a)')'(3X,10(',spaces,'x,i6))'
        write(frmt2,'(a,i2,a,i2,a)')'(5X,10(2x,f',ent,'.',dec,'))'

        do ilower=1,n,ncol
            iupper = min(ilower + ncol - 1,n) ! Last column to print in this cycle
            num = iupper - ilower + 1
            write(unit,frmt1) (j,j=ilower,iupper) ! Write column numbers
            do j=ilower,iupper
                B(j-ilower+1) = A(j)
            enddo
            do j=1,num ! Change small values with 0.d0
                if(abs(B(j)) < thresh) B(j) = 0d0
            enddo
            write(unit,frmt2) (B(j),j=1,num) ! Print columns for this row
        enddo
    end subroutine arrout

    !***********************
    !*  Print the N array A (integers)
    !*  -in
    !*  unit : Unit to write
    !*  n    : Dimension of the array
    !*  A    : Array to print
    !*  ent  : Total size of float
    !*  dec  : Number of decimal places
    !***********************
    subroutine iarrout(unit,n,A)
        implicit none

        integer,parameter             :: ncol = 5

        integer,intent(in)            :: unit,n,A(n)

        integer                       :: ilower,iupper,num,j,B(ncol)
        character(100)                :: frmt

        ! Set the format
        write(frmt,'(a,i2,a,i2,a)')'(5X,10(2x,i8))'

        do ilower=1,n,ncol
            iupper = min(ilower + ncol - 1,n) ! Last column to print in this cycle
            num = iupper - ilower + 1
            do j=ilower,iupper
                B(j-ilower+1) = A(j)
            enddo
            write(unit,frmt) (B(j),j=1,num) ! Print columns for this row
        enddo
    end subroutine iarrout

    !***********************
    !*  Create header
    !*  -in
    !*  unit  :  Unit to write
    !*  nom   :  Name of program
    !*  input :  Input name
    !*  output:  Output name
    !***********************
    subroutine writehead (unit, nom, input, output)
        implicit none
        integer,intent(in)        :: unit
        character*(*),intent(in)  :: nom,input,output

        write(unit,*) '--------------------------------------------------------------'
        write(unit,*) 'Fortran program: ',nom
        write(unit,*) ' '
        write(unit,*) 'Input: ',input
        write(unit,*) 'Output: ',output
        write(unit,*) '--------------------------------------------------------------'
        write(unit,*) ' '
    end subroutine writehead

    !***********************
    !*  Create separator
    !*  -in
    !*  unit :   Unit to write
    !*  i    :   Whether to write an empty line before and after the separator
    !*       0:  No empty lines
    !*       1:  Empty line before
    !*       2:  Empty line after
    !*       3:  Empty line before and after
    !*  j    :   Type of separator
    !*       0:  Big separator
    !*       1:  Medium separator 
    !*       2+: Small separator
    !***********************
    subroutine writesep(unit,i,j)
        implicit none
        integer,intent(in) :: unit,i,j

        if (i==1.or.i==3) then
            write(unit,*) ' '
        end if
        if (j==0)then
            write(unit,*) '--------------------------------------------------------------'
        elseif (j==1)then
            write(unit,*) '-----------------------------'
        else
            write(unit,*) '--------------'
        end if
        if (i==2.or.i==3) then
            write(unit,*) ' '
        end if
    end subroutine writesep
    
    !***********************
    !*  Create footer
    !*  -in
    !*  unit : Unit to write
    !***********************
    subroutine writeend(unit)
        integer,intent(in) :: unit

        write(unit,*) ' '
        write(unit,*) '--------------------------------------------------------------'
        write(unit,*) '---------------------  Program finished  ---------------------'
        write(unit,*) '--------------------------------------------------------------'
    end subroutine writeend

end module print_f