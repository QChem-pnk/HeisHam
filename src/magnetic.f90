program Heisenberg

    use utils
    use main_sub
    use print_f

    implicit none

    integer                   :: n_sites, n_coupl, unit, ms_spaces
    integer,allocatable       :: arr_dets(:)
    real(kind=8), allocatable :: couplings(:),f_eigenv(:,:),final_m(:,:)
    character(:), allocatable :: input, output
    real(kind=8)              :: minv,smin,secminv,ssecmin,gap,start_t,end_t, elapsed_t
    logical                   :: verb

    call cpu_time(start_t) ! Register start time

    ! Parse arguments in the script call
    call parse_args(input, output)

    ! Open file 
    if (output/='None')then
        unit=70
        open(unit,file=output)
    else
        unit=6
    end if

    ! Call subroutine to read input
    call read_input(input,verb,n_sites,n_coupl,couplings)

    ! Write header of the program
    call writehead(unit,'Heisenberg',input,output)

    call writesep(unit,0,0) ! Separator

    ! Print info
    write(unit,'(x,a30,2x,i4)')'Number of sites =',n_sites
    write(unit,'(x,a30,2x,i4)')'Number of couplings =',n_coupl

    call writesep(unit,0,1) ! Separator

    ! Print couplings matrix
    write(unit,*)'Coupling matrix ='
    call print_couplings(unit,couplings,n_sites)

    ! Iterate over al ms spaces and get eigenvalues
    call iterate(unit,verb,n_sites,couplings,f_eigenv,ms_spaces,arr_dets)

    ! Assign S to the eigenvalues
    call assign_ms(unit,verb,f_eigenv,ms_spaces,arr_dets,final_m)

    call writesep(unit,3,0) ! Separator

    ! Search for ground state and first excited state
    call two_minimum(unit,final_m,minv,smin,secminv,ssecmin,gap)

    call writesep(unit,3,0) ! Separator

    call cpu_time(end_t) ! Register finish time

    elapsed_t = end_t - start_t ! Calculate elapsed time
    
    write(unit,'(x,a30,f10.5,a)')'Time elapsed =',elapsed_t,' s' ! Print time elapsed

    call writeend(unit) ! End of program

    call writeresume(input,output,n_sites,n_coupl,couplings,ms_spaces,arr_dets,final_m, &
                   & minv,smin,secminv,ssecmin,gap,elapsed_t) ! Write resume
    
    ! Open file 
    if (output/='None')then
        write(6,'(a,a,a)')'Results saved in ',output,'.'
        close(unit)
    end if

end program Heisenberg