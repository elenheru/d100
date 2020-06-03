subroutine test_cell_partition
    use array_parameters_mod
    use interaction_mod
    implicit none
    integer i
    !logical test_if_cell_is_close,lll
    !INTEGER, DIMENSION(3,42) :: LIST
    print*, "TEST_CELL_PARTITION zapuhxen"

    !real(8) R_curr(1,cnt_q)
    !integer an_by_cn(3,atoms_max_array)
    !integer cn_by_an(0:neibors_max,atoms_max_array)

    !proveritq zadany li koordinaty atomov
    if (atoms__in_total .le. 0) stop "gde atomy? kolicxestvo atomov slisxkom malo."


!        CALL LIST_OF_42_NEAREST_CELLS_FCC(0,0,0,LIST)
!        DO I=1,42
!            PRINT*,LIST(1:3,I)
!        ENDDO
        !LLL=TEST_IF_CELL_IS_CLOSE(0,0,1)
        !PRINT"(A,L2)", "---", LLL
        !PRINT"(A,3(I7.4,1X))", "(/0,7,2 /)+(/-1,-4,-3/) = ", (/0,7,2 /)+(/-1,-4,-3/)
        !STOP ""

    do i=8000,30500,125

        !vybratq cutoff v cikle
        cutoff_param=cutoff*1d-4*i

        !zapolnitq massiv an_by_cn
        call test_shift_lattice
        call fill_an_by_cn_4d

        call wo_xyz_rhcells_4d
        !STOP "filling an by cn test"
    enddo


endsubroutine test_cell_partition

subroutine test_cell_neibors_4d
    use array_parameters_mod
    use interaction_mod
    implicit none
    integer i

    print*, "TEST_CELL_NEIBORS_4D zapuhxen"
    if (atoms__in_total .le. 0) stop "gde atomy? kolicxestvo atomov slisxkom malo."

    do i=8000,30500,125

        cutoff_param=cutoff*1d-4*i
        call test_shift_lattice
        call fill_an_by_cn_4d
        call test_lists_comparison_3d4d

    enddo

    contains
        subroutine test_lists_comparison_3d4d

            integer :: getrandom_int_fromto_inclusive
            integer, dimension(3) :: tc_fcc
            integer, dimension(4) :: tc_pc4d

            integer, dimension(4) :: nc_pc4d
            integer, dimension(3) :: nc_fcc
            integer, dimension(4) :: rc_pc4d

            integer, dimension(3,42) :: list_via3d
            integer, dimension(4,42) :: list_via4d
            integer :: i


            interface
                function fcc_to_pc4d(fcc_cell) result(pc4d_cell)
                implicit none
                integer, dimension(3),intent(in)  :: fcc_cell
                integer, dimension(4) :: pc4d_cell
                endfunction fcc_to_pc4d
            endinterface


            tc_fcc(1) = getrandom_int_fromto_inclusive(-4,4)
            tc_fcc(2) = getrandom_int_fromto_inclusive(-4,4)
            tc_fcc(3) = getrandom_int_fromto_inclusive(-4,4)
            tc_pc4d = fcc_to_pc4d(tc_fcc)


            print"(A,4(1x,I3.1),A)", "proverka spiskov dlqa jacxejki",tc_pc4d, " pc4d"
            print"(A,4(1x,I3.1),A)", "proverka spiskov dlqa jacxejki",tc_fcc,-9, "  FCC"

            call list_of_42_nearest_cells_4dpc(&
            tc_pc4d(1),tc_pc4d(2),tc_pc4d(3),tc_pc4d(4),list_via4d)

            call list_of_42_nearest_cells_fcc(&
            tc_fcc(1),tc_fcc(2),tc_fcc(3),list_via3d)

            do i=1,42

                nc_fcc = list_via3d(1:3,i)
                nc_pc4d = fcc_to_pc4d(nc_fcc)
                rc_pc4d = list_via4d(1:4,i)

                if(&
                    ( nc_pc4d(1) .ne. rc_pc4d(1) ) .or. &
                    ( nc_pc4d(2) .ne. rc_pc4d(2) ) .or. &
                    ( nc_pc4d(3) .ne. rc_pc4d(3) ) .or. &
                    ( nc_pc4d(4) .ne. rc_pc4d(4) ) &
                ) then
                    print*," ne sovpali sosedi "
                    print*,tc_fcc
                    print*,tc_pc4d
                endif
            enddo

        endsubroutine test_lists_comparison_3d4d

endsubroutine test_cell_neibors_4d

subroutine test_cell_gap
    use array_parameters_mod
    use interaction_mod
    implicit none
    integer i
    !logical test_if_cell_is_close,lll
    !INTEGER, DIMENSION(3,42) :: LIST
    print*, "TEST_CELL_GAP zapuhxen"


    !proveritq zadany li koordinaty atomov
    if (atoms__in_total .le. 0) stop "gde atomy? kolicxestvo atomov slisxkom malo."

    do i=9000,11000,27
        print*, i, " eto scxqotcxik"
        cutoff_param=cutoff*1d-4*i
        call test_shift_lattice
        call fill_an_by_cn_4d
        call fill_cn_by_an_4d
        call find_too_close_pairs
        call wo_xyz_rhcells_4d
        !STOP "filling an by cn test"
    enddo

    contains

        subroutine find_too_close_pairs
            use array_parameters_mod
            use interaction_mod
            use positions_mod
            implicit none
            integer ia1,ia2,i_l
            real(8) r_rel(3)
            real(8) dist
            integer, dimension(4,42) :: list
            integer, dimension(4) :: cell__1,cell_nb
            logical :: is_cell_present
            do ia1 = 2,atoms__in_total
                do ia2 = 1,ia1-1

                    r_rel(1) = R_curr(1,ia1) - R_curr(1,ia2)
                    r_rel(2) = R_curr(2,ia1) - R_curr(2,ia2)
                    r_rel(3) = R_curr(3,ia1) - R_curr(3,ia2)
                    dist = norm2(r_rel)

                    if (dist .gt. cutoff_param) cycle

                    cell__1=an_by_cn_4d(1:4,ia1)
                    cell_nb=an_by_cn_4d(1:4,ia2)

                    call list_of_42_nearest_cells_4dpc(&
                    cell__1(1),cell__1(2),cell__1(3),cell__1(4),list)

                    is_cell_present = .false. .or. (&
                        (cell_nb(1) .eq. cell_nb(1) ) .and.  &
                        (cell_nb(2) .eq. cell_nb(2) ) .and.  &
                        (cell_nb(3) .eq. cell_nb(3) ) .and.  &
                        (cell_nb(4) .eq. cell_nb(4) ) &
                    )

                    do i_l=1,42
                        is_cell_present = is_cell_present .or. (&
                            (cell_nb(1) .eq. list(1,i_l) ) .and.  &
                            (cell_nb(2) .eq. list(2,i_l) ) .and.  &
                            (cell_nb(3) .eq. list(3,i_l) ) .and.  &
                            (cell_nb(4) .eq. list(4,i_l) ) &
                        )
                    enddo
                    if (is_cell_present) then
                        cycle
                    else
                        print*,"Atoms are too close ",ia1,ia2
                        print*,R_curr(:,ia1), cell__1
                        print*,R_curr(:,ia2), cell_nb
                        print*," list is "
                        do i_l=1,42
                            print*, list(1:4,i_l)
                        enddo

                        stop " too bad "
                    endif

                enddo
            enddo
        print*, "Ladno, ni odin atom ne usxel…"
        endsubroutine find_too_close_pairs

endsubroutine test_cell_gap

subroutine test_shift_lattice
    use positions_mod
    implicit none
    integer i_atoms
    real(8) :: rs(3)
    call random_number(rs)
    rs=(2d0*rs-(/1d0,1d0,1d0/))*7d-1
    print*,"Resxotka celikom sdvinuta na slucxajnyj vektor: ",rs
    do i_atoms=1,atoms__in_total
        R_curr(1:3,i_atoms)=R_curr(1:3,i_atoms)+rs(1:3)
    enddo
    !R_perf=R_curr
endsubroutine test_shift_lattice

subroutine test_factor_phi
    use ackland2003potential_fe_mod
    implicit none
    !real(8) pw_fefe_an

    print*,"pw_fefe_an(1d0-epsilon(1d0))  ",pw_fefe_an(1d0-epsilon(1d0))
    print*,"pw_fefe_an(1d0+epsilon(1d0))  ",pw_fefe_an(1d0+epsilon(1d0))
    print*,pw_fefe_an(1d0+epsilon(1d0))/pw_fefe_an(1d0-epsilon(1d0))

endsubroutine test_factor_phi

subroutine test_po_dependencies
    use ackland2003potential_fe_mod
    implicit none
!    real(8) pw_fefe_an
!    real(8) ed_fefe_an
!    real(8) mf_fe_an
    real(8) pw_fefe_la
    real(8) ed_fefe_la
    real(8) mf_fe_la

    real(8) pw_fefe_dif
    real(8) ed_fefe_dif
    real(8) mf_fe_dif

    real(8) pw_fefe_difsum
    real(8) ed_fefe_difsum
    real(8) mf_fe_difsum
    real(8) argmax,arg

    integer i,intpwmax,intedmax,intmfmax

    pw_fefe_difsum = 0d0
    ed_fefe_difsum = 0d0
    mf_fe_difsum = 0d0

    print*, "test_po_dependencies ZAPUHXEN "
    call initiate_energetic_parameters
    open (1001, file = "pw_fefe_an_la.txt" )
    open (1002, file = "ed_fefe_an_la.txt" )
    open (1003, file = "mf_fe___an_la.txt" )

    !print*,isnan(pw_fefe_difsum),isnan(1d0),(pw_fefe_difsum)
    argmax = 6d0
    intpwmax = nint(argmax*1d4)
    do i=100,intpwmax+100
        arg=i*1d-4
        pw_fefe_dif = ( pw_fefe_an(arg)-pw_fefe_la(arg) )
        !pw_fefe_dif = ( pw_fefe_an(arg)-pw_fefe_la(arg) ) / ( pw_fefe_an(arg)+pw_fefe_la(arg) ) * 2d0
        !if(isnan(pw_fefe_dif-pw_fefe_dif)) pw_fefe_dif = 0d0
        pw_fefe_difsum = pw_fefe_difsum + pw_fefe_dif
        write(1001,601) arg,pw_fefe_an(arg),pw_fefe_la(arg),pw_fefe_dif
    enddo
    close(1001)

    argmax = 6d0
    intedmax = nint(argmax*1d4)
    do i=1,intedmax
        arg=i*1d-4
        ed_fefe_dif = (ed_fefe_an(arg)-ed_fefe_la(arg))
        !ED_FEFE_DIF = (ED_FEFE_AN(ARG)-ED_FEFE_LA(ARG))/(ED_FEFE_AN(ARG)+ED_FEFE_LA(ARG))*2D0
        !IF(ISNAN(ED_FEFE_DIF-ED_FEFE_DIF)) ED_FEFE_DIF = 0D0
        ed_fefe_difsum = ed_fefe_difsum + ed_fefe_dif
        write(1002,601) arg,ed_fefe_an(arg),ed_fefe_la(arg),ed_fefe_dif
    enddo
    close(1001)

    argmax = 60d0
    intmfmax = nint(argmax*1d2)
    do i=1,intmfmax
        arg=i*2d-2
        mf_fe_dif = (mf_fe_an(arg)-mf_fe_la(arg))
        !MF_FE_DIF = (MF_FE_AN(ARG)-MF_FE_LA(ARG))/(MF_FE_AN(ARG)+MF_FE_LA(ARG))*2D0
        !IF(ISNAN(MF_FE_DIF-MF_FE_DIF)) MF_FE_DIF = 0D0
        mf_fe_difsum = mf_fe_difsum + mf_fe_dif
        write(1003,601) arg,mf_fe_an(arg),mf_fe_la(arg),mf_fe_dif
    enddo
    close(1001)
    print 602,"Sumof ",intpwmax," parts and total average error for pairwise ion ion   ",pw_fefe_difsum,pw_fefe_difsum/intpwmax
    print 602,"Sumof ",intedmax," parts and total average error for electronic density ",ed_fefe_difsum,ed_fefe_difsum/intedmax
    print 602,"Sumof ",intmfmax," parts and total average error for embedding functon  ",mf_fe_difsum,mf_fe_difsum/intmfmax

    !PRINT*,ISNAN(PW_FEFE_DIFSUM),ISNAN(1D0),(PW_FEFE_DIFSUM)
    !IF(ISNAN(PW_FEFE_DIFSUM)) STOP "NAN NANA"

601 format(SP,4(ES12.5e2,1x))
602 format(SP,A,I6.2,A,2(ES12.5e2,1x))
endsubroutine test_po_dependencies

subroutine test_find_neibors_comparison
    use lists_assortiment_mod
    USE INTERACTION_MOD
    use array_parameters_mod
    implicit none

    integer i,omp_get_thread_num,rep,repmax
    INTEGER(8) :: ct1,ct2,ct3,ct4,ct5
    INTEGER(8) :: ct6,ct7,ct8,ct9,ct0

    print*, "test_find_neibors_comparison ZAPUHXEN "
    print*, "sravnim skorostq polucxenija sosedej dlqa topornogo algoritma i dlqa jacxeecxnogo "
    !STOP "pocxemu malo sosedej"
    call fill_an_by_cn_4d
    call fill_cn_by_an_4d

    repmax=40

    CALL SYSTEM_CLOCK(ct1)
        !DO REP=1,REPMAX
    do i = 1,atoms__in_total
        call get_direct_list(i)
    enddo
        !ENDDO
    CALL SYSTEM_CLOCK(ct2)

        DO REP=1,REPMAX
    do i = 1,atoms__in_total
        call get_direct_list_test_exhaustive(i)
        !print*,i
    enddo
        ENDDO
    CALL SYSTEM_CLOCK(ct3)

        DO REP=1,REPMAX
    !$OMP PARALLEL
    print*, omp_get_thread_num()
    if(omp_get_thread_num() .eq. 0) then
        do i = 1 , atoms__in_total/4
        !call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 1) then
        do i = 1 +   atoms__in_total/4 , 2*atoms__in_total/4
        !call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 2) then
        do i = 1 + 2*atoms__in_total/4 , 3*atoms__in_total/4
        !call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 3) then
        do i = 1 + 3*atoms__in_total/4 ,   atoms__in_total
        !call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    !$OMP END PARALLEL
        ENDDO
    CALL SYSTEM_CLOCK(ct4)
        DO REP=1,REPMAX
    !$OMP PARALLEL PRIVATE(I)
    print*, omp_get_thread_num()
    if(omp_get_thread_num() .eq. 0) then
        do i = 1 , atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 1) then
        do i = 1 +   atoms__in_total/4 , 2*atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 2) then
        do i = 1 + 2*atoms__in_total/4 , 3*atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    if(omp_get_thread_num() .eq. 3) then
        do i = 1 + 3*atoms__in_total/4 ,   atoms__in_total
        call get_direct_list_test_exhaustive(i)
        !call get_direct_list(i)
        enddo
    endif
    !$OMP END PARALLEL
        ENDDO
    CALL SYSTEM_CLOCK(ct5)

        DO REP=1,REPMAX
    !$OMP SECTIONS PRIVATE(I)
    print*, omp_get_thread_num()
    !$OMP SECTION
    do i = 1 , atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
    enddo
    !$OMP SECTION
    do i = 1 +   atoms__in_total/4 , 2*atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
    enddo
    !$OMP SECTION
    do i = 1 + 2*atoms__in_total/4 , 3*atoms__in_total/4
        call get_direct_list_test_exhaustive(i)
    enddo
    !$OMP SECTION
    do i = 1 + 3*atoms__in_total/4 ,   atoms__in_total
        call get_direct_list_test_exhaustive(i)
    enddo
    !$OMP END SECTIONS
        ENDDO

    CALL SYSTEM_CLOCK(ct6)

    write(*,*) "ticks get_direct_list            ",(ct2-ct1)*1d-6," ms"
    write(*,*) "ticks get_direct_list exhaustive ",(ct3-ct2)*1d-6," ms"
    write(*,*) "ticks get_direct_list exhaustive ",(ct4-ct3)*1d-6," ms (parallel) nothing to do"
    write(*,*) "ticks get_direct_list exhaustive ",(ct5-ct4)*1d-6," ms (parallel)"
    write(*,*) "ticks get_direct_list exhaustive ",(ct6-ct5)*1d-6," ms (parallel sections)"
    write(*,*) "ticks sum ",(ct3-ct1)*1d-6," ms"

endsubroutine test_find_neibors_comparison

subroutine test_find_neibors_comparison_2
    use lists_assortiment_mod
    USE INTERACTION_MOD
    use array_parameters_mod
    implicit none

    integer i,omp_get_thread_num,rep,repmax
    INTEGER(8) :: ct1,ct2,ct3,ct4,ct5
    INTEGER(8) :: ct6,ct7,ct8,ct9,ct0

    print*, "test_find_neibors_comparison ZAPUHXEN "
    print*, "sravnim skorostq polucxenija sosedej dlqa topornogo algoritma i dlqa jacxeecxnogo "
    !STOP "pocxemu malo sosedej"
    call fill_an_by_cn_4d
    call fill_cn_by_an_4d

    repmax=40


    CALL SYSTEM_CLOCK(ct2)
        DO REP=1,REPMAX
    do i = 1,atoms__in_total
        call get_direct_list_test_exhaustive(i)
    enddo
        ENDDO
    CALL SYSTEM_CLOCK(ct3)
        DO REP=1,REPMAX
    do i = 1,atoms__in_total
        call get_direct_list_test_exhaustive_linrej(i)
    enddo
        ENDDO
    CALL SYSTEM_CLOCK(ct4)



    print*, "ticks get_direct_list exhaustive ",(ct3-ct2)*1d-6," ms. only distance"
    print*, "ticks get_direct_list exhaustive ",(ct4-ct3)*1d-6," ms. linear rejections"


endsubroutine test_find_neibors_comparison_2

subroutine test_big_array
    implicit none
    real(8), allocatable :: bigarray(:)
    real :: start, finish
    integer i,j
    integer(8) k
    call cpu_time(start)
    do i=1,nint(2d9),nint(1d4)
        !do j=1,nint(2d9)
        allocate(bigarray(i))
        !print*, sizeof(bigarray), real(i)*8
        deallocate(bigarray)
        !enddo
    enddo
    call cpu_time(finish)
    print*, sizeof(bigarray), "time is ", -start + finish
    !call sleep(3)
    return
endsubroutine test_big_array

subroutine test_system_clock
  INTEGER :: count1, count_rate, count_max
  INTEGER(8) :: count2, count_rate2, count_max2,count3
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  WRITE(*,*) count1, count_rate, count_max
  CALL SYSTEM_CLOCK(count2, count_rate2, count_max2)
  WRITE(*,*) count2, count_rate2, count_max2
  CALL SYSTEM_CLOCK(count3, count_rate2, count_max2)
  WRITE(*,*) count3, count_rate2, count_max2
endsubroutine
