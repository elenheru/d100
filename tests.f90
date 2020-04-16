subroutine test_cell_partition
    use array_parameters_mod
    use interaction_mod

    integer i
    real(8) cutoff_ideal
    print*, "TEST_CELL_PARTITION zapuhxen"

    !real(8) R_curr(1,cnt_q)
    !integer an_by_cn(3,atoms_max_array)
    !integer cn_by_an(0:neibors_max,atoms_max_array)

    !proveritq zadany li koordinaty atomov
    if (atoms__in_total .le. 0) stop "gde atomy? kolicxestvo atomov slisxkom malo."

    cutoff_ideal=cutoff

    do i=7000,13000
        !vybratq cutoff v cikle
        cutoff=cutoff_ideal*1d-4*i

        !zapolnitq massiv an_by_cn

    enddo

    cutoff=cutoff_ideal
endsubroutine test_cell_partition
