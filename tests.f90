subroutine test_cell_partition
    use array_parameters_mod
    use interaction_mod

    integer i

    print*, "TEST_CELL_PARTITION zapuhxen"

    !real(8) R_curr(1,cnt_q)
    !integer an_by_cn(3,atoms_max_array)
    !integer cn_by_an(0:neibors_max,atoms_max_array)

    !proveritq zadany li koordinaty atomov
    if (atoms__in_total .le. 0) stop "gde atomy? kolicxestvo atomov slisxkom malo."



    do i=19000,195000,500

        !vybratq cutoff v cikle
        cutoff_param=cutoff*1d-4*i

        !zapolnitq massiv an_by_cn
        call fill_an_by_cn

        call wo_xyz_rhcells
        !STOP "filling an by cn test"
    enddo


endsubroutine test_cell_partition
