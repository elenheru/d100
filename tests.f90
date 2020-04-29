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

    do i=8000,145000,250

        !vybratq cutoff v cikle
        cutoff_param=cutoff*25d-5*i

        !zapolnitq massiv an_by_cn
        call fill_an_by_cn_4d

        call wo_xyz_rhcells_4d
        !STOP "filling an by cn test"
    enddo


endsubroutine test_cell_partition

subroutine test_cell_gap
    use array_parameters_mod
    use interaction_mod
    implicit none
    integer i
    !logical test_if_cell_is_close,lll
    !INTEGER, DIMENSION(3,42) :: LIST
    print*, "TEST_CELL_GAP zapuhxen"

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

    do i=9000,195000,500

        !vybratq cutoff v cikle
        cutoff_param=cutoff*1d-4*i

        !zapolnitq massiv an_by_cn
        call fill_an_by_cn
        call fill_cn_by_an
        call check_escaped_atom
        call wo_xyz_rhcells_4d
        !STOP "filling an by cn test"
    enddo

endsubroutine test_cell_gap

subroutine po_distances_fcc_vs_pc
    integer :: ixy,iyz,izx,w
    real :: dist
    w=1
    do ixy=-w,w
        do iyz=-w,w
            do izx=-w,w
    !print "(3(A,I2.1),A,$)", " ixy=",ixy," ;iyz=",iyz,";izx=",izx," "
    print "(3(A,I2.1),A,$)", "[[",ixy,",",iyz,",",izx,"]]_FCC"
    print "(3(A,I2.1),A,$)", "  ;[[",ixy+izx,",",iyz+ixy,",",izx+iyz,"]]_PC"
    dist=sqrt(real((ixy+izx)**2+(iyz+ixy)**2+(izx+iyz)**2))
    print "(A,F8.4)", "   ; Distance is " // repeat(" ",nint(6*dist)) ,dist
            enddo
        enddo
    enddo
endsubroutine po_distances_fcc_vs_pc
