
subroutine      wo_xyz_rhcells_4d !#196
    use positions_mod
    use wo_counters_mod
    use interaction_mod
    use chemical_elements_names_mod
    USE phys_parameters_mod, ONLY:A0
    implicit none
    !logical test_if_cell_is_close_4d
    integer i_atoms,cn
    character (LEN = 25) filename_trajectory
    xyz__WOcounter=xyz__WOcounter+1
    !stop "wo xyz_rhcells ustarela, ispolqzujetsqa 4d versija "
    write(filename_trajectory,'(A,I5.5,A)') "rhdch_cell_test_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)

    write(104,*) atoms__in_total
    write(104,*) "image of relaxation # ",xyz__wocounter," rhombic dodecahedron edge is ", cutoff_param

    do i_atoms=1,atoms__in_total
        !cn = 1 + mod( abs((an_by_cn_4d(1,i_atoms) + 5*an_by_cn_4d(2,i_atoms) + 13*an_by_cn_4d(3,i_atoms))) , 52)
        cn = an_by_cn_4d(4,i_atoms) + 1
!        if (.not.test_if_cell_is_close_4d(&
!                an_by_cn_4d(1,i_atoms),&
!                an_by_cn_4d(2,i_atoms),&
!                an_by_cn_4d(3,i_atoms),&
!                an_by_cn_4d(4,i_atoms)  )) cn=7
!        if ( test_if_cell_is_close_4d(&
!             an_by_cn_4d(1,i_atoms),&
!             an_by_cn_4d(2,i_atoms),&
!             an_by_cn_4d(3,i_atoms),&
!             an_by_cn_4d(4,i_atoms)  )) cn=1
!        if ( ( an_by_cn_4d(1,i_atoms) .eq. 0) .and. &
!             ( an_by_cn_4d(2,i_atoms) .eq. 0) .and. &
!             ( an_by_cn_4d(3,i_atoms) .eq. 0) .and. &
!             ( an_by_cn_4d(4,i_atoms) .eq. 0) ) cn=8

        write(104,196) elements_names(2*cn-1:2*cn), r_curr(1:3,i_atoms)
    enddo


    if(mod(xyz__WOcounter,25).eq.0)&
        print'(A,I5.5,A)'," Writing out XYZ picture # ",xyz__WOcounter," . "

    close(104)
196 format(A2,3(2x,ES11.4))
endsubroutine   wo_xyz_rhcells_4d



subroutine po_distances_fcc_vs_pc
    implicit none
    integer :: ixy,iyz,izx,w
    real :: dist
    stop " testovaja subroutine po_distances_fcc_vs_pc ne dolzxna zapuskatqsqa, eto neozxidanno"
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

subroutine po_wo_parameters
    use array_parameters_mod
    use interaction_mod
    use phys_parameters_mod
    implicit none
    open(1090, file = "parameters.txt")

    print*," parametry takovy (vyvedeno na ekran i v fajl parameters.txt) : "
    print*,"x_layers = ",x_layers
    print*,"y_layers = ",y_layers
    print*,"z_layers = ",z_layers
    print*,"atoms_max_array = ",atoms_max_array
    print*,"cells_xrange = ",cells_xrange
    print*,"cells_yrange = ",cells_yrange
    print*,"cells_zrange = ",cells_zrange

    print*,"cutoff = ",cutoff
    print*,"incell_atoms_exceed = ",incell_atoms_exceed
    print*,"incell_atoms_max = ",incell_atoms_max
    print*,"a0 = ",a0
    print*,"aepsilon = ",aepsilon
    print*,"poisson = ",poisson
    print*,"burgers = ",burgers



    write(1090,*) " parametry takovy (vyvedeno na ekran i v fajl parameters.txt) : "
    write(1090,*) "x_layers = ",x_layers
    write(1090,*) "y_layers = ",y_layers
    write(1090,*) "z_layers = ",z_layers
    write(1090,*) "atoms_max_array = ",atoms_max_array
    write(1090,*) "cells_xrange = ",cells_xrange
    write(1090,*) "cells_yrange = ",cells_yrange
    write(1090,*) "cells_zrange = ",cells_zrange

    write(1090,*) "cutoff = ",cutoff
    write(1090,*) "incell_atoms_exceed = ",incell_atoms_exceed
    write(1090,*) "incell_atoms_max = ",incell_atoms_max
    write(1090,*) "a0 = ",a0
    write(1090,*) "aepsilon = ",aepsilon
    write(1090,*) "poisson = ",poisson
    write(1090,*) "burgers = ",burgers

!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
!    print*,
    close(1090)
endsubroutine po_wo_parameters


