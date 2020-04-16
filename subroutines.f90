subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use phys_parameters_mod
    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q

    a_half  = a0*5d-1
    cnt_q = 0
    do cnt_x=-x_layers,x_layers
        if (cnt_q + 3 .ge. atoms_max_array) stop "massiv perepolnen : slisxkom mnogo atomov"
        do cnt_y=-y_layers,y_layers
            do cnt_z=-z_layers,z_layers
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0
                R_curr(2,cnt_q) = cnt_y * a0
                R_curr(3,cnt_q) = cnt_z * a0
                if(cnt_x.eq.x_layers) cycle
                if(cnt_y.eq.y_layers) cycle
                if(cnt_z.eq.z_layers) cycle
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0 + a_half
                R_curr(2,cnt_q) = cnt_y * a0 + a_half
                R_curr(3,cnt_q) = cnt_z * a0 + a_half
            enddo
        enddo
    enddo
    atoms__in_total = cnt_q
    R_perf=R_curr
    write(*,107) " Atomy OCK resxetki razmehxeny. Vsego atomov: ", cnt_q," ; "
    write(*,106) " razmery jacxejki (slojqov) :",(x_layers*2+1),"na",(y_layers*2+1),"na",(z_layers*2+1),&
        ";  parametr resxetki a_0 = ", a0, " angstrem "
106 format  (3(A,1x,I3.3,1x),A,F8.4,A)
107 format  (A,I10.3,A,$)
endsubroutine   spawn_bcc_rectangular_100

subroutine      wo_xyz_rhcells !#196
    use positions_mod
    use wo_counters_mod
    use interaction_mod
    USE phys_parameters_mod, ONLY:A0
    integer i_atoms,S,w
    character (LEN = 25) filename_trajectory
    xyz__WOcounter=xyz__WOcounter+1
    write(filename_trajectory,'(A,I5.5,A)'),"rhdch_cell_test_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)
        S=0
        DO I_ATOMS=FIRST_RELAXABLE,LAST__RELAXABLE
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            S=S+1
        ENDDO
        W=0
        DO I_ATOMS=FIRST______WALL,LAST_______WALL
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            W=W+1
        ENDDO

        WRITE(104,*),S+W
        WRITE(104,*),"IMAGE OF RELAXATION # ",XYZ__WOCOUNTER," RHOMBIC DODECAHEDRON EDGE IS ", cutoff
        DO I_ATOMS=FIRST_RELAXABLE,LAST__RELAXABLE
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            WRITE(104,196)  "Fe ", R_CURR(1:3,I_ATOMS)
        ENDDO
        DO I_ATOMS=FIRST______WALL,LAST_______WALL
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            WRITE(104,196)  " W ", R_CURR(1:3,I_ATOMS)
        ENDDO

    if(mod(xyz__WOcounter,25).eq.0)&
        print'(A,I5.5,A)'," Writing out XYZ picture # ",xyz__WOcounter," . "

    close(104)
196 format(A2,3(3x,ES11.4))
endsubroutine   wo_xyz_rhcells

subroutine fill_an_to_cn

endsubroutine fill_an_to_cn

subroutine fill_cn_to_an

endsubroutine fill_cn_to_an

