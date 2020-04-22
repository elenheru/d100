subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use phys_parameters_mod

    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q

    print*,cells_range," is cells range"

    a_half  = a0*5d-1
    cnt_q = 0

    do cnt_x=-x_layers,x_layers
        do cnt_y=-y_layers,y_layers
            do cnt_z=-z_layers,z_layers
                if (cnt_q + 1 .ge. atoms_max_array) stop "massiv perepolnen : slisxkom mnogo atomov"
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0 + aepsilon
                R_curr(2,cnt_q) = cnt_y * a0 + aepsilon
                R_curr(3,cnt_q) = cnt_z * a0 + aepsilon
                if(cnt_x.eq.x_layers) cycle
                if(cnt_y.eq.y_layers) cycle
                if(cnt_z.eq.z_layers) cycle
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0 + a_half + aepsilon
                R_curr(2,cnt_q) = cnt_y * a0 + a_half + aepsilon
                R_curr(3,cnt_q) = cnt_z * a0 + a_half + aepsilon
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
    use chemical_elements_names_mod
    USE phys_parameters_mod, ONLY:A0
    integer i_atoms,S,w,cn
    character (LEN = 25) filename_trajectory
    xyz__WOcounter=xyz__WOcounter+1
    write(filename_trajectory,'(A,I5.5,A)'),"rhdch_cell_test_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)

    write(104,*),atoms__in_total
    write(104,*),"image of relaxation # ",xyz__wocounter," rhombic dodecahedron edge is ", cutoff_param

    do i_atoms=1,atoms__in_total
        cn = 1 + mod( (an_by_cn(1,i_atoms) + 10*an_by_cn(2,i_atoms) + 100*an_by_cn(3,i_atoms))**2 , 25)

        write(104,196), elements_names(2*cn-1:2*cn), r_curr(1:3,i_atoms)
    enddo


    if(mod(xyz__WOcounter,25).eq.0)&
        print'(A,I5.5,A)'," Writing out XYZ picture # ",xyz__WOcounter," . "

    close(104)
196 format(A2,3(2x,ES11.4))
endsubroutine   wo_xyz_rhcells

subroutine fill_an_by_cn
    use positions_mod
    use interaction_mod
    use phys_parameters_mod

    integer i_atoms,s,cn_raw(3),cell_test_fcc(3)
    integer poscell(1:2,1:2,1:3)
    !alpha_z1 is poscell(1,1,1)  alpha_z2 is poscell(1,2,1)
    !beta__z1 is poscell(2,1,1)  beta__z2 is poscell(2,2,1)

    real(8) insphere_diameter,scale_factor,r_at(3),r_rel(3),r_neib(3),cellpos(3)

    cell_test_fcc(1)= 3
    cell_test_fcc(2)= 0
    cell_test_fcc(3)= 0
!    cell_test_fcc(1)= 1
!    cell_test_fcc(2)=-2
!    cell_test_fcc(3)= 3

    scale_factor=5d-1*sqrt(5d-1)/(sqrt(2d0/3d0)*cutoff_param)

    print"(A,F9.5,A,3I4.2)","Zapolnqajem massiv an_by_cn. Radius obrezanija = ",cutoff_param," angstrem",cell_test_fcc

    DO I_ATOMS=1,ATOMS__IN_TOTAL

        r_at(1:3)=R_curr(1:3,i_atoms)

        poscell(1,1,1)=nint(scale_factor*( r_at(1) +r_at(2) ) )
        poscell(2,1,1)=nint(scale_factor*(-r_at(1) +r_at(2) ) )

        poscell(1,2,1)=nint(scale_factor*( r_at(1) +r_at(2) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,1)=nint(scale_factor*(-r_at(1) +r_at(2) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )


        poscell(1,1,2)=nint(scale_factor*( r_at(2) +r_at(3) ) )
        poscell(2,1,2)=nint(scale_factor*(-r_at(2) +r_at(3) ) )

        poscell(1,2,2)=nint(scale_factor*( r_at(2) +r_at(3) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,2)=nint(scale_factor*(-r_at(2) +r_at(3) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )


        poscell(1,1,3)=nint(scale_factor*( r_at(3) +r_at(1) ) )
        poscell(2,1,3)=nint(scale_factor*(-r_at(3) +r_at(1) ) )

        poscell(1,2,3)=nint(scale_factor*( r_at(3) +r_at(1) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,3)=nint(scale_factor*(-r_at(3) +r_at(1) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )

        !CELLPOS(1)= 0D0 ; CELLPOS(2)=-SQRT(2D0)*(SQRT(2D0/3D0)*CUTOFF_PARAM) ; CELLPOS(3)= SQRT(2D0)*(SQRT(2D0/3D0)*CUTOFF_PARAM)

        CELLPOS(1)= (SQRT(2D0)*(SQRT(2D0/3D0)*CUTOFF_PARAM))*(cell_test_fcc(2)+cell_test_fcc(3))
        CELLPOS(2)= (SQRT(2D0)*(SQRT(2D0/3D0)*CUTOFF_PARAM))*(cell_test_fcc(3)+cell_test_fcc(1))
        CELLPOS(3)= (SQRT(2D0)*(SQRT(2D0/3D0)*CUTOFF_PARAM))*(cell_test_fcc(1)+cell_test_fcc(2))

        IF( (NORM2(R_AT(1:3)-CELLPOS(1:3)) .LT. (SQRT(2D0/3D0)*CUTOFF_PARAM)) ) THEN

            PRINT*,""
            PRINT 260,"Atom nomer ",I_ATOMS, " nahoditsqa v tocxke ",r_at," vozle ",CELLPOS,&
                ", [[",cell_test_fcc,"]]"," v GCK predstavlenii. Imejet takije koordinaty :"

            PRINT 261,"   Alpha Z1 =",poscell(1,1,1)
            PRINT 261," ; Beta  Z1 =",poscell(2,1,1)
            PRINT 261," ; Alpha Z2 =",poscell(1,2,1)
            PRINT 262," ; Beta  Z2 =",poscell(2,2,1)

            PRINT 261,"   Alpha X1 =",poscell(1,1,2)
            PRINT 261," ; Beta  X1 =",poscell(2,1,2)
            PRINT 261," ; Alpha X2 =",poscell(1,2,2)
            PRINT 262," ; Beta  X2 =",poscell(2,2,2)

            PRINT 261,"   Alpha Y1 =",poscell(1,1,3)
            PRINT 261," ; Beta  Y1 =",poscell(2,1,3)
            PRINT 261," ; Alpha Y2 =",poscell(1,2,3)
            PRINT 262," ; Beta  Y2 =",poscell(2,2,3)

        ENDIF

    enddo

    do i_atoms=1,atoms__in_total

        r_at(1:3)=R_curr(1:3,i_atoms)

!        r_at(1:3)=R_curr(1:3,i_atoms)*scale_factor
!
!        cn_raw(1)=nint( 5d-1*( r_at(1)+r_at(2)-r_at(3)) )
!        cn_raw(2)=nint( 5d-1*( r_at(1)-r_at(2)+r_at(3)) )
!        cn_raw(3)=nint( 5d-1*(-r_at(1)+r_at(2)+r_at(3)) )
!
!        cellpos(1:3)=(1d0/scale_factor)*an_by_cn(1:3,i_atoms)
!        r_rel(1:3)=r_at(1:3)-cellpos(1:3)
!
!

        poscell(1,1,1)=nint(scale_factor*( r_at(1) +r_at(2) ) )
        poscell(2,1,1)=nint(scale_factor*(-r_at(1) +r_at(2) ) )

        poscell(1,2,1)=nint(scale_factor*( r_at(1) +r_at(2) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,1)=nint(scale_factor*(-r_at(1) +r_at(2) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )


        poscell(1,1,2)=nint(scale_factor*( r_at(2) +r_at(3) ) )
        poscell(2,1,2)=nint(scale_factor*(-r_at(2) +r_at(3) ) )

        poscell(1,2,2)=nint(scale_factor*( r_at(2) +r_at(3) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,2)=nint(scale_factor*(-r_at(2) +r_at(3) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )


        poscell(1,1,3)=nint(scale_factor*( r_at(3) +r_at(1) ) )
        poscell(2,1,3)=nint(scale_factor*(-r_at(3) +r_at(1) ) )

        poscell(1,2,3)=nint(scale_factor*( r_at(3) +r_at(1) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )
        poscell(2,2,3)=nint(scale_factor*(-r_at(3) +r_at(1) - sqrt(2d0)*(sqrt(2d0/3d0)*cutoff_param) ) )



        !an_by_cn(1,i_atoms)=cn_raw(1)
!
!        IF( (norm2(r_at(1:3)) .GT. 0d0) .AND. (norm2(r_at(1:3)) .LT. 3d0) ) then
!
!            PRINT*,""
!            PRINT 250,"Atom nomer ",I_ATOMS, " nahoditsqa v tocxke ",r_at," i imejet takije koordinaty :"
!
!            PRINT 251,"   Alpha Z1 =",poscell(1,1,1)
!            PRINT 251," ; Beta  Z1 =",poscell(2,1,1)
!            PRINT 251," ; Alpha Z2 =",poscell(1,2,1)
!            PRINT 252," ; Beta  Z2 =",poscell(2,2,1)
!
!            PRINT 251,"   Alpha X1 =",poscell(1,1,2)
!            PRINT 251," ; Beta  X1 =",poscell(2,1,2)
!            PRINT 251," ; Alpha X2 =",poscell(1,2,2)
!            PRINT 252," ; Beta  X2 =",poscell(2,2,2)
!
!            PRINT 251,"   Alpha Y1 =",poscell(1,1,3)
!            PRINT 251," ; Beta  Y1 =",poscell(2,1,3)
!            PRINT 251," ; Alpha Y2 =",poscell(1,2,3)
!            PRINT 252," ; Beta  Y2 =",poscell(2,2,3)
!
!        ENDIF

        CELLPOS(1)=0D0 ; CELLPOS(2)=0D0 ; CELLPOS(3)=0D0

        IF( (NORM2(R_AT(1:3)-CELLPOS(1:3)) .LT. (SQRT(2D0/3D0)*CUTOFF_PARAM)) ) THEN

            PRINT*,""
            PRINT 260,"Atom nomer ",I_ATOMS, " nahoditsqa v tocxke ",r_at," vozle ",CELLPOS,&
                ", [[",0,0,0,"]]"," v GCK predstavlenii. Imejet takije koordinaty :"

            PRINT 261,"   Alpha Z1 =",poscell(1,1,1)
            PRINT 261," ; Beta  Z1 =",poscell(2,1,1)
            PRINT 261," ; Alpha Z2 =",poscell(1,2,1)
            PRINT 262," ; Beta  Z2 =",poscell(2,2,1)

            PRINT 261,"   Alpha X1 =",poscell(1,1,2)
            PRINT 261," ; Beta  X1 =",poscell(2,1,2)
            PRINT 261," ; Alpha X2 =",poscell(1,2,2)
            PRINT 262," ; Beta  X2 =",poscell(2,2,2)

            PRINT 261,"   Alpha Y1 =",poscell(1,1,3)
            PRINT 261," ; Beta  Y1 =",poscell(2,1,3)
            PRINT 261," ; Alpha Y2 =",poscell(1,2,3)
            PRINT 262," ; Beta  Y2 =",poscell(2,2,3)

        ENDIF

    enddo

250 format (A,I8.4,A,3(F9.4,1x),A)
251 format (A,1x,I4.2,$)
252 format (A,1x,I4.2)

260 format (A,I8.4,A,3(F9.4,1x),A,3(F8.3,1x),A,3(I5.4,1x),A,A,3(F8.3,1x),A)
261 format (A,1x,I4.2,$)
262 format (A,1x,I4.2)


!vycxislitq v kakom rombe otnositelqno

endsubroutine fill_an_by_cn

subroutine fill_cn_by_an

endsubroutine fill_cn_by_an

