subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use phys_parameters_mod
    implicit none
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
                R_curr(1,cnt_q) = cnt_x * a0 + aepsilon*2d0
                R_curr(2,cnt_q) = cnt_y * a0 + aepsilon*3d0
                R_curr(3,cnt_q) = cnt_z * a0 + aepsilon*5d0
                if(cnt_x.eq.x_layers) cycle
                if(cnt_y.eq.y_layers) cycle
                if(cnt_z.eq.z_layers) cycle
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0 + a_half + aepsilon*2d0
                R_curr(2,cnt_q) = cnt_y * a0 + a_half + aepsilon*3d0
                R_curr(3,cnt_q) = cnt_z * a0 + a_half + aepsilon*5d0
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
    implicit none
    integer i_atoms,S,w,cn
    character (LEN = 25) filename_trajectory
    xyz__WOcounter=xyz__WOcounter+1
    write(filename_trajectory,'(A,I5.5,A)'),"rhdch_cell_test_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)

    write(104,*),atoms__in_total
    write(104,*),"image of relaxation # ",xyz__wocounter," rhombic dodecahedron edge is ", cutoff_param

    do i_atoms=1,atoms__in_total
        cn = 1 + mod( sign(1,(an_by_cn(1,i_atoms) + an_by_cn(2,i_atoms) + an_by_cn(3,i_atoms))) , 12)
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
    implicit none
    integer i_atoms,small,cn_raw(3),cn_add(3)
    real(8) scale_factor,r_at(3),r_rel(3),r_rel_abs(3),cellpos(3)
    real(8) hc,hc_reci

    logical(1) :: xp(6),xn(6),yp(6),yn(6),zp(6),zn(6),our(6)
    !      x+y>0  ; -x+y>0  ; y+z>0  ; -y+z>0  ; z+x>0  ; -z+x>0
    xp=(/ .true.  , .false. , .true.  , .true.  , .true.  , .true.  /)
    xn=(/ .false. , .true.  , .true.  , .true.  , .false. , .false. /)

    yp=(/ .true.  , .true.  , .true.  , .false. , .true.  , .true.  /)
    yn=(/ .false. , .false. , .false. , .true.  , .true.  , .true.  /)

    zp=(/ .true.  , .true.  , .true.  , .true.  , .true.  , .false. /)
    zn=(/ .true.  , .true.  , .false. , .false. , .false. , .true.  /)

    print"(A,F9.5,A,3I4.2)","Zapolnqajem massiv an_by_cn. Radius obrezanija = ",cutoff_param," angstrem"

    hc=sqrt(16d0/27d0)*(cutoff_param)
    hc_reci=1d0/hc

    do i_atoms=1,atoms__in_total

        r_at(1:3)=R_curr(1:3,i_atoms)

        cn_raw(1:3)=nint( hc_reci*3d0*(r_at(1:3)) ) !delim na menqsxije kubiki

        cellpos(1:3)=(hc/3d0)*cn_raw(1:3)
        if ( mod( sum(cn_raw) ,2 ) .eq. 0 ) then
            !eto kubik vnutri jacxejki, nuzxno razlozxitq koordinaty po bazisu

            an_by_cn(1,i_atoms) =(-cn_raw(2)-cn_raw(3)+cn_raw(1) )/2
            an_by_cn(2,i_atoms) =(-cn_raw(3)-cn_raw(1)+cn_raw(2) )/2
            an_by_cn(3,i_atoms) =(-cn_raw(1)-cn_raw(2)+cn_raw(3) )/2

            cycle
        else
            !eto kubik dlqa razdelenija. nuzxno scxitatq logiku
            r_rel  =  r_at - cellpos

            our(1) =  r_rel(1) .gt. -r_rel(2)
            our(2) =  r_rel(2) .gt.  r_rel(1)
            our(3) =  r_rel(2) .gt. -r_rel(3)
            our(4) =  r_rel(3) .gt.  r_rel(2)
            our(5) =  r_rel(3) .gt. -r_rel(1)
            our(6) =  r_rel(1) .gt.  r_rel(3)
!                    ( our(1) .eqv. xp(1) ).and. &
!                    ( our(2) .eqv. xp(2) ).and. &
!                    ( our(3) .eqv. xp(3) ).and. &
!                    ( our(4) .eqv. xp(4) ).and. &
!                    ( our(5) .eqv. xp(5) ).and. &
!                    ( our(6) .eqv. xp(6) ) &
            if (&
                    ( our(1) .eqv. xp(1) ).and. &
                    ( our(2) .eqv. xp(2) ).and. &
                    ( our(5) .eqv. xp(5) ).and. &
                    ( our(6) .eqv. xp(6) ) &
                ) then
                cn_raw(1)=cn_raw(1)+1
            elseif (&
                    ( our(1) .eqv. xn(1) ).and. &
                    ( our(2) .eqv. xn(2) ).and. &
                    ( our(5) .eqv. xn(5) ).and. &
                    ( our(6) .eqv. xn(6) ) &
                ) then
                cn_raw(1)=cn_raw(1)-1
            elseif (&
                    ( our(1) .eqv. yp(1) ).and. &
                    ( our(2) .eqv. yp(2) ).and. &
                    ( our(3) .eqv. yp(3) ).and. &
                    ( our(4) .eqv. yp(4) ) &
                ) then
                cn_raw(2)=cn_raw(2)+1
            elseif (&
                    ( our(1) .eqv. yn(1) ).and. &
                    ( our(2) .eqv. yn(2) ).and. &
                    ( our(3) .eqv. yn(3) ).and. &
                    ( our(4) .eqv. yn(4) ) &
                ) then
                cn_raw(2)=cn_raw(2)-1
            elseif (&
                    ( our(3) .eqv. zp(3) ).and. &
                    ( our(4) .eqv. zp(4) ).and. &
                    ( our(5) .eqv. zp(5) ).and. &
                    ( our(6) .eqv. zp(6) ) &
                ) then
                cn_raw(3)=cn_raw(3)+1
            elseif (&
                    ( our(3) .eqv. zn(3) ).and. &
                    ( our(4) .eqv. zn(4) ).and. &
                    ( our(5) .eqv. zn(5) ).and. &
                    ( our(6) .eqv. zn(6) ) &
                ) then
                cn_raw(3)=cn_raw(3)-1
            else
                STOP "tablica ne dolzxna proverqatsqa dlqa nerelevantnoj pary koordinat"
                print*,"Tablica istinnosti ne sovpala ni s kem. Atom ",i_atoms
                print "(1x,6L2)", our

            endif
            an_by_cn(1,i_atoms) =(-cn_raw(2)-cn_raw(3)+cn_raw(1) )/2
            an_by_cn(2,i_atoms) =(-cn_raw(3)-cn_raw(1)+cn_raw(2) )/2
            an_by_cn(3,i_atoms) =(-cn_raw(1)-cn_raw(2)+cn_raw(3) )/2
        endif


!
!        cellpos(1) = ( 1d0 / scale_factor ) * (cn_raw(2) + cn_raw(3))
!        cellpos(2) = ( 1d0 / scale_factor ) * (cn_raw(3) + cn_raw(1))
!        cellpos(3) = ( 1d0 / scale_factor ) * (cn_raw(1) + cn_raw(2))
!
!        r_rel(1:3)=r_at(1:3)-cellpos(1:3)
!        r_rel_abs(1:3) = abs(r_rel(1:3))
!        small=minloc(r_rel_abs,DIM=1)
!        r_rel_abs(small)=0d0
!        if(&
!            cn_raw(1).eq.0 .and. &
!            cn_raw(2).eq.0 .and. &
!            cn_raw(3).eq.1 ) then
!!            print 263,"Atom nomer ",i_atoms," nahoditsqa v ", r_at ," vozle ",cellpos,&
!!                " v FCC predstavlenii [[",cn_raw,"]].",&
!!                " Otnositelqnyje koordinaty ",r_rel," po modulqu ",r_rel_abs,"; rasstojanije ",norm2(r_rel)
!
!        endif
!        if(sum(r_rel_abs).lt.insphere_radius*sqrt(2d0)) cycle !mb sqrt(2)
!
!        cn_add(1)=nint( sign( 1d0,r_rel(1) ) )
!        cn_add(2)=nint( sign( 1d0,r_rel(2) ) )
!        cn_add(3)=nint( sign( 1d0,r_rel(3) ) )
!        cn_add(small)=0
!
!        an_by_cn(1:3,i_atoms)=cn_raw(1:3)!+cn_add(1:3)

    enddo

!260 format (A,I8.4,A,3(F9.4,1x),A,3(F8.3,1x),A,3(I5.4,1x),A,A,3(F8.3,1x),A)
!261 format (A,1x,I4.2,$)
!262 format (A,1x,I4.2)
263 format (/,A,I8.4,A,3(F9.4,1x),A,3(F8.3,1x),/,A,3(I5.4,1x),A,/,A,3(F8.3,1x),A,3(F8.3,1x),A,F8.4)

endsubroutine fill_an_by_cn

subroutine fill_cn_by_an

endsubroutine fill_cn_by_an

