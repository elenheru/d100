subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use phys_parameters_mod
    implicit none
    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q

    print*,cells_xrange,cells_yrange,cells_zrange," eto x y z cells range"

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

subroutine      wo_xyz_rhcells_4d !#196
    use positions_mod
    use wo_counters_mod
    use interaction_mod
    use chemical_elements_names_mod
    USE phys_parameters_mod, ONLY:A0
    implicit none
    logical test_if_cell_is_close_4d
    integer i_atoms,s,w,cn
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

subroutine fill_an_by_cn_4d
    use positions_mod
    use interaction_mod
    use phys_parameters_mod
    implicit none
    integer i_atoms,small,cn_raw(3),cn_add(3)
    real(8) scale_factor,r_at(3),r_rel(3),r_rel_abs(3),cellpos(3)
    real(8) hc,hc_reci

    logical(1) :: our(6)
    logical(1), parameter, dimension(6) :: &
    xp=(/ .true.  , .false. , .true.  , .true.  , .true.  , .true.  /),&
    xn=(/ .false. , .true.  , .true.  , .true.  , .false. , .false. /),&
    yp=(/ .true.  , .true.  , .true.  , .false. , .true.  , .true.  /),&
    yn=(/ .false. , .false. , .false. , .true.  , .true.  , .true.  /),&
    zp=(/ .true.  , .true.  , .true.  , .true.  , .true.  , .false. /),&
    zn=(/ .true.  , .true.  , .false. , .false. , .false. , .true.  /)
    !      x+y>0  ; -x+y>0  ;  y+z>0  ; -y+z>0  ;  z+x>0  ; -z+x>0

    print"(A,F9.5,A,3I4.2)","Zapolnqajem 4d massiv an_by_cn. Radius obrezanija = ",cutoff_param," angstrem"

    ! mozxno perepakovatq massiv anbycn tak cxtoby on byl 4mernym.
    ! togda ne ponadobitsqa izbytok razmera massiva iz-za afinnosti FCC bazisa

!    hc=sqrt(16d0/27d0)*(cutoff_param) !jesli granq ravna cutoff
    !hc=sqrt(2d0/9d0)*(cutoff_param) !jesli diametr vpisannoj sfery raven cutoff, 3h
    hc=sqrt(4d0/3d0)*(cutoff_param) !jesli diametr vpisannoj sfery raven cutoff, 2h
    hc_reci=1d0/hc

    do i_atoms=1,atoms__in_total

        r_at(1:3)=R_curr(1:3,i_atoms)

        !cn_raw(1:3)=nint( hc_reci*3d0*(r_at(1:3)) ) !delim na menqsxije kubiki
        cn_raw(1:3)=nint( hc_reci*2d0*(r_at(1:3)) ) !delim na menqsxije kubiki

        !cellpos(1:3)=(hc/3d0)*cn_raw(1:3) ! ne 3h a 2h
        cellpos(1:3)=(hc/2d0)*cn_raw(1:3)
        if ( mod( sum(cn_raw) ,2 ) .eq. 0 ) then
            !eto kubik vnutri jacxejki, nuzxno razlozxitq koordinaty po bazisu

            !fcc basis
!            an_by_cn(1,i_atoms) =-(-cn_raw(2)-cn_raw(3)+cn_raw(1) )/2
!            an_by_cn(2,i_atoms) =-(-cn_raw(3)-cn_raw(1)+cn_raw(2) )/2
!            an_by_cn(3,i_atoms) =-(-cn_raw(1)-cn_raw(2)+cn_raw(3) )/2
            !fcc basis

            !pc with additions basis
            an_by_cn_4d(1,i_atoms) = ( cn_raw(1)-modulo( cn_raw(1),2 ) )/2
            an_by_cn_4d(2,i_atoms) = ( cn_raw(2)-modulo( cn_raw(2),2 ) )/2
            an_by_cn_4d(3,i_atoms) = ( cn_raw(3)-modulo( cn_raw(3),2 ) )/2
            if (     modulo(cn_raw(1),2) .eq. 1 .and. modulo(cn_raw(2),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 1
            elseif ( modulo(cn_raw(2),2) .eq. 1 .and. modulo(cn_raw(3),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 2
            elseif ( modulo(cn_raw(3),2) .eq. 1 .and. modulo(cn_raw(1),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 3
            else
                an_by_cn_4d(4,i_atoms) = 0
            endif
            !pc with additions basis
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
                print*,"Tablica istinnosti ne sovpala ni s kem. Atom ",i_atoms,r_at(1:3)
                print "(1x,6L2)", our
                STOP "tablica ne dolzxna proverqatsqa dlqa nerelevantnoj pary koordinat"
            endif
            !fcc basis
!            an_by_cn(1,i_atoms) =(-cn_raw(2)-cn_raw(3)+cn_raw(1) )/2
!            an_by_cn(2,i_atoms) =(-cn_raw(3)-cn_raw(1)+cn_raw(2) )/2
!            an_by_cn(3,i_atoms) =(-cn_raw(1)-cn_raw(2)+cn_raw(3) )/2
            !fcc basis

            !pc with additions basis
            an_by_cn_4d(1,i_atoms) = ( cn_raw(1)-modulo( cn_raw(1),2 ) )/2
            an_by_cn_4d(2,i_atoms) = ( cn_raw(2)-modulo( cn_raw(2),2 ) )/2
            an_by_cn_4d(3,i_atoms) = ( cn_raw(3)-modulo( cn_raw(3),2 ) )/2
            if (     modulo(cn_raw(1),2) .eq. 1 .and. modulo(cn_raw(2),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 1
            elseif ( modulo(cn_raw(2),2) .eq. 1 .and. modulo(cn_raw(3),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 2
            elseif ( modulo(cn_raw(3),2) .eq. 1 .and. modulo(cn_raw(1),2) .eq. 1 ) then
                an_by_cn_4d(4,i_atoms) = 3
            else
                an_by_cn_4d(4,i_atoms) = 0
            endif
            !pc with additions basis
        endif

    enddo

!263 format (/,A,I8.4,A,3(F9.4,1x),A,3(F8.3,1x),/,A,3(I5.4,1x),A,/,A,3(F8.3,1x),A,3(F8.3,1x),A,F8.4)

endsubroutine fill_an_by_cn_4d

subroutine fill_cn_by_an_4d
    use positions_mod
    use interaction_mod
    use phys_parameters_mod

    implicit none
    integer i_atoms,na
    integer anxy,anyz,anzx,anic

    print"(A,F9.5,A,3I4.2)","Zapolnqajem massiv cn_by_an. Radius obrezanija = ",cutoff_param," angstrem"
    cn_by_an_4d=0
    na=0

    do i_atoms=1,atoms__in_total

        if(na+1 .ge. incell_atoms_max) stop "massiv cn_by_an perepolnen"
        anxy=an_by_cn_4d(1,i_atoms)
        anyz=an_by_cn_4d(2,i_atoms)
        anzx=an_by_cn_4d(3,i_atoms)
        anic=an_by_cn_4d(4,i_atoms)

        na = cn_by_an_4d( 0,anxy,anyz,anzx,anic)
            !PRINT*, "NA, IATOMS", NA,I_ATOMS," AN_BY_CN = ",AN_BY_CN(1:3,I_ATOMS)
            !CALL SLEEP(1)
        cn_by_an_4d( 0,anxy,anyz,anzx,anic ) = na + 1
        cn_by_an_4d( na+1,anxy,anyz,anzx,anic ) = i_atoms
            !PRINT*, "CN_BY_AN (0)",CN_BY_AN( 0,ANXY,ANYZ,ANZX ), "CN_BY_AN (NA)",CN_BY_AN( NA,ANXY,ANYZ,ANZX ),I_ATOMS
            !PRINT*,""
    enddo

    print*,"Zapolnenije cn_by_an zaversxeno. Vsego atomov ", sum(cn_by_an_4d(0,:,:,:,:)), " = ",atoms__in_total

endsubroutine fill_cn_by_an_4d

subroutine list_of_42_nearest_cells_4dpc(ixy,iyz,izx,icell,list)
    implicit none
    integer, intent(in) :: ixy,iyz,izx,icell
    integer, intent(out),dimension(4,42) :: list
    integer, dimension(4) :: origin
    interface
        function fcc_to_pc4d(fcc_cell) result(pc4d_cell)
        implicit none
        integer, dimension(3) :: fcc_cell,cn_raw
        integer, dimension(4) :: pc4d_cell
        endfunction fcc_to_pc4d
    endinterface
    origin(1)=ixy
    origin(2)=iyz
    origin(3)=izx
    origin(4)=0!icell

    if(icell .eq. 0) then
        list(1:4, 1)=origin(1:4)+(/  0,  1, -1,  3/)!_FCC  ;[[ 1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4, 2)=origin(1:4)+(/ -1, -1,  0,  3/)!_FCC  ;[[-1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 3)=origin(1:4)+(/ -1, -1,  0,  2/)!_FCC  ;[[-2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4, 4)=origin(1:4)+(/ -1, -1,  1,  1/)!_FCC  ;[[-1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4, 5)=origin(1:4)+(/ -1, -1, -1,  3/)!_FCC  ;[[-1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4, 6)=origin(1:4)+(/  0, -1,  0,  0/)!_FCC  ;[[ 0,-2, 0]]_PC ; Distance is     2.0000
        list(1:4, 7)=origin(1:4)+(/  0, -1,  0,  3/)!_FCC  ;[[ 1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 8)=origin(1:4)+(/ -1, -1, -1,  2/)!_FCC  ;[[-2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4, 9)=origin(1:4)+(/ -1, -1,  0,  1/)!_FCC  ;[[-1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,10)=origin(1:4)+(/  0, -1,  0,  2/)!_FCC  ;[[ 0,-1, 1]]_PC ; Distance is 1.4142
        list(1:4,11)=origin(1:4)+(/  0, -1,  1,  1/)!_FCC  ;[[ 1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4,12)=origin(1:4)+(/ -1,  0,  0,  0/)!_FCC  ;[[-2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,13)=origin(1:4)+(/ -1,  0,  0,  3/)!_FCC  ;[[-1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,14)=origin(1:4)+(/  0,  0,  1,  0/)!_FCC  ;[[ 0, 0, 2]]_PC ; Distance is     2.0000
        list(1:4,15)=origin(1:4)+(/ -1,  0,  0,  2/)!_FCC  ;[[-2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,16)=origin(1:4)+(/ -1,  0,  1,  1/)!_FCC  ;[[-1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,17)=origin(1:4)+(/  0, -1, -1,  3/)!_FCC  ;[[ 1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4,18)=origin(1:4)+(/ -1, -1, -1,  1/)!_FCC  ;[[-1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,19)=origin(1:4)+(/  0, -1, -1,  2/)!_FCC  ;[[ 0,-1,-1]]_PC ; Distance is 1.4142
        list(1:4,20)=origin(1:4)+(/  0, -1,  0,  1/)!_FCC  ;[[ 1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,21)=origin(1:4)+(/  1, -1,  0,  2/)!_FCC  ;[[ 2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4,22)=origin(1:4)+(/ -1,  0, -1,  3/)!_FCC  ;[[-1, 0,-1]]_PC ; Distance is 1.4142 !list(1:4,1)=origin(1:4)+fcc_to_pc4d( (/ 0, 0, 0/) )!_FCC  ;[[ 0, 0, 0]]_PC ; Distance is0.0000
        list(1:4,23)=origin(1:4)+(/  0,  0,  0,  3/)!_FCC  ;[[ 1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,24)=origin(1:4)+(/ -1,  0, -1,  2/)!_FCC  ;[[-2, 1,-1]]_PC ; Distance is        2.4495
        list(1:4,25)=origin(1:4)+(/ -1,  0,  0,  1/)!_FCC  ;[[-1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,26)=origin(1:4)+(/  0,  0,  0,  2/)!_FCC  ;[[ 0, 1, 1]]_PC ; Distance is 1.4142
        list(1:4,27)=origin(1:4)+(/  0,  0,  1,  1/)!_FCC  ;[[ 1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,28)=origin(1:4)+(/ -1,  1,  0,  3/)!_FCC  ;[[-1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,29)=origin(1:4)+(/  0, -1, -1,  1/)!_FCC  ;[[ 1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,30)=origin(1:4)+(/  1, -1, -1,  2/)!_FCC  ;[[ 2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4,31)=origin(1:4)+(/  0,  0, -1,  0/)!_FCC  ;[[ 0, 0,-2]]_PC ; Distance is     2.0000
        list(1:4,32)=origin(1:4)+(/  0,  0, -1,  3/)!_FCC  ;[[ 1, 0,-1]]_PC ; Distance is 1.4142
        list(1:4,33)=origin(1:4)+(/  1,  0,  0,  0/)!_FCC  ;[[ 2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,34)=origin(1:4)+(/ -1,  0, -1,  1/)!_FCC  ;[[-1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,35)=origin(1:4)+(/  0,  0, -1,  2/)!_FCC  ;[[ 0, 1,-1]]_PC ; Distance is 1.4142
        list(1:4,36)=origin(1:4)+(/  0,  0,  0,  1/)!_FCC  ;[[ 1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,37)=origin(1:4)+(/  1,  0,  0,  2/)!_FCC  ;[[ 2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,38)=origin(1:4)+(/ -1,  1, -1,  3/)!_FCC  ;[[-1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4,39)=origin(1:4)+(/  0,  1,  0,  0/)!_FCC  ;[[ 0, 2, 0]]_PC ; Distance is     2.0000
        list(1:4,40)=origin(1:4)+(/  0,  1,  0,  3/)!_FCC  ;[[ 1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,41)=origin(1:4)+(/  0,  0, -1,  1/)!_FCC  ;[[ 1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,42)=origin(1:4)+(/  1,  0, -1,  2/)!_FCC  ;[[ 2, 1,-1]]_PC ; Distance is        2.4495
    elseif(icell .eq. 1) then
        list(1:4, 1)=origin(1:4)+(/  1,  1, -1,  2/)!_FCC  ;[[ 1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4, 2)=origin(1:4)+(/  0, -1,  0,  2/)!_FCC  ;[[-1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 3)=origin(1:4)+(/ -1,  0,  0,  3/)!_FCC  ;[[-2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4, 4)=origin(1:4)+(/  0,  0,  1,  0/)!_FCC  ;[[-1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4, 5)=origin(1:4)+(/  0, -1, -1,  2/)!_FCC  ;[[-1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4, 6)=origin(1:4)+(/  0, -1,  0,  1/)!_FCC  ;[[ 0,-2, 0]]_PC ; Distance is     2.0000
        list(1:4, 7)=origin(1:4)+(/  1, -1,  0,  2/)!_FCC  ;[[ 1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 8)=origin(1:4)+(/ -1,  0, -1,  3/)!_FCC  ;[[-2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4, 9)=origin(1:4)+(/  0,  0,  0,  0/)!_FCC  ;[[-1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,10)=origin(1:4)+(/  0,  0,  0,  3/)!_FCC  ;[[ 0,-1, 1]]_PC ; Distance is 1.4142
        list(1:4,11)=origin(1:4)+(/  1,  0,  1,  0/)!_FCC  ;[[ 1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4,12)=origin(1:4)+(/ -1,  0,  0,  1/)!_FCC  ;[[-2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,13)=origin(1:4)+(/  0,  0,  0,  2/)!_FCC  ;[[-1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,14)=origin(1:4)+(/  0,  0,  1,  1/)!_FCC  ;[[ 0, 0, 2]]_PC ; Distance is     2.0000
        list(1:4,15)=origin(1:4)+(/ -1,  1,  0,  3/)!_FCC  ;[[-2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,16)=origin(1:4)+(/  0,  1,  1,  0/)!_FCC  ;[[-1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,17)=origin(1:4)+(/  1, -1, -1,  2/)!_FCC  ;[[ 1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4,18)=origin(1:4)+(/  0,  0, -1,  0/)!_FCC  ;[[-1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,19)=origin(1:4)+(/  0,  0, -1,  3/)!_FCC  ;[[ 0,-1,-1]]_PC ; Distance is 1.4142
        list(1:4,20)=origin(1:4)+(/  1,  0,  0,  0/)!_FCC  ;[[ 1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,21)=origin(1:4)+(/  1,  0,  0,  3/)!_FCC  ;[[ 2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4,22)=origin(1:4)+(/  0,  0, -1,  2/)!_FCC  ;[[-1, 0,-1]]_PC ; Distance is 1.4142!list(1:4,1)=origin(1:4)+fcc_to_pc4d( (/ 0, 0, 0/) )!_FCC  ;[[ 0, 0, 0]]_PC ; Distance is0.0000
        list(1:4,23)=origin(1:4)+(/  1,  0,  0,  2/)!_FCC  ;[[ 1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,24)=origin(1:4)+(/ -1,  1, -1,  3/)!_FCC  ;[[-2, 1,-1]]_PC ; Distance is        2.4495
        list(1:4,25)=origin(1:4)+(/  0,  1,  0,  0/)!_FCC  ;[[-1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,26)=origin(1:4)+(/  0,  1,  0,  3/)!_FCC  ;[[ 0, 1, 1]]_PC ; Distance is 1.4142
        list(1:4,27)=origin(1:4)+(/  1,  1,  1,  0/)!_FCC  ;[[ 1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,28)=origin(1:4)+(/  0,  1,  0,  2/)!_FCC  ;[[-1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,29)=origin(1:4)+(/  1,  0, -1,  0/)!_FCC  ;[[ 1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,30)=origin(1:4)+(/  1,  0, -1,  3/)!_FCC  ;[[ 2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4,31)=origin(1:4)+(/  0,  0, -1,  1/)!_FCC  ;[[ 0, 0,-2]]_PC ; Distance is     2.0000
        list(1:4,32)=origin(1:4)+(/  1,  0, -1,  2/)!_FCC  ;[[ 1, 0,-1]]_PC ; Distance is 1.4142
        list(1:4,33)=origin(1:4)+(/  1,  0,  0,  1/)!_FCC  ;[[ 2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,34)=origin(1:4)+(/  0,  1, -1,  0/)!_FCC  ;[[-1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,35)=origin(1:4)+(/  0,  1, -1,  3/)!_FCC  ;[[ 0, 1,-1]]_PC ; Distance is 1.4142
        list(1:4,36)=origin(1:4)+(/  1,  1,  0,  0/)!_FCC  ;[[ 1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,37)=origin(1:4)+(/  1,  1,  0,  3/)!_FCC  ;[[ 2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,38)=origin(1:4)+(/  0,  1, -1,  2/)!_FCC  ;[[-1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4,39)=origin(1:4)+(/  0,  1,  0,  1/)!_FCC  ;[[ 0, 2, 0]]_PC ; Distance is     2.0000
        list(1:4,40)=origin(1:4)+(/  1,  1,  0,  2/)!_FCC  ;[[ 1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,41)=origin(1:4)+(/  1,  1, -1,  0/)!_FCC  ;[[ 1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,42)=origin(1:4)+(/  1,  1, -1,  3/)!_FCC  ;[[ 2, 1,-1]]_PC ; Distance is        2.4495
    elseif(icell .eq. 2) then
        list(1:4, 1)=origin(1:4)+(/  0,  1,  0,  1/)!_FCC  ;[[ 1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4, 2)=origin(1:4)+(/ -1, -1,  1,  1/)!_FCC  ;[[-1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 3)=origin(1:4)+(/ -1,  0,  1,  0/)!_FCC  ;[[-2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4, 4)=origin(1:4)+(/ -1,  0,  1,  3/)!_FCC  ;[[-1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4, 5)=origin(1:4)+(/ -1, -1,  0,  1/)!_FCC  ;[[-1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4, 6)=origin(1:4)+(/  0, -1,  0,  2/)!_FCC  ;[[ 0,-2, 0]]_PC ; Distance is     2.0000
        list(1:4, 7)=origin(1:4)+(/  0, -1,  1,  1/)!_FCC  ;[[ 1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 8)=origin(1:4)+(/ -1,  0,  0,  0/)!_FCC  ;[[-2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4, 9)=origin(1:4)+(/ -1,  0,  0,  3/)!_FCC  ;[[-1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,10)=origin(1:4)+(/  0,  0,  1,  0/)!_FCC  ;[[ 0,-1, 1]]_PC ; Distance is 1.4142
        list(1:4,11)=origin(1:4)+(/  0,  0,  1,  3/)!_FCC  ;[[ 1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4,12)=origin(1:4)+(/ -1,  0,  0,  2/)!_FCC  ;[[-2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,13)=origin(1:4)+(/ -1,  0,  1,  1/)!_FCC  ;[[-1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,14)=origin(1:4)+(/  0,  0,  1,  2/)!_FCC  ;[[ 0, 0, 2]]_PC ; Distance is     2.0000
        list(1:4,15)=origin(1:4)+(/ -1,  1,  1,  0/)!_FCC  ;[[-2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,16)=origin(1:4)+(/ -1,  1,  1,  3/)!_FCC  ;[[-1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,17)=origin(1:4)+(/  0, -1,  0,  1/)!_FCC  ;[[ 1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4,18)=origin(1:4)+(/ -1,  0, -1,  3/)!_FCC  ;[[-1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,19)=origin(1:4)+(/  0,  0,  0,  0/)!_FCC  ;[[ 0,-1,-1]]_PC ; Distance is 1.4142
        list(1:4,20)=origin(1:4)+(/  0,  0,  0,  3/)!_FCC  ;[[ 1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,21)=origin(1:4)+(/  1,  0,  1,  0/)!_FCC  ;[[ 2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4,22)=origin(1:4)+(/ -1,  0,  0,  1/)!_FCC  ;[[-1, 0,-1]]_PC ; Distance is 1.4142!list(1:4,1)=origin(1:4)+fcc_to_pc4d( (/ 0, 0, 0/) )!_FCC  ;[[ 0, 0, 0]]_PC ; Distance is0.0000
        list(1:4,23)=origin(1:4)+(/  0,  0,  1,  1/)!_FCC  ;[[ 1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,24)=origin(1:4)+(/ -1,  1,  0,  0/)!_FCC  ;[[-2, 1,-1]]_PC ; Distance is        2.4495
        list(1:4,25)=origin(1:4)+(/ -1,  1,  0,  3/)!_FCC  ;[[-1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,26)=origin(1:4)+(/  0,  1,  1,  0/)!_FCC  ;[[ 0, 1, 1]]_PC ; Distance is 1.4142
        list(1:4,27)=origin(1:4)+(/  0,  1,  1,  3/)!_FCC  ;[[ 1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,28)=origin(1:4)+(/ -1,  1,  1,  1/)!_FCC  ;[[-1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,29)=origin(1:4)+(/  0,  0, -1,  3/)!_FCC  ;[[ 1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,30)=origin(1:4)+(/  1,  0,  0,  0/)!_FCC  ;[[ 2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4,31)=origin(1:4)+(/  0,  0, -1,  2/)!_FCC  ;[[ 0, 0,-2]]_PC ; Distance is     2.0000
        list(1:4,32)=origin(1:4)+(/  0,  0,  0,  1/)!_FCC  ;[[ 1, 0,-1]]_PC ; Distance is 1.4142
        list(1:4,33)=origin(1:4)+(/  1,  0,  0,  2/)!_FCC  ;[[ 2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,34)=origin(1:4)+(/ -1,  1, -1,  3/)!_FCC  ;[[-1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,35)=origin(1:4)+(/  0,  1,  0,  0/)!_FCC  ;[[ 0, 1,-1]]_PC ; Distance is 1.4142
        list(1:4,36)=origin(1:4)+(/  0,  1,  0,  3/)!_FCC  ;[[ 1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,37)=origin(1:4)+(/  1,  1,  1,  0/)!_FCC  ;[[ 2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,38)=origin(1:4)+(/ -1,  1,  0,  1/)!_FCC  ;[[-1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4,39)=origin(1:4)+(/  0,  1,  0,  2/)!_FCC  ;[[ 0, 2, 0]]_PC ; Distance is     2.0000
        list(1:4,40)=origin(1:4)+(/  0,  1,  1,  1/)!_FCC  ;[[ 1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,41)=origin(1:4)+(/  0,  1, -1,  3/)!_FCC  ;[[ 1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,42)=origin(1:4)+(/  1,  1,  0,  0/)!_FCC  ;[[ 2, 1,-1]]_PC ; Distance is        2.4495
    elseif(icell .eq. 3) then
        list(1:4, 1)=origin(1:4)+(/  1,  1,  0,  0/)!_FCC  ;[[ 1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4, 2)=origin(1:4)+(/  0, -1,  1,  0/)!_FCC  ;[[-1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 3)=origin(1:4)+(/ -1, -1,  1,  1/)!_FCC  ;[[-2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4, 4)=origin(1:4)+(/  0, -1,  1,  2/)!_FCC  ;[[-1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4, 5)=origin(1:4)+(/  0, -1,  0,  0/)!_FCC  ;[[-1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4, 6)=origin(1:4)+(/  0, -1,  0,  3/)!_FCC  ;[[ 0,-2, 0]]_PC ; Distance is     2.0000
        list(1:4, 7)=origin(1:4)+(/  1, -1,  1,  0/)!_FCC  ;[[ 1,-2, 1]]_PC ; Distance is        2.4495
        list(1:4, 8)=origin(1:4)+(/ -1, -1,  0,  1/)!_FCC  ;[[-2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4, 9)=origin(1:4)+(/  0, -1,  0,  2/)!_FCC  ;[[-1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,10)=origin(1:4)+(/  0, -1,  1,  1/)!_FCC  ;[[ 0,-1, 1]]_PC ; Distance is 1.4142
        list(1:4,11)=origin(1:4)+(/  1, -1,  1,  2/)!_FCC  ;[[ 1,-1, 2]]_PC ; Distance is        2.4495
        list(1:4,12)=origin(1:4)+(/ -1,  0,  0,  3/)!_FCC  ;[[-2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,13)=origin(1:4)+(/  0,  0,  1,  0/)!_FCC  ;[[-1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,14)=origin(1:4)+(/  0,  0,  1,  3/)!_FCC  ;[[ 0, 0, 2]]_PC ; Distance is     2.0000
        list(1:4,15)=origin(1:4)+(/ -1,  0,  1,  1/)!_FCC  ;[[-2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,16)=origin(1:4)+(/  0,  0,  1,  2/)!_FCC  ;[[-1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,17)=origin(1:4)+(/  1, -1,  0,  0/)!_FCC  ;[[ 1,-2,-1]]_PC ; Distance is        2.4495
        list(1:4,18)=origin(1:4)+(/  0, -1, -1,  2/)!_FCC  ;[[-1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,19)=origin(1:4)+(/  0, -1,  0,  1/)!_FCC  ;[[ 0,-1,-1]]_PC ; Distance is 1.4142
        list(1:4,20)=origin(1:4)+(/  1, -1,  0,  2/)!_FCC  ;[[ 1,-1, 0]]_PC ; Distance is 1.4142
        list(1:4,21)=origin(1:4)+(/  1, -1,  1,  1/)!_FCC  ;[[ 2,-1, 1]]_PC ; Distance is        2.4495
        list(1:4,22)=origin(1:4)+(/  0,  0,  0,  0/)!_FCC  ;[[-1, 0,-1]]_PC ; Distance is 1.4142!list(1:4,1)=origin(1:4)+fcc_to_pc4d( (/ 0, 0, 0/) )!_FCC  ;[[ 0, 0, 0]]_PC ; Distance is0.0000
        list(1:4,23)=origin(1:4)+(/  1,  0,  1,  0/)!_FCC  ;[[ 1, 0, 1]]_PC ; Distance is 1.4142
        list(1:4,24)=origin(1:4)+(/ -1,  0,  0,  1/)!_FCC  ;[[-2, 1,-1]]_PC ; Distance is        2.4495
        list(1:4,25)=origin(1:4)+(/  0,  0,  0,  2/)!_FCC  ;[[-1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,26)=origin(1:4)+(/  0,  0,  1,  1/)!_FCC  ;[[ 0, 1, 1]]_PC ; Distance is 1.4142
        list(1:4,27)=origin(1:4)+(/  1,  0,  1,  2/)!_FCC  ;[[ 1, 1, 2]]_PC ; Distance is        2.4495
        list(1:4,28)=origin(1:4)+(/  0,  1,  1,  0/)!_FCC  ;[[-1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,29)=origin(1:4)+(/  1, -1, -1,  2/)!_FCC  ;[[ 1,-1,-2]]_PC ; Distance is        2.4495
        list(1:4,30)=origin(1:4)+(/  1, -1,  0,  1/)!_FCC  ;[[ 2,-1,-1]]_PC ; Distance is        2.4495
        list(1:4,31)=origin(1:4)+(/  0,  0, -1,  3/)!_FCC  ;[[ 0, 0,-2]]_PC ; Distance is     2.0000
        list(1:4,32)=origin(1:4)+(/  1,  0,  0,  0/)!_FCC  ;[[ 1, 0,-1]]_PC ; Distance is 1.4142
        list(1:4,33)=origin(1:4)+(/  1,  0,  0,  3/)!_FCC  ;[[ 2, 0, 0]]_PC ; Distance is     2.0000
        list(1:4,34)=origin(1:4)+(/  0,  0, -1,  2/)!_FCC  ;[[-1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,35)=origin(1:4)+(/  0,  0,  0,  1/)!_FCC  ;[[ 0, 1,-1]]_PC ; Distance is 1.4142
        list(1:4,36)=origin(1:4)+(/  1,  0,  0,  2/)!_FCC  ;[[ 1, 1, 0]]_PC ; Distance is 1.4142
        list(1:4,37)=origin(1:4)+(/  1,  0,  1,  1/)!_FCC  ;[[ 2, 1, 1]]_PC ; Distance is        2.4495
        list(1:4,38)=origin(1:4)+(/  0,  1,  0,  0/)!_FCC  ;[[-1, 2,-1]]_PC ; Distance is        2.4495
        list(1:4,39)=origin(1:4)+(/  0,  1,  0,  3/)!_FCC  ;[[ 0, 2, 0]]_PC ; Distance is     2.0000
        list(1:4,40)=origin(1:4)+(/  1,  1,  1,  0/)!_FCC  ;[[ 1, 2, 1]]_PC ; Distance is        2.4495
        list(1:4,41)=origin(1:4)+(/  1,  0, -1,  2/)!_FCC  ;[[ 1, 1,-2]]_PC ; Distance is        2.4495
        list(1:4,42)=origin(1:4)+(/  1,  0,  0,  1/)!_FCC  ;[[ 2, 1,-1]]_PC ; Distance is        2.4495
    endif

endsubroutine list_of_42_nearest_cells_4dpc

subroutine    check_escaped_atom
    use positions_mod
    use interaction_mod
    use phys_parameters_mod

    implicit none
    integer :: ixy,iyz,izx,icell
    integer :: ixyc,iyzc,izxc,icellc
    integer, dimension(4,42) :: list
    integer cn,anout_cl,aninn_cl,anout_gl,aninn_gl,perc
    logical :: is_cell_close
    real(8) :: dist

    !hotim uznatq jestq li takoj atom, cxto jego net v blizkih jacxejkah,
    !no on blizxe cxem katoff
    !nuzxno uznatq kakije jacxejki blizki
    !dalqsxe perebiratq atomy vnutri dalqokih jacxejek
    !vycxislqatq rasstojanija do nih, i sravnivatq s katoff
    ixyc=0
    iyzc=0
    izxc=0
    icellc=0
    call list_of_42_nearest_cells_4dpc(ixyc,iyzc,izxc,icellc,list)
    perc=100!3
    !print*,-42,ixy,iyz,izx,"+++"
    do ixy=-cells_xrange,cells_xrange
        !print"(A,(I4.3,1x))", "Proverka jacxejek ixy = ",ixy
        do iyz=-cells_yrange,cells_yrange
            do izx=-cells_zrange,cells_zrange
            do icell = 1,4
                is_cell_close=.false.
                do cn=1,42
                    if(&
                    (list(1,cn) .eq. ixy) .and. &
                    (list(2,cn) .eq. iyz) .and. &
                    (list(3,cn) .eq. izx) .and. &
                    (list(4,cn) .eq. icell)  ) then
                        is_cell_close=.true.
                    endif
                enddo
                if (is_cell_close) cycle
                !print"(A,3(I4.3,1x))", "Proverka jacxejki ",ixy,iyz,izx,icell
                do anout_cl=1,cn_by_an_4d(0,ixy,iyz,izx,icell)
                    anout_gl=cn_by_an_4d(anout_cl,ixy,iyz,izx,icell)
                    do aninn_cl=1,cn_by_an_4d(0,ixyc,iyzc,izxc,icellc)
                        aninn_gl=cn_by_an_4d(aninn_cl,ixyc,iyzc,izxc,icellc)
                        dist=norm2( R_curr(1:3,anout_gl) - R_curr(1:3,aninn_gl) )
                        if(dist .le. cutoff_param*1d-2*perc) then

                            print 269,&
                            "Atom #",aninn_gl," iz jacxejki ",ixyc,iyzc,izxc,icellc,"_PC4d, koordinaty ",R_curr(1:3,aninn_gl),&
                            " i atom #",anout_gl," iz jacxejki ",ixy,iyz,izx,icell,"_PC4d, koordinaty ",R_curr(1:3,anout_gl),&
                            "blizxe cxem katoff ",perc," %, hotqa eto ne dolzxno bytq. dist = ",dist," cutoff_param = ",cutoff_param
                            !call sleep(1)
                        endif
                    enddo
                enddo
            enddo
            enddo
        enddo
    enddo
269 format( 2( A,I7.2,A,4(I4.1,1x),A,3(F7.3,1x) ),/,A,I3.3,A,F9.5,A,F10.5 )

endsubroutine check_escaped_atom
