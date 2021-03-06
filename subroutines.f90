subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use phys_parameters_mod
    implicit none
    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q

    !print*,cells_xrange,cells_yrange,cells_zrange," eto x y z cells range"

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

subroutine fill_an_by_cn_4d
    use interaction_mod
    implicit none
    integer i_atoms

    print"(A,F9.5,A,3I4.2)"," Zapolnqajem 4d massiv an_by_cn. Radius obrezanija = ",cutoff_param," angstrem"

    do i_atoms=1,atoms__in_total
        call fill_an_by_cn_4d_peratomic(i_atoms)
        if ( an_by_cn_4d(1,i_atoms)**2 .gt. cells_xrange**2 ) stop " atom vysxel za dopustimoje kolicxestvo jacxeek po osi X"
        if ( an_by_cn_4d(2,i_atoms)**2 .gt. cells_yrange**2 ) stop " atom vysxel za dopustimoje kolicxestvo jacxeek po osi Y"
        if ( an_by_cn_4d(3,i_atoms)**2 .gt. cells_zrange**2 ) stop " atom vysxel za dopustimoje kolicxestvo jacxeek po osi Z"
    enddo

endsubroutine fill_an_by_cn_4d

subroutine fill_an_by_cn_4d_peratomic(i_atoms)
    use positions_mod
    use interaction_mod
    use phys_parameters_mod
    implicit none
    integer, intent(in) :: i_atoms
    integer cn_raw(3)
    real(8) scale_factor,r_at(3),r_rel(3),r_rel_abs(3),cellpos(3)
    real(8) hc,hc_reci

    logical(1) :: our(6)
!    logical(1), parameter, dimension(6) :: &
!    xp=(/ .true.  , .false. , .true.  , .true.  , .true.  , .true.  /),&
!    xn=(/ .false. , .true.  , .true.  , .true.  , .false. , .false. /),&
!    yp=(/ .true.  , .true.  , .true.  , .false. , .true.  , .true.  /),&
!    yn=(/ .false. , .false. , .false. , .true.  , .true.  , .true.  /),&
!    zp=(/ .true.  , .true.  , .true.  , .true.  , .true.  , .false. /),&
!    zn=(/ .true.  , .true.  , .false. , .false. , .false. , .true.  /)
!    !      x+y>0  ; -x+y>0  ;  y+z>0  ; -y+z>0  ;  z+x>0  ; -z+x>0

!    hc=sqrt(16d0/27d0)*(cutoff_param) !jesli granq ravna cutoff
    hc=sqrt(4d0/3d0)*(cutoff_param) !jesli diametr vpisannoj sfery raven cutoff, 2h
    hc_reci=1d0/hc

    r_at(1:3)=R_curr(1:3,i_atoms)

    cn_raw(1:3)=nint( hc_reci*2d0*(r_at(1:3)) ) !delim na menqsxije kubiki

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
        !cycle
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

endsubroutine fill_an_by_cn_4d_peratomic

subroutine fill_cn_by_an_4d
    use positions_mod
    use interaction_mod
    use phys_parameters_mod

    implicit none
    integer i_atoms,na
    integer anxy,anyz,anzx,anic

    print"(A,F9.5,A,3I4.2)"," Zapolnqajem massiv cn_by_an_4d. Radius obrezanija = ",cutoff_param," angstrem"
    cn_by_an_4d=0
    na=0

    do i_atoms=1,atoms__in_total

        if(na+1 .ge. incell_atoms_max) stop "massiv cn_by_an perepolnen"
        anxy=an_by_cn_4d(1,i_atoms)
        anyz=an_by_cn_4d(2,i_atoms)
        anzx=an_by_cn_4d(3,i_atoms)
        anic=an_by_cn_4d(4,i_atoms)

        na = cn_by_an_4d( 0,anxy,anyz,anzx,anic)
        cn_by_an_4d( 0,anxy,anyz,anzx,anic ) = na + 1
        cn_by_an_4d( na+1,anxy,anyz,anzx,anic ) = i_atoms
    enddo

    print*,"Zapolnenije cn_by_an_4d zaversxeno. Vsego atomov ", sum(cn_by_an_4d(0,:,:,:,:)), " = ",atoms__in_total

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

subroutine get_direct_list(i_a)
    use positions_mod
    use interaction_mod
    use phys_parameters_mod
    use lists_assortiment_mod
    implicit none

    integer, intent(in) :: i_a !atom, sosedi kotorogo trebyjetsqa najti
    !integer, dimension(4,42) :: list_cells
    integer ixy,iyz,izx,icl
    integer jxy,jyz,jzx,jcl
    integer i_aan,i_acn
    integer q_aan,n_aan, j
    !real(8) :: shorting = 1d0!3d0/3d0

    !stop "obrabotatq slucxaji vyhoda za granicu massiva jacxeek pri ispolqzovanii spiska "
    !neobhodimo vyzyvatq dlqa aktualqnogo

    ixy = an_by_cn_4d(1,i_a)
    iyz = an_by_cn_4d(2,i_a)
    izx = an_by_cn_4d(3,i_a)
    icl = an_by_cn_4d(4,i_a)

    list_atoms_qnt = 0
    list_atoms = 0
    call list_of_42_nearest_cells_4dpc(ixy,iyz,izx,icl,list_cells)

    q_aan = cn_by_an_4d(0,ixy,iyz,izx,icl)

    do i_aan = 1,q_aan
        n_aan = cn_by_an_4d(i_aan,ixy,iyz,izx,icl)
        if (n_aan .eq. i_a) cycle
        !inside the cube
!        if ( abs( R_curr(1,n_aan) - R_curr(1,i_a) ) .gt. cutoff ) cycle
!        if ( abs( R_curr(2,n_aan) - R_curr(2,i_a) ) .gt. cutoff ) cycle
!        if ( abs( R_curr(3,n_aan) - R_curr(3,i_a) ) .gt. cutoff ) cycle
!        !!inside the octahedron
        if ( abs(   +( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                    +( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                    +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                ) .gt. cutoff*sqrt(3d0) ) cycle
        if ( abs(   -( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                    +( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                    +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                ) .gt. cutoff*sqrt(3d0) ) cycle
        if ( abs(   +( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                    -( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                    +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                ) .gt. cutoff*sqrt(3d0) ) cycle
        if ( abs(   -( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                    -( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                    +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                ) .gt. cutoff*sqrt(3d0) ) cycle
        !!mozxet bytq lucxsxe ne kub s oktaedrom( t.o. usecxonnyj oktaedr) a rombododekaedr
        !IF ( NORM2( R_CURR(1:3,N_AAN) - R_CURR(1:3,I_A) ) .GT. CUTOFF ) CYCLE
        list_atoms_qnt = list_atoms_qnt + 1
        list_atoms(list_atoms_qnt) = n_aan
    enddo


    do i_acn=1,42

        jxy = list_cells(1,i_acn)
        if ( jxy**2 .gt. cells_xrange**2 ) cycle
        jyz = list_cells(2,i_acn)
        if ( jyz**2 .gt. cells_yrange**2 ) cycle
        jzx = list_cells(3,i_acn)
        if ( jzx**2 .gt. cells_zrange**2 ) cycle
        jcl = list_cells(4,i_acn)

        q_aan = cn_by_an_4d(0,jxy,jyz,jzx,jcl)

        !PRINT*, " V JACXEJKE ", JXY,JYZ,JZX,JCL, " lezxit ", Q_AAN, " ATOMOV "
        !CALL SLEEP(1)

        do i_aan = 1,q_aan

            n_aan = cn_by_an_4d(i_aan,jxy,jyz,jzx,jcl)

            !inside the cube
!            if ( abs( R_curr(1,n_aan) - R_curr(1,i_a) ) .gt. cutoff ) cycle
!            if ( abs( R_curr(2,n_aan) - R_curr(2,i_a) ) .gt. cutoff ) cycle
!            if ( abs( R_curr(3,n_aan) - R_curr(3,i_a) ) .gt. cutoff ) cycle
!            !inside the octahedron
            if ( abs(   +( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                        +( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                        +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                    ) .gt. cutoff*sqrt(3d0) ) cycle
            if ( abs(   -( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                        +( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                        +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                    ) .gt. cutoff*sqrt(3d0) ) cycle
            if ( abs(   +( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                        -( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                        +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                    ) .gt. cutoff*sqrt(3d0) ) cycle
            if ( abs(   -( R_curr(1,n_aan) - R_curr(1,i_a) ) &
                        -( R_curr(2,n_aan) - R_curr(2,i_a) ) &
                        +( R_curr(3,n_aan) - R_curr(3,i_a) ) &
                    ) .gt. cutoff*sqrt(3d0) ) cycle
            !mozxet bytq lucxsxe ne kub s oktaedrom( t.o. usecxonnyj oktaedr) a rombododekaedr

            !IF ( NORM2( R_CURR(1:3,N_AAN) - R_CURR(1:3,I_A) ) .GT. CUTOFF ) CYCLE

            list_atoms_qnt = list_atoms_qnt + 1
            list_atoms(list_atoms_qnt) = n_aan

        enddo
    enddo
    !CALL SLEEP(3)
    !print*, "atom #",i_a," imeet ",list_atoms_qnt," sosedej "
!
!    if(i_a .eq. 25001) then
!        OPEN (204, FILE = "25000_rhdc.txt")
!        DO J=1,LIST_ATOMS_QNT
!            N_AAN=LIST_ATOMS(J)
!            PRINT*,"ATOM ",I_A," SOSED DLQA ATOMA#",N_AAN, ". RASSTOJANIJE ",NORM2( R_CURR(1:3,N_AAN) - R_CURR(1:3,I_A) )
!            WRITE(204,*),"ATOM ",I_A," SOSED DLQA ATOMA#",N_AAN, &
!                ". RASSTOJANIJE ",NORM2( R_CURR(1:3,N_AAN) - R_CURR(1:3,I_A) )
!        ENDDO
!        CLOSE(204)
!    endif
endsubroutine get_direct_list

subroutine initiate_energetic_parameters
    use energetic_linappr_mod
    use ackland2003potential_fe_mod
    implicit none

    integer i
    real(8) r, rho
!    real(8) pw_fefe_an
!    real(8) ed_fefe_an
!    real(8) mf_fe_an

!    open (90, file = 'pw.txt')
!    open (91, file = 'ed.txt')
!    open (92, file = 'mf.txt')

    do i = 1,pw_steps
        r   = i * pw_step
        pw_pot_val(i) = pw_fefe_an(r)
!        pairwise_dr1(i) = pw_analytic_deriv1(r)
!        pairwise_dr2(i) = pw_analytic_deriv2(r)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            pw_ad1_val(i-1)=(pw_pot_val(i)-pw_pot_val(i-1))*pw_recistep
        endif
                !write(90, *) r, pairwise_pot(i)
    enddo
    pw_ad1_val(pw_steps) = pw_ad1_val(pw_steps-1)

    do i = 1,ed_steps
        r   = i * ed_step
        ed_pot_val(i) = ed_fefe_an(r)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            ed_ad1_val(i-1)=(ed_pot_val(i)-ed_pot_val(i-1))*ed_recistep
        endif
                !write(91, *) r, elecdens_pot(i)
    enddo
    ed_ad1_val(ed_steps) = ed_ad1_val(ed_steps-1)

    do i = 1,mf_steps
        rho = i * mf_step
        mf_pot_val(i) = mf_fe_an(rho)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            mf_ad1_val(i-1)=(mf_pot_val(i)-mf_pot_val(i-1))*mf_recistep
        endif
                !write(92, *) r, embefunc_pot(i)
    enddo
    mf_ad1_val(mf_steps) = mf_ad1_val(mf_steps-1)


    !close(90);close(91);close(92)
    write(*,*) "Priblizxonnyje znacxenija potencialov polucxeny iz analiticxeskih. OK"

endsubroutine initiate_energetic_parameters

subroutine  get_verlet_list_long(i_a)
    !use array_parameters_mod
    use positions_mod
    use interaction_mod, only : cutoff
    use verlet_lists_mod
    implicit none
    integer, intent(in) :: i_a !atom, sosedi kotorogo trebyjetsqa najti
    integer j,list_atoms_qnt_loc,omp_get_thread_num
    real(8) delx,dely,delz,delc,dist2


        integer :: l_vl_len! skolqko atomov v spiske
        integer :: l_vl_len_ext! skolqko atomov v spiske
        integer :: l_verlet_list(neibors_max)
        integer :: l_verlet_list_long(neibors_max*12) !inside double cutoff
        real(8) :: l_distan_list(neibors_max)
        !real(8) :: l_alm_ed_list(neibors_max) !elektronnyje plotnosti dlqa sosedej, krome vklada centralqnogo atoma
        !real(8) :: l_eldens_list(neibors_max) !elektronnyje plotnosti dlqa sosedej

    atom_n = i_a
    delc = cutoff*cutoff

    l_vl_len = 0
    l_vl_len_ext = 0

    do j=1,atoms__in_total

        delx = r_curr(1,j) - r_curr(1,i_a)
        if ( delx .gt.  cutoff*2) cycle
        if ( delx .lt. -cutoff*2) cycle
        dely = r_curr(2,j) - r_curr(2,i_a)
        if ( dely .gt.  cutoff*2) cycle
        if ( dely .lt. -cutoff*2) cycle
        delz = r_curr(3,j) - r_curr(3,i_a)
        if ( delz .gt.  cutoff*2) cycle
        if ( delz .lt. -cutoff*2) cycle

        dist2 = delx*delx + dely*dely + delz*delz
        if ( dist2 .gt. delc*4 ) cycle

        l_vl_len_ext = l_vl_len_ext + 1
        !PRINT*,J," J-",L_VL_LEN_EXT
        l_verlet_list_long(l_vl_len_ext) = j

        if ( delx .gt.  cutoff) cycle
        if ( delx .lt. -cutoff) cycle
        if ( dely .gt.  cutoff) cycle
        if ( dely .lt. -cutoff) cycle
        if ( delz .gt.  cutoff) cycle
        if ( delz .lt. -cutoff) cycle
        if ( dist2 .gt. delc ) cycle

        if (j .eq. i_a) cycle
        l_vl_len = l_vl_len + 1
        !PRINT*,J," J+",L_VL_LEN
        l_verlet_list(l_vl_len) = j
        l_distan_list(l_vl_len) = sqrt(dist2)
    enddo

    !PRINT*,VL_LEN_EXT,VL_LEN,I_A," i_a"
    !CALL SLEEP(1)
    vl_len = l_vl_len
    vl_len_ext = l_vl_len_ext
    distan_list = l_distan_list
    verlet_list = l_verlet_list
    verlet_list_long = l_verlet_list_long

endsubroutine get_verlet_list_long

subroutine  get_verlet_list_short(i_a)
    !use array_parameters_mod
    use positions_mod
    use interaction_mod, only : cutoff
    use verlet_lists_mod
    implicit none
    integer, intent(in) :: i_a !atom, sosedi kotorogo trebyjetsqa najti
    integer j
    real(8) delx,dely,delz,delc,dist2

    atom_n = i_a
    delc = cutoff*cutoff

    vl_len = 0

    do j=1,atoms__in_total
        delx = r_curr(1,j) - r_curr(1,i_a)
        if ( delx .gt.  cutoff) cycle
        if ( delx .lt. -cutoff) cycle
        dely = r_curr(2,j) - r_curr(2,i_a)
        if ( dely .gt.  cutoff) cycle
        if ( dely .lt. -cutoff) cycle
        delz = r_curr(3,j) - r_curr(3,i_a)
        if ( delz .gt.  cutoff) cycle
        if ( delz .lt. -cutoff) cycle

        dist2 = delx*delx + dely*dely + delz*delz

        if ( dist2 .gt. delc ) cycle
        if (j .eq. i_a) cycle

        vl_len = vl_len + 1
        !PRINT*,J," J+",VL_LEN," dist ",dist2
        verlet_list(vl_len) = j
        distan_list(vl_len) = sqrt(dist2)
    enddo

    !PRINT*,VL_LEN_EXT,VL_LEN,I_A," i_a"
    !CALL SLEEP(1)

endsubroutine get_verlet_list_short

subroutine part_and_sort_three_zones_100
    use positions_mod
    use phys_parameters_mod
    implicit none
    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q
    integer i_atp, i_atc

    !print*,cells_xrange,cells_yrange,cells_zrange," eto x y z cells range"

    i_atc = 0

    do i_atp = 1,atoms__in_total
        if ( ( norm2( R_perf(1:2,i_atp) ) .lt. x_layers*a0 - wall_thickness ) .and. &
             ( abs(   R_perf(3,  i_atp) ) .lt. z_layers*a0 - wall_thickness ) ) then
            if (atoms__in_total .le. i_atc) stop "error 1"
            i_atc = i_atc + 1
            R_curr(1:3,i_atc) = R_perf(1:3,i_atp)
        endif
    enddo

    last_inbody = i_atc

    do i_atp = 1,atoms__in_total
        if ( ( norm2( R_perf(1:2,i_atp) ) .lt. x_layers*a0 - wall_thickness ) .and. &
             ( abs(   R_perf(3,  i_atp) ) .gt. z_layers*a0 - wall_thickness ) ) then
            if (atoms__in_total .le. i_atc) stop "error 2"
            i_atc = i_atc + 1
            R_curr(1:3,i_atc) = R_perf(1:3,i_atp)
        endif
    enddo

    last_inmirr = i_atc

    do i_atp = 1,atoms__in_total
        if ( norm2( R_perf(1:2,i_atp) ) .gt. x_layers*a0 - wall_thickness )  then
            if (atoms__in_total .le. i_atc) stop "error 3"
            i_atc = i_atc + 1
            R_curr(1:3,i_atc) = R_perf(1:3,i_atp)
        endif
    enddo

    if (atoms__in_total .gt. i_atc) stop "poterqalisq atomy v processe razdelenija na zony"
    R_perf = R_curr

    print 107, " Sistema razbita na zony. Vsego atomov: ", i_atc," ; "
    print 106, " Osnovnaja jacxejka : ",last_inbody," atomov." // &
        " Dublirujuhxije sloi : ",last_inmirr-last_inbody," atomov." // &
        " Uprugaja zona ",1+i_atc-last_inmirr,&
        "; tolsxina uprugoj zony (stenki) >= ", wall_thickness, " angstrem "

106 format  (3(A,1x,I7.3,1x),A,F8.4,A)
107 format  (A,I10.3,A,$)
endsubroutine part_and_sort_three_zones_100


subroutine find_symmetric_pairs_100
    use positions_mod
    use symmetry_pairs_mod
    use phys_parameters_mod
    implicit none
    !real(8) a_half
    !integer cnt_x, cnt_y, cnt_z, cnt_q
    integer i_atb, i_atp

    bodymarknumber = -1

    print*,"Ihxem simmetricxnyje pary."

    do i_atp = 1 + last_inbody, last_inmirr
        do i_atb = 1, last_inbody
            if ( ( norm2( R_perf(1:2,i_atp) - R_perf(1:2,i_atb) ) .lt. 1d-2 ) .and. &
                 ( R_perf(3,  i_atb) .lt. a0*6d-1 ) .and. &
                 ( R_perf(3,  i_atb) .gt.   -1d-2 ) ) then
                bodymarknumber(i_atp) = i_atb
            endif
        enddo
        if (bodymarknumber(i_atp) .eq. -1) then
            print*, "cant find pair for atom ", i_atp
            print*, R_perf(1:3,i_atp)
            stop "cant find pair"
        endif
    enddo

    print *, " Najdeny vse pary "
endsubroutine find_symmetric_pairs_100


subroutine  atom_energy(i_a,e_a)
    use verlet_lists_mod
    USE ACKLAND2003POTENTIAL_FE_MOD
    implicit none
    integer, intent(in) :: i_a
    real(8), intent(inout) :: e_a
    real(8) :: eld,iip,dist
    real(8) :: ed_fefe_la
    real(8) :: pw_fefe_la
    real(8) :: mf_fe_la
    integer :: n_n,n_a

    call get_verlet_list_short(i_a)

    eld = 0
    iip = 0
    do n_n = 1,vl_len
        n_a = verlet_list(n_n)
        dist = distan_list(n_n)
        !PRINT*, N_A, IIP, ELD, "A_E",DIST
        eld = eld + ed_fefe_la(dist)
        iip = iip + pw_fefe_la(dist)
!        ELD = ELD + ED_FEFE_AN(DIST)
!        IIP = IIP + PW_FEFE_AN(DIST)
    enddo
    e_a = iip + mf_fe_la(eld)
!    e_a = iip + mf_fe_an(eld)
endsubroutine
