real(8) function distance_inunits_from_origin_fcc(ixy,iyz,izx)

    implicit none
    integer, intent(in) :: ixy,iyz,izx
    real(8) :: dist

    dist = 5d-1*sqrt(0d0+&
        ( ixy-iyz+izx)**2+&
        ( ixy+iyz-izx)**2+&
        (-ixy+iyz+izx)**2 )
    distance_inunits_from_origin_fcc = dist

endfunction distance_inunits_from_origin_fcc

logical function test_if_cell_is_close_4d(ixy,iyz,izx,iic)

    implicit none
    integer, intent(in) :: ixy,iyz,izx,iic
    integer, dimension(4,42) :: list
    integer cn
    logical :: lr

    lr=.true.
    if ( (ixy .eq. 0) .and. ( iyz .eq. 0) .and. (izx .eq. 0) .and. (iic .eq. 0) ) then
        test_if_cell_is_close_4d= .true.
        return
    endif
    call list_of_42_nearest_cells_4dpc(0,0,0,0,list)
    lr=.false.
    !print*,-42,ixy,iyz,izx,"+++"
    do cn=1,42
        !print*, cn,list(1:3,cn)
        lr = lr .or. (&
        list(1,cn) .eq. ixy .and. &
        list(2,cn) .eq. iyz .and. &
        list(3,cn) .eq. izx .and. &
        list(4,cn) .eq. iic )
        if( lr .eqv. .true.) exit
    enddo
    test_if_cell_is_close_4d=lr
endfunction test_if_cell_is_close_4d

function fcc_to_pc4d(fcc_cell) result(pc4d_cell)
    implicit none
    integer, dimension(3) :: fcc_cell,cn_raw
    integer, dimension(4) :: pc4d_cell

    !cnraw is pc
    cn_raw(1) = (fcc_cell(3) + fcc_cell(1) - fcc_cell(2))*2
    cn_raw(2) = (fcc_cell(1) + fcc_cell(2) - fcc_cell(3))*2
    cn_raw(3) = (fcc_cell(2) + fcc_cell(3) - fcc_cell(1))*2
    !pc with additions basis
    pc4d_cell(1) = ( cn_raw(1)-modulo( cn_raw(1),2 ) )/2
    pc4d_cell(2) = ( cn_raw(2)-modulo( cn_raw(2),2 ) )/2
    pc4d_cell(3) = ( cn_raw(3)-modulo( cn_raw(3),2 ) )/2
    if (     modulo(cn_raw(1),2) .eq. 1 .and. modulo(cn_raw(2),2) .eq. 1 ) then
        pc4d_cell(4) = 1
    elseif ( modulo(cn_raw(2),2) .eq. 1 .and. modulo(cn_raw(3),2) .eq. 1 ) then
        pc4d_cell(4) = 2
    elseif ( modulo(cn_raw(3),2) .eq. 1 .and. modulo(cn_raw(1),2) .eq. 1 ) then
        pc4d_cell(4) = 3
    else
        pc4d_cell(4) = 0
    endif
    !pc with additions basis

endfunction fcc_to_pc4d


!    print*,(/ 2, 0,-1/),"->",fcc_to_pc4d( (/ 2, 0,-1/) )!_FCC  ;[[ 1, 2,-1]]_PC   ; Distance is                  2.4495
!    print*,(/-2, 0, 1/),"->",fcc_to_pc4d( (/-2, 0, 1/) )!_FCC  ;[[-1,-2, 1]]_PC   ; Distance is                  2.4495
!    print*,(/-2, 1, 0/),"->",fcc_to_pc4d( (/-2, 1, 0/) )!_FCC  ;[[-2,-1, 1]]_PC   ; Distance is                  2.4495
!    print*,(/-2, 1, 1/),"->",fcc_to_pc4d( (/-2, 1, 1/) )!_FCC  ;[[-1,-1, 2]]_PC   ; Distance is                  2.4495
!    print*,(/-1,-1, 0/),"->",fcc_to_pc4d( (/-1,-1, 0/) )!_FCC  ;[[-1,-2,-1]]_PC   ; Distance is                  2.4495
!    print*,(/-1,-1, 1/),"->",fcc_to_pc4d( (/-1,-1, 1/) )!_FCC  ;[[ 0,-2, 0]]_PC   ; Distance is               2.0000
!    print*,(/-1,-1, 2/),"->",fcc_to_pc4d( (/-1,-1, 2/) )!_FCC  ;[[ 1,-2, 1]]_PC   ; Distance is                  2.4495
!    print*,(/-1, 0,-1/),"->",fcc_to_pc4d( (/-1, 0,-1/) )!_FCC  ;[[-2,-1,-1]]_PC   ; Distance is                  2.4495
!    print*,(/-1, 0, 0/),"->",fcc_to_pc4d( (/-1, 0, 0/) )!_FCC  ;[[-1,-1, 0]]_PC   ; Distance is           1.4142
!    print*,(/-1, 0, 1/),"->",fcc_to_pc4d( (/-1, 0, 1/) )!_FCC  ;[[ 0,-1, 1]]_PC   ; Distance is           1.4142
!    print*,(/-1, 0, 2/),"->",fcc_to_pc4d( (/-1, 0, 2/) )!_FCC  ;[[ 1,-1, 2]]_PC   ; Distance is                  2.4495
!    print*,(/-1, 1,-1/),"->",fcc_to_pc4d( (/-1, 1,-1/) )!_FCC  ;[[-2, 0, 0]]_PC   ; Distance is               2.0000
!    print*,(/-1, 1, 0/),"->",fcc_to_pc4d( (/-1, 1, 0/) )!_FCC  ;[[-1, 0, 1]]_PC   ; Distance is           1.4142
!    print*,(/-1, 1, 1/),"->",fcc_to_pc4d( (/-1, 1, 1/) )!_FCC  ;[[ 0, 0, 2]]_PC   ; Distance is               2.0000
!    print*,(/-1, 2,-1/),"->",fcc_to_pc4d( (/-1, 2,-1/) )!_FCC  ;[[-2, 1, 1]]_PC   ; Distance is                  2.4495
!    print*,(/-1, 2, 0/),"->",fcc_to_pc4d( (/-1, 2, 0/) )!_FCC  ;[[-1, 1, 2]]_PC   ; Distance is                  2.4495
!    print*,(/ 0,-2, 1/),"->",fcc_to_pc4d( (/ 0,-2, 1/) )!_FCC  ;[[ 1,-2,-1]]_PC   ; Distance is                  2.4495
!    print*,(/ 0,-1,-1/),"->",fcc_to_pc4d( (/ 0,-1,-1/) )!_FCC  ;[[-1,-1,-2]]_PC   ; Distance is                  2.4495
!    print*,(/ 0,-1, 0/),"->",fcc_to_pc4d( (/ 0,-1, 0/) )!_FCC  ;[[ 0,-1,-1]]_PC   ; Distance is           1.4142
!    print*,(/ 0,-1, 1/),"->",fcc_to_pc4d( (/ 0,-1, 1/) )!_FCC  ;[[ 1,-1, 0]]_PC   ; Distance is           1.4142
!    print*,(/ 0,-1, 2/),"->",fcc_to_pc4d( (/ 0,-1, 2/) )!_FCC  ;[[ 2,-1, 1]]_PC   ; Distance is                  2.4495
!    print*,(/ 0, 0,-1/),"->",fcc_to_pc4d( (/ 0, 0,-1/) )!_FCC  ;[[-1, 0,-1]]_PC   ; Distance is           1.4142
!    print*,(/ 0, 0, 0/),"->",fcc_to_pc4d( (/ 0, 0, 0/) )!_FCC  ;[[ 0, 0, 0]]_PC   ; Distance is   0.0000
!    print*,(/ 0, 0, 1/),"->",fcc_to_pc4d( (/ 0, 0, 1/) )!_FCC  ;[[ 1, 0, 1]]_PC   ; Distance is           1.4142
!    print*,(/ 0, 1,-2/),"->",fcc_to_pc4d( (/ 0, 1,-2/) )!_FCC  ;[[-2, 1,-1]]_PC   ; Distance is                  2.4495
!    print*,(/ 0, 1,-1/),"->",fcc_to_pc4d( (/ 0, 1,-1/) )!_FCC  ;[[-1, 1, 0]]_PC   ; Distance is           1.4142
!    print*,(/ 0, 1, 0/),"->",fcc_to_pc4d( (/ 0, 1, 0/) )!_FCC  ;[[ 0, 1, 1]]_PC   ; Distance is           1.4142
!    print*,(/ 0, 1, 1/),"->",fcc_to_pc4d( (/ 0, 1, 1/) )!_FCC  ;[[ 1, 1, 2]]_PC   ; Distance is                  2.4495
!    print*,(/ 0, 2,-1/),"->",fcc_to_pc4d( (/ 0, 2,-1/) )!_FCC  ;[[-1, 2, 1]]_PC   ; Distance is                  2.4495
!    print*,(/ 1,-2, 0/),"->",fcc_to_pc4d( (/ 1,-2, 0/) )!_FCC  ;[[ 1,-1,-2]]_PC   ; Distance is                  2.4495
!    print*,(/ 1,-2, 1/),"->",fcc_to_pc4d( (/ 1,-2, 1/) )!_FCC  ;[[ 2,-1,-1]]_PC   ; Distance is                  2.4495
!    print*,(/ 1,-1,-1/),"->",fcc_to_pc4d( (/ 1,-1,-1/) )!_FCC  ;[[ 0, 0,-2]]_PC   ; Distance is               2.0000
!    print*,(/ 1,-1, 0/),"->",fcc_to_pc4d( (/ 1,-1, 0/) )!_FCC  ;[[ 1, 0,-1]]_PC   ; Distance is           1.4142
!    print*,(/ 1,-1, 1/),"->",fcc_to_pc4d( (/ 1,-1, 1/) )!_FCC  ;[[ 2, 0, 0]]_PC   ; Distance is               2.0000
!    print*,(/ 1, 0,-2/),"->",fcc_to_pc4d( (/ 1, 0,-2/) )!_FCC  ;[[-1, 1,-2]]_PC   ; Distance is                  2.4495
!    print*,(/ 1, 0,-1/),"->",fcc_to_pc4d( (/ 1, 0,-1/) )!_FCC  ;[[ 0, 1,-1]]_PC   ; Distance is           1.4142
!    print*,(/ 1, 0, 0/),"->",fcc_to_pc4d( (/ 1, 0, 0/) )!_FCC  ;[[ 1, 1, 0]]_PC   ; Distance is           1.4142
!    print*,(/ 1, 0, 1/),"->",fcc_to_pc4d( (/ 1, 0, 1/) )!_FCC  ;[[ 2, 1, 1]]_PC   ; Distance is                  2.4495
!    print*,(/ 1, 1,-2/),"->",fcc_to_pc4d( (/ 1, 1,-2/) )!_FCC  ;[[-1, 2,-1]]_PC   ; Distance is                  2.4495
!    print*,(/ 1, 1,-1/),"->",fcc_to_pc4d( (/ 1, 1,-1/) )!_FCC  ;[[ 0, 2, 0]]_PC   ; Distance is               2.0000
!    print*,(/ 1, 1, 0/),"->",fcc_to_pc4d( (/ 1, 1, 0/) )!_FCC  ;[[ 1, 2, 1]]_PC   ; Distance is                  2.4495
!    print*,(/ 2,-1,-1/),"->",fcc_to_pc4d( (/ 2,-1,-1/) )!_FCC  ;[[ 1, 1,-2]]_PC   ; Distance is                  2.4495
!    print*,(/ 2,-1, 0/),"->",fcc_to_pc4d( (/ 2,-1, 0/) )!_FCC  ;[[ 2, 1,-1]]_PC   ; Distance is                  2.4495

!  02   00  -01 (/ 0, 1,-1, 3/)
! -02   00   01 (/-1,-1, 0, 3/)
! -02   01   00 (/-1,-1, 0, 2/)
! -02   01   01 (/-1,-1, 1, 1/)
! -01  -01   00 (/-1,-1,-1, 3/)
! -01  -01   01 (/ 0,-1, 0, 0/)
! -01  -01   02 (/ 0,-1, 0, 3/)
! -01   00  -01 (/-1,-1,-1, 2/)
! -01   00   00 (/-1,-1, 0, 1/)
! -01   00   01 (/ 0,-1, 0, 2/)
! -01   00   02 (/ 0,-1, 1, 1/)
! -01   01  -01 (/-1, 0, 0, 0/)
! -01   01   00 (/-1, 0, 0, 3/)
! -01   01   01 (/ 0, 0, 1, 0/)
! -01   02  -01 (/-1, 0, 0, 2/)
! -01   02   00 (/-1, 0, 1, 1/)
!  00  -02   01 (/ 0,-1,-1, 3/)
!  00  -01  -01 (/-1,-1,-1, 1/)
!  00  -01   00 (/ 0,-1,-1, 2/)
!  00  -01   01 (/ 0,-1, 0, 1/)
!  00  -01   02 (/ 1,-1, 0, 2/)
!  00   00  -01 (/-1, 0,-1, 3/)
!  00   00   00 (/ 0, 0, 0, 0/)
!  00   00   01 (/ 0, 0, 0, 3/)
!  00   01  -02 (/-1, 0,-1, 2/)
!  00   01  -01 (/-1, 0, 0, 1/)
!  00   01   00 (/ 0, 0, 0, 2/)
!  00   01   01 (/ 0, 0, 1, 1/)
!  00   02  -01 (/-1, 1, 0, 3/)
!  01  -02   00 (/ 0,-1,-1, 1/)
!  01  -02   01 (/ 1,-1,-1, 2/)
!  01  -01  -01 (/ 0, 0,-1, 0/)
!  01  -01   00 (/ 0, 0,-1, 3/)
!  01  -01   01 (/ 1, 0, 0, 0/)
!  01   00  -02 (/-1, 0,-1, 1/)
!  01   00  -01 (/ 0, 0,-1, 2/)
!  01   00   00 (/ 0, 0, 0, 1/)
!  01   00   01 (/ 1, 0, 0, 2/)
!  01   01  -02 (/-1, 1,-1, 3/)
!  01   01  -01 (/ 0, 1, 0, 0/)
!  01   01   00 (/ 0, 1, 0, 3/)
!  02  -01  -01 (/ 0, 0,-1, 1/)
!  02  -01   00 (/ 1, 0,-1, 2/)
