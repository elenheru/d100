pure real(8) function biersack_ziegler(x)
    implicit none !ackland 2003
    real(8), intent(in) :: x
    biersack_ziegler = &
    ( 0.18180d0 * exp( -3.2000d0 * x ) ) + &
    ( 0.50990d0 * exp( -0.9423d0 * x ) ) + &
    ( 0.28020d0 * exp( -0.4029d0 * x ) ) + &
    ( 0.02817d0 * exp( -0.2016d0 * x ) )
endfunction biersack_ziegler

pure real(8) function hvs(x)
    implicit none
    real(8), intent(in) :: x
    if(x.gt.0d0) then
        hvs=1d0
    else
        hvs=0d0
    endif
endfunction hvs

elemental real(8) function pw_fefe_an(r)
    use energetic_parameters_mod
    implicit none
    real(8), intent(in) :: r
    real(8) res
    integer i

    interface
        pure real(8) function biersack_ziegler(x)
            real(8), intent(in) :: x
        endfunction biersack_ziegler
        pure real(8) function hvs(x)
            real(8), intent(in) :: x
        endfunction hvs
    endinterface

    res = 0d0

    if(r .le. r1_fe) then
        res = electron_charge*electron_charge * z_fe*z_fe * biersack_ziegler( r/r_s_fe )/r
    elseif (r .lt. r1_fe) then
        res = exp(    &
        b_thr_fe(0) +       &
        b_thr_fe(1)*r +     &
        b_thr_fe(2)*r*r +   &
        b_thr_fe(3)*r*r*r )
    else
        do i=2,14
            res = res + a_thr_fe(i) * hvs(r_thr_fe(i)-r) * (r_thr_fe(i)-r)**3
        enddo
    endif

    pw_fefe_an = res
endfunction pw_fefe_an

elemental real(8) function ed_fefe_an(r)
    use energetic_parameters_mod
    implicit none
    real(8), intent(in) :: r
    real(8) res
    integer i
    interface
        pure real(8) function hvs(x)
            real(8), intent(in) :: x
        endfunction hvs
    endinterface
    res = 0d0
    do i=1,3
        res = res + a_thr_fe_ed(i) * hvs(r_thr_fe_ed(i)-r) * (r_thr_fe_ed(i)-r)**3
    enddo
    ed_fefe_an = res
endfunction ed_fefe_an

elemental real(8) function mf_fe_an(rho)
    use energetic_parameters_mod
    implicit none
    real(8), intent(in) :: rho
    mf_fe_an = -sqrt(rho) + a_fe_mf * rho * rho
endfunction mf_fe_an

elemental real(8) function pw_fefe_la(r)
    use energetic_linappr_mod
    implicit none
    real(8), intent(in) :: r
    integer intdist
    intdist = floor(r * pw_recistep)
    if( intdist .gt. pw_steps ) then
        pw_fefe_la = 0d0
        return
    elseif( intdist .lt. 1 ) then
        pw_fefe_la = pw_pot_val(1) + pw_ad1_val(1) * (r - intdist * pw_step)
        return
    endif
    pw_fefe_la = pw_pot_val(intdist) + pw_ad1_val(intdist) * (r - intdist * pw_step)
endfunction pw_fefe_la

elemental real(8) function ed_fefe_la(r)
    use energetic_linappr_mod
    implicit none
    real(8), intent(in) :: r
    integer intdist
    intdist = floor(r * ed_recistep)
    if( intdist .gt. ed_steps ) then
        ed_fefe_la = 0d0
        return
    elseif( intdist .lt. 1 ) then
        return
        ed_fefe_la = ed_pot_val(1) + ed_ad1_val(1) * (r - intdist * ed_step)
    endif
    ed_fefe_la = ed_pot_val(intdist) + ed_ad1_val(intdist) * (r - intdist * ed_step)
endfunction ed_fefe_la

elemental real(8) function mf_fe_la(rho)
    use energetic_linappr_mod
    implicit none
    real(8), intent(in) :: rho
    interface
        elemental real(8) function mf_fe_an(rho)
            real(8), intent(in) :: rho
        endfunction mf_fe_an
    endinterface
    integer intrho
    intrho = floor(rho * mf_recistep)
    if( intrho .gt. mf_steps ) then
        mf_fe_la = mf_fe_an(rho)
        return
    elseif( intrho .le. 10 ) then !threshold is negotiable
        mf_fe_la = mf_fe_an(rho)
        return
    endif
    mf_fe_la = mf_pot_val(intrho) + mf_ad1_val(intrho) * (rho - intrho * mf_step)
endfunction mf_fe_la

