module array_parameters_mod
    save
    integer, parameter  ::  x_layers= 67,y_layers= 67,z_layers= 6
    integer, parameter  ::  atoms_max_array = nint(&
        (2*x_layers+1)*(2*y_layers+1)*(2*z_layers+1)*21d-1)
    integer             ::  atoms__in_total
    integer, parameter  ::  cells_xrange=5*x_layers/3
    integer, parameter  ::  cells_yrange=5*y_layers/3
    integer, parameter  ::  cells_zrange=5*z_layers/3

endmodule array_parameters_mod

module      positions_mod
    use array_parameters_mod
    save
    real(8), dimension(1:3,atoms_max_array) ::  R_perf,R_curr
    !vse koordinaty v angstremah
endmodule   positions_mod

module phys_parameters_mod
    !use comp_parameters_mod
    save
    real(8), parameter  ::  a0 = 2.8600d0
    real(8), parameter  ::  aepsilon = 2.8600d-7 !-7 !to avoid match with cells centers positions
    real(8), parameter  ::  poisson=369d-3
    real(8), parameter  ::  burgers=-a0*50d-2*2d0!dislocation is doubled
    !real(8), parameter  ::  core_sign= 1d0
    !core is compressed if sign=1 or decomressed if sign=-1
endmodule phys_parameters_mod

module      interaction_mod
    use array_parameters_mod
    use phys_parameters_mod
    save
    real(8), parameter ::   cutoff = 6.0d0 !sphere where the neighbor is
    real(8)            ::   cutoff_param = cutoff
    real(8), parameter ::   incell_atoms_exceed = 2650d-3*2
    integer, parameter ::   incell_atoms_max = &
        nint( (16d0*sqrt(3d0)/9d0*(cutoff*sqrt(375d-3))**3)/(5d-1*(a0**3))*incell_atoms_exceed )
!    integer an_by_cn(3,atoms_max_array)
!    integer cn_by_an(0:incell_atoms_max,-cells_xrange:cells_xrange,-cells_yrange:cells_yrange,-cells_zrange:cells_zrange)
    integer an_by_cn_4d(4,atoms_max_array)
    integer cn_by_an_4d(0:incell_atoms_max,-cells_xrange:cells_xrange,-cells_yrange:cells_yrange,-cells_zrange:cells_zrange,0:3)
    !real(8) distan_list(  neibors_max,atoms_max_array)

    logical(1), parameter, dimension(6) :: &
    xp=(/ .true.  , .false. , .true.  , .true.  , .true.  , .true.  /),&
    xn=(/ .false. , .true.  , .true.  , .true.  , .false. , .false. /),&
    yp=(/ .true.  , .true.  , .true.  , .false. , .true.  , .true.  /),&
    yn=(/ .false. , .false. , .false. , .true.  , .true.  , .true.  /),&
    zp=(/ .true.  , .true.  , .true.  , .true.  , .true.  , .false. /),&
    zn=(/ .true.  , .true.  , .false. , .false. , .false. , .true.  /)
    !      x+y>0  ; -x+y>0  ;  y+z>0  ; -y+z>0  ;  z+x>0  ; -z+x>0

endmodule   interaction_mod

module chemical_elements_names_mod
    save
    character(LEN=236), parameter :: elements_names = &
    "H He"   // &
    "LiBeB C N O F Ne"   // &
    "NaMgAlSiP S ClAr"   // &
    "K Ca"   // &
    "ScTiV CrMnFeCoNiCuZn"   // &
    "GaGeAsSeBrKr" // &
    "RbSr" // &
    "Y ZrNbMoTcRuRhPdAgCd"   // &
    "InSnSbTeI Xe"   // &
    "CsBa" // &
    "LaCePrNdPmSmEuGdTbDyHoErTmYbLu"   // &
    "HfTaW ReOsIrPtAuHg"   // &
    "TlPbBiPoAtRn" // &
    "FrRa" // &
    "AcThPaU NpPuAmCmBkCfEsFmMdNoLr"   // &
    "RfDbSgBhHsMtDsRgCn"   // &
    "NhFlMcLvTsOg"

endmodule chemical_elements_names_mod

module wo_counters_mod
    save
    integer :: xyz__WOcounter=0
endmodule wo_counters_mod

module energetic_parameters_mod
    save
    real(8), parameter :: r_bohr = 0.52917721067d0
    !real(8), parameter :: electron_charge = 1d0!although ackland gives no value in 2003 article, 1d0 sews phi smoothly ? no
    real(8), parameter :: electron_charge = 1.644429636713d0 !this sews values, but not derivatives
    !real(8), parameter :: electron_charge = sqrt(1.4399764*1d1) !value in sqrt( eV*Å ) from wiki
    real(8), parameter :: z_fe = 26d0
    real(8), parameter :: r_s_fe = 0.88534d0*sqrt(5d-1)*r_bohr/(z_fe**(1d0/3d0))
    real(8), dimension(0:3), parameter :: &
        b_thr_fe = (/&
        6.4265260576348d0,&
        1.7900488524286d0,&
        -4.5108316729807d0,&
        1.0866199373306d0 /)
    real(8), dimension(14), parameter :: &
        a_thr_fe = (/&
        0d0,&
        -24.028204854115d0,&
        11.300691696477d0,&
        5.3144495820462d0,&
        -4.6659532856049d0,&
        5.9637758529194d0,&
        1.7710262006061d0,&
        0.85913830768731d0,&
        -2.1845362968261d0,&
        2.6424377007466d0,&
        -1.0358345370208d0,&
        0.33548264951582d0,&
        -0.046448582149334d0,&
        -0.0070294963048689d0/)
    real(8), dimension(14), parameter :: &
        r_thr_fe = (/&
        2.1d0,&
        2.2d0,&
        2.3d0,&
        2.4d0,&
        2.5d0,&
        2.6d0,&
        2.7d0,&
        2.8d0,&
        3.0d0,&
        3.3d0,&
        3.7d0,&
        4.2d0,&
        4.7d0,&
        5.3d0/)
    real(8), parameter :: r1_fe = 1d0,r2_fe = 2d0

    real(8), dimension(3), parameter :: &
        r_thr_fe_ed = (/&
        2.4d0,&
        3.2d0,&
        4.2d0/)
    real(8), dimension(3), parameter :: &
        a_thr_fe_ed = (/&
        11.686859407970d0,&
        -0.014710740098830d0,&
        0.47193527075943d0/)

    real(8), parameter :: a_fe_mf = -0.35387096579929d-3

endmodule energetic_parameters_mod

module energetic_linappr_mod
    implicit none
    save
    real(8), parameter ::   pw_domain = 5.35d0!from 0.0 A to 5.35 A. cutoff otherwise
    integer, parameter ::   pw_steps = 50000

    real(8), parameter ::   ed_domain = 4.25d0!from 0.0 A to 4.25 A. cutoff otherwise
    integer, parameter ::   ed_steps = 50000

    real(8), parameter ::   mf_domain = 650d0!from 0 to 650 charge units. analytical otherwise
    integer, parameter ::   mf_steps = 50000
    !for function   itself we only need its actual value
    !for first  derivative we both need its actual value and it's average value on step interval
    !for second derivative we only need it's average value on step interval
    real(8), parameter ::   pw_step = pw_domain/real(pw_steps), pw_recistep = 1d0/pw_step
    real(8), dimension(pw_steps) :: pw_pot_val, pw_ad1_val
    !pairwise potential, reciprocal argument step, first derivative average
    real(8), parameter ::   ed_step = ed_domain/real(ed_steps), ed_recistep = 1d0/ed_step
    real(8), dimension(ed_steps) :: ed_pot_val, ed_ad1_val
    !electronic density, reciprocal argument step, first derivative average
    real(8), parameter ::   mf_step = mf_domain/real(mf_steps), mf_recistep = 1d0/mf_step
    real(8), dimension(mf_steps) :: mf_pot_val, mf_ad1_val
    !embedding function, reciprocal argument step, first derivative average
endmodule energetic_linappr_mod

!module      potentials_fefe_linear_mod
!
!    !contains
!!    elemental real(8) function pairwise_deriv1(r)
!!    !pairwise potential using  initiated massive
!!        !use energetic_parameters_mod
!!        implicit none
!!        real(8),intent(in)::r;real(8) dist
!!        integer intdistance
!!        pairwise_deriv1= 0d0
!!        dist = r*pairwise_recistep
!!        intdistance = floor(dist)
!!        if (intdistance .ge. pot_steps) return
!!        if (intdistance .le. 0        ) return
!!        pairwise_deriv1 =  pairwise_dr1(intdistance) + pairwise_ad2(intdistance) * &
!!                    (r - pairwise_step*intdistance)
!!    endfunction pairwise_deriv1
!
!!    elemental real(8) function edensity_deriv1(r)
!!        !electronic density using  initiated massive
!!        !use energetic_parameters_mod
!!        implicit none
!!        real(8),intent(in)::r;real(8) dist
!!        integer intdistance
!!        edensity_deriv1 = 0d0
!!            !RETURN!to check if III Newton's law is correct here
!!            !surprisingly it is incorrect.
!!        dist = r*elecdens_recistep
!!        intdistance = floor(dist)
!!        if (intdistance .ge. pot_steps) return
!!        if (intdistance .le. 0        ) return
!!        edensity_deriv1 =  elecdens_dr1(intdistance) + elecdens_ad2(intdistance) * &
!!                    (r - elecdens_step*intdistance)
!!    endfunction edensity_deriv1
!
!!    elemental real(8) function embeddin_deriv1(rho)
!!        !embedding function using  initiated massive
!!        !of course argument is not really a distance
!!        !use energetic_parameters_mod
!!        implicit none
!!        real(8),intent(in)::rho;real(8) dist
!!        integer intdistance
!!        embeddin_deriv1 = 0d0
!!        if(rho .lt. 0d0) return
!!        dist = rho*embefunc_recistep
!!        intdistance = floor(dist)
!!        if ( (intdistance .ge. pot_steps) .or. (intdistance .lt. 10) ) then
!!            embeddin_deriv1 = -1d9 !for zero rho we do not want enormous derivative
!!            if(rho .gt. 1d-7) embeddin_deriv1 = -dsqrt(25d-2/rho)
!!            return  !square root is very fast changing when close to 0.0
!!        endif
!!        embeddin_deriv1 =  embefunc_dr1(intdistance) + embefunc_ad2(intdistance) * &
!!                    (rho - embefunc_step*intdistance)
!!    endfunction embeddin_deriv1
!
!endmodule   potentials_fefe_linear_mod

module lists_assortiment_mod
    use interaction_mod
    save
    integer, dimension(incell_atoms_max*12) :: list_atoms !atom, sosedi kotorogo trebyjetsqa najtimodule
    integer :: list_atoms_qnt
    integer, dimension(4,42) :: list_cells !
endmodule lists_assortiment_mod
