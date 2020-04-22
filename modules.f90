module array_parameters_mod
    save
    integer, parameter  ::  x_layers= 12,y_layers= 20,z_layers= 12
    integer, parameter  ::  atoms_max_array = &
        (2*x_layers+1)*(2*y_layers+1)*(2*z_layers+1)*3
    integer             ::  atoms__in_total
    integer, parameter  ::  cells_range=nint(15d-1*atoms_max_array**(1d0/3d0))

endmodule array_parameters_mod

module      interaction_mod
    use array_parameters_mod
    save
    real(8), parameter ::   cutoff = 3.5d0*2.5 !sphere where the neighbor is
    real(8)            ::   cutoff_param = cutoff
    integer, parameter ::   neibors_max = 150
    integer an_by_cn(3,atoms_max_array)
    integer cn_by_an(0:neibors_max,-cells_range:cells_range,-cells_range:cells_range,-cells_range:cells_range)
    !real(8) distan_list(  neibors_max,atoms_max_array)
endmodule   interaction_mod

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
    real(8), parameter  ::  aepsilon = 2.8600d-7 !to avoid match with cells centers positions
    real(8), parameter  ::  poisson=369d-3
    real(8), parameter  ::  burgers=-a0*50d-2*2d0!dislocation is doubled
    real(8), parameter  ::  core_sign= 1d0
    !core is compressed if sign=1 or decomressed if sign=-1
endmodule phys_parameters_mod

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
