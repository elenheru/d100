PROGRAM main
    use chemical_elements_names_mod
    use interaction_mod, only : incell_atoms_max
    implicit none
    !integer i,ri
    !real(8) bi,ri
    !integer getrandom_int_fromto_inclusive


    call random_seed()

    print *, "Eta programma budet vycxislqatq strukturu jadra dislokacii …"
    print *, "s vektorom bqurgersa po napravleniiju 100 v OCK zxeleze."
    print *, ""

    print *, "Sejcxas programma budet vyvoditq seriju fajlov v formate xyz dlqa ovito …"
    print *, "cxtoby protestirovatq algoritm razbijenija rascxqotnoj jacxejki na rombododekaedry."
    print *, ""

    print *, "TODO : V procedure polucxajuhxej spiski atomov"
    print *, ""

    !call test_quadratic_minimizer
    !call test_randomwalk_minimizer
    !STOP " TEST"
    call po_wo_parameters
    call spawn_bcc_rectangular_100
    call initiate_energetic_parameters
    !CALL TEST_CELL_PARTITION
    !CALL TEST_CELL_NEIBORS_4D
    !CALL TEST_CELL_GAP
    !CALL TEST_PO_DEPENDENCIES
    !call test_big_array
    !CALL TEST_FIND_NEIBORS_COMPARISON
    !CALL TEST_FIND_NEIBORS_COMPARISON_2
    !call test_system_clock
    call part_and_sort_three_zones_100
    !call find_symmetric_pairs_100
    CALL TEST_ENERGY_ATOM_ORIGIN

    call wo_xyz_whole_system
ENDPROGRAM

!
!модуль списков

!массив an_to_cn(1:N,1:3)
!сопоставляет номеру атома три индекса ячейки в которой лежит атом
!an_to_cn(ia,kс) хранит k-ый индекс ячейки

!массив cn_to_an(-300:300,-300:300,-300:300,0:100)
!сопоставляет тройке индексов ячейки количество атомов в ней и номера атомов
!cn_to_an(kx,ky,kz,0) хранит количество атомов в ячейке с индексами kx ky kz
!cn_to_an(kx,ky,kz,i) хранит номер первого атома в ячейке с индексами kx ky kz

!конец модуля списков

!fill_an_to_cn
!процедура заполнения an_to_cn
!нужно разложить радиус-вектор атома по базису гцк решётки
!дальше выполнить nint от каждой координаты в этом базисе
!может быть стоит

!fill_cn_to_an
!процедура заполнения cn_to_an
!можно проходить цикл по атомам для каждой ячейки и смотреть геометрические условия
!условие в том, чтобы атом был ближе к центру ячейки чем к центрам любой из соседних ячеек
!но это будет много холостой работы - большинство атомов не будут лежать в произвольно взятой ячейке
!лучше пройти по всем атомам, и записывать в соответствующую ему ячейку (cn_to_an) нужную информацию

!процедура get_direct_list(i,)
!получения списка соседей для атома для расчета энергии или силы

