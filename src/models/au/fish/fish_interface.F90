#include "fabm_driver.h"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
module fish_interface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    use fabm_types
    use fabm_standard_variables
    implicit none
    private
    
    !===============================================================================================================
    ! derived types:
    !===============================================================================================================
    ! type 1: get_pvar. Loops through space and saves pelagic variables as diagnostic variables that can be accessed by horizontal fish_cohort
    type,extends(type_base_model),public :: type_get_pvar
       !variables
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_pvar(:,:)                          ! for saving local variable values, holds a number of local values equal to number of depth cells x variables
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_thickness(:)                       ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_temperature(:)                       ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       
       type (type_state_variable_id), allocatable                    :: id_targets_in(:)                            ! the target variables (e.g. zooplankton,light,oxygen) (NOTE: NOT the same as in fish_rate_distributor!)
       type (type_dependency_id)                                     :: id_thickness                                ! cell thickness; this will be needed for calculating/allocating fish predation
       type (type_dependency_id)                                     :: id_uTm
       !parameters
       integer                                                       :: nlev, nvars                                 ! number of layers in host model setup, and number of variables to get 
    contains
       procedure :: initialize => get_pvar_initialize
       procedure :: get_light  => get_pvar_do ! currently only get_light supports non-local action in depth
    end type
    !-------------------------------------------------------------------------------------------------------
    
    ! type 2: fish_rate_distributer. Loops through space and aplies local rates calculated by fish_cohort (and saved as horizontal diagnostics) to local variables
    type,extends(type_base_model),public  :: type_fish_rate_dist
       type (type_state_variable_id)                                 :: id_DPOMpoolW, id_PPOMpoolW, id_NPOMpoolW
       type (type_state_variable_id)                                 :: id_NH4poolW, id_PO4poolW, id_DDOMpoolW
       type (type_state_variable_id)                                 :: id_PDOMpoolW,id_NDOMpoolW
       type (type_diagnostic_variable_id)                            :: id_BOT_FEED

       type (type_state_variable_id), allocatable                    :: id_targets_preyDW(:), id_targets_preyN(:)   ! prey variable ids that should absorp the sources-sinks (NOTE: NOT the same as in get_pvar!)
       type (type_state_variable_id), allocatable                    :: id_targets_preyP(:)
       type (type_horizontal_dependency_id), allocatable             :: id_pas_rates(:,:)                           ! total depth integrated rates of change for the target passive variables, as calculated by fish_cohort
       type (type_horizontal_dependency_id), allocatable             :: id_prey_loss_DW(:,:), id_prey_loss_N(:,:)   ! depth specific prey loss rates (as g m**-2)
       type (type_horizontal_dependency_id), allocatable             :: id_prey_loss_P(:,:)                         ! depth specific prey loss rates (as g m**-2)
       
       integer                                                       :: nlev, nprey                           ! number of layers in host model setup, number of passive and prey variables to get, and number of cohorts 
    contains
       procedure :: initialize => fish_rate_dist_initialize
       procedure :: do         => fish_rate_dist_do ! currently only get_light supports non-local action in depth
    end type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
contains
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! get_pvar code
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine get_pvar_initialize(self,configunit)
        class (type_get_pvar),intent(inout),target :: self
        integer,                    intent(in)     :: configunit
        
        integer                                    :: i_z, i_v
        character(len=64)   :: index_z, index_v

        
        call self%get_parameter(self%nlev,'nlev', '-', 'Number of depth levels in water column (same as in host model)')       ! number of depth levels in simulation
        call self%get_parameter(self%nvars,'nvars', '-', 'Number of variables to get')                                         ! number of depth levels in simulation
        
        call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)                                        ! current cell thickness
        call self%register_dependency(self%id_uTm,standard_variables%temperature)                                        ! current cell thickness
        
        ! this should be allocated:
        allocate(self%id_local_pvar(self%nlev,self%nvars))
        allocate(self%id_local_thickness(self%nlev))
        allocate(self%id_targets_in(self%nvars))
        do i_z=1,self%nlev
            write(index_z,'(i0)') i_z
            call self%register_diagnostic_variable(self%id_local_thickness(i_z),'dz_at_z'//trim(index_z),'','local cell thickness at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none)                                        ! current cell thickness
            call self%register_diagnostic_variable(self%id_local_temperature(i_z),'uTm_at_z'//trim(index_z),'','local cell temperature at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none)                                        ! current cell thickness
            do i_v=1,self%nvars
                 write(index_v,'(i0)') i_v
                 if (i_z==1) then
                      call self%register_state_dependency(self%id_targets_in(i_v), 'target'//trim(index_v), '', 'saved variable nr '//trim(index_v))                          ! the target variable (e.g zooplankton)
                 endif
                 call self%register_diagnostic_variable(self%id_local_pvar(i_z,i_v),'pvar'//trim(index_v)//'_at_z'//trim(index_z),'','local value of var'//trim(index_v)//' at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none)
            end do
        end do
    end subroutine get_pvar_initialize  
    !===============================================================================================================
    ! do part
    !===============================================================================================================
    subroutine get_pvar_do(self,_ARGUMENTS_VERTICAL_)
        class (type_get_pvar),intent(in) :: self
        _DECLARE_ARGUMENTS_VERTICAL_

        real(rk) :: local,thickness
        integer  :: i_z, i_v

        i_z=1
        
        _VERTICAL_LOOP_BEGIN_
        do i_v=1,self%nvars
           _GET_(self%id_targets_in(i_v),local)
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_pvar(i_z,i_v),local)
        end do
        _GET_(self%id_thickness,thickness)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_thickness(i_z),thickness)
        
        i_z=i_z+1
        _VERTICAL_LOOP_END_
    end subroutine get_pvar_do
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! fish_rate_distributor code
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    subroutine fish_rate_dist_initialize(self,configunit)
        class (type_fish_rate_dist),intent(inout),target :: self
        integer,                                     intent(in)           :: configunit

        integer                                    :: i_z, i_pas, i_prey
        character(len=64)   :: index_z, index_pas, index_prey
        
        call self%get_parameter(self%nlev,'nlev', '-', 'Number of depth levels in water column (same as in host model)')       ! number of depth levels in simulation
        call self%get_parameter(self%nprey,'nprey', '-', 'Number of target prey variables')
        call self%register_diagnostic_variable(self%id_BOT_FEED,'bottom_activity','g m-3 ','bottom activity', output=output_none)

        allocate(self%id_targets_preyDW(self%nprey))
        allocate(self%id_targets_preyN(self%nprey))
        allocate(self%id_targets_preyP(self%nprey))
        allocate(self%id_pas_rates(self%nlev,9))
        allocate(self%id_prey_loss_DW(self%nlev,self%nprey))
        allocate(self%id_prey_loss_N(self%nlev,self%nprey))
        allocate(self%id_prey_loss_P(self%nlev,self%nprey))
        do i_z=1,self%nlev
           write(index_z,'(i0)') i_z 
           do i_prey=1,self%nprey
                write(index_prey,'(i0)') i_prey
                if (i_z==1) then
                    call self%register_state_dependency(self%id_targets_preyDW(i_prey), 'prey_target'//trim(index_prey)//'_DW', '-', 'prey DW target variable number '//trim(index_prey))
                    call self%register_state_dependency(self%id_targets_preyN(i_prey), 'prey_target'//trim(index_prey)//'_N', '-', 'prey N target variable number '//trim(index_prey))
                    call self%register_state_dependency(self%id_targets_preyP(i_prey), 'prey_target'//trim(index_prey)//'_P', '-', 'prey P target variable number '//trim(index_prey))
                endif
                call self%register_horizontal_dependency(self%id_prey_loss_DW(i_z,i_prey), 'DW_loss_rate_of_prey'//trim(index_prey)//'_at_z'//trim(index_z), '-', 'DW loss rate (m**-2) of prey number '//trim(index_prey)//'_at_z'//trim(index_z))
                call self%register_horizontal_dependency(self%id_prey_loss_N(i_z,i_prey), 'N_loss_rate_of_prey'//trim(index_prey)//'_at_z'//trim(index_z), '-', 'N loss rate (m**-2) of prey number '//trim(index_prey)//'_at_z'//trim(index_z))
                call self%register_horizontal_dependency(self%id_prey_loss_P(i_z,i_prey), 'P_loss_rate_of_prey'//trim(index_prey)//'_at_z'//trim(index_z), '-', 'P loss rate (m**-2) of prey number '//trim(index_prey)//'_at_z'//trim(index_z))
            end do
            do i_pas=1,9
                write(index_pas,'(i0)') i_pas
                call self%register_horizontal_dependency(self%id_pas_rates(i_z,i_pas), 'rate_of_passive_target'//trim(index_pas)//'_at_z'//trim(index_z), '-', 'rate of change of passive target variable number '//trim(index_pas)//' at depth '//trim(index_z))
            end do
        end do        
    end subroutine fish_rate_dist_initialize

    !===============================================================================================================
    ! do part
    !===============================================================================================================
    subroutine fish_rate_dist_do(self,_ARGUMENTS_DO_)
        class (type_fish_rate_dist),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        real(rk)                :: lossDW, lossN, lossP
        real(rk), dimension(9)  :: rate
        integer                 :: i_pas, i_prey, i_z
        
        i_z=1
        
        _LOOP_BEGIN_
           ! First update passive target variables.
           do i_pas=1,9
              _GET_HORIZONTAL_(self%id_pas_rates(i_z,i_pas),rate(i_pas))
           end do
           _SET_ODE_(self%id_NH4poolW, rate(1))
           _SET_ODE_(self%id_PO4poolW, rate(2))
           _SET_ODE_(self%id_DPOMpoolW, rate(3))
           _SET_ODE_(self%id_DDOMpoolW, rate(4))
           _SET_ODE_(self%id_PPOMpoolW, rate(5))
           _SET_ODE_(self%id_PDOMpoolW, rate(6))
           _SET_ODE_(self%id_NPOMpoolW, rate(7))
           _SET_ODE_(self%id_NDOMpoolW, rate(8))
           _SET_DIAGNOSTIC_(self%id_BOT_FEED, rate(9))
           ! then update prey variables
           do i_prey=1,self%nprey
              _GET_HORIZONTAL_(self%id_prey_loss_DW(i_z,i_prey),lossDW)
              _GET_HORIZONTAL_(self%id_prey_loss_N(i_z,i_prey),lossN)
              _GET_HORIZONTAL_(self%id_prey_loss_P(i_z,i_prey),lossP)
              
              _SET_ODE_(self%id_targets_preyDW(i_prey),lossDW)
              _SET_ODE_(self%id_targets_preyN(i_prey),lossN)
              _SET_ODE_(self%id_targets_preyP(i_prey),lossP)
           end do
           
           i_z=i_z+1
        _LOOP_END_
    end subroutine fish_rate_dist_do
end module