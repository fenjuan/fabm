#include "fabm_driver.h"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
module fish_interface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    use fabm_types
    use fabm_standard_variables
    use fabm_particle
    use fabm_expressions
    use fabm_builtin_models
    implicit none
    private
    
    !===============================================================================================================
    ! derived types:
    !===============================================================================================================
    ! type 1: get_pvar. Loops through space and saves pelagic variables as diagnostic variables that can be accessed by horizontal fish_cohort
    type,extends(type_particle_model),public :: type_get_pvar
       !variables
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_pvar(:,:)                          ! for saving local variable values, holds a number of local values equal to number of depth cells x variables
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_thickness(:)                       ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_temperature(:)                     ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_area(:)                            ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_volume(:)                          ! for saving local cell thicknesses, holds a number of local values equal to number of depth cells
       
       type (type_state_variable_id), allocatable                    :: id_targets_in(:)                            ! the target variables (e.g. zooplankton,light,oxygen) (NOTE: NOT the same as in fish_rate_distributor!)
       type (type_dependency_id)                                     :: id_thickness                                ! cell thickness; this will be needed for calculating/allocating fish predation
       type (type_dependency_id)                                     :: id_uTm
       type (type_dependency_id)                                     :: id_volume
       type (type_dependency_id)                                     :: id_area
       !parameters
       integer                                                       :: nlev, nvars                                 ! number of layers in host model setup, and number of variables to get 
    contains
       procedure :: initialize => get_pvar_initialize
       procedure :: get_light  => get_pvar_do ! currently only get_light supports non-local action in depth
    end type
    !-------------------------------------------------------------------------------------------------------
    
    ! type 2: fish_rate_distributer. Loops through space and aplies local rates calculated by fish_cohort (and saved as horizontal diagnostics) to local variables
    type,extends(type_particle_model),public  :: type_fish_rate_dist
       type (type_state_variable_id)                                 :: id_DPOMpoolW, id_PPOMpoolW, id_NPOMpoolW
       type (type_state_variable_id)                                 :: id_NH4poolW, id_PO4poolW, id_DDOMpoolW
       type (type_state_variable_id)                                 :: id_PDOMpoolW,id_NDOMpoolW
       type (type_diagnostic_variable_id)                 :: id_BOT_FEED
       !type (type_horizontal_diagnostic_variable_id)                 :: id_BOT_FEED

       type (type_model_id)                                          :: id_input
       type (type_state_variable_id), allocatable                    :: id_targets_preyDW(:), id_targets_preyN(:)   ! prey variable ids that should absorp the sources-sinks (NOTE: NOT the same as in get_pvar!)
       type (type_state_variable_id), allocatable                    :: id_targets_preyP(:)
       type (type_horizontal_dependency_id), allocatable             :: id_pas_rates(:,:)                           ! total depth integrated rates of change for the target passive variables, as calculated by fish_cohort
       type (type_horizontal_dependency_id), allocatable             :: id_prey_loss_DW(:,:), id_prey_loss_N(:,:)   ! depth specific prey loss rates (as g m**-2)
       type (type_horizontal_dependency_id), allocatable             :: id_prey_loss_P(:,:)                         ! depth specific prey loss rates (as g m**-2)
       
       integer                                                       :: nlev, nprey                                 ! number of layers in host model setup, number of passive and prey variables to get, and number of cohorts
       
    contains
       procedure :: initialize => fish_rate_dist_initialize
       !procedure :: get_light  => fish_rate_dist_do ! currently only get_light supports non-local action in depth
       procedure :: do  => fish_rate_dist_do
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
        call self%register_dependency(self%id_volume,standard_variables%cell_volume)
        call self%register_dependency(self%id_area,standard_variables%cell_area)
        
        ! this should be allocated:
        allocate(self%id_local_pvar(self%nlev,self%nvars))
        allocate(self%id_local_thickness(self%nlev))
        allocate(self%id_local_area(self%nlev))
        allocate(self%id_local_volume(self%nlev))
        allocate(self%id_local_temperature(self%nlev))
        allocate(self%id_targets_in(self%nvars))
        do i_z=1,self%nlev
            write(index_z,'(i0)') i_z
            call self%register_diagnostic_variable(self%id_local_thickness(i_z),'dz_z'//trim(index_z),'','local cell thickness at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none,source=source_do_column,missing_value=0.0_rk)                                        ! current cell thickness
            call self%register_diagnostic_variable(self%id_local_area(i_z),'Af_z'//trim(index_z),'','local cell area at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none,source=source_do_column,missing_value=0.0_rk)                                        ! current cell thickness
            call self%register_diagnostic_variable(self%id_local_volume(i_z),'Vc_z'//trim(index_z),'','local cell volume at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none,source=source_do_column,missing_value=0.0_rk)                                        ! current cell thickness
            
            call self%register_diagnostic_variable(self%id_local_temperature(i_z),'uTm_z'//trim(index_z),'','local cell temperature at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none,source=source_do_column,missing_value=0.0_rk)                                        ! current cell thickness
            do i_v=1,self%nvars
                 write(index_v,'(i0)') i_v
                 if (i_z==1) then
                      call self%register_state_dependency(self%id_targets_in(i_v), 'target'//trim(index_v), '', 'saved variable nr '//trim(index_v))                          ! the target variable (e.g zooplankton)
                 endif
                 call self%register_diagnostic_variable(self%id_local_pvar(i_z,i_v),'pvar'//trim(index_v)//'_at_z'//trim(index_z),'','local value of var'//trim(index_v)//' at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,domain=domain_surface,output=output_none,source=source_do_column,missing_value=0.0_rk)
            end do
        end do
    end subroutine get_pvar_initialize  
    !===============================================================================================================
    ! do part
    !===============================================================================================================
    subroutine get_pvar_do(self,_ARGUMENTS_VERTICAL_)
        class (type_get_pvar),intent(in) :: self
        _DECLARE_ARGUMENTS_VERTICAL_

        real(rk)                                    :: Vc, Af
        real(rk), dimension(self%nlev)              :: thickness, uTm
        real(rk), dimension(self%nlev,self%nvars)   :: locals
        integer  :: i_z, i_v

        i_z=0
        
        _VERTICAL_LOOP_BEGIN_
            i_z = i_z +1
            
            do i_v=1,self%nvars
               _GET_(self%id_targets_in(i_v),locals(i_z,i_v))
            end do
            _GET_(self%id_thickness,thickness(i_z))
            _GET_(self%id_uTm,uTm(i_z))
            _GET_(self%id_volume,Vc)
            _GET_(self%id_area,Af)
!            write (*,*) CA
        _VERTICAL_LOOP_END_
!        stop
        ! calculate (relative) benthic area?
        ! calculate relative volume?
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_thickness,thickness)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_area,Af)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_volume,Vc)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_temperature,uTm)
        do i_v=1,self%nvars
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_pvar(:,i_v),locals(:,i_v))
        end do           
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
        
        call self%register_diagnostic_variable(self%id_BOT_FEED,'bottom_activity','g m-2','bottom activity', & ! not sure about act_as_state_variable
                          act_as_state_variable=.true.,output=output_none,missing_value=0.0_rk)
        
        call self%register_state_dependency(self%id_DPOMpoolW,    'POM_DW_pool_water',     'g m-3', 'POM DW pool in water')
        call self%register_state_dependency(self%id_NPOMpoolW,    'POM_N_pool_water',      'g m-3', 'POM N pool in water')
        call self%register_state_dependency(self%id_PPOMpoolW,    'POM_P_pool_water',      'g m-3', 'POM P pool in water')
        call self%register_state_dependency(self%id_NH4poolW,     'NH4_pool_water',        'g m-3', 'NH4 pool in water')
        call self%register_state_dependency(self%id_PO4poolW,     'PO4_pool_water',        'g m-3', 'PO4 pool in water')
        call self%register_state_dependency(self%id_DDOMpoolW,    'DOM_DW_pool_water',     'g m-3', 'DOM DW in water')
        call self%register_state_dependency(self%id_NDOMpoolW,    'DOM_N_pool_water',      'g m-3', 'DOM N in water')
        call self%register_state_dependency(self%id_PDOMpoolW,    'DOM_P_pool_water',      'g m-3', 'DOM P in water')
        
        allocate(self%id_targets_preyDW(self%nprey))
        allocate(self%id_targets_preyN(self%nprey))
        allocate(self%id_targets_preyP(self%nprey))
        allocate(self%id_pas_rates(self%nlev,9))
        allocate(self%id_prey_loss_DW(self%nlev,self%nprey))
        allocate(self%id_prey_loss_N(self%nlev,self%nprey))
        allocate(self%id_prey_loss_P(self%nlev,self%nprey))
        
        call self%register_model_dependency(self%id_input,'fish')
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
                
                                    
                call self%request_coupling_to_model(self%id_prey_loss_DW(i_z,i_prey),self%id_input,'PL_DW_z'//trim(index_z))
                call self%request_coupling_to_model(self%id_prey_loss_N(i_z,i_prey),self%id_input,'PL_N_z'//trim(index_z))
                call self%request_coupling_to_model(self%id_prey_loss_P(i_z,i_prey),self%id_input,'PL_P_z'//trim(index_z))

            end do
            do i_pas=1,9
                write(index_pas,'(i0)') i_pas
                call self%register_horizontal_dependency(self%id_pas_rates(i_z,i_pas), 'rate_of_passive_target'//trim(index_pas)//'_at_z'//trim(index_z), '-', 'rate of change of passive target variable number '//trim(index_pas)//' at depth '//trim(index_z))
                call self%request_coupling_to_model(self%id_pas_rates(i_z,i_pas),self%id_input,'PR'//trim(index_pas)//'_z'//trim(index_z))
            end do
        end do        
    end subroutine fish_rate_dist_initialize

    !===============================================================================================================
    ! do part
    !===============================================================================================================
    subroutine fish_rate_dist_do(self,_ARGUMENTS_DO_)
        class (type_fish_rate_dist),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        real(rk), dimension(self%nlev,self%nprey)  :: lossDW, lossN, lossP
        real(rk), dimension(self%nlev,9)           :: rate
        integer                                    :: i_pas, i_prey, i_z
        
        i_z=0
        
        !!! GET HORIZONTAL DIAGNOSTICS BEFORE LOOP !!!
        do i_pas=1,9
              _GET_HORIZONTAL_(self%id_pas_rates(:,i_pas),rate(:,i_pas))
        end do
        do i_prey=1,self%nprey
              _GET_HORIZONTAL_(self%id_prey_loss_DW(:,i_prey),lossDW(:,i_prey))
              _GET_HORIZONTAL_(self%id_prey_loss_N(:,i_prey),lossN(:,i_prey))
              _GET_HORIZONTAL_(self%id_prey_loss_P(:,i_prey),lossP(:,i_prey))
        end do
        !write (*,*) lossDW(1,1)
        _LOOP_BEGIN_
            i_z=i_z+1
           ! First update passive target variables.
           _SET_ODE_(self%id_NH4poolW, rate(i_z,1))
           _SET_ODE_(self%id_PO4poolW, rate(i_z,2))
           _SET_ODE_(self%id_DPOMpoolW, rate(i_z,3))
           _SET_ODE_(self%id_DDOMpoolW, rate(i_z,4))
           _SET_ODE_(self%id_PPOMpoolW, rate(i_z,5))
           _SET_ODE_(self%id_PDOMpoolW, rate(i_z,6))
           _SET_ODE_(self%id_NPOMpoolW, rate(i_z,7))
           _SET_ODE_(self%id_NDOMpoolW, rate(i_z,8))
           _SET_DIAGNOSTIC_(self%id_BOT_FEED, rate(i_z,9))
           !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_BOT_FEED, rate(i_z,9))
           
           ! then update prey variables
           do i_prey=1,self%nprey
                _SET_ODE_(self%id_targets_preyDW(i_prey),lossDW(i_z,i_prey))
                _SET_ODE_(self%id_targets_preyN(i_prey),lossN(i_z,i_prey))
                _SET_ODE_(self%id_targets_preyP(i_prey),lossP(i_z,i_prey))
           end do
           
        _LOOP_END_
    end subroutine fish_rate_dist_do
end module