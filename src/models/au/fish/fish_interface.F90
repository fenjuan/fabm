#include "fabm_driver.h"
   
module fish_interface

   use fabm_types
   use fabm_standard_variables

   implicit none

   private
    
   ! type_depth_integral: module that computes weighted depth integral of a specified variable.
   ! Weights are taken from an externally computed variable [id_target below],
   ! which typically is a diagnostic computed by a model of type type_vertical_distribution.
   ! The variable to be integrated is a dependency [id_target below] that must be coupled at run time.
   ! The resultant integral is a diagnsotic that acts like a state variable, that is to say,
   ! other models can use it as if it were a state avriable, and provide sources and sinks,
   ! and thsi model will then automagically distribute those sinks and sources again over their
   ! original depth-explicit source variable, using the appropriate weights.
   type,extends(type_base_model),public :: type_get_pvar
      !variables
      type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_pvar(:,:)         ! needs flexibility, to hold a number of local values equal to number of depth cells x variables
      type (type_horizontal_diagnostic_variable_id), allocatable    :: id_local_thickness(:)         ! needs flexibility, to hold a number of local values equal to number of depth cells
      
      type (type_state_variable_id), allocatable                    :: id_targets(:)                ! the target variables (e.g. zooplankton,light,oxygen)
      type (type_dependency_id)                                     :: id_thickness          ! this will be needed for calculating/allocating fish predation
      !parameters
      integer                                                       :: nlev, nvars         ! number of layers in host model setup, and number of variables to get 
   contains
      procedure :: initialize => get_pvar_initialize
      procedure :: get_light  => get_pvar_do ! currently only get_light supports non-local action in depth
   end type
   
   type,extends(type_base_model) :: type_fish_rate_distributor
      type (type_state_variable_id), allocatable                    :: id_targets(:)     ! Depth-explicit variables that should absorp the sources-sinks
      type (type_horizontal_dependency_id), allocatable             :: id_totals(:)              ! total depth integrated rates of change for the target variables, as calculated by fish_cohort
      type (type_horizontal_dependency_id), allocatable             :: id_fish_dist(:)    ! Weights for the vertical distribution of the sinks and sources, as calculated by fish_cohort
      type (type_horizontal_dependency_id), allocatable             :: id_weight_order(:)    ! how many times the weight should factor in on the final rates of the target variables
      
      integer                                                       :: nlev, nvars         ! number of layers in host model setup, and number of variables to get 
   contains
      procedure :: initialize => fish_rate_distributor_initialize
      procedure :: do         => fish_rate_distributor_do ! currently only get_light supports non-local action in depth
   end type
   
    contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! get_pvar code
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine get_pvar_initialize(self,configunit)
      class (type_get_pvar),intent(inout),target :: self
      integer,                    intent(in)     :: configunit
      
      integer                                    :: i_z, i_v    
      
      call self%get_parameter(self%nlev,'nlev', '-', 'Number of depth levels in water column (same as in host model)')       ! number of depth levels in simulation
      call self%get_parameter(self%nvars,'nvars', '-', 'Number of variables to get')                                         ! number of depth levels in simulation
      
      call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)                                        ! current cell thickness
      
      ! this should be allocated:
      allocate(self%id_local_pvar(self%nlev,self%nvars))
      allocate(self%id_local_thickness(self%nlev))
      allocate(self%id_targets(self%nvars))
      do i_z=1,self%nlev
          write(index_z,'(i0)') i_z
          call self%register_diagnostic_variable(self%id_local_thickness(i_z),'dz_at_z'//trim(index_z),'','local cell thickness at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                        act_as_state_variable=.true.,domain=domain_surface,output=output_none)                                        ! current cell thickness
          do i_v=1,self%nvars
               write(index_v,'(i0)') i_v
               if (i_z==1)
                    call self%register_state_dependency(self%id_targets(i_v), 'target'//trim(index_v), '', 'saved variable nr '//trim(index_v))                          ! the target variable (e.g zooplankton)
               endif
               call self%register_diagnostic_variable(self%id_local_pvar(i_z,i_v),'pvar'//trim(index_v)//'_at_z'//trim(index_z),'','local value of var'//trim(index_v)//' at depth interval '//trim(index_z), & ! not sure about act_as_state_variable
                        act_as_state_variable=.true.,domain=domain_surface,output=output_none)
          end do
   end do
   
  end subroutine get_pvar_initialize  
    
    !======================================================================================================================
    ! do part
    !======================================================================================================================
   subroutine get_pvar_do(self,_ARGUMENTS_VERTICAL_)
      class (type_get_pvar),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: local,thickness
      integer  :: i_z

      i_z=1
      
      _VERTICAL_LOOP_BEGIN_
      do i_v=1,self%nvars
         _GET_(self%id_targets(i_v),local)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_pvar(i_z,i_v),local)
      end do
         _GET_(self%id_thickness,thickness)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local_thickness(i_z),thickness)
         
         i_z=i_z+1
         _VERTICAL_LOOP_END_
   end subroutine get_pvar_do
   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! fish_rate_distributor code
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   subroutine type_fish_rate_distributor_initialize(self,configunit)
      class (type_fish_rate_distributor),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit

      call self%get_parameter(self%nlev,'nlev', '-', 'Number of depth levels in water column (same as in host model)')       ! number of depth levels in simulation
      call self%get_parameter(self%nvars,'nvars', '-', 'Number of target variables')       ! number of depth levels in simulation

´     allocate(self%id_fish_dist(self%nlev))
      allocate(self%id_targets(self%nvars))
      allocate(self%id_totals(self%nvars))
      allocate(self%id_weight_order(self%nvars))
      do i_z=1,self%nlev
         write(index,'(i0)') i_z 
         call self%register_horizontal_dependency(self%id_fish_dist(i_z),'fish_dist_at_z'//trim(index),'-','fish distribution weight at depth interval '//trim(index))                              ! the calculated weights from fish_cohort
      end do
      do i_v=1,self%nvars
         write(index,'i0') i_v
         call self%register_state_dependency(self%id_targets(i_v), 'target'//trim(index), '-', 'target variable number '//trim(index))
         call self%register_horizontal_dependency(self%id_totals(i_v), 'target'//trim(index)//'_total_rate', '-', 'total rate of change of target variable number '//trim(index))
         call self%register_horizontal_dependency(self%id_weight_order(i_v), 'weight_exponent_of_target'//trim(index), '-', 'exponent for applying weights to target '//trim(index))
      end do

    end subroutine type_fish_rate_distributor_initialize

    !======================================================================================================================
    ! do part
    !======================================================================================================================
   subroutine type_fish_rate_distributor_do(self,_ARGUMENTS_DO_)
      class (type_fish_rate_distributor),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)  :: total,weight,weight_exp
      integer   :: i_v, i_z
    
      i_z=1
      
      _LOOP_BEGIN_
         ! First compute relative rate of change of depth-integrated target variable.
         _GET_HORIZONTAL_(self%id_fish_dist(i_z),weight)
         do i_v=1,self%nvars
            _GET_HORIZONTAL_(self%id_totals(i_v),total)
            _GET_HORIZONTAL_(self%id_weight_order(i_v),weight_exp)
            
            _SET_ODE_(self%id_targets(i_v),total*weight**weight_exp)
         end do
         i_z=i_z+1
      _LOOP_END_
   end subroutine type_fish_rate_distributor_do

end module