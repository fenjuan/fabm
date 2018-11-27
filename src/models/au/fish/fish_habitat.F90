#include "fabm_driver.h"
   
module fish_habitat

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
   ! and this model will then automagically distribute those sinks and sources again over their
   ! original depth-explicit source variable, using the appropriate weights.

!  food_integral module: added food source together, which will be used in distribution module
   type,extends(type_base_model),public :: type_depth_integral
      type (type_horizontal_diagnostic_variable_id) :: id_integral, id_habitdepth,id_foodsource
      type (type_state_variable_id)                 :: id_target
      type (type_dependency_id)                     :: id_thickness, id_depth
      type (type_dependency_id)                     :: id_temp, id_light
      
      integer :: env_depend  ! environmentl dependent type
      real(rk)   :: min_light,max_light, min_temp, max_temp
   contains
      procedure :: initialize => depth_integral_initialize
      procedure :: get_light  => depth_integral_do ! currently only get_light supports non-local action in depth
   end type
   
   ! type_vertical_distribution: module that computes a "vertical distribution" defined in terms
   ! of weights applied to different depth levels when computing a vertical integral (e.g., of prey).
   ! In this simple model, the vertical distribution is defined by an absolute upper depth limit (z_top),
   ! an absolute lower depth limit (z_bot), a weight to apply at the top of the domain (w_top) and a weight to
   ! apply at the bottom of the domain (w_bot). Linear interpolation is used to determine weights in between
   ! these levels. In the special case where both weights equal 1, a vertical integral computed using
   ! this presence distribution will in fact be the real vertical integral between z_top and z_bot.
   type,extends(type_base_model),public :: type_vertical_distribution
      type (type_diagnostic_variable_id) :: id_weights
      type (type_horizontal_dependency_id) :: id_integral
!   make a state dependency on food, use to calculate the distribution factor
      type (type_state_variable_id)       :: id_target
!   added environmental dependencies of light and tempeture, for defining fish habitat
      type (type_dependency_id)          :: id_depth, id_light, id_temp
!   add critial light and temperature factor for fish habitat
      integer    :: env_depend
      real(rk)   :: min_light,max_light, min_temp, max_temp
!      real(rk) :: z_top, z_bot, w_top, w_bot    ! this will be automatically calculated depending on real case
   contains
      procedure :: initialize => vertical_distribution_initialize
      procedure :: do         => vertical_distribution_do
   end type


   

   type,extends(type_base_model), public :: type_depth_integral_rate_distributor
!      type (type_bottom_state_id) :: id_predator   ! the predator, depth integrated state, works for bottom for now. should be horizontal
      type (type_horizontal_dependency_id) :: id_sms        ! Depth-integrated sources-sinks of predation variable, rates retured from predator
      type (type_state_variable_id)        :: id_target     ! Depth-explicit variable that should absorp the sources-sinks
      type (type_dependency_id)            :: id_weights    ! Weights for the vertical distribution of the sinks and sources
   contains
      procedure :: initialize => depth_integral_rate_distributor_initialize
      procedure :: do         => depth_integral_rate_distributor_do ! currently only get_light supports non-local action in depth
   end type
   

   contains
    
    subroutine depth_integral_initialize(self,configunit)
      class (type_depth_integral),intent(inout),target :: self
      integer,                    intent(in)           :: configunit

      
      call self%get_parameter(self%env_depend, 'env_depend', '[-]','environmental dependency for fish habitat: 0, ignore environmental dependency;  &
          & 1. habitat is temperature dependent; 2. habitat is light dependent; 3. habitat is both light/temperature dependent',default = 0)
      call self%get_parameter(self%min_light, 'min_light', 'W m-2','lower limit light for fish habitat',default = 0.0_rk)
      call self%get_parameter(self%max_light, 'max_light', 'W m-2','upper limit light for fish habitat', default = 300.0_rk)
      call self%get_parameter(self%min_temp, 'min_temp', 'degree_Celsius','lower limit temperature for fish habitat',default = 0.0_rk)
      call self%get_parameter(self%max_temp, 'max_temp', 'degree_Celsius','upper limit temperature for fish habitat', default = 30.0_rk)  ! default value for temperate fresh water fish

      call self%register_state_dependency(self%id_target, 'target', 'g/m-3', 'variable to depth-integrate')
!      call self%register_dependency(self%id_weights,'weights','-','weights for vertical integration')
      call self%register_diagnostic_variable(self%id_integral,'integral','g/m2','depth-integrated target variable',missing_value=0.0_rk, &
         act_as_state_variable=.true.,domain=domain_bottom)  !  this variable in the end should be calculated with lake volume, and units in g?? 
      call self%register_diagnostic_variable(self%id_habitdepth,'habitdepth','m','habitat depth',missing_value=0.0_rk, &
         act_as_state_variable=.false.,domain=domain_bottom)
      call self%register_diagnostic_variable(self%id_foodsource,'foodsource','g/m3','food in concentration',missing_value=0.0_rk, &
         act_as_state_variable=.true.,domain=domain_bottom)
      call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)
      call self%register_dependency(self%id_depth,standard_variables%depth)
      call self%register_dependency(self%id_light,standard_variables%downwelling_shortwave_flux)
      call self%register_dependency(self%id_temp,standard_variables%temperature)


    end subroutine depth_integral_initialize

    !!!!! fen: calculate the sum of the integral target variable
!!!   
   subroutine depth_integral_do(self,_ARGUMENTS_VERTICAL_)
      class (type_depth_integral),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: thickness, z,light,temp
      real(rk) :: integral,local, habitdepth

      integral = 0.0_rk
      habitdepth = 0.0_rk
      _VERTICAL_LOOP_BEGIN_
      
         _GET_(self%id_thickness,thickness)
         _GET_(self%id_depth, z)
         _GET_(self%id_light,light)
         _GET_(self%id_temp,temp)
        select case (self%env_depend)
        case(0)
             _GET_(self%id_target, local)
             integral = integral + local * thickness  !! for now tot_food is g/m2
             habitdepth = habitdepth + thickness
        case(1)
            if (temp < self%min_temp .OR. temp > self%max_temp) then
               integral = integral + 0.0_rk   ! if out of habitat, not predating on current place
            else 
               _GET_(self%id_target, local)
               integral = integral + local * thickness  !! for now tot_food is g/m2
               habitdepth = habitdepth + thickness
            end if
        case(2)
            if (light < self%min_light .OR. light > self%max_light) then
                integral = integral + 0.0_rk
            else
                _GET_(self%id_target, local)
                integral = integral + local * thickness  !! for now tot_food is g/m2
                habitdepth = habitdepth + thickness
            end if
        case(3)
            if ( temp > self%min_temp .AND. temp < self%max_temp &
                & .AND. light > self%min_light .AND. light < self%max_light) then
                _GET_(self%id_target, local)
                integral = integral + local * thickness
                habitdepth = habitdepth + thickness
            else
                integral = integral + 0.0_rk
            end if
        end select

      _VERTICAL_LOOP_END_

      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_integral,integral)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_habitdepth,habitdepth)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_foodsource,integral/habitdepth)
      
   end subroutine depth_integral_do
   
   
   subroutine vertical_distribution_initialize(self,configunit)
      class (type_vertical_distribution),intent(inout),target :: self
      integer,                           intent(in)           :: configunit

      call self%get_parameter(self%env_depend, 'env_depend', '[-]','environmental dependency for fish habitat: 0, ignore environmental dependency;  &
          & 1. habitat is temperature dependent; 2. habitat is light dependent; 3. habitat is both light/temperature dependent',default = 0)
      call self%get_parameter(self%min_light, 'min_light', 'W m-2','lower limit light for fish habitat',default = 0.0_rk)
      call self%get_parameter(self%max_light, 'max_light', 'W m-2','upper limit light for fish habitat', default = 300.0_rk)
      call self%get_parameter(self%min_temp, 'min_temp', 'degree_Celsius','lower limit temperature for fish habitat',default = 0.0_rk)
      call self%get_parameter(self%max_temp, 'max_temp', 'degree_Celsius','upper limit temperature for fish habitat', default = 30.0_rk)  ! default value for temperate fresh water fish

      call self%register_diagnostic_variable(self%id_weights,'weights','-','weights')
      call self%register_dependency(self%id_integral,'integral','g/m2','integrated target variable')
      call self%register_dependency(self%id_depth,standard_variables%depth)
      call self%register_dependency(self%id_light,standard_variables%downwelling_shortwave_flux)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      
      
   end subroutine vertical_distribution_initialize
   
   
      subroutine vertical_distribution_do(self,_ARGUMENTS_DO_)
      class (type_vertical_distribution),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: z,light, temp    ! environmental dependencies for defining fish habitat
      real(rk) :: weight, local
      real(rk) :: integral  ! this should be a dependency on intergral, which will get the total amount of food in whole water column

      _LOOP_BEGIN_
!        get environmental dependencisy
         _GET_(self%id_depth,z)
         _GET_(self%id_light,light)
         _GET_(self%id_temp,temp)
         _GET_HORIZONTAL_(self%id_integral,integral)
!       get state dependencies(usually food) which can decide distribution weights
         _GET_(self%id_target,local)
!        define fish habitat depending on light and temperature
        select case (self%env_depend)
        case(0)
            weight = local /integral
        case(1)
            if (temp < self%min_temp .OR. temp > self%max_temp) then
               weight = 0.0_rk
            else 
               weight = local /integral
!               weight = food_conc * current_water_volumn /total_foodbiomass total food biomass should be the total mass in the habitat
            end if
        case(2)
            if (light < self%min_light .OR. light > self%max_light) then
                weight = 0.0_rk
            else
                weight = local /integral
            end if
        case(3)
            if ( temp > self%min_temp .AND. temp < self%max_temp &
                & .AND. light > self%min_light .AND. light < self%max_light) then
                weight = local/integral
            else
            weight = 0.0_rk
            end if
        end select
      _SET_DIAGNOSTIC_(self%id_weights, weight)
      _LOOP_END_

      end subroutine vertical_distribution_do
      
   
! fen: distribute feedback to local trophic levels, such as zooplankton, nutrients, detritus, etc.
   subroutine depth_integral_rate_distributor_initialize(self,configunit)
      class (type_depth_integral_rate_distributor),intent(inout),target :: self
      integer,                                     intent(in)           :: configunit
      
!      class (type_depth_integral_rate_distributor),pointer :: rate_distributor

      call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
      call self%register_dependency(self%id_weights,'weights','-','weights for vertical distribution')
!      call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')   !the applied state variable
      call self%register_dependency(self%id_sms,'sms','','depth-integrated sources-sinks')  ! the applied local rates, sink-source terms
      

     ! Create a child model that receives the intended rate of the change of the depth-integrated variable,
     ! and redistributes that rate of change over the pelagic variable that the integral was originally computed from.
     ! All dependencies of the child model can be resolved by coupling to our own variables.
      
    ! fen:  this is not working for lake branch since each of the benthic layer will make the child model
!      allocate(rate_distributor)
!      call self%add_child(rate_distributor,'rate_distributor',configunit=-1)
!      call rate_distributor%request_coupling(rate_distributor%id_target,'target')
!      call rate_distributor%request_coupling(rate_distributor%id_weights,'weights')
!!      call rate_distributor%request_coupling(rate_distributor%id_integral,'result')
!      call rate_distributor%request_coupling(rate_distributor%id_sms,'result_sms_tot')
!      call self%register_dependency(self%id_target,'target','g/m3','target variable')
!      call self%register_dependency(self%id_weights,'weights','-','distrubtion weights')
!      call self%register_dependency(self%id_sms,'result_sms_tot','g/m2/s','integrated changing rate')
      
   end subroutine depth_integral_rate_distributor_initialize
   
      subroutine depth_integral_rate_distributor_do(self,_ARGUMENTS_DO_)
      class (type_depth_integral_rate_distributor),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: integrated_sms,change
      real(rk) :: local,weight

      _LOOP_BEGIN_
         ! First compute relative rate of change of depth-integrated target variable.
!         _GET_HORIZONTAL_(self%id_integral,integral)
         _GET_HORIZONTAL_(self%id_sms,integrated_sms)
         _GET_(self%id_target,local)
         _GET_(self%id_weights,weight)
         if (integrated_sms/=0.0_rk) then
            change = integrated_sms * weight
         else
            change = 0.0_rk
         end if

         ! Now distribute the same relative change across the column.

         _SET_ODE_(self%id_target,change)
      _LOOP_END_
   end subroutine depth_integral_rate_distributor_do

end module