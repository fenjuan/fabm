#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module fish_cohort
!
! !DESCRIPTION:
! Fen, 5th, oct, 2018: the FABM structure fish cohort model, that defined groups of state variables
!and 
!  Nicolas to fill out
!
!  USES:
   use fabm_types

   implicit none

   private

!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_fish_cohort
!     Variable identifiers
!     id_prey:     the component which is uptaken by predator;
!     id_predator: the predator uptaking the prey
      type (type_surface_state_variable_id), allocatable   :: id_N(:)
      type (type_surface_state_variable_id), allocatable   :: id_r_mass(:)
      type (type_surface_state_variable_id), allocatable   :: id_i_mass(:)
      
      type (type_surface_state_variable_id)   :: id_Z_N
!
!     Model parameters
      integer         :: cht_nC
      real(rk)        :: cht_ext,cht_Tmin,cht_Tmax,cht_q_J,cht_q_A
      real(rk)        :: cht_lambda_1,cht_lambda_2,cht_A_mx,cht_W_opt,cht_alpha,cht_x_mu
      real(rk)        :: cht_xi_1 ,cht_xi_2,cht_k_e,cht_rho_1,cht_rho_2,cht_q_s
      real(rk)        :: cht_s,cht_mu_b,cht_mu_0,cht_phi,cht_delta,cht_eps
      real(rk)        :: cht_beta,cht_sigma,cht_x_f,cht_w_b,cht_k_r,cht_theta_m_s
      real(rk)        :: cht_theta_m_e,cht_theta_a_s,cht_theta_a_e,cht_T_ref,cht_K,cht_r
      real(rk)        :: cht_m,cht_Q10_z,cht_R_day
      

!      real(rk), allocatable       :: wb(:) ! mass at birth
!                        
      contains           
!
! Model procedures
      procedure :: initialize
      procedure :: do_surface

   end type type_fish_cohort

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the prey and predator variables and parameters
!
! !INTERFACE:
!
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_fish_cohort), intent(inout), target :: self
   integer,                           intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   integer             :: n
   character(len=64)   :: name, longname
   


   
!EOP
!-----------------------------------------------------------------------
!BOC
!   Register model parameters

   call self%get_parameter(self%cht_nC        ,'number_of_cohorts'      , '[-]'               , 'number of cohorts'                                        , default=10       )
   call self%get_parameter(self%cht_ext       ,'extinction_abundcance'  , '[-]'               , 'extinction threshold abundance'                           , default=2E-9_rk  )
   call self%get_parameter(self%cht_Tmin      ,'min_winter_temp'        , '°C'                , 'minimum winter temperature'                               , default=1.0_rk   )
   call self%get_parameter(self%cht_Tmax      ,'max_summer_temp'        , '°C'                , 'maximum summer temperature'                               , default=22.0_rk  )
   call self%get_parameter(self%cht_q_J       ,'max_jv_con'             , '[-]'               , 'juvenile maximum condition'                               , default=0.74_rk  )
   call self%get_parameter(self%cht_q_A       ,'max_ad_con'             , '[-]'               , 'adult maximum condition'                                  , default=1.37_rk  )
   call self%get_parameter(self%cht_lambda_1  ,'length-weight_scalar'   , '(mm g^(-lambda_2)' , 'allometric length-weight scalar'                          , default=48.3_rk  )
   call self%get_parameter(self%cht_lambda_2  ,'length-weight_exponent' , '[-'                , 'allometric length-weight exponent'                        , default=0.317_rk )
   call self%get_parameter(self%cht_A_mx      ,'max_attack_rate'        , 'L d^-1'            , 'max attack rate'                                          , default=3E4_rk   )
   call self%get_parameter(self%cht_W_opt     ,'att_opt_size'           , 'g'                 , 'size for optimal attack rate on zooplankton'              , default=8.2_rk   )
   call self%get_parameter(self%cht_alpha     ,'att_rate_exp'           , '[-]'               , 'allometric exponent for attack rate on zooplankton'       , default=0.62_rk  )
   call self%get_parameter(self%cht_x_mu      ,'mort_cha_size'          , 'g'                 , 'shape parameter for size-dependent mortality'             , default=0.5_rk   )
   call self%get_parameter(self%cht_xi_1      ,'handling_time_scalar'   , 'd g^-(1+xi_2)'     , 'allometric scalar for handling time'                      , default=5.0_rk   )
   call self%get_parameter(self%cht_xi_2      ,'handling_time_exp'      , '[-]'               , 'allometric exponent for handling time'                    , default=-0.8_rk  )
   call self%get_parameter(self%cht_k_e       ,'ass_efficiency'         , '[-]'               , 'assimilation efficiency'                                  , default=0.61_rk  )
   call self%get_parameter(self%cht_rho_1     ,'metab_scalar'           , '(g^(1-rho2) d^-1'  , 'allometric scalar for basic metabolism'                   , default=0.033_rk )
   call self%get_parameter(self%cht_rho_2     ,'metab_exp'              , '[-]'               , 'allometric exponent for basic metabolism'                 , default=0.77_rk  )
   call self%get_parameter(self%cht_q_s       ,'starv_con'              , '[-]'               , 'starvation condition'                                     , default=0.2_rk   )
   call self%get_parameter(self%cht_s         ,'starv_coeff'            , 'd^-1'              , 'starvation coefficient'                                   , default=0.2_rk   )
   call self%get_parameter(self%cht_mu_b      ,'background_mort'        , 'd^-1'              , 'background mortality'                                     , default=0.0_rk   )
   call self%get_parameter(self%cht_mu_0      ,'mort_scalar_size'       , 'd^-1'              , 'scaling constant for size-dependent mortality'            , default=0.034_rk )
   call self%get_parameter(self%cht_eps       ,'pred_size_max'          , '[-]'               , 'maximum victim/cannibal size ratio'                       , default=0.2_rk   )
   call self%get_parameter(self%cht_phi       ,'pred_size_opt'          , '[-]'               , 'optimal victim/cannibal size ratio'                       , default=0.05_rk  )
   call self%get_parameter(self%cht_delta     ,'pred_size_min'          , '[-]'               , 'minimum victim/cannibal size ratio'                       , default=0.45_rk  )
   call self%get_parameter(self%cht_beta      ,'pred_max_att_scalar'    , 'L d^-1 mm^alpha'   , 'cannibalistic voracity (allometric scalar for piscivory)' , default=0.0_rk   )
   call self%get_parameter(self%cht_sigma     ,'pred_max_att_exp'       , '[-]'               , 'allometric exponent for cannibalism (piscivory)'          , default=0.6_rk   )
   call self%get_parameter(self%cht_x_f       ,'mat_irr_mass'           , 'g'                 , 'maturation irreversible mass'                             , default=4.6_rk   )
   call self%get_parameter(self%cht_w_b       ,'egg_size'               , 'g'                 , 'egg size'                                                 , default=1.8e-3_rk)
   call self%get_parameter(self%cht_k_r       ,'repro_eff'              , '[-]'               , 'gonad-offspring conversion factor'                        , default=0.5_rk   )
   call self%get_parameter(self%cht_theta_m_s ,'Q10_metab_scalar'       , '[-]'               , 'allometric scalar of metabolism Q10'                      , default=2.0_rk   )
   call self%get_parameter(self%cht_theta_m_e ,'Q10_metab_exp'          , '[-]'               , 'allometric exponent of metabolism Q10'                    , default=0.073_rk )
   call self%get_parameter(self%cht_theta_a_s ,'Q10_att_scalar'         , '[-]'               , 'allometric scalar of metabolism Q10'                      , default=2.8_rk   )
   call self%get_parameter(self%cht_theta_a_e ,'Q10_att_exp'            , '[-]'               , 'allometric exponent of metabolism Q10 '                   , default=0.072_rk )
   call self%get_parameter(self%cht_T_ref     ,'temp_ref_fish'          , '°C'                , 'reference temperature for perch '                         , default=20.0_rk  )
   call self%get_parameter(self%cht_K         ,'carr_food'              , '[-]'               , 'resource carrying capacity'                               , default=100.0_rk )
   call self%get_parameter(self%cht_r         ,'food_growth'            , '[-]'               , 'resource growth rate'                                     , default=0.1_rk   )
   call self%get_parameter(self%cht_m         ,'food_ind_weight'        , 'g'                 , 'weight of resource individual'                            , default=3e-5_rk  )
   call self%get_parameter(self%cht_Q10_z     ,'Q10_food'               , '[-]'               , 'Q10 value of zooplankton population growth'               , default=1.799_rk )
   
   call self%get_parameter(self%cht_R_day     ,'repro_day'              , '[-]'               , 'reproduction day of year for fish '                       , default=152.0_rk )


!  allocate state variable pointers, according to user-defined number of cohorts(cht_nC)
   allocate(self%id_N(self%cht_nC))
   allocate(self%id_r_mass(self%cht_nC))
   allocate(self%id_i_mass(self%cht_nC))
!   allocate(self%wb(self%Nclasses))


   ! get parameters for specific class
!   do n=1,self%Nclasses
!      write (name,"(I02.2,A3)") n-1,"/wb"
!      call self%get_parameter(self%wb(n),name,'d-1','growth rate of prey',default=1.00_rk,scale_factor=1._rk)
!   end do
!
!write(*,*) self%wb
!stop

   ! setup fish cohorts 
   do n=1,self%cht_nC
      write (name,"(A1,I02.2)") "N", n-1
      write (longname, "(A9, I02.2)") "abundance", n-1
      call self%register_state_variable(self%id_N(n),name, '-', longname,initial_value= 0._rk, minimum=0.0_rk)
      write (name,"(A6,I02.2)") "r_mass", n-1
      write (longname, "(A15, I02.2)") "reversible mass", n-1
      call self%register_state_variable(self%id_r_mass(n),name,'kg/individual',longname,initial_value= 0._rk,minimum=0.0_rk)
      write (name,"(A6,I02.2)") "i_mass", n-1
      write (longname, "(A17, I02.2)") "irreversible mass", n-1
      call self%register_state_variable(self%id_i_mass(n),name,'kg/individual',longname,initial_value=0._rk,minimum=0.0_rk)
   end do
   
! seup zooplankton
call self%register_state_variable(self%id_Z_N,'prey','ind./L^-1','food concentration',initial_value=self%cht_K,minimum=0.0_rk)
   
#if 0
!   call self%get_parameter(self%b,'b','d-1','growth rate of prey',   default=1.00_rk,scale_factor=d_per_s)
!   call self%get_parameter(self%p,'p','d-1','impact of predation',   default=0.05_rk,scale_factor=d_per_s)
!   call self%get_parameter(self%r,'r','d-1','growth efficiency rate',default=0.02_rk,scale_factor=d_per_s)
!   call self%get_parameter(self%d,'d','d-1','death rate',            default=0.50_rk,scale_factor=d_per_s)
!
!  Register state variables
!   call self%register_state_variable(self%id_prey,'prey','mmol/m**3','nutrient',     &
!                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)
!   call self%register_state_variable(self%id_predator,'predator','mmol/m**3','phytoplankton',     &
!                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)
!
!!  Register conserved quantities
!   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_prey)
!   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_predator)
#endif

   return
   end subroutine initialize

!-----------------------------------------------------------------------
! !IROUTINE:the type bound precedure: do(),right hand sides of prey and predator model
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_fish_cohort),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
!  state variables
   real(rk)   :: N,r_mass,i_mass,Z_N
   real(rk)   :: d_N, d_rmass, d_imass, d_ZN
   real(rk)    :: g, f,tot_r
   integer    :: i
!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

! get zooplankton (food) concentration:
_GET_HORIZONTAL_(self%id_Z_N,Z_N)

! now do fish:
if (Day >= self%cht_R_day .AND. Day <= self%cht_R_day + 1.0_rk) then  ! determine if reproduction is due
    !-----------------------------------------------------------------------
    !  fish reproduction
    !-----------------------------------------------------------------------
    do i=1,self%cht_nC
      _GET_HORIZONTAL_(self%id_N(i),N)
      _GET_HORIZONTAL_(self%id_r_mass(i),r_mass)
      _GET_HORIZONTAL_(self%id_i_mass(i),i_mass)
      ! do something more useful here

      !      write(*,*) N,r_mass,i_mass
      _SET_SURFACE_ODE_(self%id_N,g)
      _SET_SURFACE_ODE_(self%id_r_mass,tot_r)
      _SET_SURFACE_ODE_(self%id_i_mass,f)
   end do
else                                                                  ! reproduction is not due
    !-----------------------------------------------------------------------
    !  fish dynamics
    !-----------------------------------------------------------------------
    
    !      do yearly dynamics:
    do i=1,self%cht_nC
      _GET_HORIZONTAL_(self%id_N(i),N)
      _GET_HORIZONTAL_(self%id_r_mass(i),r_mass)
      _GET_HORIZONTAL_(self%id_i_mass(i),i_mass)
      ! do something more useful here
      
!      write(*,*) N,r_mass,i_mass
      _SET_SURFACE_ODE_(self%id_N,g)
      _SET_SURFACE_ODE_(self%id_r_mass,tot_r)
      _SET_SURFACE_ODE_(self%id_i_mass,f)
   end do
endif
   
!   do i=1,self%cht_nC
!      _GET_HORIZONTAL_(self%id_N(i),N)
!      _GET_HORIZONTAL_(self%id_r_mass(i),r_mass)
!      _GET_HORIZONTAL_(self%id_i_mass(i),i_mass)
!      ! do something more useful here
!      g = 0.0001_rk*self%cht_r
!      f = 0.00001_rk*self%cht_mu_b
!      
!      tot_r = g*r_mass - f*r_mass
!!      write(*,*) N,r_mass,i_mass
!      _SET_SURFACE_ODE_(self%id_N,g)
!      _SET_SURFACE_ODE_(self%id_r_mass,tot_r)
!      _SET_SURFACE_ODE_(self%id_i_mass,f)
!   end do


   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

!-----------------------------------------------------------------------

end module fish_cohort

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
