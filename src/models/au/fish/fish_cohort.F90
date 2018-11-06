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
      
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm
      type (type_global_dependency_id)         :: id_Day
!
!     Model parameters
      integer         :: cht_nC,cht_nc_init
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
   call self%get_parameter(self%cht_nc_init   ,'initial_cohorts'        , '[-]'               , 'number of cohorts present at simulation start '           , default=1 )


!  allocate state variable pointers, according to user-defined number of cohorts(cht_nC)
   allocate(self%id_N(self%cht_nC+self%cht_nc_init))
   allocate(self%id_r_mass(self%cht_nC+self%cht_nc_init))
   allocate(self%id_i_mass(self%cht_nC+self%cht_nc_init))
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
   do n=1,self%cht_nC+self%cht_nc_init
      write (name,"(A1,I02.2)") "N", n-1
      write (longname, "(A9, I02.2)") "abundance", n-1
      call self%register_state_variable(self%id_N(n),name, '-', longname,initial_value= 0._rk, minimum=0.0_rk)
      write (name,"(A6,I02.2)") "r_mass", n-1
      write (longname, "(A15, I02.2)") "reversible mass", n-1
      call self%register_state_variable(self%id_r_mass(n),name,'kg/individual',longname,initial_value= self%cht_q_J*self%cht_w_b/(1+self%cht_q_J),minimum=0.0_rk)
      write (name,"(A6,I02.2)") "i_mass", n-1
      write (longname, "(A17, I02.2)") "irreversible mass", n-1
      call self%register_state_variable(self%id_i_mass(n),name,'kg/individual',longname,initial_value=self%cht_w_b/(1+self%cht_q_J),minimum=0.0_rk)
   end do
   
! setup zooplankton
call self%register_state_variable(self%id_Z_N,'prey','ind./L^-1','food concentration',initial_value=self%cht_K,minimum=0.0_rk)

!  register environmental dependencies
call self%register_dependency(self%id_uTm,    standard_variables%temperature)  ! eventually, this shoud be replaced
call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)

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
!  LOCAL VARIABLES:
!  state variables
   real(rk)                     :: N,r_mass,i_mass,Z_N
!  derivatives
   real(rk)                     :: d_N, d_rmass, d_imass, d_ZN
!  local variables
   real(rk),save                :: F_tot=0.0_rk
   real(rk), allocatable, save  :: g_mass(:) !                                                                                                                                                              <==================<<<
!                                           
!  carriers for environment dependencies and logical switches
   real(rk)                     :: uTm,Day
   integer                      :: i, ix_repro !                                                                                                                                                            <==================<<<
   integer                      :: tcheck=0
   integer,save                 :: year=0

!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1. retrieve variable values, update trackers & variables:
!===============================================================================================
    _GET_HORIZONTAL_(self%id_Z_N,Z_N)                                                           ! get zooplankton (food) concentration
    
    !  retrieve environmental dependencies:
    _GET_(self%id_uTm,uTm)                                                                      ! get temperature. eventually, this shoud be replaced
    _GET_GLOBAL_(self%id_Day,Day)                                                               ! get current julian day
    
    ! check & update year - used for selecting new cohorts to start. To be removed if
    !                       cohort implementation changes or if functionality to get 
    !                       simulation year identification is added to FABM
    if (day==1 .AND. tcheck == 0) then
        year=year+1
        tcheck=1
    elseif (day==1 .AND. tcheck == 1) then
        year=year
    else
        tcheck=0
    endif
    
    ! allocate dimensions for saved cohort-specific variables
    if (.not. allocated(g_mass)) allocate(g_mass(self%cht_nC+self%cht_nc_init)) !                                                                                                                           <==================<<<
             
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2. now do fish:
!===============================================================================================
    do i=1,self%cht_nC ! loop over the cohorts
        _GET_HORIZONTAL_(self%id_N(i),N)                                                        ! get cohort abundance
        _GET_HORIZONTAL_(self%id_r_mass(i),r_mass)                                              ! get cohort reversible mass (reserves+gonads)
        _GET_HORIZONTAL_(self%id_i_mass(i),i_mass)                                              ! get cohort irreversible mass (structural tissues, organs)

        if (N>self%cht_ext) then ! if cohort is not empty or extinct
    !===========================================================================================
    ! 2.1 fish reproduction
    !-------------------------------------------------------------------------------------------
            if (Day >= self%cht_R_day .AND. Day <= self%cht_R_day + 1.0_rk .AND. ix_repro==0) then ! if reproduction is due                                                                                 <==================<<<
                if (i_mass >= self%cht_x_f .AND. r_mass > self%cht_q_J*i_mass) then ! cohort is mature and non-starving
                    g_mass(i) = (r_mass-self%cht_q_J*i_mass)                                    ! calculate total loss of r_mass due to reproduction (gonad mass)
                    F_tot = F_tot + g_mass(i)*N*self%cht_k_r/self%cht_w_b                       ! add offspring from cohort i to total produced offspring                                                   <==================<<<
                else ! cohort is not ready to reproduce
                    g_mass(i) = 0.0_rk                                                          ! no change in r_mass from reproduction for cohort i
                endif   
                ix_repro = int(i/self%cht_nc)                                                   ! only calculate reproduction once per reproductive event; stop when all cohorts have been cycled through   <==================<<<
            else ! reproduction is not due
                ix_repro=0                                                                      ! reset reproductive switch                                                                                 <==================<<<
                g_mass(i) = 0.0_rk                                                              ! no change in r_mass from reproduction for cohort i
                F_tot = 0.0_rk                                                                  ! no reproduction the rest of the year
            endif
    !===========================================================================================
    ! 2.2 fish dynamics
    !-------------------------------------------------------------------------------------------
            ! calculate size-specific rates and parameters
            ! calculate temperature correction rates
            ! do predation on other fish (looping over other fish)
            ! calculate total encounter rates and final growth rates
            ! calculate final mortality rates
            ! allocate growth to reversible/irreversible mass
            d_N = 0.0_rk
            d_imass = 0.0_rk
            d_rmass = 0.0_rk
        else ! cohort is empty
            ! rates of change of epmpty cohorts are zero
            d_N = 0.0_rk
            d_imass = 0.0_rk
            d_rmass = 0.0_rk     
        endif
    !===========================================================================================
    ! 2.3 update fish rates of change
    !-------------------------------------------------------------------------------------------
        _SET_SURFACE_ODE_(self%id_N(i),d_N)
        _SET_SURFACE_ODE_(self%id_r_mass(i),d_rmass)
        _SET_SURFACE_ODE_(self%id_i_mass(i),d_imass)
    end do
    !============================================================================================
    ! 2.4 special case: new cohort
    !--------------------------------------------------------------------------------------------
    if (Day >= self%cht_R_day .AND. Day <= self%cht_R_day + 1.0_rk) then                     ! reproduction is due                                                                                           <==================<<<
        _SET_SURFACE_ODE_(self%id_N(year+self%cht_nc_init),F_tot)                            ! start new cohort with the produced offspring                                                                  <==================<<<
    endif
   
    
    ! do zooplankton
    ! d_ZN = 
     _SET_SURFACE_ODE_(self%id_Z_N,d_ZN)
     
    
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
