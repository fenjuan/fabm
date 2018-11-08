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
   call self%get_parameter(self%cht_nc_init   ,'initial_cohorts'        , '[-]'               , 'number of cohorts present at simulation start '           , default=1        )


!  allocate state variable pointers, according to user-defined number of cohorts(cht_nC)
   allocate(self%id_N(self%cht_nC+self%cht_nc_init))
   allocate(self%id_r_mass(self%cht_nC+self%cht_nc_init))
   allocate(self%id_i_mass(self%cht_nC+self%cht_nc_init))

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
!  state variables -------------------------------------------------------------------------------------------------
   real(rk), DIMENSION(self%cht_nC+self%cht_nc_init)        :: N,       r_mass,   i_mass
   real(rk)                                                 :: Z_N
!  derivatives -----------------------------------------------------------------------------------------------------
   real(rk), dimension(self%cht_nC+self%cht_nc_init)        :: d_N,     d_rmass,  d_imass
   real(rk)                                                 :: d_ZN
!  local variables which carry over between time steps -------------------------------------------------------------                          
   real(rk), allocatable, save                              :: g_mass(:)                                            ! g_mass should be saved between time step in another way before multiple fish instances are run in FABM, e.g. as a derived variable
!  local cohort vector variables -----------------------------------------------------------------------------------
   real(rk), dimension(self%cht_nC+self%cht_nc_init)        :: lengths, st_mass,  H,        A_z,    mu_size,    E_m         
   real(rk), dimension(self%cht_nC+self%cht_nc_init)        :: eta_z,   Q10_f_m,  Q10_f_a,  rt_f_m, rt_f_a,     eta
   real(rk), dimension(self%cht_nC+self%cht_nc_init)        :: mu_c,    E_a,      Ing,      E_g,    mu_s,       mu
! local cohort matrix variables ------------------------------------------------------------------------------------
   real(rk), dimension(self%cht_nC+self%cht_nc_init,    &
       self%cht_nC+self%cht_nc_init)                        :: A_c,     L_c,      L_v,      eta_c
! local single within-loop variables -------------------------------------------------------------------------------
   real(rk)                                                 :: F_tot
!  carriers for environment dependencies and logical switches ------------------------------------------------------
   real(rk)                                                 :: uTm,     Day
   integer                                                  :: i,       j,        ix_repro, nC_fin
   integer                                                  :: tcheck=0
   integer, save                                            :: year=0

!EOP
!-------------------------------------------------------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1. retrieve variable values, update trackers & variables:
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _GET_HORIZONTAL_(self%id_Z_N,Z_N)                                                                           ! get zooplankton (food) concentration
        _GET_HORIZONTAL_(self%id_N,N)                                                                               ! get cohort abundance
        _GET_HORIZONTAL_(self%id_r_mass,r_mass)                                                                     ! get cohort reversible mass (reserves+gonads)
        _GET_HORIZONTAL_(self%id_i_mass,i_mass)                                                                     ! get cohort irreversible mass (structural tissues, organs)
    
        !  retrieve environmental dependencies:
        _GET_(self%id_uTm,uTm)                                                                                      ! get temperature. eventually, this shoud be replaced
        _GET_GLOBAL_(self%id_Day,Day)                                                                               ! get current julian day
        
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
        nC_fin = self%cht_nC+self%cht_nc_init                                                                       ! calculate total number of cohorts
        if (.not. allocated(g_mass)) allocate(g_mass(nC_fin))     
                 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2. do fish:
    !===============================================================================================================
        where (N<self%cht_ext) ! find cohorts that are extinct
            N = 0.0_rk                                                                                              ! treat extinct cohorts as empty
            r_mass = 0.0_rk                                                                                         ! treat extinct cohorts as empty
            i_mass = 0.0_rk                                                                                         ! treat extinct cohorts as empty
        endwhere
        !===========================================================================================================
        ! 2.1 calculate size-specific rates and parameters
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        st_mass = i_mass*(1+self%cht_q_J)                                                                           ! calculate fish standard mass
        lengths = self%cht_lambda_1*st_mass**self%cht_lambda_2                                                      ! calculate fish lengths
        H = self%cht_xi_1*st_mass**self%cht_xi_2                                                                    ! calculate size-specific handling times (digestion)
        A_z = self%cht_A_mx*((st_mass/self%cht_W_opt)*exp(1-st_mass/self%cht_W_opt))**self%cht_alpha                ! calculate size-specific attack rate on zooplankton
        mu_size = self%cht_mu_0*exp(-i_mass/self%cht_x_mu)                                                          ! calculate size-spcific mortality
        E_m = self%cht_rho_1*(i_mass+r_mass)**self%cht_rho_2                                                        ! calculate size-specific respiration
        eta_z = A_z*Z_N*self%cht_m                                                                                  ! calculate encounter rate with zooplankton
        !===========================================================================================================
        ! 2.2 calculate temperature correction factors
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                
        Q10_f_m = self%cht_theta_m_s*(r_mass+i_mass)**self%cht_theta_m_e                                            ! calculate size-specific Q10 for fish metabolism
        Q10_f_a = self%cht_theta_a_s*(r_mass+i_mass)**self%cht_theta_a_e                                            ! calculate size-specific Q10 for fish feeding
        rt_f_m = Q10_f_m**((uTm-self%cht_T_ref)/10)                                                                 ! calculate temperature correction factor for metabolism
        rt_f_a = Q10_f_a**((uTm-self%cht_T_ref)/10)                                                                 ! calculate temperature correction factor for feeding     
        !===========================================================================================================
        ! 2.3 calculate cannibalistic (eventually all piscivorous) interactions
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        L_c = spread(lengths, DIM=1, NCOPIES=nC_fin)                                                                ! generate predator interaction matrix
        L_v = spread(lengths, DIM=2, NCOPIES=nC_fin)                                                                ! generate victim interaction matrix
        where (L_v>(L_c*self%cht_delta) .AND. L_v<=(L_c*self%cht_phi)) ! where victims are within the lower size range of predators
            A_c = (self%cht_beta*L_c**self%cht_sigma)*(L_v-L_c*self%cht_delta)/((self%cht_phi-self%cht_delta)*L_c)  ! calculate piscivorous attack rates
        elsewhere (L_v>(L_c*self%cht_phi) .AND. L_v<(L_c*self%cht_eps)) ! where victims are within the upper size range of predators
            A_c = (self%cht_beta*L_c**self%cht_sigma)*(self%cht_eps*L_c-L_v)/((self%cht_eps-self%cht_phi)*L_c)      ! calculate piscivorous attack rates
        elsewhere ! where victims are outside the predation window of predators
            A_c = 0.0_rk                                                                                            ! no predation
        endwhere
        !===========================================================================================================
        ! 2.4 scale rates to temperature
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                
        H = H*(1/rt_f_a)                                                                                            ! temperature correct handling times
        where (isnan(H)) H = 0                                                                                      ! get rid of NANs
        A_z = A_z*rt_f_a                                                                                            ! temperature correct attack rate on zooplankton
        E_m = E_m*rt_f_m                                                                                            ! temperature correct basal metabolism
        A_c = A_c*spread(rt_f_a, DIM=1,  NCOPIES=nC_fin)                                                            ! temperature correct piscivorous attack rates according to attackers temperature correction factor
        !===========================================================================================================
        ! 2.5 fish reproduction
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        if (Day >= self%cht_R_day .AND. Day <= self%cht_R_day + 1.0_rk .AND. ix_repro==0) then ! first time step of reproduction event
            where (i_mass >= self%cht_x_f .AND. r_mass > self%cht_q_J*i_mass) ! where cohorts are mature and non-starving
                g_mass = (r_mass-self%cht_q_J*i_mass)                                                               ! calculate total loss of r_mass due to reproduction (gonad mass), which is equivalent to the daily loss rate when reproduction lasts one day
                A_z = 0.0_rk                                                                                        ! reproducing individuals are not feeding
            elsewhere ! cohort is not ready to reproduce
                g_mass = 0.0_rk                                                                                     ! no change in r_mass from reproduction for immature or starving cohorts
            endwhere
            where  (spread(i_mass, DIM=1,  NCOPIES=nC_fin) >= self%cht_x_f .AND.    &
                    spread(r_mass, DIM=1,  NCOPIES=nC_fin) >  self%cht_q_J*spread(i_mass, DIM=1,  NCOPIES=nC_fin)) ! where cohorts are mature and non-starving (matrix edition)
                    A_c = 0.0_rk                                                                        ! reproducing individuals are not feeding
            endwhere
            A_z(year+self%cht_nc_init) = 0.0_rk                                                                     ! new cohort only starts feeding after the reproductive event
            A_c(1:nC_fin,year+self%cht_nc_init) = 0.0_rk                                                            ! new cohort only starts feeding after the reproductive event
            ix_repro = 1                                                                                            ! only calculate reproduction at the start of the reproductive event
        elseif (Day < self%cht_R_day .OR. Day > self%cht_R_day + 1.0_rk) then ! rest of the year
            ix_repro=0                                                                                              ! reset reproductive switch                                                                              
            g_mass = 0.0_rk                                                                                         ! no reproduction outside reproductive event
        endif
        F_tot = sum(g_mass*N*self%cht_k_r/self%cht_w_b)                                                             ! calculate total number of produced offspring (zero when g_mass is zero)
        !===========================================================================================================
        ! 2.6 calculate final vital rates
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        eta_c = A_c*spread((r_mass+i_mass)*N, DIM=2, NCOPIES=nC_fin)                                                ! calculate encounter rates between piscivores  and fish prey
        eta = eta_z+sum(eta_c,1)                                                                                    ! calculate total encounter rate of piscivores with prey
        mu_c = sum(A_c*spread(N/(1+H*eta), DIM=1, NCOPIES=nC_fin),2)                                                ! calculate piscivorous mortality on cohorts
        Ing = eta/(1+H*eta)                                                                                         ! calculate total food intake rate
        E_a = self%cht_k_e*Ing                                                                                      ! calculate assimilation
        E_g = E_a-E_m                                                                                               ! final growth (or degrowth)
        where (r_mass<self%cht_q_s*i_mass) ! find which cohorts are starving
            mu_s = (self%cht_s*(self%cht_q_s*i_mass/r_mass-1))                                                      ! calculate starvation mortality
        elsewhere
            mu_s = 0.0_rk                                                                                           ! otherwise no starvation mortality
        endwhere
        where (mu_s>10) ! find where starvation mortality is very high
            mu_s=10                                                                                                 ! too high a starvation mortality might produce problems with integration
        endwhere
        mu = self%cht_mu_b + mu_size + mu_s + mu_c                                                                  ! total mortality
        !===========================================================================================================
        ! 2.7 set ODEs
        !-----------------------------------------------------------------------------------------------------------
        d_N = -mu*N                                                                                                 ! rate of change of cohort abundances (mortality)
        where (E_g>=0 .AND. i_mass>=self%cht_x_f) ! mature cohorts with positive growth
            d_imass = (1/((1+self%cht_q_A)*self%cht_q_A))*(r_mass/i_mass)*E_g                                       ! calculate part of growth allocated to structure
            d_rmass = (1-(1/((1+self%cht_q_A)*self%cht_q_A))*(r_mass/i_mass))*E_g - g_mass                          ! calculate part of growth allocated to reserves and gonads
        elsewhere (E_g>=0 .AND. i_mass<self%cht_x_f) ! juvenile cohorts with positive growth
            d_imass = (1/((1+self%cht_q_J)*self%cht_q_J))*(r_mass/i_mass)*E_g                                       ! calculate part of growth allocated to structure
            d_rmass = (1-(1/((1+self%cht_q_J)*self%cht_q_J))*(r_mass/i_mass))*E_g - g_mass                          ! calculate part of growth allocated to reserves
        elsewhere (E_g<0) ! cohorts with negative growth
            d_imass = 0.0_rk                                                                                        ! no change in structural tissues when growth is negative
            d_rmass = E_g - g_mass                                                                                  ! loss of reserves due to negative growth
        endwhere
        
        _SET_SURFACE_ODE_(self%id_N,d_N)                                                                            ! update zooplankton ODE
        _SET_SURFACE_ODE_(self%id_i_mass,d_imass)                                                                   ! update zooplankton ODE
        _SET_SURFACE_ODE_(self%id_r_mass,d_rmass)                                                                   ! update zooplankton ODE
        !===========================================================================================================
        ! 2.8 start new cohort
        !-----------------------------------------------------------------------------------------------------------
        _SET_SURFACE_ODE_(self%id_N(year+self%cht_nc_init),F_tot)                                                   ! start new cohort with the produced offspring

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. do zooplankton
    !===============================================================================================================
        d_ZN = self%cht_r*self%cht_Q10_z**((uTm-self%cht_T_ref)/10) * (self%cht_K-Z_N)-Z_N*sum((A_z*N)/(1+H*eta))   ! set zooplankton derivative (K limited growth with growth rate r - minus predation)
        
        _SET_SURFACE_ODE_(self%id_Z_N,d_ZN)                                                                         ! update zooplankton ODE
         
        
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

    
    
    
    
    
    
        
!        if (0) then ! temp code exclusion
!        !===========================================================================================================
!        ! 2.4 do individual cohorts
!        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!        do i=1,nC_fin ! loop over the cohorts
!            if (N(i)>self%cht_ext) then ! if cohort is not empty or extinct
!            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!            ! 2.4.1 fish reproduction
!            !-------------------------------------------------------------------------------------------------------
!                if (Day >= self%cht_R_day .AND. Day <= self%cht_R_day + 1.0_rk .AND. ix_repro==0) then ! if reproduction is due
!                    if (i_mass(i) >= self%cht_x_f .AND. r_mass(i) > self%cht_q_J*i_mass(i)) then ! cohort is mature and non-starving
!                        g_mass(i) = (r_mass(i)-self%cht_q_J*i_mass(i))                                              ! calculate total loss of r_mass due to reproduction (gonad mass)
!                        F_tot = F_tot + g_mass(i)*N(i)*self%cht_k_r/self%cht_w_b                                    ! add offspring from cohort i to total produced offspring
!                    else ! cohort is not ready to reproduce
!                        g_mass(i) = 0.0_rk                                                                          ! no change in r_mass from reproduction for cohort i
!                    endif
!                    ix_repro = int(i/nC_fin)                                                                        ! only calculate reproduction once per reproductive event; stop when all cohorts have been cycled through
!                else ! reproduction is not due
!                    ix_repro=0                                                                                      ! reset reproductive switch                                                                              
!                    g_mass(i) = 0.0_rk                                                                              ! no change in r_mass from reproduction for cohort i
!                    F_tot = 0.0_rk                                                                                  ! no reproduction the rest of the year
!                endif
!            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!            ! 2.4.2 fish dynamics
!            !-------------------------------------------------------------------------------------------------------
!                ! 2.4.2.1 calculate cannibalistic interactions
!                !---------------------------------------------------------------------------------------------------
!                do j=1,nC_fin ! loop over the cohorts                                                                         THE FOLLOWING LOOP (AND POSSIBLY ITS PARENT DO-LOOP) SHOULD BE REWRITTEN TO ELEMENT-WISE MATRIX OPERATIONS
!                    
!                    ! first calculate predation rates by the current cohort on the other cohorts
!                    if (lengths(j)>(lengths(i)*self%cht_delta) .AND. lengths(j)<=(lengths(i)*self%cht_phi)) then ! if cohort j is within the lower size range for predation by cohort i
!                        A_vict(j) = (self%cht_beta*lengths(i)**self%cht_sigma)*     &
!                            (lengths(j)-lengths(i)*self%cht_delta)/((self%cht_phi-self%cht_delta)*lengths(i))       ! calculate attack rate by cohort i on cohort j
!                    elseif (lengths(j)>(lengths(i)*self%cht_phi) .AND. lengths(j)<(lengths(i)*self%cht_eps)) then ! if cohort j is within the upper size range for predation by cohort i
!                        A_vict(j) = (self%cht_beta*lengths(i)**self%cht_sigma)*     &
!                            (self%cht_eps*lengths(i)-lengths(j))/((self%cht_eps-self%cht_phi)*lengths(i))           ! calculate attack rate by cohort i on cohort j
!                    else ! if cohort j is outside the predation window of cohort i
!                        A_vict(j)=0.0                                                                               ! no predation by cohort i on cohort j
!                    endif
!                    
!                    ! then calculate predation on the current cohort by the other cohorts
!                    if (lengths(i)>(lengths(j)*self%cht_delta) .AND. lengths(i)<=(lengths(j)*self%cht_phi)) then ! if cohort i is within the lower size range for predation of cohort j
!                        A_can(j) = (self%cht_beta*lengths(j)**self%cht_sigma)*      &
!                            (lengths(i)-lengths(j)*self%cht_delta)/((self%cht_phi-self%cht_delta)*lengths(j))       ! calculate attack rate by cohort j on cohort i
!                    elseif (lengths(i)>(lengths(j)*self%cht_phi) .AND. lengths(i)<(lengths(j)*self%cht_eps)) then ! if cohort i is within the upper size range for predation of cohort j
!                        A_can(j) = (self%cht_beta*lengths(j)**self%cht_sigma)*      &
!                            (self%cht_eps*lengths(j)-lengths(i))/((self%cht_eps-self%cht_phi)*lengths(j))           ! calculate attack rate by cohort j on cohort i
!                    else ! if cohort i is outside the predation window of cohort j
!                        A_can(j)=0.0                                                                                ! no predation by cohort j on cohort i
!                    endif
!                end do
!                !---------------------------------------------------------------------------------------------------
!                ! 2.4.2.2 scale piscivory rates to temperature
!                !---------------------------------------------------------------------------------------------------
!                A_vict = A_vict*rt_fish_att(i)                                                                      ! temperature correct piscivorous attack rates of cohort i
!                A_can = A_can*rt_fish_att                                                                           ! temperature correct piscivorous attack rates of other cohorts on cohort i
!                !---------------------------------------------------------------------------------------------------
!                ! 2.4.2.3 calculate final vital rates
!                !---------------------------------------------------------------------------------------------------
!                eta_c = sum(A_c*((r_mass+i_mass)*N))                                                                ! calculate total piscivorous encounter rate of cohort i
!                eta = eta_z+eta_c                                                                                   ! calculate total encounter rate
!                
!                mu_c = sum(A_can*N/(1+H*eta))                                                                       ! calculate cannibalistic mortality on cohort i
!                Ing = eta./(1+H.*eta)                                                                               ! total food intake rate
!                E_a = self%cht_k_e*Ing                                                                                     ! assimilation
!                E_g = E_a-E_m                                                                                       ! scope for growth
!                starvindex = CWR<self%cht_q_s*CWI                                                                        ! find which cohorts are starving
!                mu_s = starvindex.*(self%cht_s.*(self%cht_q_s.*CWI./CWR-1))                                                    ! starvation mortality
!                mu_s(isnan(mu_s))=0                                                                                 ! get rid of NaN's from empty cohorts
!                mu_s(mu_s>10)=10                                                                                    ! too high a starvation mortality gives matlab problems with integration
!                mu = self%cht_mu_b + mu_size + mu_s + mu_c                                                               ! total mortality
!                !calculate final mortality rates
!                ! allocate growth to reversible/irreversible mass
!                d_N = 0.0_rk
!                d_imass = 0.0_rk
!                d_rmass = 0.0_rk
!            else ! cohort is empty
!                ! rates of change of epmpty cohorts are zero
!                d_N = 0.0_rk
!                d_imass = 0.0_rk
!                d_rmass = 0.0_rk     
!            endif
!            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!            ! 2.4.3 update individual cohort rates of change
!            !-------------------------------------------------------------------------------------------------------
!            _SET_SURFACE_ODE_(self%id_N(i),d_N)
!            _SET_SURFACE_ODE_(self%id_r_mass(i),d_rmass)
!            _SET_SURFACE_ODE_(self%id_i_mass(i),d_imass)
!        end do
!        endif
