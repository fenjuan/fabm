#include "fabm_driver.h"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
module fish_cohort
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !DESCRIPTION:
    ! Fen, 5th, oct, 2018: the FABM structure fish cohort model, that defined groups of state variables
    !and 
    ! Nicolas to fill out
    
    !USES:
    use fabm_types
    implicit none
    private
    
    !PUBLIC DERIVED TYPE:
    type,extends(type_base_model),public :: type_fish_cohort
       ! variable identifiers---------------------------------------------------------------------------------------
       type (type_surface_state_variable_id), allocatable    ::  id_N(:)
       type (type_surface_state_variable_id), allocatable    ::  id_r_mass(:)
       type (type_surface_state_variable_id), allocatable    ::  id_i_mass(:)
       type (type_surface_state_variable_id)                 ::  id_Z_N
       
       type (type_surface_state_variable_id), allocatable    ::  id_N_mass(:)
       type (type_surface_state_variable_id), allocatable    ::  id_P_mass(:)

       ! environmental dependencies---------------------------------------------------------------------------------
       type (type_dependency_id)                             ::  id_uTm
       type (type_global_dependency_id)                      ::  id_Day
       
       ! temprorary placeholder parameters--------------------------------------------------------------------------
       real(rk)     ::  cht_Tmin,       cht_Tmax,       cht_K,          cht_r,          cht_m,          cht_Q10_z
       real(rk)     ::  cht_PDZoo,      cht_NDZoo
       ! model setup parameters-------------------------------------------------------------------------------------
       integer      ::  cht_nC,         cht_nc_init
       real(rk)     ::  cht_ext
       ! core model parameters--------------------------------------------------------------------------------------
       real(rk)     ::  cht_q_J,        cht_q_A,        cht_theta_m_e,  cht_theta_a_s,  cht_theta_a_e,  cht_T_ref
       real(rk)     ::  cht_lambda_1,   cht_lambda_2,   cht_A_mx,       cht_W_opt,      cht_alpha,      cht_x_mu
       real(rk)     ::  cht_xi_1,       cht_xi_2,       cht_k_e,        cht_rho_1,      cht_rho_2,      cht_q_s
       real(rk)     ::  cht_s,          cht_mu_b,       cht_mu_0,       cht_phi,        cht_delta,      cht_eps
       real(rk)     ::  cht_beta,       cht_sigma,      cht_x_f,        cht_w_b,        cht_theta_m_s,  cht_k_r      
       real(rk)     ::  cht_R_day
       ! stochiometric parameters-----------------------------------------------------------------------------------
       real(rk)     ::  cht_PD_Ref,     cht_ND_Ref,     cht_P_Bone,     cht_D_Bone,     cht_Diss_Eges,  cht_DOMW
       real(rk)     ::  cht_Diss_Mort
       
       contains                                                                                                  
       ! model procedures-------------------------------------------------------------------------------------------
       procedure :: initialize
       procedure :: do_surface
       
    end type type_fish_cohort
    
    ! private data members(API0.92)
    real(rk),parameter :: secs_pr_day=86400.0_rk
    real(rk),parameter :: L_pr_m3 = 1000.0_rk
    
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
contains
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !INTERFACE:
    subroutine initialize(self,configunit)
    
    !INPUT PARAMETERS:
    class (type_fish_cohort), intent(inout), target  :: self
    integer, intent(in)                              :: configunit
    
    !LOCAL VARIABLES:
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
    integer             :: n
    character(len=64)   :: name, longname
    
    !REGISTER MODEL PARAMETERS:
    ! temporary placeholder parameters--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    call self%get_parameter(self%cht_Tmin      ,'min_winter_temp'        , '°C'                 , 'minimum winter temperature'                               , default=1.0_rk                                                    )
    call self%get_parameter(self%cht_Tmax      ,'max_summer_temp'        , '°C'                 , 'maximum summer temperature'                               , default=22.0_rk                                                   )
    call self%get_parameter(self%cht_K         ,'carr_food'              , '[L^-1]'             , 'resource carrying capacity'                               , default=100.0_rk    ,    scale_factor=L_pr_m3                     )   
    call self%get_parameter(self%cht_r         ,'food_growth'            , '[d^-1]'             , 'resource growth rate'                                     , default=0.1_rk      ,    scale_factor=1.0_rk/secs_pr_day          )   
    call self%get_parameter(self%cht_m         ,'food_ind_weight'        , 'gDW'                , 'weight of resource individual'                            , default=7.8e-6_rk                                                 ) ! OW = 3e-5 , NW = 7.8e-6 CHECK
    call self%get_parameter(self%cht_Q10_z     ,'Q10_food'               , '[-]'                , 'Q10 value of zooplankton population growth'               , default=1.799_rk                                                  )
    call self%get_parameter(self%cht_PDZoo     ,'Z_PD_ratio'             , '[-]'                , 'Phosphorous to Dry Weight ratio'                          , default=0.01_rk                                                   )
    call self%get_parameter(self%cht_NDZoo     ,'Z_ND_ratio'             , '[-]'                , 'Nitrogen to Dry Weight ratio'                             , default=0.07_rk                                                   )
   ! model setup parameters------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    call self%get_parameter(self%cht_nC        ,'number_of_cohorts'      , '[-]'                , 'number of cohorts'                                        , default=10                                                        )
    call self%get_parameter(self%cht_nc_init   ,'initial_cohorts'        , '[-]'                , 'number of cohorts present at simulation start '           , default=1                                                         )
    call self%get_parameter(self%cht_ext       ,'extinction_abundance'   , '[L^-1]'             , 'extinction threshold abundance'                           , default=1E-15_rk    ,     scale_factor=L_pr_m3                    )
    ! core model parameters-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    call self%get_parameter(self%cht_q_J       ,'max_jv_con'             , '[-]'                , 'juvenile maximum condition'                               , default=0.74_rk                                                   )
    call self%get_parameter(self%cht_q_A       ,'max_ad_con'             , '[-]'                , 'adult maximum condition'                                  , default=1.37_rk                                                   )
    call self%get_parameter(self%cht_lambda_1  ,'length-weight_scalar'   , '(mm gDW^(-lambda_2)', 'allometric length-weight scalar'                          , default=74.0288_rk                                                ) ! OW = 48.3 , NW = 74.0288 (original lambda_1=48.3 * 10**(lambda_2*log10(1/0.26=dw:ww)) )) CHECK
    call self%get_parameter(self%cht_lambda_2  ,'length-weight_exponent' , '[-]'                , 'allometric length-weight exponent'                        , default=0.317_rk                                                  )   
    call self%get_parameter(self%cht_A_mx      ,'max_attack_rate'        , 'L d^-1'             , 'max attack rate'                                          , default=3E4_rk      ,    scale_factor=1.0_rk/(L_pr_m3*secs_pr_day))   
    call self%get_parameter(self%cht_W_opt     ,'att_opt_size'           , 'gDW'                , 'size for optimal attack rate on zooplankton (WW)'         , default=2.1320_rk                                                 ) ! OW = 8.2 , NW = 2.1320 (see conversion for lambda_1) CHECK
    call self%get_parameter(self%cht_alpha     ,'att_rate_exp'           , '[-]'                , 'allometric exponent for attack rate on zooplankton'       , default=0.62_rk                                                   )   
    call self%get_parameter(self%cht_x_mu      ,'mort_cha_size'          , 'gDW'                , 'shape parameter for size-dependent mortality'             , default=0.26_rk                                                   ) ! OW = 1.0 , NW = 0.26 CHECK
    call self%get_parameter(self%cht_xi_1      ,'handling_time_scalar'   , 'd gDW^-(1+xi_2)'    , 'allometric scalar for handling time'                      , default=6.5462_rk   ,    scale_factor=secs_pr_day                 ) ! OW = 5.0 , NW = 6.5462 (changed using dw:wW ratios for fish, then divided by zooplankton dw:ww ratio (assumed identical to the one of fish)) CHECK
    call self%get_parameter(self%cht_xi_2      ,'handling_time_exp'      , '[-]'                , 'allometric exponent for handling time'                    , default=-0.8_rk                                                   )   
    call self%get_parameter(self%cht_k_e       ,'ass_efficiency'         , '[-]'                , 'assimilation efficiency'                                  , default=0.61_rk                                                   )   
    call self%get_parameter(self%cht_rho_1     ,'metab_scalar'           , '(gDW^(1-rho_2) d^-1', 'allometric scalar for basic metabolism'                   , default=0.0242_rk   ,    scale_factor=1.0_rk/secs_pr_day          ) ! OW = 0.033 , NW = 0.0242 (changed using dw:wW ratios for fish, then multiplied by the dw:ww ratio to convert ww loss into dw loss) CHECK
    call self%get_parameter(self%cht_rho_2     ,'metab_exp'              , '[-]'                , 'allometric exponent for basic metabolism'                 , default=0.77_rk                                                   )   
    call self%get_parameter(self%cht_q_s       ,'starv_con'              , '[-]'                , 'starvation condition'                                     , default=0.2_rk                                                    )   
    call self%get_parameter(self%cht_s         ,'starv_coeff'            , 'd^-1'               , 'starvation coefficient'                                   , default=0.2_rk      ,    scale_factor=1.0_rk/secs_pr_day          )   
    call self%get_parameter(self%cht_mu_b      ,'background_mort'        , 'd^-1'               , 'background mortality'                                     , default=0.0025_rk   ,    scale_factor=1.0_rk/secs_pr_day          )   
    call self%get_parameter(self%cht_mu_0      ,'mort_scalar_size'       , 'd^-1'               , 'scaling constant for size-dependent mortality'            , default=0.034_rk    ,    scale_factor=1.0_rk/secs_pr_day          )   
    call self%get_parameter(self%cht_eps       ,'pred_size_max'          , '[-]'                , 'maximum victim/cannibal size ratio'                       , default=0.45_rk                                                   )   
    call self%get_parameter(self%cht_phi       ,'pred_size_opt'          , '[-]'                , 'optimal victim/cannibal size ratio'                       , default=0.2_rk                                                    )   
    call self%get_parameter(self%cht_delta     ,'pred_size_min'          , '[-]'                , 'minimum victim/cannibal size ratio'                       , default=0.05_rk                                                   )   
    call self%get_parameter(self%cht_beta      ,'pred_max_att_scalar'    , 'L d^-1 mm^sigma'    , 'cannibalistic voracity (allometric scalar for piscivory)' , default=0.0_rk      ,    scale_factor=1.0_rk/(L_pr_m3*secs_pr_day))   
    call self%get_parameter(self%cht_sigma     ,'pred_max_att_exp'       , '[-]'                , 'allometric exponent for cannibalism (piscivory)'          , default=0.6_rk                                                    )   
    call self%get_parameter(self%cht_x_f       ,'mat_irr_mass'           , 'gDW'                , 'maturation irreversible mass'                             , default=1.1960_rk                                                 ) ! OW = 4.6 , NW = 1.1960 CHECK
    call self%get_parameter(self%cht_w_b       ,'egg_size'               , 'gDW'                , 'egg size'                                                 , default=4.68e-4_rk                                                ) ! OW = 1.8e-3 , NW = 4.68e-4 CHECK
    call self%get_parameter(self%cht_k_r       ,'repro_eff'              , '[-]'                , 'gonad-offspring conversion factor'                        , default=0.5_rk                                                    )   
    call self%get_parameter(self%cht_theta_m_s ,'Q10_metab_scalar'       , '[gDW^-theta_m_e]'   , 'allometric scalar of metabolism Q10'                      , default=2.2067_rk                                                 ) ! OW = 2.0 , NW = 2.2067 CHECK
    call self%get_parameter(self%cht_theta_m_e ,'Q10_metab_exp'          , '[-]'                , 'allometric exponent of metabolism Q10'                    , default=0.073_rk                                                  )   
    call self%get_parameter(self%cht_theta_a_s ,'Q10_att_scalar'         , '[gDW^-theta_a_e]'   , 'allometric scalar of attack Q10'                          , default=3.085176_rk                                               ) ! OW = 2.8 , NW = 3.085176 CHECK
    call self%get_parameter(self%cht_theta_a_e ,'Q10_att_exp'            , '[-]'                , 'allometric exponent of attack Q10 '                       , default=0.072_rk                                                  )   
    call self%get_parameter(self%cht_T_ref     ,'temp_ref_fish'          , '°C'                 , 'reference temperature for perch '                         , default=20.0_rk                                                   )   
    call self%get_parameter(self%cht_R_day     ,'repro_day'              , '[d]'                , 'reproduction day of year for fish '                       , default=152.0_rk                                                  )
    ! stochiometric parameters----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    call self%get_parameter(self%cht_PD_Ref    ,'P_DW_ratio_ref'         , '[-]'                , 'reference P/DW ratio of fish'                             , default=0.022_rk                                                  )
    call self%get_parameter(self%cht_ND_Ref    ,'N_DW_ratio_ref'         , '[-]'                , 'reference N/DW ratio of fish'                             , default=0.1_rk                                                    )
    call self%get_parameter(self%cht_D_Bone    ,'DW_bone_fraction'       , '[-]'                , 'fraction of fish DW fixed in bones and scales'            , default=0.35_rk                                                   )
    call self%get_parameter(self%cht_P_Bone    ,'P_bone_faction'         , '[-]'                , 'fraction of fish P fixed in bones and scales'             , default=0.5_rk                                                    )
    call self%get_parameter(self%cht_Diss_Eges ,'soluble_egestion'       , '[-]'                , 'soluble nutrient fraction of fish egested food'           , default=0.25_rk                                                   )
    call self%get_parameter(self%cht_Diss_Mort ,'soluble_carcass'        , '[-]'                , 'sol. nut. fract. of dead fish (excl bones and scales)'    , default=0.1_rk                                                    )
    call self%get_parameter(self%cht_DOMW      ,'DOMW_fish'              , '[-]'                , 'dissolved organic matter fraction from fish'              , default=0.5_rk                                                    )
    
    !STATE VARIABLE POINTERS:
    ! setup zooplankton
    call self%register_state_variable(self%id_Z_N,'prey','[m^-3]','food concentration',initial_value=self%cht_K,minimum=0.0_rk)

    ! set up fish cohorts 
    allocate(self%id_N(self%cht_nC+self%cht_nc_init))
    allocate(self%id_r_mass(self%cht_nC+self%cht_nc_init))
    allocate(self%id_i_mass(self%cht_nC+self%cht_nc_init))
    
    allocate(self%id_N_mass(self%cht_nC+self%cht_nc_init))
    allocate(self%id_P_mass(self%cht_nC+self%cht_nc_init))

    
    do n=1,self%cht_nC+self%cht_nc_init
       write (name,"(A1,I02.2)") "N", n-1
       write (longname, "(A9, I02.2)") "abundance", n-1
       call self%register_state_variable(self%id_N(n),name, '[m^-2]', longname,initial_value= 0._rk, minimum=0.0_rk)
       write (name,"(A6,I02.2)") "r_mass", n-1
       write (longname, "(A18, I02.2)") "reversible mass DW", n-1
       call self%register_state_variable(self%id_r_mass(n),name,'[gDW]',longname,initial_value= self%cht_q_J*self%cht_w_b/(1+self%cht_q_J), minimum=0.0_rk)
       write (name,"(A6,I02.2)") "i_mass", n-1
       write (longname, "(A20, I02.2)") "irreversible mass DW", n-1
       call self%register_state_variable(self%id_i_mass(n),name,'[gDW]',longname,initial_value=self%cht_w_b/(1+self%cht_q_J), minimum=0.0_rk)
        
       write (name,"(A6,I02.2)") "N_mass", n-1
       write (longname, "(A17, I02.2)") "individual mass n", n-1
       call self%register_state_variable(self%id_N_mass(n),name,'[gN]',longname,initial_value= self%cht_w_b*self%cht_ND_Ref, minimum=0.0_rk)
       write (name,"(A6,I02.2)") "P_mass", n-1
       write (longname, "(A17, I02.2)") "individual mass P", n-1
       call self%register_state_variable(self%id_P_mass(n),name,'[gP]',longname,initial_value= self%cht_w_b*self%cht_PD_Ref, minimum=0.0_rk)
    end do
    
    !  register environmental dependencies
    call self%register_dependency(self%id_uTm,    standard_variables%temperature)  ! eventually, this shoud be replaced
    call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
    
    return
    end subroutine initialize
    
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !INTERFACE:
    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    
    !INPUT PARAMETERS:
    class (type_fish_cohort),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_SURFACE_
    
    !LOCAL VARIABLES:
    ! state variables ----------------------------------------------------------------------------------------------
    real(rk), DIMENSION(self%cht_nC+self%cht_nc_init)   :: N,           r_mass,     i_mass,     N_mass,     P_mass
    real(rk)                                            :: Z_N
    ! derivatives --------------------------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: d_N,         d_rmass,    d_imass,    d_Nmass,    d_Pmass
    real(rk)                                            :: d_ZN
    ! local variables which carry over between time steps ----------------------------------------------------------                          
    real(rk), allocatable, save                         :: g_mass(:)                                                ! g_mass should be saved between time step in another way before multiple fish instances are run in FABM, e.g. as a derived variable
    ! local cohort vector variables --------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: lengths,     st_mass,    H,          A_z,        mu_size
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: eta_z,       Q10_f_m,    Q10_f_a,    rt_f_m,     rt_f_a
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: mu_c,        E_a,        Ing,        E_g,        mu_s
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: mu,          E_m,        eta
    ! local cohort matrix variables --------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nC+self%cht_nc_init)                   :: A_c,         L_c,        L_v,        eta_c
    ! local cohort single within-loop variables --------------------------------------------------------------------
    real(rk)                                            :: F_tot
    ! local stoichiometric variables--------------------------------------------------------------------------------
    real(rk)                                            :: NH4_tot,     PO4_tot,    DW_tot,     DW_POM,     DW_DOM
    real(rk)                                            :: OP_tot,      OP_POM,     OP_DOM,     ON_POM,     ON_DOM      
    real(rk)                                            :: ON_tot
    ! local stoichiometric vector variables-------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: PD,          ND,         E_m_corr,   Ing_zoo,    ke_P_zoo
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: P_ass_zoo,   ke_P_pis,   P_ass,      ke_N_zoo,   N_ass_zoo
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: ke_N_pis,    N_ass,      P_exc,      N_exc,      D_ege
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: Ing_P,       P_ege,      Ing_N,      N_ege,      NH4_ege
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: NH4_mor,     PO4_ege,    PO4_mor,    DW_mor,     OP_mor   
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: OP_ege,      ON_ege,     ON_mor
    ! local stoichiometric matrix variables ------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nC+self%cht_nc_init)                   :: Ing_pis, P_ass_pis,  N_ass_pis
    ! carriers for environment dependencies and logical switches ---------------------------------------------------
    real(rk)                                            :: uTm,     Day
    integer                                             :: i,       j,        ix_repro, nC_fin
    integer                                             :: tcheck=0
    integer, save                                       :: year=0
    
    !---------------------------------------------------------------------------------------------------------------
    ! Enter spatial loops (if any)
    _HORIZONTAL_LOOP_BEGIN_
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1. retrieve variable values, update trackers & variables:
    !===============================================================================================================
        ! 1.1 get cohort state variables
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _GET_HORIZONTAL_(self%id_Z_N,Z_N)                                                                           ! get zooplankton (food) concentration
        
        _GET_HORIZONTAL_(self%id_N,N)                                                                               ! get cohort abundance
        _GET_HORIZONTAL_(self%id_r_mass,r_mass)                                                                     ! get cohort reversible mass (reserves+gonads)
        _GET_HORIZONTAL_(self%id_i_mass,i_mass)                                                                     ! get cohort irreversible mass (structural tissues, organs)
        _GET_HORIZONTAL_(self%id_N_mass,N_mass)                                                                     ! get cohort total Nitrogen mass (structure+reserves+gonads)
        _GET_HORIZONTAL_(self%id_P_mass,P_mass)                                                                     ! get cohort total Phosphorous mass (structure+reserves+gonads)
        !===========================================================================================================
        ! 1.2 retrieve environmental dependencies and update local trackers & switches
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
            P_mass = 0.0_rk                                                                                         ! treat extinct cohorts as empty
            N_mass = 0.0_rk                                                                                         ! treat extinct cohorts as empty
        endwhere
        !===========================================================================================================
        ! 2.1 calculate size-specific rates and parameters
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        st_mass = i_mass*(1.0_rk+self%cht_q_J)                                                                      ! calculate fish standard mass
        lengths = self%cht_lambda_1*st_mass**self%cht_lambda_2                                                      ! calculate fish lengths
        where (N>0.0_rk) ! non-empty cohorts
            H = self%cht_xi_1*st_mass**self%cht_xi_2                                                                ! calculate size-specific handling times (digestion)
        elsewhere ! empty cohorts
            H = 0.0_rk
        endwhere
        A_z = self%cht_A_mx*((st_mass/self%cht_W_opt)*exp(1.0_rk-st_mass/self%cht_W_opt))**self%cht_alpha           ! calculate size-specific attack rate on zooplankton
        mu_size = self%cht_mu_0*exp(-i_mass/self%cht_x_mu)                                                          ! calculate size-spcific mortality
        E_m = self%cht_rho_1*(i_mass+r_mass)**self%cht_rho_2                                                        ! calculate size-specific respiration
        !===========================================================================================================
        ! 2.2 calculate temperature correction factors
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        where (N>0.0_rk) ! only do calculations on the cohorts that are non-empty
            Q10_f_m = self%cht_theta_m_s*(r_mass+i_mass)**self%cht_theta_m_e                                        ! calculate size-specific Q10 for fish metabolism
            Q10_f_a = self%cht_theta_a_s*(r_mass+i_mass)**self%cht_theta_a_e                                        ! calculate size-specific Q10 for fish feeding
            rt_f_m = Q10_f_m**((uTm-self%cht_T_ref)/10.0_rk)                                                        ! calculate temperature correction factor for metabolism
            rt_f_a = Q10_f_a**((uTm-self%cht_T_ref)/10.0_rk)                                                        ! calculate temperature correction factor for feeding
        elsewhere ! empty cohorts
            rt_f_m = 0.0_rk                                                                                         ! avoid getting inf in empty cohorts
            rt_f_a = 0.0_rk                                                                                         ! avoid getting inf in empty cohorts
        endwhere
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
        where (isnan(H)) H = 0.0_rk ! get rid of NANs
        A_z = A_z*rt_f_a                                                                                            ! temperature correct attack rate on zooplankton
        E_m = E_m*rt_f_m                                                                                            ! temperature correct basal metabolism
        A_c = A_c*spread(rt_f_a, DIM=1,  NCOPIES=nC_fin)                                                            ! temperature correct piscivorous attack rates according to attackers temperature correction factor
        !===========================================================================================================
        ! 2.5 stochiometric correction of respiration
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                
        where (N>0.0_rk) ! non-empty cohorts
            PD = P_mass / (r_mass*i_mass)                                                                         ! calculate current cohort P/D ratios
            ND = N_mass / (r_mass*i_mass)                                                                         ! calculate current cohort N/D ratios
        elsewhere ! empty cohorts
            PD = 1.0_rk                                                                                             ! prevents INF and NAN results later
            ND = 1.0_rk                                                                                             ! prevents INF and NAN results later
        endwhere
        E_m_corr = E_m*max(self%cht_PD_Ref/PD,self%cht_ND_Ref/ND)                                                   ! correct respiration to maintain stochiometric equilibrium
        !===========================================================================================================
        ! 2.6 fish reproduction
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        if (Day >= self%cht_R_day - 1.0_rk .AND. Day < self%cht_R_day .AND. ix_repro==0) then ! first time step of reproduction event
            where (i_mass >= self%cht_x_f .AND. r_mass > self%cht_q_J*i_mass) ! where cohorts are mature and non-starving
                g_mass = (r_mass-self%cht_q_J*i_mass)/secs_pr_day                                                   ! calculate total loss of r_mass due to reproduction (gonad mass), which is equivalent to the daily loss rate when reproduction lasts one day
                A_z = 0.0_rk                                                                                        ! reproducing individuals are not feeding
            elsewhere ! cohort is not ready to reproduce
                g_mass = 0.0_rk                                                                                     ! no change in r_mass from reproduction for immature or starving cohorts
            endwhere
            where  (spread(i_mass, DIM=1,  NCOPIES=nC_fin) >= self%cht_x_f .AND.    &
                    spread(r_mass, DIM=1,  NCOPIES=nC_fin) >  self%cht_q_J*spread(i_mass, DIM=1,  NCOPIES=nC_fin)) ! where cohorts are mature and non-starving (matrix edition)
                    A_c = 0.0_rk                                                                                    ! reproducing individuals are not feeding
            endwhere
            A_z(year+self%cht_nc_init) = 0.0_rk                                                                     ! new cohort only starts feeding after the reproductive event
            A_c(1:nC_fin,year+self%cht_nc_init) = 0.0_rk                                                            ! new cohort only starts feeding after the reproductive event
            ix_repro = 1                                                                                            ! only calculate reproduction at the start of the reproductive event
        elseif (Day < self%cht_R_day - 1.0_rk .OR. Day >= self%cht_R_day) then ! rest of the year
            ix_repro=0                                                                                              ! reset reproductive switch                                                                              
            g_mass = 0.0_rk                                                                                         ! no reproduction outside reproductive event
        endif
        F_tot = sum(g_mass*N*self%cht_k_r/self%cht_w_b)                                                             ! calculate total number of produced offspring (zero when g_mass is zero)
        !===========================================================================================================
        ! 2.7 calculate final vital rates
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        eta_z = A_z*Z_N*self%cht_m                                                                                  ! calculate encounter rate with zooplankton
        eta_c = A_c*spread((r_mass+i_mass)*N, DIM=2, NCOPIES=nC_fin)                                                ! calculate encounter rates between piscivores  and fish prey
        eta = eta_z+sum(eta_c,1)                                                                                    ! calculate total encounter rate of piscivores with prey
        mu_c = sum(A_c*spread(N/(1.0_rk+H*eta), DIM=1, NCOPIES=nC_fin),2)                                           ! calculate piscivorous mortality on cohorts
        Ing = eta/(1.0_rk+H*eta)                                                                                    ! calculate total food intake rate
        E_a = self%cht_k_e*Ing                                                                                      ! calculate assimilation
        E_g = E_a-E_m_corr                                                                                          ! final growth (or degrowth)
        where (r_mass<self%cht_q_s*i_mass) ! find which cohorts are starving
            mu_s = (self%cht_s*(self%cht_q_s*i_mass/r_mass-1.0_rk))                                                 ! calculate starvation mortality
        elsewhere
            mu_s = 0.0_rk                                                                                           ! otherwise no starvation mortality
        endwhere
        where (mu_s>10.0_rk/secs_pr_day) ! find where starvation mortality is very high
            mu_s=10.0_rk/secs_pr_day                                                                                ! too high a starvation mortality might produce problems with integration
        endwhere
        mu = self%cht_mu_b + mu_size + mu_s + mu_c                                                                  ! total mortality
        !===========================================================================================================
        ! 2.8 set ODEs
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        d_N = -mu*N                                                                                                 ! rate of change of cohort abundances (mortality)
        where (E_g>=0.0_rk .AND. i_mass>=self%cht_x_f .AND. N>0.0_rk) ! mature non-empty cohorts with positive growth
            d_imass = (1.0_rk/((1.0_rk+self%cht_q_A)*self%cht_q_A))*(r_mass/i_mass)*E_g                             ! calculate part of growth allocated to structure
            d_rmass = (1.0_rk-(1.0_rk/((1.0_rk+self%cht_q_A)*self%cht_q_A))*(r_mass/i_mass))*E_g - g_mass           ! calculate part of growth allocated to reserves and gonads
        elsewhere (E_g>=0.0_rk .AND. i_mass<self%cht_x_f .AND. N>0.0_rk) ! juvenile non-empty cohorts with positive growth
            d_imass = (1.0_rk/((1.0_rk+self%cht_q_J)*self%cht_q_J))*(r_mass/i_mass)*E_g                             ! calculate part of growth allocated to structure
            d_rmass = (1.0_rk-(1.0_rk/((1.0_rk+self%cht_q_J)*self%cht_q_J))*(r_mass/i_mass))*E_g - g_mass           ! calculate part of growth allocated to reserves
        elsewhere (E_g<0.0_rk) ! cohorts with negative growth
            d_imass = 0.0_rk                                                                                        ! no change in structural tissues when growth is negative
            d_rmass = E_g - g_mass                                                                                  ! loss of reserves due to negative growth
        elsewhere (N==0.0_rk) ! empty cohorts
            d_imass = 0.0_rk                                                                                        ! rate of change is zero for empty cohorts
            d_rmass = 0.0_rk                                                                                        ! rate of change is zero for empty cohorts
        endwhere

        _SET_SURFACE_ODE_(self%id_N,d_N)                                                                            ! update abundance ODE
        _SET_SURFACE_ODE_(self%id_i_mass,d_imass)                                                                   ! update irreversible mass ODE
        _SET_SURFACE_ODE_(self%id_r_mass,d_rmass)                                                                   ! update reversible mass ODE
        !===========================================================================================================
        ! 2.9 start new cohort
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _SET_SURFACE_ODE_(self%id_N(year+self%cht_nc_init),F_tot)                                                   ! start new cohort with the produced offspring
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. do zooplankton - temporary section
    !===============================================================================================================
        d_ZN = self%cht_r*self%cht_Q10_z**((uTm-self%cht_T_ref)/10) * (self%cht_K-Z_N)      &
            - Z_N*sum((A_z*N)/(1.0_rk+H*eta))                                                                       ! set zooplankton derivative (K limited growth with growth rate r - minus predation)
        
        _SET_SURFACE_ODE_(self%id_Z_N,d_ZN)                                                                         ! update zooplankton ODE
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. do stochiometry & and update external processes
    !===============================================================================================================
        Ing_zoo = eta_z/(1.0_rk+H*eta)                                                                              ! calculate zooplankton ingestion rates
        Ing_pis = eta_c/spread(1.0_rk+H*eta, DIM=1, NCOPIES=nC_fin)                                                 ! calculate piscivorous ingestion rates
        !===========================================================================================================
        ! 4.1 assimilation
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.1.1 phosphorous
            !-------------------------------------------------------------------------------------------------------
            ke_P_zoo = min(1.0_rk,self%cht_PD_Ref/self%cht_PDZoo * self%cht_k_e)                                    ! calculate P assimilation efficiencies of zooplanktivory
            P_ass_zoo = ke_P_zoo * Ing_zoo * self%cht_PDZoo                                                         ! calculate P assimilation rates from zooplanktivory
            
            ke_P_pis = min(1.0_rk,self%cht_PD_Ref/PD * self%cht_k_e)                                                ! calculate P assimilation efficiencies of piscivory
            P_ass_pis = spread(ke_P_pis, DIM=2, NCOPIES=nC_fin) * Ing_pis * spread(PD, DIM=2, NCOPIES=nC_fin)       ! calculate P assimilation rates from piscivory

            P_ass = P_ass_zoo+sum(P_ass_pis,1)                                                                      ! calculate total P assimilation rates
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.1.2 nitrogen
            !-------------------------------------------------------------------------------------------------------
            ke_N_zoo = min(1.0_rk,self%cht_ND_Ref/self%cht_NDZoo * self%cht_k_e)                                    ! calculate N assimilation efficiencies of zooplanktivory
            N_ass_zoo = ke_N_zoo * Ing_zoo * self%cht_NDZoo                                                         ! calculate N assimilation rates from zooplanktivory
            
            ke_N_pis = min(1.0_rk,self%cht_ND_Ref/ND * self%cht_k_e)                                                ! calculate N assimilation efficiencies of piscivory
            N_ass_pis = spread(ke_N_pis, DIM=2, NCOPIES=nC_fin) * Ing_pis * spread(ND, DIM=2, NCOPIES=nC_fin)       ! calculate N assimilation rates from piscivory

            N_ass = N_ass_zoo+sum(N_ass_pis,1)                                                                      ! calculate total N assimilation rates
        !===========================================================================================================
        ! 4.2 excretion
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.2.1 phosphorous
            !-------------------------------------------------------------------------------------------------------
            P_exc = (PD/self%cht_PD_Ref) * E_m                                                                      ! calculate rates of P excretion
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.2.2 nitrogen
            !-------------------------------------------------------------------------------------------------------
            N_exc = (ND/self%cht_ND_Ref) * E_m                                                                      ! calculate rates of N excretion
        !===========================================================================================================
        ! 4.3 egestion
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.1 dry weight
            !-------------------------------------------------------------------------------------------------------
            D_ege = Ing - E_a                                                                                       ! calculate egestion rate
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.2 phosphorous
            !-------------------------------------------------------------------------------------------------------
            Ing_P = Ing_zoo * self%cht_PDZoo + sum(Ing_pis * spread(PD, DIM=2, NCOPIES=nC_fin),1)                   ! calculate total P ingestion rates
            P_ege = Ing_P - P_ass                                                                                   ! calculate egested P
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.3 nitrogen
            !-------------------------------------------------------------------------------------------------------
            Ing_N = Ing_zoo * self%cht_NDZoo + sum(Ing_pis * spread(ND, DIM=2, NCOPIES=nC_fin),1)                   ! calculate total P ingestion rates
            N_ege = Ing_N - N_ass                                                                                   ! calculate egested P
        !===========================================================================================================
        ! 4.4 set nutrient ODEs for fish
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        d_Pmass = (P_ass-P_exc) - g_mass*self%cht_PD_Ref                                                            ! calculate rates of change in P content for individual fish
        d_Nmass = (N_ass-N_exc) - g_mass*self%cht_ND_Ref                                                            ! calculate rates of change in N content for individual fish
        
        _SET_SURFACE_ODE_(self%id_P_mass,d_Pmass)                                                                   ! update P_mass ODE
        _SET_SURFACE_ODE_(self%id_N_mass,d_Nmass)                                                                   ! update N_mass ODE
        !===========================================================================================================
        ! 4.5 update external processes
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.1 NH4
            !-------------------------------------------------------------------------------------------------------
            NH4_ege = self%cht_Diss_Eges * N_ege                                                                    ! calculate NH4 flux from egestion
            NH4_mor = self%cht_Diss_Mort * (N_mass-self%cht_D_Bone*self%cht_ND_Ref*st_mass)*mu                      ! calculate per capita NH4 flux from mortality - bone mass calculated from standard mass, since can't be starved
            
            NH4_tot = sum((N_exc + NH4_ege + NH4_mor)*N)                                                            ! calculate total population NH4 flux to the water column
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.2 PO4
            !-------------------------------------------------------------------------------------------------------
            PO4_ege = self%cht_Diss_Eges * P_ege                                                                    ! calculate PO4 flux from egestion
            PO4_mor = self%cht_Diss_Mort * (P_mass-self%cht_P_Bone*self%cht_PD_Ref*st_mass)*mu                      ! calculate per capita PO4 flux from mortality - bone mass calculated from standard mass, since can't be starved
            
            PO4_tot = sum((P_exc + PO4_ege + PO4_mor)*N)                                                            ! calculate total population PO4 flux to the water column
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.3 organic matter DW
            !-------------------------------------------------------------------------------------------------------
            DW_mor = ((r_mass+i_mass)-self%cht_D_Bone*st_mass)*mu                                                   ! calculate per capita DW flux from mortality - bone mass calculated from standard mass, since can't be starved
            
            DW_tot = sum((DW_mor + D_ege)*N)                                                                        ! calculate total DW flux to the water column
            DW_POM =DW_tot*(1.0_rk-self%cht_DOMW)                                                                   ! particulate fraction of DW_tot
            DW_DOM = DW_tot*self%cht_DOMW                                                                           ! dissolved fraction of DW_tot
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.4 organic matter P
            !-------------------------------------------------------------------------------------------------------
            OP_mor = (1.0_rk-self%cht_Diss_Mort)*(P_mass-self%cht_P_Bone*self%cht_PD_Ref*st_mass)*mu                ! calculate per capita organic P flux from mortality - bone mass calculated from standard mass, since can't be starved
            OP_ege = P_ege*(1.0_rk-self%cht_Diss_Eges)                                                              ! calculate organic P flux from egestion
            
            OP_tot = sum((OP_mor + OP_ege)*N)                                                                       ! calculate total organic P flux to the water column
            OP_POM = OP_tot*(1.0_rk-self%cht_DOMW)                                                                  ! particulate fraction of OP_tot
            OP_DOM = OP_tot*self%cht_DOMW                                                                           ! dissolved fraction of OP_tot
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.4 organic matter N
            !-------------------------------------------------------------------------------------------------------
            ON_mor = (1.0_rk-self%cht_Diss_Mort)*(N_mass-self%cht_D_Bone*self%cht_ND_Ref*st_mass)*mu                ! calculate per capita organic N flux from mortality - bone mass calculated from standard mass, since can't be starved
            ON_ege = N_ege*(1.0_rk-self%cht_Diss_Eges)                                                              ! calculate organic N flux from egestion
            
            ON_tot = sum((ON_mor + ON_ege)*N)                                                                       ! calculate total organic N flux to the water column
            ON_POM = ON_tot*(1.0_rk-self%cht_DOMW)                                                                  ! particulate fraction of ON_tot
            ON_DOM = ON_tot*self%cht_DOMW                                                                           ! dissolved fraction of ON_tot
        
    _HORIZONTAL_LOOP_END_

    end subroutine do_surface

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
end module fish_cohort
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------