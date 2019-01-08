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
       type (type_surface_state_variable_id), allocatable             ::  id_N(:)
       type (type_surface_state_variable_id), allocatable             ::  id_r_mass(:)
       type (type_surface_state_variable_id), allocatable             ::  id_i_mass(:)
                                                                      
       type (type_surface_state_variable_id), allocatable             ::  id_N_mass(:)
       type (type_surface_state_variable_id), allocatable             ::  id_P_mass(:)
              
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_g_mass_out(:)
       type (type_horizontal_dependency_id), allocatable              ::  id_g_mass_in(:)
       
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_dist_out(:,:)
       type (type_horizontal_dependency_id), allocatable              ::  id_dist_in(:,:)
       
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_prey_loss_DW(:)
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_prey_loss_P(:)
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_prey_loss_N(:)
       type (type_horizontal_diagnostic_variable_id), allocatable     ::  id_pas_rates(:,:)
 
       ! environmental dependencies---------------------------------------------------------------------------------
       type (type_horizontal_dependency_id), allocatable     ::  id_ZD(:)
       type (type_horizontal_dependency_id), allocatable     ::  id_ZN(:)
       type (type_horizontal_dependency_id), allocatable     ::  id_ZP(:)
       type (type_horizontal_dependency_id), allocatable     ::  id_dz(:)
       type (type_horizontal_dependency_id), allocatable     ::  id_uTm(:)
       
       type (type_global_dependency_id)                      ::  id_Day
       
       ! temprorary placeholder parameters--------------------------------------------------------------------------
       real(rk)     ::  cht_m,          cht_Q10_z
       real(rk)     ::  cht_PDZoo,      cht_NDZoo
       ! model setup parameters-------------------------------------------------------------------------------------
       integer      ::  cht_nC,         cht_nc_init,    cht_nlev
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
    integer             :: i_nc, i_z, i_pas
    character(len=64)   :: index_nc, index_z, index_pas
    
    !REGISTER MODEL PARAMETERS:
    ! model setup parameters------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    call self%get_parameter(self%cht_nC        ,'number_of_cohorts'      , '[-]'                , 'number of cohorts'                                        , default=10                                                        )
    call self%get_parameter(self%cht_nc_init   ,'initial_cohorts'        , '[-]'                , 'number of cohorts present at simulation start '           , default=1                                                         )
    call self%get_parameter(self%cht_ext       ,'extinction_abundance'   , '[L^-1]'             , 'extinction threshold abundance'                           , default=1E-15_rk    ,     scale_factor=L_pr_m3                    )
    call self%get_parameter(self%cht_nlev      ,'number_of_depth_lvls'   , '[-]'                , 'number of depth levels'                                                                                                       )
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
    ! setup depth-specific variables
    allocate(self%id_ZD(self%cht_nlev))
    allocate(self%id_ZN(self%cht_nlev))
    allocate(self%id_ZP(self%cht_nlev))
    allocate(self%id_dz(self%cht_nlev))
    allocate(self%id_uTm(self%cht_nlev))
    allocate(self%id_pas_rates(self%cht_nlev, 9))
    allocate(self%id_prey_loss_DW(self%cht_nlev))
    allocate(self%id_prey_loss_P(self%cht_nlev))
    allocate(self%id_prey_loss_N(self%cht_nlev))
    
    do i_z=1,self%cht_nlev
        write(index_z,'(i0)') i_z
        call self%register_horizontal_dependency(self%id_ZD(i_z),'ZD_at_z'//trim(index_z),'gDW m^-3','zooplankton dry weight concentration at depth interval '//trim(index_z))
        call self%register_horizontal_dependency(self%id_ZN(i_z),'ZN_at_z'//trim(index_z),'gN m^-3','zooplankton nitrogen concentration at depth interval '//trim(index_z))
        call self%register_horizontal_dependency(self%id_ZP(i_z),'ZP_at_z'//trim(index_z),'gP m^-3','zooplankton phosphorous concentration at depth interval '//trim(index_z))
        
        call self%register_horizontal_dependency(self%id_dz(i_z),'dz_at_z'//trim(index_z),'m','layer depth at depth interval '//trim(index_z))
        call self%register_horizontal_dependency(self%id_uTm(i_z),'uTm_at_z'//trim(index_z),'*C','temperature at depth interval '//trim(index_z))
        
        call self%register_diagnostic_variable(self%id_prey_loss_DW(i_z),'prey_lossDW_z'//trim(index_z),'-','prey DW loss rate at depth '//trim(index_z),act_as_state_variable=.true.,domain=domain_surface,output=output_none)
        call self%register_diagnostic_variable(self%id_prey_loss_P(i_z),'prey_lossP_z'//trim(index_z),'-','prey P loss rate at depth '//trim(index_z),act_as_state_variable=.true.,domain=domain_surface,output=output_none)
        call self%register_diagnostic_variable(self%id_prey_loss_N(i_z),'prey_lossN_z'//trim(index_z),'-','prey N loss rate at depth '//trim(index_z),act_as_state_variable=.true.,domain=domain_surface,output=output_none)
        
        do i_pas=1,9
            write(index_pas,'(i0)') i_pas
            call self%register_diagnostic_variable(self%id_pas_rates(i_z,i_pas),'pas_rate_var_'//trim(index_pas)//'_z'//trim(index_z),'-','rate of change of passive variable '//trim(index_pas)//' at depth '//trim(index_z),act_as_state_variable=.true.,domain=domain_surface,output=output_none)
        end do
    end do
    
    ! set up cohort-specific variables
    allocate(self%id_N(self%cht_nC+self%cht_nc_init))
    allocate(self%id_r_mass(self%cht_nC+self%cht_nc_init))
    allocate(self%id_i_mass(self%cht_nC+self%cht_nc_init))
    
    allocate(self%id_N_mass(self%cht_nC+self%cht_nc_init))
    allocate(self%id_P_mass(self%cht_nC+self%cht_nc_init))
    
    allocate(self%id_g_mass_in(self%cht_nC+self%cht_nc_init))
    allocate(self%id_g_mass_out(self%cht_nC+self%cht_nc_init))
    allocate(self%id_dist_in(self%cht_nC+self%cht_nc_init,self%cht_nlev))
    allocate(self%id_dist_out(self%cht_nC+self%cht_nc_init,self%cht_nlev))
    
    do i_nc=1,self%cht_nC+self%cht_nc_init
       write (index_nc,'(i0)'), i_nc
       call self%register_state_variable(self%id_N(i_nc),'N'//trim(index_nc), '[m^-2]','abundance of cohort number '//trim(index_nc),initial_value= 0._rk, minimum=0.0_rk)
       call self%register_state_variable(self%id_r_mass(i_nc),'r_mass'//trim(index_nc),'[gDW]','reversible mass of cohort number '//trim(index_nc),initial_value= self%cht_q_J*self%cht_w_b/(1+self%cht_q_J), minimum=0.0_rk)
       call self%register_state_variable(self%id_i_mass(i_nc),'i_mass'//trim(index_nc),'[gDW]','irrevarsible mass of cohort number '//trim(index_nc),initial_value=self%cht_w_b/(1+self%cht_q_J), minimum=0.0_rk)
        
       call self%register_state_variable(self%id_N_mass(i_nc),'N_mass'//trim(index_nc),'[gN]','nitrogen mass of cohort number '//trim(index_nc),initial_value= self%cht_w_b*self%cht_ND_Ref, minimum=0.0_rk)
       call self%register_state_variable(self%id_P_mass(i_nc),'P_mass'//trim(index_nc),'[gP]','phospherous mass of cohort number '//trim(index_nc),initial_value= self%cht_w_b*self%cht_PD_Ref, minimum=0.0_rk)
       
       call self%register_diagnostic_variable(self%id_g_mass_out(i_nc),'reproductive_effort'//trim(index_nc),'gDW','reproductive effort of cohort number '//trim(index_nc),act_as_state_variable=.true.,domain=domain_surface,output=output_none)
       call self%register_horizontal_dependency(self%id_g_mass_in(i_nc),'reproductive_effort'//trim(index_nc),'gDW','reproductive effort of cohort number '//trim(index_nc))
       do i_z=1,self%cht_nlev
           write(index_z,'(i0)') i_z
           call self%register_diagnostic_variable(self%id_dist_out(i_nc,i_z),'dist_cohort_'//trim(index_nc)//'_at_z'//trim(index_z),'-','proportion time spent at depth '//trim(index_z)//' by cohort '//trim(index_nc),act_as_state_variable=.true.,domain=domain_surface)
           call self%register_horizontal_dependency(self%id_dist_in(i_nc,i_z),'dist_cohort_'//trim(index_nc)//'_at_z'//trim(index_z),'-','proportion time spent at depth '//trim(index_z)//' by cohort '//trim(index_nc))
       end do
    end do
    
    !  register environmental dependencies
    call self%register_dependency(self%id_Day,standard_variables%number_of_days_since_start_of_the_year)
    
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
    ! state variables & diagnostics --------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: N,           r_mass,     i_mass,     N_mass,     P_mass
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: g_mass
    real(rk), dimension(self%cht_nlev)                  :: ZD,          ZN,         ZP
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: dist
    ! derivatives --------------------------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: d_N,         d_rmass,    d_imass,    d_Nmass,    d_Pmass
    real(rk), dimension(self%cht_nlev)                  :: d_ZD,        d_ZN,       d_ZP
    ! local cohort vector variables --------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: lengths,     st_mass,    H,          A_z,        mu_size
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: Q10_m,       Q10_a      
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: mu_c,        E_a,        Ing,        E_g,        mu_s
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: mu,          E_m                  
    ! local cohort matrix variables --------------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nC+self%cht_nc_init)                   :: A_c,         L_c,        L_v
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: Q10_m_M,     Q10_a_M,    uTm_M,      rt_m,       rt_a
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: H_M,         A_z_M,      E_m_P,      I_P,        p_i
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: I_r,         eta_z,      eta,        mu_c_M,     E_m_M
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nC+self%cht_nc_init, self%cht_nlev)    :: eta_c,       A_c_M
    ! local cohort single within-loop variables --------------------------------------------------------------------
    real(rk)                                            :: F_tot
    ! local stoichiometric vector variables-------------------------------------------------------------------------
    real(rk), dimension(self%cht_nlev)                  :: ke_P_zoo,    ke_N_zoo,   NH4_tot,    PO4_tot,    DW_tot
    real(rk), dimension(self%cht_nlev)                  :: DW_POM,      DW_DOM,     OP_tot,     OP_POM,     OP_DOM
    real(rk), dimension(self%cht_nlev)                  :: ON_tot,      ON_POM,     ON_DOM,     PDZoo,      NDZoo
    real(rk), dimension(self%cht_nlev)                  :: BOT_FEED
    real(rk), dimension(self%cht_nC+self%cht_nc_init)   :: ke_P_pis,    ke_N_pis,   PD,         ND,         E_m_corr
    ! local stoichiometric matrix variables ------------------------------------------------------------------------
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: I_r_zoo,     I_r_N_zoo,  N_ass_zoo,  N_ass,      I_r_P_zoo
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: P_ass_zoo,   P_ass,      P_exc,      N_exc,      D_ege
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: P_ege,       N_ege,      NH4_ege,    NH4_mor,    PO4_ege
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: PO4_mor,     DW_mor,     OP_mor,     OP_ege,     ON_mor
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nlev)                                  :: ON_ege
    real(rk), dimension(self%cht_nC+self%cht_nc_init, &
        self%cht_nC+self%cht_nc_init, self%cht_nlev)    :: I_r_pis,     I_r_N_pis,  N_ass_pis,  I_r_P_pis,  P_ass_pis
    ! carriers for environment dependencies and logical switches ---------------------------------------------------
    real(rk), dimension(self%cht_nlev)                  :: dz,          uTm
    real(rk)                                            :: Day
    integer                                             :: i,           j,          ix_repro,   nC_fin
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
        _GET_HORIZONTAL_(self%id_N,N)                                                                               ! get cohort abundance
        _GET_HORIZONTAL_(self%id_r_mass,r_mass)                                                                     ! get cohort reversible mass (reserves+gonads)
        _GET_HORIZONTAL_(self%id_i_mass,i_mass)                                                                     ! get cohort irreversible mass (structural tissues, organs)
        _GET_HORIZONTAL_(self%id_N_mass,N_mass)                                                                     ! get cohort total Nitrogen mass (structure+reserves+gonads)
        _GET_HORIZONTAL_(self%id_P_mass,P_mass)                                                                     ! get cohort total Phosphorous mass (structure+reserves+gonads)
        _GET_HORIZONTAL_(self%id_g_mass_in,g_mass)                                                                  ! get cohort reproductive effort (gonads)
        !===========================================================================================================
        ! 1.2 get food concentrations & calculate food nutrient ratios
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _GET_HORIZONTAL_(self%id_ZD,ZD)                                                                             ! get zooplankton (food) DW concentrations
        _GET_HORIZONTAL_(self%id_ZN,ZN)                                                                             ! get zooplankton (food) N concentrations
        _GET_HORIZONTAL_(self%id_ZP,ZP)                                                                             ! get zooplankton (food) P concentrations
        
        where (ZD>0)
            PDZoo = ZP/ZD                                                                                           ! depth-specific zooplankton PW:DW ratio
            NDZoo = ZN/ZD                                                                                           ! depth-specific zooplankton NW:DW ratio
        elsewhere
            PDZoo = 0.0_rk
            NDZoo = 0.0_rk
        endwhere   
        !===========================================================================================================
        ! 1.3 retrieve environmental dependencies and update local trackers & switches
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _GET_GLOBAL_(self%id_Day,Day)                                                                               ! get current julian day
        _GET_HORIZONTAL_(self%id_dz,dz)                                                                             ! get layer depths
        _GET_HORIZONTAL_(self%id_uTm,uTm)                                                                           ! get layer temperatures
        do i=1,nC_fin
            _GET_HORIZONTAL_(self%id_dist_in(i,:),dist(i,:))                                                        ! get current fidh distributions (dimensions: nC_fin x nlev)
        end do
        ! check & update year - used for selecting new cohorts to start. To be removed if cohort implementation
        !                       changes or if functionality to get simulation year identification is added to FABM
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
            Q10_m = self%cht_theta_m_s*(r_mass+i_mass)**self%cht_theta_m_e                                          ! calculate size-specific Q10 for fish metabolism
            Q10_a = self%cht_theta_a_s*(r_mass+i_mass)**self%cht_theta_a_e                                          ! calculate size-specific Q10 for fish feeding
        elsewhere
            Q10_m = 1.0_rk
            Q10_a = 1.0_rk
        endwhere
        Q10_m_M = spread(Q10_m, DIM=2, NCOPIES=self%cht_nlev)                                                       ! copy Q10s to nC_fin x nlev matrix
        Q10_a_M = spread(Q10_a, DIM=2, NCOPIES=self%cht_nlev)                                                       ! copy Q10s to nC_fin x nlev matrix
        uTm_M = spread(uTm, DIM=1, NCOPIES=nC_fin)                                                                  ! copy temperatures to nC_fin x nlev matrix
        rt_m = Q10_m_M**((uTm_M-self%cht_T_ref)/10.0_rk)                                                            ! calculate temperature correction factor for metabolism for each cohort and layer
        rt_a = Q10_a_M**((uTm_M-self%cht_T_ref)/10.0_rk)                                                            ! calculate temperature correction factor for each cohort and layer
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
        ! 2.4 fish reproduction
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
        ! 2.5 calculate depth- and temperature-specific rates
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        A_z_M = spread(A_z, DIM=2, NCOPIES=self%cht_nlev)*rt_a                                                      ! depth-specific temperature corrected attack rates (dimensions: nC_fin x nlev)
        H_M = spread(H, DIM=2, NCOPIES=self%cht_nlev)*(1/rt_a)                                                      ! depth-specific temperature corrected handling times (dimensions: nC_fin x nlev)
        where (isnan(H)) H = 0.0_rk                                                                                 ! get rid of NANs
        A_c_M = spread(A_c, DIM=3, NCOPIES=self%cht_nlev)*spread(rt_a, DIM=1, NCOPIES=nC_fin)                       ! depth-specific temperature corrected piscivorous attack rate matrix (dimensions: nC_fin x nC_fin x nlev)
        E_m_P = spread(E_m, DIM=2, NCOPIES=self%cht_nlev)*rt_m                                                      ! depth-specific temperature corrected basal metabolism (dimensions: nC_fin x nlev)
        !===========================================================================================================
        ! 2.6  determine spatial distributions & calculate depth-specific rates
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        eta_z = A_z_M*spread(ZD, DIM=1, NCOPIES=nC_fin)                                                             ! depth-specific encounter rates with zooplankton prey (dimensions: nC_fin x nlev)
        eta_c = A_c_M*spread(dist*spread((r_mass+i_mass)*N, DIM=2, NCOPIES=self%cht_nlev)/ &
            spread(dz, DIM=1, NCOPIES=nC_fin), DIM=2, NCOPIES=nC_fin)   ! depth-specific encounter rates with fish prey (dimensions: nC_fin x nC_fin x nlev)
        eta = eta_z+sum(eta_c,1)                                                                                    ! depth-specific total encounter rates with prey (dimensions: nC_fin x nlev)
        
        I_P = eta/(1.0_rk+H_M*eta)                                                                                  ! layer-specific potential ingestion rates from feeding on zooplankton, corrected with layer depth (dimensions: nC_fin x nlev)
        where (spread(sum(I_P,2), DIM=2, NCOPIES=self%cht_nlev)>0) ! avoid creating NAN's where cohort are empty or there is no food in the watercolumn
            p_i = (I_P*spread(dz, DIM=1, NCOPIES=nC_fin))/spread(sum(I_P*spread(dz, DIM=1, NCOPIES=nC_fin),2), &
                DIM=2, NCOPIES=self%cht_nlev)                                                                       ! calculate distribution weights (dimensions: nC_fin x nlev)
        elsewhere
            p_i = spread(dz/sum(dz), DIM=1, NCOPIES=nC_fin)                                                         ! in case cohort is not empty, but ingestion is zero everywhere, fish distribute evenly
        endwhere
               
        I_r= p_i*I_P                                                                                                ! depth-specific realized ingestion rates (dimensions: nC_fin x nlev)
        Ing = sum(I_r, 2)                                                                                           ! total realized ingestion rates (dimensions: nC_fin)
        E_m_M = E_m_P*p_i                                                                                           ! depth-specific realized metabolism (dimensions: nC_fi)
        
        mu_c_M = sum(A_c_M*spread(dist*(spread(N, DIM=2, NCOPIES=self%cht_nlev)/spread(dz, DIM=1, NCOPIES=nC_fin)) &
            /(1.0_rk+H_M*eta), DIM=1, NCOPIES=nC_fin), 2)*p_i                                                       ! depth-specific mortality from piscivory (dimensions: nC_fin x nlev)
        mu_c = sum(mu_c_M, 2)                                                                                       ! total piscivorous mortality (dimensions: nC_fin)
        !===========================================================================================================
        ! 2.7 stochiometric correction of respiration
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        where (N>0.0_rk) ! non-empty cohorts
            PD = P_mass / (r_mass*i_mass)                                                                           ! calculate current cohort P/D ratios
            ND = N_mass / (r_mass*i_mass)                                                                           ! calculate current cohort N/D ratios
        elsewhere ! empty cohorts
            PD = 1.0_rk                                                                                             ! prevents INF and NAN results later
            ND = 1.0_rk                                                                                             ! prevents INF and NAN results later
        endwhere
        E_m_corr = sum(E_m_M,2)*max(self%cht_PD_Ref/PD,self%cht_ND_Ref/ND)                                          ! correct respiration to maintain stochiometric equilibrium (dimensions: nC_fin)
        !===========================================================================================================
        ! 2.8 calculate final vital rates
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
        ! 2.9 set ODEs
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
        
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_g_mass_out,g_mass)                                                      ! update reproductive effort diagnostic
        !===========================================================================================================
        ! 2.10 start new cohort
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        _SET_SURFACE_ODE_(self%id_N(year+self%cht_nc_init),F_tot)                                                   ! start new cohort with the produced offspring
        
        p_i(year+self%cht_nc_init,:)=dz/sum(dz)                                                                     ! new cohorts are initially evenly distibuted
        do i=1,nC_fin
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dist_out(i,:),p_i(i,:))                                             ! update fish distributions
        end do
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. calculate and update zooplankton losses
    !===============================================================================================================
        d_ZD = - sum(spread(N,DIM=2, NCOPIES=self%cht_nlev)*p_i*eta_z/(1.0_rk+H_M*eta), 1)/dz                       ! set zooplankton DW derivative
        d_ZP = d_ZD*PDZoo                                                                                           ! set zooplankton PW derivative
        d_ZN = d_ZD*NDZoo                                                                                           ! set zooplankton NW derivative
        
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_prey_loss_DW,d_ZD)                                                      ! update zooplankton DW loss rate
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_prey_loss_P,d_ZP)                                                       ! update zooplankton P loss rate
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_prey_loss_N,d_ZN)                                                       ! update zooplankton N loss rate
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. do stochiometry & and update external abiotic variables
    !===============================================================================================================
        I_r_zoo = (eta_z/(1.0_rk+H_M*eta))*p_i                                                                      ! calculate depth-spcific total zooplankton ingestion rates (dimensions: nC_fin x nlev)
        I_r_pis = (eta_c/(1.0_rk+spread(H_M*eta,DIM=1, NCOPIES=nC_fin)))*spread(p_i, DIM=1, NCOPIES=nC_fin)         ! calculate depth- and cohort-spcefic piscivorous ingestion rates (dimensions: nC_fin x nC_fin x nlev)
        !===========================================================================================================
        ! 4.1 assimilation
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.1.1 phosphorous
            !-------------------------------------------------------------------------------------------------------
            ke_P_zoo = min(1.0_rk,self%cht_PD_Ref/PDZoo * self%cht_k_e)                                             ! calculate P assimilation efficiencies of zooplanktivory (dimensions: nlev)
            I_r_P_zoo = spread(PDZoo, DIM=1, NCOPIES=nC_fin) * I_r_zoo                                              ! calculate P ingestion rates from zooplanktivory (dimensions: nC_fin x nlev)
            P_ass_zoo = spread(ke_P_zoo, DIM=1, NCOPIES=nC_fin) * I_r_P_zoo                                         ! calculate P assimilation rates from zooplanktivory (dimensions: nC_fin x nlev)
            
            ke_P_pis = min(1.0_rk,self%cht_PD_Ref/PD * self%cht_k_e)                                                ! calculate P assimilation efficiencies of piscivory (dimensions: nC_fin)
            I_r_P_pis = spread(spread(PD, DIM=2, NCOPIES=nC_fin), DIM=3, NCOPIES=self%cht_nlev) * I_r_pis           ! calculate P ingestion rates from piscivory (dimensions: nC_fin x nC_fin x nlev)
            P_ass_pis = spread(spread(ke_P_pis, DIM=2, NCOPIES=nC_fin), DIM=3, NCOPIES=self%cht_nlev) * I_r_P_pis   ! calculate P assimilation rates from piscivory (dimensions: nC_fin x nC_fin x nlev)

            P_ass = P_ass_zoo+sum(P_ass_pis,1)                                                                      ! calculate total P assimilation rates (dimensions: nC_fin x nlev)
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.1.2 nitrogen
            !-------------------------------------------------------------------------------------------------------
            ke_N_zoo = min(1.0_rk,self%cht_ND_Ref/NDZoo * self%cht_k_e)                                             ! calculate N assimilation efficiencies of zooplanktivory (dimensions: nlev)
            I_r_N_zoo = spread(NDZoo, DIM=1, NCOPIES=nC_fin) * I_r_zoo                                              ! calculate P ingestion rates from zooplanktivory (dimensions: nC_fin x nlev)
            N_ass_zoo = spread(ke_N_zoo, DIM=1, NCOPIES=nC_fin) * I_r_N_zoo                                         ! calculate N assimilation rates from zooplanktivory (dimensions: nC_fin x nlev)
            
            ke_N_pis = min(1.0_rk,self%cht_ND_Ref/ND * self%cht_k_e)                                                ! calculate N assimilation efficiencies of piscivory (dimensions: nC_fin)
            I_r_N_pis = spread(spread(ND, DIM=2, NCOPIES=nC_fin), DIM=3, NCOPIES=self%cht_nlev) * I_r_pis           ! calculate P ingestion rates from piscivory (dimensions: nC_fin x nC_fin x nlev)
            N_ass_pis = spread(spread(ke_N_pis, DIM=2, NCOPIES=nC_fin), DIM=3, NCOPIES=self%cht_nlev) * I_r_N_pis   ! calculate N assimilation rates from piscivory (dimensions: nC_fin x nC_fin x nlev)

            N_ass = N_ass_zoo+sum(N_ass_pis,1)                                                                      ! calculate total N assimilation rates (dimensions: nC_fin x nlev)
        !===========================================================================================================
        ! 4.2 excretion
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.2.1 phosphorous
            !-------------------------------------------------------------------------------------------------------
            P_exc = spread(PD/self%cht_PD_Ref, DIM=2, NCOPIES=self%cht_nlev) * E_m_M                                ! calculate rates of P excretion (dimensions: nC_fin x nlev)
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.2.2 nitrogen
            !-------------------------------------------------------------------------------------------------------
            N_exc = spread(ND/self%cht_ND_Ref, DIM=2, NCOPIES=self%cht_nlev) * E_m_M                                ! calculate rates of N excretion (dimensions: nC_fin x nlev)
        !===========================================================================================================
        ! 4.3 egestion
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.1 dry weight
            !-------------------------------------------------------------------------------------------------------
            D_ege = I_r*(1-self%cht_k_e)                                                                            ! calculate egestion rate (dimensions: nC_fin x nlev)
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.2 phosphorous
            !-------------------------------------------------------------------------------------------------------
            P_ege = I_r_P_zoo + sum(I_r_P_pis,1) - P_ass                                                            ! calculate egested P (dimensions: nC_fin x nlev)
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.3.3 nitrogen
            !-------------------------------------------------------------------------------------------------------
            N_ege = I_r_N_zoo + sum(I_r_N_pis,1) - N_ass                                                            ! calculate egested P (dimensions: nC_fin x nlev)
        !===========================================================================================================
        ! 4.4 set nutrient ODEs for fish
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        d_Pmass = sum(P_ass-P_exc,2) - g_mass*self%cht_PD_Ref                                                       ! calculate rates of change in P content for individual fish (dimensions: nC_fin)
        d_Nmass = sum(N_ass-N_exc,2) - g_mass*self%cht_ND_Ref                                                       ! calculate rates of change in N content for individual fish (dimensions: nC_fin)
        
        _SET_SURFACE_ODE_(self%id_P_mass,d_Pmass)                                                                   ! update P_mass ODE
        _SET_SURFACE_ODE_(self%id_N_mass,d_Nmass)                                                                   ! update N_mass ODE
        !===========================================================================================================
        ! 4.5 calculate and update contributions to external processes
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.1 NH4
            !-------------------------------------------------------------------------------------------------------
            NH4_ege = self%cht_Diss_Eges * N_ege                                                                    ! calculate depth-specific total NH4 flux from egestion (dimensions: nC_fin x nlev)
            NH4_mor = spread(self%cht_Diss_Mort * (N_mass-self%cht_D_Bone*self%cht_ND_Ref*st_mass)*mu, &
                DIM=2, NCOPIES=self%cht_nlev)                                                                       ! calculate per capita NH4 flux from mortality - bone mass calculated from standard mass, since can't be starved (dimensions: nC_fin x nlev)
            
            NH4_tot = sum((N_exc + NH4_ege + NH4_mor)*spread(N, DIM=2, NCOPIES=self%cht_nlev)*p_i,1)                ! calculate total population NH4 flux to the water column (dimensions: nlev)                                                - final export!
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.2 PO4
            !-------------------------------------------------------------------------------------------------------
            PO4_ege = self%cht_Diss_Eges * P_ege                                                                    ! calculate PO4 flux from egestion (dimensions: nC_fin x nlev)
            PO4_mor = spread(self%cht_Diss_Mort * (P_mass-self%cht_P_Bone*self%cht_PD_Ref*st_mass)*mu, &
                DIM=2, NCOPIES=self%cht_nlev)                                                                       ! calculate per capita PO4 flux from mortality - bone mass calculated from standard mass, since can't be starved (dimensions: nC_fin x nlev)
            
            PO4_tot = sum((P_exc + PO4_ege + PO4_mor)*spread(N, DIM=2, NCOPIES=self%cht_nlev)*p_i,1)                ! calculate total population PO4 flux to the water column (dimensions: nlev)                                                - final export!
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.3 organic matter DW
            !-------------------------------------------------------------------------------------------------------
            DW_mor = spread(((r_mass+i_mass)-self%cht_D_Bone*st_mass)*mu, DIM=2, NCOPIES=self%cht_nlev)             ! calculate per capita DW flux from mortality - bone mass calculated from standard mass, since can't be starved (dimensions: nC_fin x nlev)
            
            DW_tot = sum((DW_mor + D_ege)*spread(N, DIM=2, NCOPIES=self%cht_nlev)*p_i, 1)                           ! calculate total DW flux to the water column (dimensions: nlev)
            DW_POM = DW_tot*(1.0_rk-self%cht_DOMW)                                                                  ! particulate fraction of DW_tot (dimensions: nlev)                                                                         - final export!
            DW_DOM = DW_tot*self%cht_DOMW                                                                           ! dissolved fraction of DW_tot (dimensions: nlev)                                                                           - final export!
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.4 organic matter P
            !-------------------------------------------------------------------------------------------------------
            OP_mor = spread((1.0_rk-self%cht_Diss_Mort)*(P_mass-self%cht_P_Bone*self%cht_PD_Ref*st_mass)*mu, &
                DIM=2, NCOPIES=self%cht_nlev)                                                                       ! calculate per capita organic P flux from mortality - bone mass calculated from standard mass, since can't be starved (dimensions: nC_fin x nlev)
            OP_ege = P_ege*(1.0_rk-self%cht_Diss_Eges)                                                              ! calculate organic P flux from egestion (dimensions: nC_fin x nlev)
            
            OP_tot = sum((OP_mor + OP_ege)*spread(N, DIM=2, NCOPIES=self%cht_nlev)*p_i, 1)                          ! calculate total organic P flux to the water column (dimensions: nlev)
            OP_POM = OP_tot*(1.0_rk-self%cht_DOMW)                                                                  ! particulate fraction of OP_tot (dimensions: nlev)                                                                         - final export!
            OP_DOM = OP_tot*self%cht_DOMW                                                                           ! dissolved fraction of OP_tot (dimensions: nlev)                                                                           - final export!
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.4 organic matter N
            !-------------------------------------------------------------------------------------------------------
            ON_mor = spread((1.0_rk-self%cht_Diss_Mort)*(N_mass-self%cht_D_Bone*self%cht_ND_Ref*st_mass)*mu, &      
                DIM=2, NCOPIES=self%cht_nlev)                                                                       ! calculate per capita organic N flux from mortality - bone mass calculated from standard mass, since can't be starved (dimensions: nC_fin x nlev)
            ON_ege = N_ege*(1.0_rk-self%cht_Diss_Eges)                                                              ! calculate organic N flux from egestion (dimensions: nC_fin x nlev)
            
            ON_tot = sum((ON_mor + ON_ege)*spread(N, DIM=2, NCOPIES=self%cht_nlev)*p_i, 1)                          ! calculate total organic N flux to the water column (dimensions: nlev)
            ON_POM = ON_tot*(1.0_rk-self%cht_DOMW)                                                                  ! particulate fraction of ON_tot (dimensions: nlev)                                                                         - final export!
            ON_DOM = ON_tot*self%cht_DOMW                                                                           ! dissolved fraction of ON_tot (dimensions: nlev)                                                                           - final export!
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            ! 4.5.5 oupdate output
            !-------------------------------------------------------------------------------------------------------
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,1), NH4_tot)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,2), PO4_tot)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,3), DW_POM)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,4), DW_DOM)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,5), OP_POM)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,6), OP_DOM)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,7), ON_POM)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,8), ON_DOM)
            
            BOT_FEED=0.0_rk                                                                                         ! replace with biomass of bottomfeeding fish (dist*DW*%bottomfeeding)
            
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pas_rates(:,9), BOT_FEED)
    _HORIZONTAL_LOOP_END_

    end subroutine do_surface

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
end module fish_cohort
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------