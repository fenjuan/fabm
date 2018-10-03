#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module fish_weight_length
!
! !DESCRIPTION:
!  This is fish module, to calculating the body weight and lenght of individual fish/fish group
!  the bioenergetic model is adopted from Bernard A. Megrey, Kenneth A. Rose, Robert A. Klumb, Douglas E. Hay, Francisco E. Werner, David L. Eslinger, S. Lan Smith,
!  A bioenergetics-based population dynamics model of Pacific herring (Clupea harengus pallasi) coupled to a lower trophic level nutrient–phytoplankton–zooplankton model:
!  Description, calibration, and sensitivity analysis,Ecological Modelling,Volume 202, Issues 1–2,2007,Pages 144-164
   
   
!  fen: Sep. 24th, 2018, step 1, simple fish weight model
   
   
   
!  USES:
   use fabm_types

   implicit none

   private

!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_fish_weight_length
!     Variable identifiers
!     id_w:  fish weight
!     id_l:  fish length
      type (type_state_variable_id)   :: id_w
!
!     Model parameters
      real(rk)        :: con_r, res_r,dig_r, ege_r, exc_r, egg_r
      real(rk)        :: cal_f, cal_z
!
      contains
!
! Model procedures
      procedure :: initialize
      procedure :: do

      end type type_fish_weight_length

   real(rk), parameter :: secs_pr_day = 86400_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise 
!
! !INTERFACE:
!
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_fish_weight_length), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Register model parameters
!   phisological rates

   call self%get_parameter(self%con_r,'consumption','g/g/d',      'consumption rate',default=0.1_rk,scale_factor=1.0/ secs_pr_day)
   call self%get_parameter(self%res_r,'respiration','g/g/d',      'repiration rate',default=0.01_rk,scale_factor=1.0/ secs_pr_day)
   call self%get_parameter(self%dig_r,'digestion','g/g/d',      'digestion rate',default=0.01_rk,scale_factor=1.0/ secs_pr_day)
   call self%get_parameter(self%ege_r,'egestion','g/g/d',      'egestion rate',default=0.01_rk,scale_factor=1.0/ secs_pr_day)
   call self%get_parameter(self%exc_r,'excretion','g/g/d',      ' excretion rate',default=0.01_rk,scale_factor=1.0/ secs_pr_day)
   call self%get_parameter(self%egg_r,'reproduction','g/g/d',      'reproduction rate',default=0.01_rk,scale_factor=1.0/ secs_pr_day)
!   enegy density (fixed for each group/individual, see equation 2.2.2.1)
   call self%get_parameter(self%cal_f,'fish_energy_density','J/g fish',      'fish energy density',default=5750.0_rk)
   call self%get_parameter(self%cal_z,'zoo_energy_density','J/g zooplankton', 'zooplankton energy density',default=2580.0_rk)

!  Register state variables
   call self%register_state_variable(self%id_w,'weight','g','fish body weight',1._rk,minimum=0.0_rk,no_river_dilution=.TRUE.)

!  Register conserved quantities


   return
   end subroutine initialize

!-----------------------------------------------------------------------
! !IROUTINE:the type bound precedure: do(),right hand sides of prey and predator model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_fish_weight_length),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
!  state variables
   real(rk)   :: weight

   real(rk)   :: d_weight
!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _LOOP_BEGIN_

!  Retrieve current state variables values

   _GET_(self%id_w,weight)

   d_weight =((self%con_r -self%res_r - self%dig_r -self%ege_r -self%exc_r) * self%cal_z/self%cal_f - self%egg_r) * weight

!  Set temporal derivatives,  subtracted by the secs_pr_day
   _SET_ODE_(self%id_w,d_weight)

!  Export diagnostic variables
!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module fish_weight_length

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
