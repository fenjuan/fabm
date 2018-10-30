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
      type (type_state_variable_id), allocatable   :: id_N(:)
      type (type_state_variable_id), allocatable   :: id_r_mass(:)
      type (type_state_variable_id), allocatable   :: id_i_mass(:)
!
!     Model parameters
      integer         :: Nclasses ! number of classes
      real(rk), allocatable       :: wb(:) ! mass at birth
!
      contains
!
! Model procedures
      procedure :: initialize
      procedure :: do

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

   call self%get_parameter(self%Nclasses,'Nclasses','','Number of classes',default=1)

   allocate(self%wb(self%Nclasses))

   allocate(self%id_N(self%Nclasses))
   allocate(self%id_r_mass(self%Nclasses))
   allocate(self%id_i_mass(self%Nclasses))

   ! get parameters for specific class
   do n=1,self%Nclasses
      write (name,"(I02.2,A3)") n-1,"/wb"
      call self%get_parameter(self%wb(n),name,'d-1','growth rate of prey',default=1.00_rk,scale_factor=1._rk)
   end do

write(*,*) self%wb
!stop

   do n=1,self%Nclasses
      write (name,"(A1,I02.2)") "N", n-1
      write (longname, "(A9, I02.2)") "abundance", n-1
      call self%register_state_variable(self%id_N(n),name, '-', longname,     &
                                       0._rk,minimum=0.0_rk,no_river_dilution=.TRUE.)
      write (name,"(A6,I02.2)") "r_mass", n-1
      write (longname, "(A15, I02.2)") "reversible mass", n-1
      call self%register_state_variable(self%id_r_mass(n),name,'kg/individual',longname,     &
                                       0._rk,minimum=0.0_rk,no_river_dilution=.TRUE.)
      write (name,"(A6,I02.2)") "i_mass", n-1
      write (longname, "(A17, I02.2)") "irreversible mass", n-1
      call self%register_state_variable(self%id_i_mass(n),name,'kg/individual',longname,     &
                                       0._rk,minimum=0.0_rk,no_river_dilution=.TRUE.)
   end do

#if 0
   call self%get_parameter(self%b,'b','d-1','growth rate of prey',   default=1.00_rk,scale_factor=d_per_s)
   call self%get_parameter(self%p,'p','d-1','impact of predation',   default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r,'r','d-1','growth efficiency rate',default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%d,'d','d-1','death rate',            default=0.50_rk,scale_factor=d_per_s)

!  Register state variables
   call self%register_state_variable(self%id_prey,'prey','mmol/m**3','nutrient',     &
                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_predator,'predator','mmol/m**3','phytoplankton',     &
                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)

!  Register conserved quantities
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_prey)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_predator)
#endif

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
   class (type_fish_cohort),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
!  state variables
   real(rk)   :: N,r_mass,i_mass
   real(rk)   :: d_N, d_rmass, d_imass
   integer    :: i
!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _LOOP_BEGIN_

   do i=1,self%Nclasses
      _GET_(self%id_N(i),N)
      _GET_(self%id_r_mass(i),r_mass)
      _GET_(self%id_i_mass(i),i_mass)
      ! do something more useful here
!      write(*,*) N,r_mass,i_mass
!      _SET_ODE_(self%id_N,g)
!      _SET_ODE_(self%id_r_mass,f)
!      _SET_ODE_(self%id_i_mass,f)
   end do

#if 0
!  Retrieve current state variables values
!  prey density - pelagic
   _GET_(self%id_prey,prey)
!  predator density - pelagic
   _GET_(self%id_predator,predator)

!  Use Lotka-Volterra model
!  Calculate change of predator
   g = (self%r*prey-self%d)*predator
!  Calculate change of prey
   f = (self%b-self%p*predator)*prey
!
!  Set temporal derivatives,  subtracted by the secs_pr_day
   _SET_ODE_(self%id_predator,g)
   _SET_ODE_(self%id_prey,f)
!  Export diagnostic variables
!  Leave spatial loops (if any)
#endif
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module fish_cohort

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
