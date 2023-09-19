module Reaction_Sandbox_DWP_class
! Distributed Waste Package 
! equations (considering same units):
!
! flux_to_mobile = diffusion_rate * ( conc_waste_package - conc_mobile ) 
! waste_decay = decay_rate * conc_waste_package
! conc_mobile += flux_to_mobile
! conc_waste_package -= flux_to_mobile + decay_halflife
!
! - different scaling has to be applied to the mobaile and immobile contribution in order to have mass conservation (without decay)
! - both `diffusion_rate` and `decay rate` could be space dependent using auxiliary minerals

#include "petsc/finclude/petscsys.h"
  use petscsys

! 1. Change all references to "DWP" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.,
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_DWP_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: mobile_name, immobile_name, diffusion_scaling_mineral_name, decay_scaling_mineral_name, &
    waste_mineral_name
    PetscInt :: species_id, species_im_id, diffusion_scaling_mineral_id, decay_scaling_mineral_id, &
    waste_mineral_id
    PetscReal :: diffusion_rate, decay_rate, waste_rate
    PetscInt :: auxiliary_offset
    
  contains
    procedure, public :: ReadInput => DWPRead
    procedure, public :: Setup => DWPSetup
    procedure, public :: Evaluate => DWPEvaluate
    procedure, public :: AuxiliaryPlotVariables => DWPAuxiliaryPlotVariables
    procedure, public :: UpdateKineticState => DWPUpdateKineticState
    procedure, public :: Destroy => DWPDestroy
  end type reaction_sandbox_DWP_type

  public :: DWPCreate

contains

! ************************************************************************** !

function DWPCreate()
  !
  ! Allocates DWP reaction object.
  !

  implicit none

  class(reaction_sandbox_DWP_type), pointer :: DWPCreate

! 4. Add code to allocate the object, initialize all variables to zero and
!    nullify all pointers. E.g.,
  allocate(DWPCreate)
  DWPCreate%auxiliary_offset = 0
  DWPCreate%mobile_name = ''
  DWPCreate%immobile_name = ''
  DWPCreate%waste_mineral_name = ''
  DWPCreate%species_id = 0
  !DWPCreate%rate_immobile = 0.d0
  DWPCreate%diffusion_rate = 0.d0
  DWPCreate%waste_rate = 0.d0
  nullify(DWPCreate%next)

end function DWPCreate

! ************************************************************************** !

subroutine DWPRead(this,input,option)
  !
  ! Reads input deck for DWP reaction parameters (if any)
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_DWP_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'CHEMISTRY,REACTION_SANDBOX,DISTRIBUTED_WASTE_PACKAGE'
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', trim(error_string))
    call StringToUpper(word)

    select case(trim(word))

      ! DWP Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !     # begin user-defined input
      !     DWP
      !       DWP_INTEGER 1
      !       DWP_INTEGER_ARRAY 2 3 4
      !     END
      !     # end user defined input
      !   END
      !   ...
      ! END

! 5. Add a case statement for reading variables.
      case('MOBILE_NAME')
        call InputReadWord(input,option,this%mobile_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'mobile_name', trim(error_string))
      case('IMMOBILE_NAME')
        call InputReadWord(input,option,this%immobile_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'IMMOBILE_NAME', trim(error_string))
! 8. Repeat for other variables.
!       case('RATE_IMMOBILE')
!         ! Read the double precision rate constant.
!         call InputReadDouble(input,option,this%rate_immobile)
!         ! Note the use of character variable 'word' instead of 'RATE_CONSTANT'
!         ! in the error message, as they are identical.
!         call InputErrorMsg(input,option,word, &
!                            'CHEMISTRY,REACTION_SANDBOX,DWP')
!         ! Read the optional units and convert to internal
!         ! units of 1/s.
!         internal_units = 'unitless/sec'
!         call InputReadAndConvertUnits(input,this%rate_immobile, &
!                                 internal_units,'CHEMISTRY,REACTION_SANDBOX,&
!                                 &DWP,RATE_IMMOBILE',option)

      case('DIFFUSION_SCALE_MINERAL')
        call InputReadWord(input,option,this%diffusion_scaling_mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, trim(error_string))
      case('DECAY_SCALE_MINERAL')
        call InputReadWord(input,option,this%decay_scaling_mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, trim(error_string))
        
        
      case('WASTE_MINERAL')
        call InputReadWord(input,option,this%waste_mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, trim(error_string))
        
      case('WASTE_RATE')
        call InputReadDouble(input,option,this%waste_rate)
        call InputErrorMsg(input,option,word, trim(error_string))
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%waste_rate, &
                                 internal_units, trim(error_string), option)
        
      case('DIFFUSION_RATE')
        call InputReadDouble(input,option,this%diffusion_rate)
        call InputErrorMsg(input,option,word, trim(error_string))
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%diffusion_rate, &
                                 internal_units, trim(error_string), option)
      case('DECAY_RATE')
        call InputReadDouble(input,option,this%decay_rate)
        call InputErrorMsg(input,option,word, trim(error_string))
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%decay_rate, &
                                 internal_units, trim(error_string), option)

    case default
        call InputKeywordUnrecognized(input,word, trim(error_string),option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine DWPRead

! ************************************************************************** !

subroutine DWPSetup(this,reaction,option)
  !
  ! Sets up the DWP reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only: GetKineticMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_DWP_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 1
  
! 9. Add code to initialize.
  this%species_id = &
    GetPrimarySpeciesIDFromName(this%mobile_name, reaction, option)
  this%species_im_id = &
    GetImmobileSpeciesIDFromName(this%immobile_name, reaction%immobile, option)  
  this%diffusion_scaling_mineral_id = &
    GetKineticMineralIDFromName(this%diffusion_scaling_mineral_name, reaction%mineral, option)
  this%decay_scaling_mineral_id = &
    GetKineticMineralIDFromName(this%decay_scaling_mineral_name, reaction%mineral, option)
  this%waste_mineral_id = &
    GetKineticMineralIDFromName(this%waste_mineral_name, reaction%mineral, option)
    

end subroutine DWPSetup

! ************************************************************************** !

subroutine DWPEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Reaction_Mineral_Aux_module

  implicit none

  class(reaction_sandbox_DWP_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscInt :: idof_immobile, idof_mobile, idof_waste
  PetscReal :: L_water, flux_to_aq, rate, decay_rate, decay_flux, waste_rate, waste_flux
  PetscInt :: iauxiliary

  iauxiliary = this%auxiliary_offset + 1

  ! Description of subroutine arguments:

  ! Residual - 1D array storing Residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entries in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analytical derivatives will
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN defined within the NUMERICAL_METHODS TRANSPORT,
  !   NEWTON_SOLVER block of the input deck.  If the use of
  !   NUMERICAL_JACOBIAN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.,
  !
  !   if (compute_derivative) then
  !     option%io_buffer = 'NUMERICAL_JACOBIAN must be specified within &
  !       &the NEWTON_SOLVER block of NUMERICAL_METHODS TRANSPORT due to &
  !       &assumptions made in DWPEvaluate.'
  !     call PrintErrMsg(option)
  !   endif
  !
  ! rt_auxvar - Object holding chemistry information (e.g., concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g., saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - liquid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - liquid density [kg/m^3]
  !     global_auxvar%sat(iphase) - liquid saturation [m^3 water/m^3 pore]
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.,
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_mobile_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.

! 10. Add code for the Residual evaluation.

  ! Units of the Residual must be in moles/second.
  ! 1.d3 converts m^3 water -> L water
  !L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
  !       material_auxvar%volume*1.d3
  
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)*1.d3
  rate = this%diffusion_rate * rt_auxvar%mnrl_volfrac(this%diffusion_scaling_mineral_id)
  decay_rate = this%decay_rate * rt_auxvar%mnrl_volfrac(this%decay_scaling_mineral_id)
  waste_rate = this%waste_rate * rt_auxvar%mnrl_volfrac(this%waste_mineral_id)  !waste_rate = this%waste_rate * rt_auxvar%mnrl_area(this%waste_mineral_id)
  rt_auxvar%auxiliary_data(iauxiliary) = waste_rate
  
  idof_mobile = this%species_id
  idof_immobile = this%species_im_id + reaction%offset_immobile
  
  flux_to_aq = rate * ( rt_auxvar%immobile(this%species_im_id) - rt_auxvar%total(this%species_id,iphase) ) 
  decay_flux = decay_rate * rt_auxvar%immobile(this%species_im_id)
  
  
  !* material_auxvar%volume 
  ! Always "subtract" the contribution from the Residual.
  Residual(idof_mobile) = Residual(idof_mobile) - flux_to_aq * L_water * material_auxvar%volume ! * rt_auxvar%total(this%species_id,iphase)
  Residual(idof_immobile) = Residual(idof_immobile) + (flux_to_aq + decay_flux - waste_rate) * material_auxvar%volume

  !this%rate_D_m* &                   ! 1/s
  !                         rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
  !                         material_auxvar%volume                                 ! m3 bulk
 
  !
  
  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for the Jacobian evaluation.

    ! Always "add" the contribution to the Jacobian.
    ! Units = (mol/sec)*(kg water/mol) = kg water/sec
    
    Jacobian(idof_mobile, idof_mobile) =   Jacobian(idof_mobile, idof_mobile) + rate * L_water * material_auxvar%volume
    Jacobian(idof_mobile, idof_immobile) =   Jacobian(idof_mobile, idof_immobile) - rate * L_water * material_auxvar%volume
    Jacobian(idof_immobile, idof_mobile) =   Jacobian(idof_immobile, idof_mobile) - rate * material_auxvar%volume
    Jacobian(idof_immobile, idof_immobile) =   Jacobian(idof_immobile, idof_immobile) + (rate + decay_rate - waste_rate) * & 
      material_auxvar%volume
  
  endif

  
end subroutine DWPEvaluate

! ************************************************************************** !

subroutine DWPAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds calcite auxiliary plot variables to output list
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY

  class(reaction_sandbox_DWP_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units

  word = 'Dissolutional Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)

end subroutine DWPAuxiliaryPlotVariables

! *************************************************************************** !

subroutine DWPUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                     material_auxvar,reaction,option)
  !
  ! Updates mineral volume fraction at end converged timestep based on latest
  ! rate
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_DWP_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: waste_mnrl
  PetscReal :: delta_volfrac

  waste_mnrl = this%waste_mineral_id
  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
  !                                mol_vol (m^3 mnrl/mol mnrl)
  !delta_volfrac = this%waste_rate * rt_auxvar%mnrl_volfrac(waste_mnrl) * &
  !       option%tran_dt
  
  
  delta_volfrac = rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(waste_mnrl)* &
                  option%tran_dt
  ! m^3 mnrl/m^3 bulk
  rt_auxvar%mnrl_volfrac(waste_mnrl) = rt_auxvar%mnrl_volfrac(waste_mnrl) - &
                                  delta_volfrac
  ! zero to avoid negative volume fractions
  if (rt_auxvar%mnrl_volfrac(waste_mnrl) < 0.d0) &
    rt_auxvar%mnrl_volfrac(waste_mnrl) = 0.d0

end subroutine DWPUpdateKineticState

! ************************************************************************** !

subroutine DWPDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  implicit none

  class(reaction_sandbox_DWP_type) :: this

! 12. Add code to deallocate dynamic members of reaction_sandbox_DWP_type.

end subroutine DWPDestroy

end module Reaction_Sandbox_DWP_class
