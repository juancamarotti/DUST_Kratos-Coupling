!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2023 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!!
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!!
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!!
!! Authors:
!!          Federico Gentile
!!          Matteo Dall'Ora
!!          Alessandro Cocco
!!=========================================================================

module mod_test

use mod_param, only: &
    wp, max_char_len

use mod_hdf5_io, only: &
    h5loc, write_hdf5_attr, open_hdf5_file, close_hdf5_file, read_hdf5

use mod_handling, only: &
    error, warning, info, printout, dust_time, t_realtime, internal_error

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

use mod_basic_io, only: &
  read_real_array_from_file  

implicit none

public :: t_sim_param, sim_param, init_sim_param, create_param_test_particle

private    

!> sim_param
type t_sim_param
  !Time:
  !> Start time
  real(wp) :: t0
  !> Time step
  real(wp) :: dt
  !> Final time
  real(wp) :: tend
  !> Number of timesteps
  integer  :: n_timesteps
  !> Vector of time instants
  real(wp) , allocatable :: time_vec(:)
  !> Actual time
  real(wp) :: time
  !> Previous time
  real(wp) :: time_old
  !> ndt between 2 wake updates
  integer :: ndt_update_wake
  !> debug level
  integer :: debug_level
  !> file particles
  character(len=max_char_len) :: particles_file

  !Physical parameters:
  !> Free stream pressure
  real(wp) :: P_inf
  !> Free stream density
  real(wp) :: rho_inf
  !> Free stream velocity
  real(wp) :: u_inf(3)
  !> Reference velocity (magnitude of u_inf unless specified)
  real(wp) :: u_ref
  !> Free stream speed of sound
  real(wp) :: a_inf
  !> Free stream dynamic viscosity
  real(wp) :: mu_inf
  !> Free stream kinematic viscosity
  real(wp) :: nu_inf

  !Wake
  !> Number of wake particles
  integer :: n_wake_particles
  !> Minimum and maximum of the particles box
  real(wp) :: particles_box_min(3)
  real(wp) :: particles_box_max(3)
  !> Minimun particle magnitude allowed (suppress if lower)
  real(wp) :: mag_threshold

  !> Wake initial condition
  integer  :: part_n0 
  real(wp), allocatable :: part_pos0(:,:)
  real(wp), allocatable :: part_vort0_dir(:,:)
  real(wp), allocatable :: part_vort0_mag(:)
  real(wp), allocatable :: part_vol(:)

  !Method parameters
  !> Rankine Radius for vortices
  real(wp) :: RankineRad
  !> Vortex Radius for vortex particles
  real(wp) :: VortexRad
  !> Complete cutoff radius
  real(wp) :: CutoffRad
  !> use the vortex stretching or not
  logical :: use_vs
  !> use the divergence filtering
  logical :: use_divfilt
  !> Pedrizzetti relaxation coefficient for divergence filtering
  real(wp):: alpha_divfilt
  !> time scale of the divergence filter
  real(wp) :: filt_eta
  !> use the vorticity diffusion or not
  logical :: use_vd
  !> use turbulent viscosity or not
  logical :: use_tv
  !> integrators
  character(len=max_char_len) :: integrator
  !> use reformulated formulation rVPM (Alvarez 2023)
  logical   :: use_reformulated
    !> rVPM coefficients
    real(wp)  :: f                   
    real(wp)  :: g

  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
    !> Employing the FMM method also for panels
    logical :: use_fmm_pan
    !> Size of the Octree box
    real(wp) :: BoxLength
    !> Number of boxes in each direction
    integer :: NBox(3)
    !> Origin of the Octree system of boxes
    real(wp) :: OctreeOrigin(3)
    !> Number of Octree levels
    integer :: NOctreeLevels
    !> Minimum number of particles for each box
    integer :: MinOctreePart
    !> Multipole expansion degree
    integer :: MultipoleDegree
    !> Use dynamic levels
    logical :: use_dyn_layers
      !> Maximum number of octree levels
      integer :: NMaxOctreeLevels
      !> Time ratio that triggers the increase of levels
      real(wp) :: LeavesTimeRatio
    !> use particles redistribution
    logical :: use_pr
      !> Level at which is checked the presence of panels
      integer :: lvl_solid
      real(wp) :: part_redist_ratio

  !> Output interval
  real(wp) :: dt_out
  !> Basename
  character(len=max_char_len) :: basename
  !> Restart from file
  logical :: restart_from_file
  !> Restart file
  character(len=max_char_len) :: restart_file
  !> Reset the time after restart
  logical :: reset_time

!----------------------------------------------------------------------

contains
  procedure, pass(this) :: save_param => save_sim_param
end type t_sim_param

type(t_sim_param) :: sim_param
character(len=*), parameter     :: this_mod_name = 'mod_test'

real(wp), allocatable :: particlesMat(:,:)
!----------------------------------------------------------------------
contains

!> Subroutines for main dust, dust_pre and dust_post prms creation 
!  to avoid clutter in source files
!> Create parameter object for dust main parameters
subroutine create_param_test_particle(prms)
  type(t_parse), intent(inout) :: prms
  
  !> Define the parameters to be read
  !> Time
  call prms%CreateRealOption('tstart', "Starting time")
  call prms%CreateRealOption('tend',   "Ending time")
  call prms%CreateRealOption('dt',     "time step")
  call prms%CreateIntOption ('timesteps', "number of timesteps")
  call prms%CreateRealOption('dt_out', "output time interval")
  call prms%CreateRealOption('dt_debug_out', "debug output time interval")
  call prms%CreateIntOption ('ndt_update_wake', 'n. dt between two wake updates', '1')

  call prms%CreateStringOption('basename','output basename','./')
  call prms%CreateStringOption('basename_debug','output basename for debug','./')
  call prms%CreateLogicalOption('output_start', "output values at starting &
                                                            & iteration", 'F')
  call prms%CreateIntOption('debug_level', 'Level of debug verbosity/output', '1') 
  call prms%CreateStringOption('particles_file', 'file with particles initial condition', 'particles.dat') 
  
  !> Restart
  call prms%CreateLogicalOption('restart_from_file','restarting from file?','F')
  call prms%CreateStringOption('restart_file','restart file name')
  call prms%CreateLogicalOption('reset_time','reset the time from previous execution?','F')

  !> Parameters: reference conditions 
  call prms%CreateRealArrayOption('u_inf', "free stream velocity", '(/1.0, 0.0, 0.0/)')
  call prms%CreateRealOption('u_ref', "reference velocity")             
  call prms%CreateRealOption('P_inf', "free stream pressure", '101325')    
  call prms%CreateRealOption('rho_inf', "free stream density", '1.225')   
  call prms%CreateRealOption('a_inf', "Speed of sound", '340.0')        ! m/s   for dimensional sim
  call prms%CreateRealOption('mu_inf', "Dynamic viscosity", '0.000018') ! kg/ms

  call prms%CreateIntOption('n_wake_particles', 'number of wake particles', '100000')
  call prms%CreateRealArrayOption('particles_box_min', 'min coordinates of the &
                                  &particles bounding box', '(/-10.0, -10.0, -10.0/)')
  call prms%CreateRealArrayOption('particles_box_max', 'max coordinates of the &
                                  &particles bounding box', '(/10.0, 10.0, 10.0/)')
  call prms%CreateRealOption('mag_threshold', "Minimum particle magnitude allowed", '1.0e-9')

  !> Regularisation 
  call prms%CreateRealOption('rankine_rad', &
        "Radius of Rankine correction for vortex induction near core", '0.1')
  call prms%CreateRealOption('vortex_rad', &
        "Radius of vortex core, for particles", '0.1')
  call prms%CreateRealOption('k_vortex_rad', &
        "Radius coefficient of vortex core, for particles", '1.0') ! default is ON
  call prms%CreateRealOption('cutoff_rad', &
        "Radius of complete cutoff  for vortex induction near core", '0.001')

  !> Octree and multipole data 
  call prms%CreateLogicalOption('fmm','Employ fast multipole method?','T')
  call prms%CreateLogicalOption('fmm_panels','Employ fast multipole method &
                                &also for panels?','F')
  call prms%CreateRealOption('box_length','length of the octree box')
  call prms%CreateIntArrayOption('n_box','number of boxes in each direction')
  call prms%CreateRealArrayOption( 'octree_origin', "rigid wake velocity" )
  call prms%CreateIntOption('n_octree_levels','number of octree levels')
  call prms%CreateIntOption('min_octree_part','minimum number of octree particles')
  call prms%CreateIntOption('multipole_degree','multipole expansion degree')
  call prms%CreateLogicalOption('dyn_layers','Use dynamic layers','F')
  call prms%CreateIntOption('nmax_octree_levels','maximum number of octree levels')
  call prms%CreateRealOption('leaves_time_ratio','Ratio that triggers the &
                                            &increase of the number of levels')

    !> Models options
  call prms%CreateLogicalOption('vortstretch','Employ vortex stretching','T')
  call prms%CreateLogicalOption('vortstretch_from_elems','Employ vortex stretching&
                                & from geometry elements','F')
  call prms%CreateLogicalOption('divergence_filtering','Employ divergence filtering','T')
  call prms%CreateRealOption('alpha_divfilt','Pedrizzetti relaxation coefficient','0.3')
  call prms%CreateLogicalOption('diffusion','Employ vorticity diffusion','T')
  call prms%CreateLogicalOption('turbulent_viscosity','Employ turbulent &
                                &viscosity','F') 
  call prms%CreateLogicalOption('viscosity_effects','Simulate viscosity &
                                                                & effects','F')
  call prms%CreateLogicalOption('particles_redistribution','Employ particles &
                                                          &redistribution','F')
  call prms%CreateIntOption('octree_level_solid','Level at which the panels &
                            & are considered for particles redistribution')
  call prms%CreateRealOption('particles_redistribution_ratio','How many times &
            &a particle need to be smaller than the average of the cell to be&
            & eliminated','3.0')

  !> Reformulated formulation                                         
  call prms%CreateLogicalOption('reformulated','Employ rVPM by Alvarez','F')
  call prms%CreateRealOption('f','rVPM coefficient f','0.0')
  call prms%CreateRealOption('g','rVPM coefficient g','0.2')

  !> Integrators
  call prms%CreateStringOption('integrator', 'integrator solver: Euler or low storage RK', &
                              'Euler') 


end subroutine create_param_test_particle  

!----------------------------------------------------------------------
subroutine save_sim_param(this, loc)
  class(t_sim_param) :: this
  integer(h5loc), intent(in) :: loc

  call write_hdf5_attr(this%t0, 't0', loc)
  call write_hdf5_attr(this%dt, 'dt', loc)
  call write_hdf5_attr(this%tend, 'tend', loc)
  call write_hdf5_attr(this%P_inf, 'P_inf', loc)
  call write_hdf5_attr(this%rho_inf, 'rho_inf', loc)
  call write_hdf5_attr(this%u_inf, 'u_inf', loc)
  call write_hdf5_attr(this%u_ref, 'u_ref', loc)
  call write_hdf5_attr(this%a_inf, 'a_inf', loc)
  call write_hdf5_attr(this%mu_inf, 'mu_inf', loc)
  call write_hdf5_attr(this%n_wake_particles, 'n_wake_particles', loc)
  call write_hdf5_attr(this%particles_box_min, 'particles_box_min', loc)
  call write_hdf5_attr(this%particles_box_max, 'particles_box_max', loc)
  call write_hdf5_attr(this%mag_threshold, 'mag_threshold', loc)

  call write_hdf5_attr(this%RankineRad, 'RankineRad', loc)
  call write_hdf5_attr(this%CutoffRad, 'CutoffRad', loc)
  call write_hdf5_attr(this%use_vs, 'Vortstretch', loc)
  if(this%use_vs) then
    call write_hdf5_attr(this%use_divfilt, 'DivergenceFiltering', loc)
    call write_hdf5_attr(this%alpha_divfilt, 'AlphaDivFilt', loc)
  endif
  call write_hdf5_attr(this%use_vd, 'vortdiff', loc)
  call write_hdf5_attr(this%use_tv, 'turbvort', loc)

  call write_hdf5_attr(this%use_fmm, 'use_fmm', loc)
  if(this%use_fmm) then
    call write_hdf5_attr(this%BoxLength, 'BoxLength', loc)
    call write_hdf5_attr(this%Nbox, 'Nbox', loc)
    call write_hdf5_attr(this%OctreeOrigin, 'OctreeOrigin', loc)
    call write_hdf5_attr(this%NOctreeLevels, 'NOctreeLevels', loc)
    call write_hdf5_attr(this%MinOctreePart, 'MinOctreePart', loc)
    call write_hdf5_attr(this%MultipoleDegree, 'MultipoleDegree', loc)
  endif

  call write_hdf5_attr(this%use_reformulated, 'use_reformulated', loc)
  if(this%use_reformulated) then
    call write_hdf5_attr(this%f, 'f', loc)
    call write_hdf5_attr(this%g, 'g', loc)
  endif


  call write_hdf5_attr(this%dt_out, 'dt_out', loc)
  call write_hdf5_attr(this%basename, 'basename', loc)
  
end subroutine save_sim_param

!> Initialize all the parameters reading them from the the input file
subroutine init_sim_param(sim_param, prms, nout, output_start)
  class(t_sim_param)          :: sim_param
  type(t_parse)               :: prms
  integer, intent(inout)      :: nout
  logical, intent(inout)      :: output_start
  
  character(len=*), parameter :: this_sub_name = 'init_sim_param'
  !> Timing
  sim_param%t0                  = getreal(prms, 'tstart')
  sim_param%tend                = getreal(prms, 'tend')
  if(CountOption(prms,'dt') .gt. 0) then
    if( CountOption(prms,'timesteps') .gt. 0) then
      call error(this_mod_name, this_sub_name, 'Both number of timesteps and dt are&
                                              & set, but only one of the two can be specified')
    else
      !> get dt and compute number of timesteps
      sim_param%dt     = getreal(prms, 'dt')
      sim_param%n_timesteps = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
                              !(+1 for the zero time step)
    endif
  else
    !> get number of steps, compute dt
    sim_param%n_timesteps = getint(prms, 'timesteps')
    sim_param%dt =  (sim_param%tend-sim_param%t0)/&
                      real(sim_param%n_timesteps,wp)
    sim_param%n_timesteps = sim_param%n_timesteps + 1
                            !add one for the first step
  endif 

  sim_param%dt_out              = getreal(prms,'dt_out')
  sim_param%n_timesteps         = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
  sim_param%debug_level         = 1
  sim_param%debug_level         = getint(prms, 'debug_level')  
  
  !> file dat 
  sim_param%particles_file      = getstr(prms, 'particles_file') 
  !> Reference environment values
  sim_param%P_inf               = getreal(prms,'P_inf')
  sim_param%rho_inf             = getreal(prms,'rho_inf')
  sim_param%a_inf               = getreal(prms,'a_inf')
  sim_param%mu_inf              = getreal(prms,'mu_inf')
  sim_param%nu_inf              = sim_param%mu_inf/sim_param%rho_inf
  sim_param%u_inf               = getrealarray(prms, 'u_inf', 3)
  
  !> Check on reference velocity
  if ( countoption(prms,'u_ref') .gt. 0 ) then
    sim_param%u_ref = getreal(prms, 'u_ref')
  else
    sim_param%u_ref = norm2(sim_param%u_inf)
    if (sim_param%u_ref .le. 0.0_wp) then
      call error(this_mod_name, this_sub_name,'No reference velocity u_ref provided but &
      &zero free stream velocity. Provide a non-zero reference velocity. &
      &Stopping now before producing invalid results')
    endif
  end if
  
  !> Wake parameters
  sim_param%n_wake_particles      = getint(prms, 'n_wake_particles')
  sim_param%particles_box_min     = getrealarray(prms, 'particles_box_min',3)
  sim_param%particles_box_max     = getrealarray(prms, 'particles_box_max',3)
  sim_param%mag_threshold         = getreal(prms, 'mag_threshold')

  sim_param%basename              = getstr(prms, 'basename')
  !Integrators 
  sim_param%integrator            = getstr(prms, 'integrator')
  !> Manage restart
  sim_param%restart_from_file             = getlogical(prms,'restart_from_file')
  if (sim_param%restart_from_file) then

    sim_param%reset_time                  = getlogical(prms,'reset_time')
    sim_param%restart_file                = getstr(prms,'restart_file')
    
    !> Removing leading "./" if present to avoid issues when restarting
    if(sim_param%basename(1:2) .eq. './') sim_param%basename = sim_param%basename(3:)
    
    if(sim_param%restart_file(1:2) .eq. './') sim_param%restart_file = sim_param%restart_file(3:)
    
    !call printout('RESTART: restarting from file: '//trim(sim_param%restart_file))
    !sim_param%GeometryFile = sim_param%restart_file(1:len(trim(sim_param%restart_file))-11) //'geo.h5'

    !restarting the same simulation, advance the numbers
    if(sim_param%restart_file(1:len(trim(sim_param%restart_file))-12).eq. &
                                                trim(sim_param%basename)) then
    read(sim_param%restart_file(len(trim(sim_param%restart_file))-6:len(trim(sim_param%restart_file))-3),*) nout
      call printout('Identified restart from the same simulation, keeping the&
                    & previous output numbering')
      !> avoid rewriting the same timestep
      output_start = .false.
    endif
    if(.not. sim_param%reset_time) call load_time(sim_param%restart_file, sim_param%t0)
  endif

  !> Method parameters
  sim_param%RankineRad            = getreal(prms,    'rankine_rad')
  sim_param%VortexRad             = getreal(prms,    'vortex_rad')
  sim_param%CutoffRad             = getreal(prms,    'cutoff_rad')
  sim_param%use_vs                = getlogical(prms, 'vortstretch')
  sim_param%use_vd                = getlogical(prms, 'diffusion')
  sim_param%use_divfilt           = getlogical(prms, 'divergence_filtering')
  sim_param%use_tv                = getlogical(prms, 'turbulent_viscosity')
  sim_param%alpha_divfilt         = getreal(prms,    'alpha_divfilt')

  !> Reformulated formulation (Alvarez rVPM 2023)
  sim_param%use_reformulated      = getlogical(prms, 'reformulated')

  if(sim_param%use_reformulated) then
    sim_param%f                   = getreal(prms, 'f')
    sim_param%g                   = getreal(prms, 'g')
  endif

  !> Octree and FMM parameters
  sim_param%use_fmm                       = getlogical(prms, 'fmm')

  if(sim_param%use_fmm) then
    sim_param%use_fmm_pan                 = getlogical(prms, 'fmm_panels')
    sim_param%BoxLength                   = getreal(prms, 'box_length')
    sim_param%NBox                        = getintarray(prms, 'n_box',3)
    sim_param%OctreeOrigin                = getrealarray(prms, 'octree_origin',3)
    sim_param%NOctreeLevels               = getint(prms, 'n_octree_levels')
    sim_param%MinOctreePart               = getint(prms, 'min_octree_part')
    sim_param%MultipoleDegree             = getint(prms,'multipole_degree')
    sim_param%use_dyn_layers              = getlogical(prms,'dyn_layers')

    if(sim_param%use_dyn_layers) then
      sim_param%NMaxOctreeLevels          = getint(prms, 'nmax_octree_levels')
      sim_param%LeavesTimeRatio           = getreal(prms, 'leaves_time_ratio')
    else
      sim_param%NMaxOctreeLevels          = sim_param%NOctreeLevels
    endif

    sim_param%use_pr                      = getlogical(prms, 'particles_redistribution')

    if(sim_param%use_pr) then
      sim_param%part_redist_ratio         = getreal(prms,'particles_redistribution_ratio')
      if ( countoption(prms,'octree_level_solid') .gt. 0 ) then
        sim_param%lvl_solid               = getint(prms, 'octree_level_solid')
      else
        sim_param%lvl_solid               = max(sim_param%NOctreeLevels-2,1)
      endif
    endif
  else
    sim_param%use_fmm_pan = .false.
  endif

  sim_param%ndt_update_wake       = 1

end subroutine init_sim_param

subroutine load_time(filename, time)
  character(len=*), intent(in) :: filename
  real(wp), intent(out)        :: time

  integer(h5loc)               :: floc

  call open_hdf5_file(filename, floc)
  call read_hdf5(time,'time',floc)
  call close_hdf5_file(floc)

end subroutine load_time

end module mod_test