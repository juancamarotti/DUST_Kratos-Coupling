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
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================

module mod_test

use mod_param, only: &
    wp, max_char_len

use mod_hdf5_io, only: &
    h5loc, write_hdf5_attr, open_hdf5_file, close_hdf5_file, read_hdf5

use mod_handling, only: &
    error, warning, info, printout, dust_time, t_realtime, internal_error

implicit none

public :: t_sim_param, sim_param, init_sim_param

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

  !> Wake initial condition
  integer  :: part_n0 
  real(wp) :: part_pos0(3)
  real(wp) :: part_vel0(3)
  real(wp) :: part_vort0_dir(3)
  integer  :: part_vort0_mag

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

  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
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
    !> Maximum number of octree levels
    integer :: NMaxOctreeLevels

  !> Output interval
  real(wp) :: dt_out
  !> Basename
  character(len=max_char_len) :: basename

!----------------------------------------------------------------------

contains
  procedure, pass(this) :: save_param => save_sim_param
end type t_sim_param

type(t_sim_param) :: sim_param

!----------------------------------------------------------------------
contains
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

  call write_hdf5_attr(this%dt_out, 'dt_out', loc)
  call write_hdf5_attr(this%basename, 'basename', loc)
  
end subroutine save_sim_param

!> Initialize all the parameters reading them from the the input file
subroutine init_sim_param(sim_param)
    class(t_sim_param)          :: sim_param 

    !> Timing
    sim_param%t0                  = 0.0_wp
    sim_param%tend                = 1.0_wp
    sim_param%dt                  = 0.001_wp
    sim_param%dt_out              = sim_param%dt
    sim_param%n_timesteps         = ceiling((sim_param%tend-sim_param%t0)/sim_param%dt) + 1
    sim_param%debug_level         = 1

    !> Reference environment values
    sim_param%P_inf               = 101325.0_wp
    sim_param%rho_inf             = 1.225_wp
    sim_param%a_inf               = 340.0_wp
    sim_param%mu_inf              = 0.000018_wp
    sim_param%nu_inf              = sim_param%mu_inf/sim_param%rho_inf
    sim_param%u_inf               = 0.0_wp
    sim_param%u_ref               = 10.0_wp  
    
    !> Wake parameters
    sim_param%n_wake_particles      = 100000
    sim_param%particles_box_min     = (/ -10.0_wp, -10.0_wp, -10.0_wp /)
    sim_param%particles_box_max     = (/ +10.0_wp, +10.0_wp, +10.0_wp /)

    !> Wake initial condition
    sim_param%part_n0               = 1
    sim_param%part_vel0             = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    sim_param%part_pos0             = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
    sim_param%part_vort0_dir        = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
    sim_param%part_vort0_mag        = 1.0_wp

    !> Names
    sim_param%basename              = 'Output'
  
    !> Method parameters
    sim_param%RankineRad            = 0.1_wp
    sim_param%VortexRad             = 0.1_wp
    sim_param%CutoffRad             = 0.001_wp
    sim_param%use_vs                = .true.
    sim_param%use_vd                = .true.
    sim_param%use_tv                = .false.
    sim_param%use_divfilt           = .true.
    sim_param%alpha_divfilt         = 0.3_wp
    !> Octree and FMM parameters
    sim_param%use_fmm               = .false.
  
    if(sim_param%use_fmm) then
      sim_param%BoxLength           = 20
      sim_param%NBox                = (/ 1, 1, 1 /)
      sim_param%OctreeOrigin        = (/ -10.0_wp, -10.0_wp, -10.0_wp /)
      sim_param%NOctreeLevels       = 6
      sim_param%MinOctreePart       = 5
      sim_param%MultipoleDegree     = 2
      sim_param%NMaxOctreeLevels    = sim_param%NOctreeLevels 
    endif

    sim_param%ndt_update_wake       = 1

end subroutine init_sim_param

end module mod_test