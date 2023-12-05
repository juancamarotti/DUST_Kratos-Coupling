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
!!          Andrea Colli
!!          Alberto Savino
!!=========================================================================


!> Module to treat the whole wake
module mod_wake_test

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_math, only: &
  cross, infinite_plate_spline, tessellate

use mod_test, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_aeroel, only: &
  c_elem, c_vort_elem, &
  t_elem_p, t_vort_elem_p

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_hdf5_io, only: &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  write_hdf5_attr, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5

use mod_octree_test, only: &
  t_octree, sort_particles, calculate_multipole, apply_multipole

use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: t_wake, initialize_wake, update_wake, &
          prepare_wake, complete_wake
private


!> Type containing wake panels information
type :: t_wake

  !! Particles data

  !> Maximum number of particles
  integer :: nmax_prt

  !> Actual number of particles
  integer :: n_prt

  !> Wake particles
  type(t_vortpart), allocatable :: wake_parts(:)

  !> Magnitude of particles vorticity
  real(wp), allocatable :: prt_ivort(:)

  !> Wake particles pointer
  type(t_vortpart_p), allocatable :: part_p(:)

  !> Bounding box
  real(wp) :: part_box_min(3), part_box_max(3)

  type(t_vort_elem_p), allocatable :: vort_p(:)

end type

!> Class to change methods from different wake implementations
type, abstract :: c_wake_mov
  contains
  procedure(i_get_vel), deferred, pass(this) :: get_vel
end type

abstract interface
  subroutine i_get_vel(this, wake, pos, vel)
    import                                :: c_wake_mov, wp, t_wake
    class(c_wake_mov)                     :: this
    type(t_wake), intent(in)              :: wake
    real(wp), intent(in)                  :: pos(3)
    real(wp), intent(out)                 :: vel(3)
  end subroutine
end interface

type, extends(c_wake_mov) :: t_free_wake
contains
  procedure, pass(this) :: get_vel => get_vel_free
end type

class(c_wake_mov), allocatable  :: wake_movement
character(len=max_char_len)     :: msg
real(t_realtime)                :: t1 , t0
character(len=*), parameter     :: this_mod_name='mod_wake_test'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the particle wake
subroutine initialize_wake(wake)
  type(t_wake), intent(out),target     :: wake
  integer                              :: ip


allocate(t_free_wake::wake_movement)

!Particles

wake%nmax_prt = sim_param%part_n0
allocate(wake%wake_parts(wake%nmax_prt))
allocate(wake%prt_ivort(wake%nmax_prt))

wake%n_prt = wake%nmax_prt
allocate(wake%part_p(wake%n_prt))


do ip = 1,wake%n_prt
   wake%wake_parts(ip)%mag      => wake%prt_ivort(ip)
   wake%wake_parts(ip)%cen      = sim_param%part_pos0(ip,:)
   wake%wake_parts(ip)%dir      = sim_param%part_vort0_dir(ip,:)
   wake%wake_parts(ip)%mag      = sim_param%part_vort0_mag(ip)
   wake%wake_parts(ip)%r_Vortex = sim_param%VortexRad
   wake%wake_parts(ip)%free     = .false.
   wake%part_p(ip)%p            => wake%wake_parts(ip)
enddo

wake%part_box_min = sim_param%particles_box_min
wake%part_box_max = sim_param%particles_box_max


allocate(wake%vort_p(wake%n_prt))
end subroutine initialize_wake

!----------------------------------------------------------------------


!----------------------------------------------------------------------

!> Prepare the wake before the timestep
!!
!! Mainly prepare all the structures for the octree
subroutine prepare_wake(wake, octree)
  type(t_wake), intent(inout), target   :: wake
  type(t_octree), intent(inout)         :: octree
  integer                               :: k, ip, ir, iw, ie, n_end_vort

  if (sim_param%use_fmm) then
    call sort_particles(wake%wake_parts, wake%n_prt, octree)
    call calculate_multipole(wake%part_p, octree)
  endif

  !==>Recreate structures and pointers, if particles are present
  if( wake%n_prt.gt.0 ) then

    !Recreate the pointer vector
    if(allocated(wake%part_p)) then 
      deallocate(wake%part_p)
    endif

    allocate(wake%part_p(wake%n_prt))

    
    k = 1
    do ip = 1, wake%n_prt
      do ir=k,wake%nmax_prt
        if(.not. wake%wake_parts(ir)%free) then
          k = ir+1
          wake%part_p(ip)%p => wake%wake_parts(ir)
          wake%vort_p(ip)%p => wake%wake_parts(ir)
          exit
        endif
      enddo
    enddo
  endif
end subroutine prepare_wake

!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels 
!  Brings them to the next time step
!  Only updates the "old" panels, ie not the first two, and existing particles
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll)
!! elements
subroutine update_wake(wake, octree)
  type(t_wake), intent(inout), target :: wake
  type(t_octree), intent(inout)       :: octree

  integer                             :: iw, ipan, ie, ip, np, iq
  integer                             :: id, ir
  real(wp)                            :: pos_p(3), vel_p(3)
  real(wp)                            :: str(3), stretch(3)
  real(wp)                            :: ru(3), rotu(3)
  real(wp)                            :: df(3), diff(3)
  real(wp)                            :: hcas_vel(3)
  real(wp), allocatable               :: point_old(:,:,:)
  real(wp), allocatable               :: points(:,:,:)
  logical                             :: increase_wake
  integer                             :: size_old
  character(len=*), parameter         :: this_sub_name='update_wake'

  !==>    Particles: evolve the position in time

  !calculate the velocities at the points
!$omp parallel do private(pos_p, vel_p, ip, iq,  stretch, diff, df, str, ru, rotu)
  do ip = 1, wake%n_prt
    wake%part_p(ip)%p%vel_old = wake%part_p(ip)%p%vel
    wake%part_p(ip)%p%stretch_old = wake%part_p(ip)%p%stretch
    wake%part_p(ip)%p%stretch = 0.0_wp
    wake%part_p(ip)%p%rotu = 0.0_wp

    !If not using the fast multipole, update particles position now
    if (.not.sim_param%use_fmm) then
      pos_p = wake%part_p(ip)%p%cen

      call wake_movement%get_vel(wake, pos_p, vel_p)

      wake%part_p(ip)%p%vel =  vel_p
      !if using vortex stretching, calculate it now
      if(sim_param%use_vs) then
        stretch = 0.0_wp
        rotu = 0.0_wp
        do iq = 1, wake%n_prt
          if (ip.ne.iq) then
            call wake%part_p(iq)%p%compute_stretch(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, str)
            ! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            stretch = stretch + str/(4.0_wp*pi)

            if(sim_param%use_divfilt) then
              call wake%part_p(iq)%p%compute_rotu(wake%part_p(ip)%p%cen, &
                    wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, ru)
              rotu = rotu + ru/(4.0_wp*pi)

            endif
          endif
        enddo
      
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + stretch
        
        if(sim_param%use_divfilt) then 
          wake%part_p(ip)%p%rotu = wake%part_p(ip)%p%rotu + rotu
        endif
      
      endif !use_vs

      !if using the vortex diffusion, calculate it now
      if(sim_param%use_vd) then
        diff = 0.0_wp

        do iq = 1, wake%n_prt

          if (ip.ne.iq) then
            call wake%part_p(iq)%p%compute_diffusion(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, &
                  wake%part_p(ip)%p%r_Vortex, df)
            diff = diff + df*sim_param%nu_inf
          endif

        enddo !iq
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + diff
      endif !use_vd

    end if ! not use_fmm

  enddo
!$omp end parallel do

  if (sim_param%use_fmm) then
    t0 = dust_time()
    call apply_multipole(wake%part_p, octree)
    t1 = dust_time()
    write(msg,'(A,F9.3,A)') 'Multipoles calculation: ' , t1 - t0,' s.'
    if(sim_param%debug_level.ge.3) call printout(msg)
    write(msg,'(A,I0)') 'Number of particles: ' , wake%n_prt
    if(sim_param%debug_level.ge.1) call printout(msg)
  endif


end subroutine update_wake

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!! Completes the updating to the next time step begun in update_wake
!! first and second row are updated to the next step and new particles are
!! created if necessary; they will appear at the save_date in the next time step
subroutine complete_wake(wake)
  type(t_wake), target, intent(inout)   :: wake

  integer                               :: p1, p2
  integer                               :: ip, iw, id, is, nprev
  real(wp)                              :: dist(3) , vel_te(3), pos_p(3)
  real(wp)                              :: dir(3), partvec(3), ave, alpha_p(3), alpha_p_n
  integer                               :: n_part
  real(wp)                              :: vel_in(3), vel_out(3), wind(3), filt_eta
  
  character(len=max_char_len)           :: msg
  character(len=*), parameter           :: this_sub_name='complete_wake'

!==> Particles: update the position and intensity in time, avoid penetration
!               and chech if remain into the boundaries
n_part = wake%n_prt
!$omp parallel do schedule(dynamic,4) private(ip,pos_p,alpha_p,alpha_p_n,vel_in,vel_out)
  do ip = 1, n_part
    if(.not. wake%part_p(ip)%p%free) then
      pos_p = wake%part_p(ip)%p%cen + wake%part_p(ip)%p%vel* &
              sim_param%dt*real(sim_param%ndt_update_wake,wp)
      if(all(pos_p .ge. wake%part_box_min) .and. &
          all(pos_p .le. wake%part_box_max)) then
        wake%part_p(ip)%p%cen = pos_p
        if(sim_param%use_vs .or. sim_param%use_vd) then

          !add filtering (Pedrizzetti Relaxation)
          if(sim_param%use_divfilt) then
            filt_eta = sim_param%alpha_divfilt/sim_param%dt
            wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch - &
              filt_eta/real(sim_param%ndt_update_wake,wp)*( wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag - &
              wake%part_p(ip)%p%rotu*wake%part_p(ip)%p%mag/norm2(wake%part_p(ip)%p%rotu))
          endif

          !Explicit Euler
          alpha_p = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag + &
                          wake%part_p(ip)%p%stretch* &
                          sim_param%dt*real(sim_param%ndt_update_wake,wp)
          alpha_p_n = norm2(alpha_p)

! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
          if(alpha_p_n .ne. 0.0_wp) &
              wake%part_p(ip)%p%dir = alpha_p/alpha_p_n
          endif
      else
        wake%part_p(ip)%p%free = .true.
!$omp atomic update
        wake%n_prt = wake%n_prt -1
!$omp end atomic
      endif
    endif
  enddo
!$omp end parallel do

end subroutine complete_wake

!----------------------------------------------------------------------

subroutine compute_vel_from_all(wake, pos, vel)
  type(t_wake), intent(in)        :: wake
  real(wp), intent(in)            :: pos(3)
  real(wp), intent(out)           :: vel(3)

  integer                         :: ie
  real(wp)                        :: v(3)

  vel = 0.0_wp

  !calculate the influence of particles
  do ie=1,size(wake%part_p)
    call wake%part_p(ie)%p%compute_vel(pos, v)
    vel = vel + v/(4*pi)
  enddo

end subroutine compute_vel_from_all

!----------------------------------------------------------------------

subroutine get_vel_free(this, wake, pos, vel)
  class(t_free_wake)                    :: this
  type(t_wake), intent(in)              :: wake
  real(wp), intent(in)                  :: pos(3)
  real(wp), intent(out)                 :: vel(3)

  call compute_vel_from_all(wake, pos, vel)

  !vel = vel + sim_param%u_inf
  vel = vel + variable_wind(pos, sim_param%time) 

end subroutine get_vel_free

!---------------------------------------------------------------------

end module mod_wake_test
