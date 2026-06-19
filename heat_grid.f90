module heat_grid

        use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_int, c_char, &
                c_null_char, c_ptr, c_null_ptr, c_f_pointer
        implicit none 
        private

        !--- Two grid arrays for double buffering ---
        real(c_double), allocatable, target, save :: grid_src(:), grid_dst(:)       ! module-level storage

        !--- Public Interface ---
        public :: allocate_grid_src, allocate_grid_dst, destroy_grids
        public :: fill_initial,  apply_boundary_conditions, get_element
        public :: solve_poisson, write_grid_binary    
        public :: solve_gauss_seidel, solve_sor

contains

        ! ----------------------------------------------------------------------------------
        ! Allocate the source grid
        ! ----------------------------------------------------------------------------------
        function allocate_grid_src(nx, ny) bind( C, name = "allocate_grid_src") result(ptr)
                integer(c_int), intent(in), value :: nx, ny
                type(c_ptr) :: ptr
                
                if(allocated(grid_src)) then
                        deallocate(grid_src)
                end if

                allocate(grid_src(nx * ny))         
                grid_src = 0.0_c_double             ! initialize to zero
                ptr = c_loc(grid_src(1))           ! get c pointer to first element
        end function allocate_grid_src

        ! ----------------------------------------------------------------------------------
        ! Allocate the destination grid
        ! ----------------------------------------------------------------------------------
        function allocate_grid_dst(nx, ny) bind( C, name = "allocate_grid_dst") result(ptr)
                integer(c_int), intent(in), value :: nx, ny
                type(c_ptr) :: ptr
                
                if(allocated(grid_dst)) then
                        deallocate(grid_dst)
                end if

                allocate(grid_dst(nx * ny))         
                grid_dst = 0.0_c_double             ! initialize to zero
                ptr = c_loc(grid_dst(1))           ! get c pointer to first element
        end function allocate_grid_dst


        ! -----------------------------------------------------------------------------------
        ! Destroy both grids
        ! -----------------------------------------------------------------------------------
        subroutine destroy_grids() bind(C, name="destroy_grids")
                if (allocated(grid_src)) deallocate(grid_src)
                if (allocated(grid_dst)) deallocate(grid_dst)
        end subroutine destroy_grids

        ! -----------------------------------------------------------------------------------
        ! Fill the grid with initial conditions (example: a Gaussian bump)
        ! ----------------------------------------------------------------------------------
        subroutine fill_initial(grid_ptr, nx, ny) bind(C, name="fill_initial")
                type(c_ptr), intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx,ny

                real(c_double), pointer:: arr(:)
                integer :: i, j, idx
                real(c_double) :: x, y, center_x, center_y, sigma

                ! Associates a Fortran pointer with the C pointer
                call c_f_pointer(grid_ptr, arr, [nx*ny])

                center_x = (nx-1)/2.0_c_double
                center_y = (ny-1)/2.0_c_double
                sigma = min(nx, ny)/6.0_c_double

                do j=1, ny
                        do i=1, nx
                                x = real(i-1, c_double)
                                y = real(j-1, c_double)
                                idx = (j-1)*nx + i
                                arr(idx) = exp(-((x-center_x)**2+(y-center_y)**2) / (2.0_c_double * sigma **2))
                        end do
                end do
        end subroutine fill_initial

        ! ===========================================================================================
        ! NEW SUBROUTINE: Apply Dirichlet boundary conditions
        ! ===========================================================================================
        subroutine apply_boundary_conditions(grid_ptr, nx, ny, top, bottom, left, right) &
                        bind(C, name = "apply_boundary_conditions")
                
                type(c_ptr), intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny
                real(c_double), intent(in), value :: top, bottom, left, right

                real(c_double), pointer :: arr(:)
                integer :: i, j, idx

                ! Turn the C pointer into a Fortran 1-D array view
                call c_f_pointer(grid_ptr, arr, [nx * ny])

                ! -----------Bottom Edge ( j = 1 ) ---------------
                j = 1
                do i = 1, nx
                        idx = (j-1)*nx + i
                        arr(idx) = bottom
                end do

                ! -----------Top Edge ( j = ny ) ---------------
                j = ny
                do i = 1, nx
                        idx = (j-1)*nx + i
                        arr(idx) = top
                end do

                ! -----------Left Edge ( i = 1 ) ---------------
                i = 1
                do j = 1, ny
                        idx = (j-1)*nx + i
                        arr(idx) = left
                end do

                ! -----------Right Edge ( i = nx ) ---------------
                i = nx
                do j = 1, ny
                        idx = (j-1)*nx + i
                        arr(idx) = right
                end do


        end subroutine apply_boundary_conditions

        ! ------------------------------------------------------------------------------------------
        ! Get a single grid element ( for testing )
        ! ------------------------------------------------------------------------------------------
        function get_element(grid_ptr, i, j, nx, ny) bind(C, name="get_element") result(val)
                type(c_ptr), intent(in), value :: grid_ptr
                integer(c_int),intent(in), value :: i, j, nx, ny
                real(c_double) :: val
                real(c_double), pointer :: arr(:)

                call c_f_pointer(grid_ptr, arr, [nx*ny])
                val = arr((j-1)*nx + i)
        end function get_element

        ! ------------------------------------------------------------------------------------------
        ! Private: perform one Jacobi sweep from src to dst, return max change 
        ! ------------------------------------------------------------------------------------------
        subroutine jacobi_step_internal(src, dst, nx, ny, max_change)

        real(c_double), intent(in)  :: src(:)
        real(c_double), intent(out) :: dst(:)
        integer,        intent(in)  :: nx, ny
        real(c_double), intent(out) :: max_change
        integer :: i, j, idx

        real(c_double) :: stencil_val

        ! Copy the entire source to destination (include boundaries)
        dst = src

        ! Update only the interior (i = 2....nx-1,  j = 2....ny-1)
        do j = 2, ny-1
                do i = 2, nx - 1
                        idx = (j-1) * nx + i
                        stencil_val = 0.25_c_double * &
                                ( src((j-1)*nx + i-1)   & ! west
                                + src((j-1)*nx + i+1)   & ! east
                                + src((j-2)*nx + i  )   & ! south
                                + src((j)  *nx + i  )  )  ! North
                        dst(idx) = stencil_val
                end do 
        end do

        ! Compute maximum absolute change (interior only)
        max_change = 0.0_c_double
        do j = 2, ny -1
                do i = 2, nx -1
                        idx = (j-1) * nx + i
                        max_change = max(max_change, abs(dst(idx) - src(idx)))
                end do 
        end do

        end subroutine jacobi_step_internal

        ! ========================================================================================
        ! Public: solve the heat equation to steady state
        ! ========================================================================================
        subroutine solve_poisson(src_ptr, dst_ptr, nx, ny, tol, max_iter, actual_iter, residual ) &
                  bind(C, name="solve_poisson")

          type(c_ptr),    intent(in), value     :: src_ptr, dst_ptr
          integer(c_int), intent(in), value     :: nx, ny, max_iter
          real(c_double), intent(in), value     :: tol
          integer(c_int), intent(out)           :: actual_iter
          real(c_double), intent(out)           :: residual

          real(c_double), pointer       :: curr(:), next(:), tmp(:)
          real(c_double), pointer       :: src_array(:)         ! will point to src memory 
          integer       :: iter
          real(c_double)        :: change

          ! Associate Fortran pointers with the C memory
          call c_f_pointer(src_ptr, src_array, [nx*ny])
          call c_f_pointer(dst_ptr, next,      [nx*ny])

          ! Start: current = source, next = destination
          curr => src_array

          do iter = 1, max_iter
                ! Perform one Jacobi sweep from curr into next
                call jacobi_step_internal(curr, next, nx, ny, change)

                if (change < tol) then
                        ! Converged: ensure the result is in curr
                        if (.not. associated(curr, src_array)) then
                                src_array = curr        ! copy result back to src
                        end if
                        residual = change
                        actual_iter = iter
                        return
                end if

                ! Swap pointers from next iteration
                tmp  => curr
                curr => next
                next => tmp
          end do

          ! Maximum iterations reached: guarantee result in src
          if (.not. associated(curr, src_array)) then
                  src_array = curr
          end if
          residual = change             ! best residual from final iteration
          actual_iter = max_iter
       
        end subroutine solve_poisson

        ! =================================================================================
        ! Write the Grid to binary files (stram access, no record markers) 
        ! =================================================================================
        subroutine write_grid_binary(grid_ptr, nx, ny, filename) & 
                        bind(C, name="write_grid_binary")
        
                type(c_ptr),    intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny
                character(kind=c_char), intent(in):: filename(*)

                real(c_double), pointer :: arr(:)
                integer :: unit, iostat
                character(len=256) :: fname
                integer :: k

                ! Convert C string to Fortran string
                fname =' '
                k = 1
                do while(filename(k) /= c_null_char .and. k <= 256)
                        fname(k:k) = filename(k)
                        k = k+1
                end do

                call c_f_pointer(grid_ptr, arr, [nx*ny])

                ! Open a new file for unformatted sequential access
                open(newunit=unit, file=trim(fname), form='unformatted', &
                        access='stream', action='write', iostat=iostat)

                if (iostat /= 0) then
                        print *, 'Error opening file:  ', trim(fname)
                        return
                end if
                
                ! Write header: nx, ny
                write(unit) nx, ny

                ! Write the whole temperature array
                write(unit) arr(1:nx*ny)

                close(unit)

        end subroutine write_grid_binary

        ! =================================================================================
        ! Gauss-Seidel subroutine
        ! Perform one Gauss-Seidel sweep on grid (in-place).
        ! Returns maximal change (before update) as residual.
        ! =================================================================================
        subroutine gauss_seidel_step(grid_ptr, nx, ny, max_change) & 
                        bind(C, name="gauss_seidel_step")
                type(c_ptr),    intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny
                real(c_double), intent(out)       :: max_change

                real(c_double), pointer :: arr(:)
                integer :: i, j, idx
                real(c_double) :: new_val, old_val

                call c_f_pointer(grid_ptr, arr, [nx*ny])

                max_change = 0.0_c_double

                do j = 2, ny-1
                        do i = 2, nx-1
                                idx = (j-1)*nx + i
                                old_val = arr(idx)

                                ! New value uses updated west ( i -1 ) and south  ( j -1 )
                                ! and old east ( i + 1 ) and north (j+1)
                                new_val = 0.25_c_double * &
                                        ( arr((j-1)*nx + i-1)  &    ! west
                                        + arr((j-1)*nx + i+1)  &    ! east
                                        + arr((j-2)*nx + i  )  &    ! south
                                        + arr(  j  *nx + i  ) )     ! north   
                                arr(idx) = new_val
                                max_change = max(max_change, abs(new_val - old_val))
                        end do
                end do
        end subroutine gauss_seidel_step


        ! =================================================================================
        ! Gauss-Seidel with Successive Over Relaxation (SOR) subroutine
        ! =================================================================================
        subroutine sor_sweep(grid_ptr, nx, ny, omega, max_change)
                type(c_ptr),    intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny
                real(c_double), intent(in)        :: omega
                real(c_double), intent(out)       :: max_change

                real(c_double), pointer :: arr(:)
                integer :: i, j, idx
                real(c_double) :: new_gs, old_val

                call c_f_pointer(grid_ptr, arr, [nx*ny])

                max_change = 0.0_c_double

                do j = 2, ny-1
                        do i = 2, nx-1
                                idx = (j-1)*nx + i
                                old_val = arr(idx)

                                ! New value uses updated west ( i -1 ) and south  ( j -1 )
                                ! and old east ( i + 1 ) and north (j+1)
                                new_gs = 0.25_c_double * &
                                        ( arr((j-1)*nx + i-1)  &    ! west
                                        + arr((j-1)*nx + i+1)  &    ! east
                                        + arr((j-2)*nx + i  )  &    ! south
                                        + arr(  j  *nx + i  ) )     ! north   
                                
                                ! SOR update: relax by omega
                                arr(idx) = old_val + omega * (new_gs - old_val)
                                max_change = max(max_change, abs(arr(idx) - old_val))
                        end do
                end do
        end subroutine sor_sweep


        ! ===================================================================================
        ! Gauss - Seidel Solver
        ! ===================================================================================
        subroutine solve_gauss_seidel(grid_ptr, nx, ny, tol, max_iter, actual_iter, residual) &
                        bind(C, name="solve_gauss_seidel")
                type(c_ptr),    intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny, max_iter
                real(c_double), intent(in), value :: tol
                integer(c_int), intent(out)   :: actual_iter
                real(c_double), intent(out)   :: residual

                integer :: iter
                real(c_double) :: change

                do iter = 1, max_iter
                        call gauss_seidel_step(grid_ptr, nx, ny, change)
                        if (change < tol) then
                                residual = change
                                actual_iter = iter
                                return
                        end if
                end do
                residual = change
                actual_iter = max_iter
        end subroutine solve_gauss_seidel


        ! ===================================================================================
        ! SOR solver 
        ! ===================================================================================
        subroutine solve_sor(grid_ptr, nx, ny, tol, max_iter, omega,  actual_iter, residual) &
                        bind(C, name="solve_sor")
                type(c_ptr),    intent(in), value :: grid_ptr
                integer(c_int), intent(in), value :: nx, ny, max_iter
                real(c_double), intent(in), value :: tol, omega
                integer(c_int), intent(out)   :: actual_iter
                real(c_double), intent(out)   :: residual

                integer :: iter
                real(c_double) :: change

                do iter = 1, max_iter
                        call sor_sweep(grid_ptr, nx, ny, omega, change)
                        if (change < tol) then
                                residual = change
                                actual_iter = iter
                                return
                        end if
                residual = change
                actual_iter = max_iter
                end do
        end subroutine solve_sor

end module heat_grid
