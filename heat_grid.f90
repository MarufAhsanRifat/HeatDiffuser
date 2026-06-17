module heat_grid

        use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_int, c_ptr, c_null_ptr, c_f_pointer
        implicit none 
        private

        real(c_double), allocatable, target, save :: grid_data(:)       ! module-level storage

        public :: allocate_grid, fill_initial, get_element, destroy_grid
        public :: apply_boundary_conditions             

contains

        ! Allocate a 1D array of size nx * ny and return a C pointer to it
        function allocate_grid(nx, ny) bind( C, name = "allocate_grid") result(grid_ptr)
                integer(c_int), intent(in), value :: nx, ny
                type(c_ptr) :: grid_ptr
                
                if(allocated(grid_data)) then
                        deallocate(grid_data)
                end if

                allocate(grid_data(nx * ny))         
                grid_data = 0.0_c_double             ! initialize to zero
                grid_ptr = c_loc(grid_data(1))       ! get c pointer to first element
                ! CRUCIAL: we must save a reference to temp so it isn't deallocated. 
                ! We'll store a pointer in a module variable later; for now, we rely
                ! on the caller never losing the pointer.
        end function allocate_grid

        ! Fill the grid with initial conditions (example: a Gaussian bump)
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
! ===========================================================================================


        ! Get a single grid element ( for testing )
        function get_element(grid_ptr, i, j, nx, ny) bind(C, name="get_element") result(val)
                type(c_ptr), intent(in), value :: grid_ptr
                integer(c_int),intent(in), value :: i, j, nx, ny
                real(c_double) :: val
                real(c_double), pointer :: arr(:)

                call c_f_pointer(grid_ptr, arr, [nx*ny])
                val = arr((j-1)*nx + i)
        end function get_element

        ! Free the allocated memory
        subroutine destroy_grid() bind(C, name="destroy_grid")
                if (allocated(grid_data)) then
                        deallocate(grid_data)
                end if
        end subroutine destroy_grid


end module heat_grid
