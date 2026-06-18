module heat_grid

        use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_int, c_ptr, c_null_ptr, c_f_pointer
        implicit none 
        private

        !--- Two grid arrays for double buffering ---
        real(c_double), allocatable, target, save :: grid_src(:), grid_dst(:)       ! module-level storage

        !--- Public Interface ---
        public :: allocate_grid_src, allocate_grid_dst, destroy_grids
        public :: fill_initial,  apply_boundary_conditions, get_element
        public :: jacobi_step       

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
        ! One Jacobi step: update interior from src to dst, return max change 
        ! ------------------------------------------------------------------------------------------
        function jacobi_step(src_ptr, dst_ptr, nx, ny) &
                       bind(C, name="jacobi_step") result(max_change)

        type(c_ptr), intent(in), value :: src_ptr, dst_ptr
        integer(c_int), intent(in), value :: nx, ny
        real(c_double) :: max_change

        real(c_double), pointer :: src(:), dst(:)
        integer :: i, j, idx
        real(c_double) :: stencil_val

        call c_f_pointer(src_ptr, src, [nx*ny])
        call c_f_pointer(dst_ptr, dst, [nx*ny])

        ! Copy the entire source to destination (include boundaries)
        dst = src

        ! Update only the interior (i = 2....nx-1,  j = 2....ny-1)
        do j = 2, ny-1
                do i = 2, nx - 1
                        idx = (j-1) * nx + i
                        stencil_val = 0.25_c_double * &
                                ( src((j-1)*nx + i-1)   & ! west
                                + src((j-1)*nx + i-1)   & ! east
                                + src((j-2)*nx + i  )   & ! south
                                + src((j)  *nx + i  )  )  ! North
                        dst(idx) =stencil_val
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

        end function jacobi_step

end module heat_grid
