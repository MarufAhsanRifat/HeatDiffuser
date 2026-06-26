module heat_grid_oop
        use, intrinsic :: iso_c_binding, only: c_loc, c_double, c_int, c_ptr, &
                c_f_pointer, c_null_char, c_null_ptr

        implicit none
        private

        ! ----- Pulic type and factory -----
        public :: grid_t

        type :: grid_t
                private
                integer(c_int) :: nx = 0, ny = 0
                real(c_double), allocatable :: data(:)
        contains
                ! Initialization and destruction
                procedure :: init               => grid_init
                procedure :: destroy            => grid_destroy

                ! Data access
                procedure :: get_element        => grid_get_element
                procedure :: get_data_ptr       => grid_get_data_ptr  ! for C++ interop

                ! Physics
                procedure :: fill_gaussian      => grid_fill_gaussian
                procedure :: set_boundaries     => grid_set_boundaries

                ! Solvers
                procedure :: solve_jacobi       => grid_solve_jacobi
                procedure :: solve_gauss_seidel => grid_solve_gauss_seidel
                procedure :: solve_sor          => grid_solve_sor

                ! I/O
                procedure :: write_binary       => grid_write_binary
        end type grid_t

        ! ----- C interop: factory and wrapper functions -----
        public :: grid_create, grid_destroy_c, grid_get_cptr
        public :: grid_fill_gaussian, grid_ser_boundaries_c
        public :: grid_get_element_c
        public :: grid_solve_jacobi_c, grid_solve_gauss_seidel, grid_solve_sor_c
        public :: grid_write_binary_c

contains

        ! ==========================================================================
        ! TYPE-BOUND PROCEDRUES
        ! ==========================================================================

        subroutine grid_init(self, nx, ny)
                class(grid_t), intent(inout) :: self
                integer(c_int), intent(in) :: nx, ny
                self%nx = nx
                self%ny = ny
                if (allocated(self%data)) deallocated(self%data)
                allocate(self%data(nx * ny))
                self%data = 0.0_c_double
        end subroutine grid_init(self, nx, ny)
                
        subroutine grid_destroy(self)
                class(grid_t), intent(inout) :: self
                if (allocated(self%data)) deallocate(self%data)
                self%nx = 0
                self%ny = 0
        end subroutine grid_destroy
        
        function grid_get_element(self, i, j) result(val)
                class(grid_t), intent(in) :: self
                integer(c_int), intent(in) :: i, j
                real(c_double) :: val
                val = self%data((j-1)* self%nx + i)
        end function grid_get_element

        function grid_get_data_ptr(self) result(ptr)
                class(grid_t), intent(in) :: self
                type(c_ptr) :: ptr
                if (allocated(self%data)) then
                        ptr = c_loc(self%data(1))
                else
                        ptr = c_null_ptr
                end if
        end function grid_get_data_ptr

        subroutine grid_fill_gaussian(self)
                class(grid_t), intent(inout) :: self
                integer :: i, j, idx
                real(c_double) :: x, y, center_x, center_y, sigma

                center_x = (self%nx - 1) / 2.0_c_double
                center_y = (self%ny - 1) / 2.0_c_double
                sigma    = min(self%nx, self%ny) / 6.0_c_double

                do j = 1, self%ny
                        do i = 1, self%nx
                                x = real(i-1, c_double)
                                y = real(j-1, c_double)
                                idx = (j-1) * self%nx + i
                                self%data(idx) = exp( -((x-center_x)**2 + (y-center_y)**2) / &
                                       ( 2.0_c_double * sigma **2) )
                        end do
                end do
        end subroutine grid_fill_gaussian

        subroutine grid_set_boundaries(self, top, bottom, left, right)
        end subroutine grid_set_boundaries

        ! ----- Private Jacovi step -----
        subroutine jacobi_step_internal(src, dst, nx, ny, max_change)
        end subroutine jacobi_step_internal

        subroutine grid_solve_jacobi(self, tol, max_iter, actual_iter, residual)
        end subroutine grid_solve_jacobi

        subroutine grid_solve_gauss_seidel(self, tol, max_iter, actual_iter, residual)
        end subroutine grid_solve_gauss_seidel


        subroutine grid_solve_sor(self, tol, max_iter, actual_iter, residual)
        end subroutine grid_solve_sor

        subroutine grid_write_binary(self, filename)
        end subroutine grid_write_binary

        ! ==============================================================================
        ! C INTEROPERABILITY WRAPPERS (bind(c) functions)
        ! These are thin wrappers that call type-bound procedures via a stored
        ! pointer. We store grid_t instance in a simple array indexed by a 
        ! handle, or we return a c_ptr to the grid_t itself
        !
        ! For simplicity, we return a c_ptr to the grid_t object, allocated on
        ! the Fortran heap. C++ holds this as void*. 
        ! ==============================================================================

        function grid_create(nx, ny) bind(C, name='grid_create') result(ptr)
        end function grid_create

        subroutine grid_destroy_c(ptr) bind(C, name="grid_destroy")
        end subroutine grid_destroy_c

        function grid_get_cptr(ptr) bind(C, name="grid_get_data_ptr" result(data_ptr)
        end function grid_get_cptr

        subroutine grid_fill_gaussian_c(ptr) bind(C, name="grid_fill_gaussian")
        end subroutine grid_fill_gaussian_c

        subroutine grid_set_boundaries_c(ptr, top, bottom, left, right) &
                        bind(C, name="grid_set_boundaries")
        end subroutine grid_set_boundaries_c

        function grid_get_element_c(ptr, i, j) bind(C, name="grid_get_element") result(val)
        end function grid_get_element_c

        subroutine grid_solve_jacobi_c(self, tol, max_iter, actual_iter, residual) &
                        bind(C, name="grid_solve_jacobi")
        end subroutine grid_solve_jacobi_c

        subroutine grid_solve_gauss_seidel_c(self, tol, max_iter, actual_iter, residual) &
                        bind(C, name="grid_solve_gauss_seidel")
        end subroutine grid_solve_gauss_seidel_c


        subroutine grid_solve_sor_c(self, tol, max_iter, actual_iter, residual) &
                        bind(C, name="grid_solve_sor")
        end subroutine grid_solve_sor_c

        subroutine grid_write_binary_c(self, filename) bind(C, name="grid_write_binary")
        end subroutine grid_write_binary_c

end module heat_grid_oop
