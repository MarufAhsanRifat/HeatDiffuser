module heat_grid_oop
        use, intrinsic :: iso_c_binding, only: c_char, c_loc, c_double, c_int, c_ptr, &
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
        public :: grid_fill_gaussian, grid_set_boundaries_c
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
                if (allocated(self%data)) deallocate(self%data)
                allocate(self%data(nx * ny))
                self%data = 0.0_c_double
        end subroutine grid_init
                
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
                class(grid_t), intent(in), target :: self
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
                class(grid_t), intent(inout) :: self
                real(c_double), intent(in), value :: top, bottom, left, right
                integer :: i, j, idx

                ! Bottom
                j = 1
                do i = 1, self%nx
                        idx = (j - 1) * self%nx + i
                        self%data(idx) = bottom
                end do 

                ! Top
                j = self%ny
                do i = 1, self%nx
                        idx = (j-1) * self%nx + i 
                        self%data(idx) = top
                end do

                ! Left
                i = 1
                do j = 1, self%ny
                        idx = ( j - 1 ) * self%nx + i
                        self%data(idx) = left
                end do         


                ! Right
                i = self%nx
                do j = 1, self%ny
                        idx = (j - 1) * self%nx + i
                        self%data(idx) = right
                end do
        end subroutine grid_set_boundaries

        ! ----- Private Jacovi step -----
        subroutine jacobi_step_internal(src, dst, nx, ny, max_change)
                real(c_double), intent(in) :: src(:)
                real(c_double), intent(inout),  allocatable :: dst(:)
                integer, intent(in), value :: nx, ny
                real(c_double), intent(out) :: max_change
                integer :: i, j, idx
                dst = src

                do j = 2, ny-1
                        do i = 2, nx-1
                                idx = (j-1) * nx + i
                                dst(idx) = 0.25_c_double * &
                                        ( src((j-1)*nx + i-1) &
                                        + src((j-1)*nx + i+1) &
                                        + src((j-2)*nx + i) &
                                        + src( j   *nx + i ) )
                        end do
                end do

                max_change = 0.0_c_double
                do j = 2, ny-1
                        do i = 2, nx-2
                                idx = (j-1)*nx +i
                                max_change = max(max_change, abs(dst(idx)-src(idx)))
                        end do
                end do
        end subroutine jacobi_step_internal

        subroutine grid_solve_jacobi(self, tol, max_iter, actual_iter, residual)
                class(grid_t), intent(inout) :: self
                real(c_double), intent(in) :: tol
                integer(c_int), intent(in) :: max_iter
                integer(c_int), intent(out):: actual_iter
                real(c_double), intent(out):: residual
                real(c_double), allocatable :: tmp(:)
                integer :: iter
                real(c_double) :: change

                allocate(tmp(self%nx * self%ny))
                change = 0.0_c_double
                do iter= 1, max_iter
                        call jacobi_step_internal(self%data, tmp, self%nx, self%ny, change)
                        self%data = tmp

                        if (change < tol) then
                                residual = change
                                actual_iter = iter
                                deallocate(tmp)
                                return
                        end if
                end do

        residual = change
        actual_iter = max_iter
        end subroutine grid_solve_jacobi

        subroutine grid_solve_gauss_seidel(self, tol, max_iter, actual_iter, residual)
                class(grid_t), intent(inout) :: self
                real(c_double), intent(in) :: tol
                integer(c_int), intent(in) :: max_iter
                integer(c_int), intent(out):: actual_iter
                real(c_double), intent(out):: residual
                real(c_double) :: change, old
                integer :: i, j, iter, idx

                do iter=1, max_iter
                        change = 0.0_c_double
                        do j = 2, self%ny - 1
                                do i =2, self%nx - 1
                                        idx = ( j - 1 ) * self%nx + i
                                        old = self%data(idx)
                                        self%data(idx) = 0.25_c_double * ( & 
                                                         self%data((j-1) * self%nx + i - 1 )  & 
                                                       + self%data((j-1) * self%nx + i + 1 )  & 
                                                       + self%data((j-2) * self%nx + i     )  & 
                                                       + self%data((j  ) * self%nx + i     )  ) 
                                         change = max(change, abs(self%data(idx)- old))
                                end do
                        end do 
                         
                if (change < tol) then
                        residual = change
                        actual_iter = iter
                        return 
                end if
                end do

        residual = change
        actual_iter = max_iter
        end subroutine grid_solve_gauss_seidel


        subroutine grid_solve_sor(self, tol, max_iter, omega, actual_iter, residual)
            class(grid_t), intent(inout) :: self
            real(c_double), intent(in) :: tol, omega
            integer(c_int), intent(in) :: max_iter
            integer(c_int), intent(out) :: actual_iter
            real(c_double), intent(out) :: residual
            integer :: i, j, idx, iter
            real(c_double) :: old_val, gs_val, change
    
            do iter = 1, max_iter
                change = 0.0_c_double
                do j = 2, self%ny-1
                    do i = 2, self%nx-1
                        idx = (j-1)*self%nx + i
                        old_val = self%data(idx)
                        gs_val = 0.25_c_double * &
                            ( self%data((j-1)*self%nx + i-1) &
                            + self%data((j-1)*self%nx + i+1) &
                            + self%data((j-2)*self%nx + i  ) &
                            + self%data( j   *self%nx + i  ) )
                        self%data(idx) = old_val + omega * (gs_val - old_val)
                        change = max(change, abs(self%data(idx) - old_val))
                    end do
                end do
                if (change < tol) then
                    residual = change
                    actual_iter = iter
                    return
                end if
            end do
            residual = change
            actual_iter = max_iter
        end subroutine grid_solve_sor
    
        subroutine grid_write_binary(self, filename)
                class(grid_t), intent (in) :: self
                character(len=*), intent(in) :: filename
                integer :: unit, iostat

                open(newunit=unit, file=trim(filename), form='unformatted',&
                        access='stream', action='write', iostat=iostat)
                if (iostat /= 0) then
                        print *, 'Error opening file: ', trim(filename)
                        return
                end if
                write(unit) self%nx, self%ny
                write(unit) self%data(1:self%nx*self%ny)
                close(unit)
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
                integer(c_int), intent(in), value :: nx, ny
                type(c_ptr) :: ptr
                type(grid_t), pointer :: g
                allocate(g)
                call g%init(nx, ny)
                ptr = c_loc(g)
        end function grid_create

        subroutine grid_destroy_c(ptr) bind(C, name="grid_destroy")
                type(c_ptr), intent(in), value :: ptr
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                call g%destroy()
                deallocate(g)
        end subroutine grid_destroy_c

        function grid_get_cptr(ptr) bind(C, name="grid_get_data_ptr") result(data_ptr)
                type(c_ptr), intent(in), value :: ptr
                type(c_ptr) :: data_ptr
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                data_ptr = g%get_data_ptr()
        end function grid_get_cptr

        subroutine grid_fill_gaussian_c(ptr) bind(C, name="grid_fill_gaussian")
                type(c_ptr), intent(in), value :: ptr
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                call g%fill_gaussian()
        end subroutine grid_fill_gaussian_c

        subroutine grid_set_boundaries_c(ptr, top, bottom, left, right) &
                        bind(C, name="grid_set_boundaries")
                type(c_ptr), intent(in), value :: ptr
                real(c_double), intent(in), value :: top, bottom, left, right
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                call g%set_boundaries(top, bottom, left, right)
        end subroutine grid_set_boundaries_c

        function grid_get_element_c(ptr, i, j) bind(C, name="grid_get_element") result(val)
                type(c_ptr), intent(in), value :: ptr
                integer(c_int), intent(in), value :: i, j
                real(c_double) :: val
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                val = g%get_element(i, j)
        end function grid_get_element_c

        subroutine grid_solve_jacobi_c(ptr, tol, max_iter, actual_iter, residual) &
                        bind(C, name="grid_solve_jacobi")
                type(c_ptr), intent(in), value :: ptr
                real(c_double), intent(in), value :: tol
                integer(c_int), intent(in), value :: max_iter
                integer(c_int), intent(out) :: actual_iter
                real(c_double), intent(out) :: residual
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                call g%solve_jacobi(tol, max_iter, actual_iter, residual)
        end subroutine grid_solve_jacobi_c

        subroutine grid_solve_gauss_seidel_c(ptr, tol, max_iter, actual_iter, residual) &
                        bind(C, name="grid_solve_gauss_seidel")
                type(c_ptr), intent(in), value :: ptr
                real(c_double), intent(in), value :: tol
                integer(c_int), intent(in), value :: max_iter
                integer(c_int), intent(out) :: actual_iter
                real(c_double), intent(out) :: residual
                type(grid_t), pointer ::  g
                call c_f_pointer(ptr, g)
                call g%solve_gauss_seidel(tol, max_iter, actual_iter, residual)
        end subroutine grid_solve_gauss_seidel_c


        subroutine grid_solve_sor_c(ptr, tol, max_iter, omega, actual_iter, residual) &
                        bind(C, name="grid_solve_sor")
                type(c_ptr), intent(in), value :: ptr
                real(c_double), intent(in), value :: tol, omega
                integer(c_int), intent(in), value :: max_iter
                integer(c_int), intent(out) :: actual_iter
                real(c_double), intent(out) :: residual
                type(grid_t), pointer :: g
                call c_f_pointer(ptr, g)
                call g%solve_sor(tol, max_iter, omega, actual_iter, residual)
        end subroutine grid_solve_sor_c

        subroutine grid_write_binary_c(ptr, filename) bind(C, name="grid_write_binary")
                type(c_ptr), intent(in), value :: ptr
                character(kind=c_char), intent(in) :: filename(*)
                type(grid_t), pointer :: g
                character(len=256) :: fname
                integer :: k

                fname = ' '
                k = 1
                do while (filename(k) /= c_null_char .and. k <= 256)
                        fname(k:k) = filename(k)
                        k = k+1
                end do 

                call c_f_pointer(ptr, g)
                call g%write_binary(trim(fname))
        end subroutine grid_write_binary_c


end module heat_grid_oop
