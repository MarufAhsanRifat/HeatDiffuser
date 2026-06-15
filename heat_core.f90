module heat_core
  implicit none
  private           ! hide everything by default

  public :: get_version, compute_square  ! expose only these

contains

  ! Fill a C string buffer with the version string      
  subroutine get_version(buffer, bufsize) bind(C, name="get_version")
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    implicit none
    integer(c_int), intent(in), value :: bufsize
    character(kind=c_char), intent(out) :: buffer(*)
    character(len=*), parameter :: str = "HeatDiffuser v0.1.0" // c_null_char
    integer :: i, n

    n = min(len(str), bufsize - 1)
    do i = 1, n
        buffer(i) = str(i:i)
    end do
    buffer(n+1) = c_null_char

  end subroutine get_version

  ! Compute square of an integer (trivial, but demonstrates return value)
  function compute_square(n) bind(C, name="compute_square")
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int), value :: n
    integer(c_int)        :: compute_square
    compute_square = n * n
  end function compute_square

end module heat_core
