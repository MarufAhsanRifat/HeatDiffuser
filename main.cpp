#include <iostream>
#include <cmath>
#include <cstring>

exter "C"{
	void* grid_create(int nx, int ny);
	void grid_destroy(void* grid);

	void grid_fill_gaussian(void* grid);
	void grid_set_boundaries(void* grid, double top, double bottom, double left, double right);

	double grid_get_element(void* grid, int i, int j);
		
	void grid_solve_jacobi(void* grid, double tol, int max_iter, int* actual iter, 
			double* residual);
	void grid_solve_gauss_seidel(void* grid, double tol, int max_iter, int* actual iter, 
			double* residual);
	void grid_solve_sor(void* grid, double tol, int max_iter, int* actual iter, 
			double* residual);

	void grid_write_binary(void* grid, char* filename);
}

int main() {
	const int nx = 21, ny = 21;
	const double tol = 1e-6;
	const int max_iter = 5000;

	// Create grid
	void* grid = grid_create(nx, ny)
		
	// Initialize
	grid_fill_gaussian(grid);
	grid_set_boundaries(grid, 1.0, 0.0, 0.0, 0.0);

	// Solve with Jacobi
	int iter;
	double resid;
	grid_solve_jacobi(grid, tol, max_iter, &iter, &resid);
	std::cout<<"Jacobi:	"<< iter <<" iterations, residual= "<< resid << '\n';

	char fname1[] = "temperature_oop_jacobi.bin";
	grid_write_binary(grid, fname1);

	// Reset and test Gauss-Seidel
	grid_fill_gaussian(grid);
	grid_set_boundaries(grid, 1.0, 0.0, 0.0, 0.0);
	grid_solve_gauss_seidel(grid, tol, max_iter, &iter, &resid);
	std::cout<<"Gauss-Seidel:	"<< iter <<" iterations, residual= "<< resid << '\n';

	char fname2[] = "temperature_oop_gs.bin";
	grid_write_binary(grid, fname2);


	// Reset and test SOR
	grid_fill_gaussian(grid);
	grid_set_boundaries(grid, 1.0, 0.0, 0.0, 0.0);
	double pi	 = 3.14159265358979323846;
	double rho	 = cos(pi / nx);
	double omega	 = 2.0 / ( 1.0 + sqrt( 1.0 - rho*rho ));
	grid_solve_sor(grid, tol, max_iter, omega, &iter, &resid);
	std::cout<<"SOR (opt omega):	"<< iter <<" iterations, residual= "<< resid << '\n';

	char fname3[] = "temperature_oop_sor.bin";
	grid_write_binary(grid, fname3);

	// Clean up
	grid_destroy(grid);
	return 0;
}
