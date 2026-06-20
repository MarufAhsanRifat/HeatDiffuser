#include <iostream>
#include <cstring>
#include <cmath>

// Declare the Fortran functions as extern "C"
extern "C" {
	void get_version(char version[], int bufsize);
	int compute_square(int n);

	void* allocate_grid_src(int nx, int ny);
	void* allocate_grid_dst(int nx, int ny);
	void destroy_grids();

	void fill_initial(void* grid, int nx, int ny);
	void apply_boundary_conditions(void* grid, int nx, int ny,
		       	               double top, double bottom,
				       double left, double right);
	double get_element(void* grid, int i, int j, int nx, int ny);

        void solve_poisson(void* src, void* dst, int nx, int ny,
		           double tol, int max_iter,
			   int* actual_iter, double* residual);

        void solve_gauss_seidel(void* src, int nx, int ny,
		           double tol, int max_iter,
			   int* actual_iter, double* residual);

        void solve_sor(void* src, int nx, int ny,
		           double tol, int max_iter, double  omega_opt,
			   int* actual_iter, double* residual);

	void write_grid_binary(void* src, int nx, int ny, char* filename);
}


void setup_initial_condition(void* src, int nx, int ny, double top, 
	       	             double bottom, double left, double right) {
	// Set up the initial field with Dirichlet boundaries
	fill_initial(src, nx, ny);
	apply_boundary_conditions(src, nx, ny, top, bottom, left, right );
}


// Print the grid with top row first ( j = ny down to 1 )
static void print_grid(void* grid, int nx, int ny) {
	for (int j=ny; j >= 1; --j){
		for(int i = 1; i <= nx; ++i){
			std::cout << get_element(grid, i, j, nx, ny) << "\t";
		}
		std::cout << '\n';
	}
}

int main() {
	// 1. Version 
	constexpr int bufsize = 32;
	char version[bufsize] = {};
	get_version(version, bufsize);
	std::cout << "Fortran version: " << version <<std::endl;

	// square test
	int val = 7;
	int sq = compute_square(val);
	std::cout<< val << " squared = " << sq << std::endl;

	// 2. Problem size and parameters 
	const int nx = 21, ny = 21;
	const double tol  = 1.0e-6;
	const int max_iter = 5000;

	// 3. Allocate both grids
	void* src = allocate_grid_src( nx, ny );
	void* dst = allocate_grid_dst( nx, ny );

	// 4. Setting up the initial conditions
	double top = 1.0; double bottom = 0.0;
	double left = 0.0; double right = 0.0;
	setup_initial_condition(src, nx, ny, top, bottom, left, right);
	std::cout << "\nInitial source grid(top row printed first):  \n";
	print_grid(src, nx, ny);

	int actual_iter = 0;
	double residual = 0.0;

	// 5.0 Solve Poisson equation with JACOBI solver
	solve_poisson(src, dst, nx, ny, tol, max_iter, &actual_iter, &residual);

	std::cout << "\n ================ JACOBI Solver finished ================== \n";
	std::cout << "  Iterations: " << actual_iter << '\n';
	std::cout << "  Residual  : " << std::scientific << residual << '\n';

	// 5.1  Display final solution (guaranteed in src)
	std::cout << "\nFinal steady-state temperature: \n";
	print_grid(src, nx, ny);
	
	// 5.2 Write the final grid to binary file
	char filename[]= "temperature_jacobi.bin";
	write_grid_binary(src, nx, ny, filename);
	std::cout << "Final grid written to temperature_jacobi.bin\n"; 

	
	// 6.0 Solve Poisson equation with Gauss-Seidel

	setup_initial_condition(src, nx, ny, top, bottom, left, right);
	solve_gauss_seidel(src, nx, ny, tol, max_iter, &actual_iter, &residual);

	std::cout << "\n ================ GAUSS-SEIDEL Solver finished ================== \n";
	std::cout << "  Iterations: " << actual_iter << '\n';
	std::cout << "  Residual  : " << std::scientific << residual << '\n';

	// 6.1  Display final solution (guaranteed in src)
	std::cout << "\nFinal steady-state temperature: \n";
	print_grid(src, nx, ny);

	
	// 6.2 Write the final grid to binary file
        char filename2[] =  "temperature_gauss_seidel.bin";
	write_grid_binary(src, nx, ny, filename2);
	std::cout << "Final grid written to temperature_gauss_seidel.bin\n"; 
	
	// 7.0 Solve Poisson equation with SOR (optimal omega for laplace)
	double pi = 3.14159265358979323846;
	double rho_jacobi = cos(pi/nx);   // approximate spectral radius for square grid
	double omega_opt = 2.0/ (1.0 + sqrt(1.0 - rho_jacobi * rho_jacobi));

	setup_initial_condition(src, nx, ny, top, bottom, left, right);
	solve_sor(src, nx, ny, tol, max_iter, omega_opt, &actual_iter, &residual);

	std::cout << "\n ================ SOR Solver finished ================== \n";
	std::cout<< " Optimal omega = " << omega_opt << '\n';
	std::cout << "  Iterations: " << actual_iter << '\n';
	std::cout << "  Residual  : " << std::scientific << residual << '\n';

	// 7.1  Display final solution (guaranteed in src)
	std::cout << "\nFinal steady-state temperature: \n";
	print_grid(src, nx, ny);
	
	// 7.2 Write the final grid to binary file
        char filename3[]  = "temperature_sor.bin";
	write_grid_binary(src, nx, ny, filename3);
	std::cout << "Final grid written to temperature_sor.bin\n"; 
	
	// 8. clean up
	destroy_grids();
	return 0;
}
