#include <iostream>
#include <cstring>

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

	void write_grid_binary(void* src, int nx, int ny, char* filename);
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

	// 4. Set up the initial field with Dirichlet boundaries
	fill_initial(src, nx, ny);
	apply_boundary_conditions(src, nx, ny, 
					1.0, 	// top
					0.0,	// bottom
					0.0,	// left
					0.0);	// right
	
	std::cout << "\nInitial source grid(top row printed first):  \n";
	print_grid(src, nx, ny);

	// 5. Solve Poisson equation
	int actual_iter = 0;
	double residual = 0.0;
	solve_poisson(src, dst, ny, ny, tol, max_iter, &actual_iter, &residual);

	std::cout << "\nSolver finished: \n";
	std::cout << "  Iterations: " << actual_iter << '\n';
	std::cout << "  Residual  : " << std::scientific << residual << '\n';

	// 6. Display final solution (guaranteed in src)
	std::cout << "\nFinal steady-state temperature: \n";
	print_grid(src, nx, ny);
	
	// 7. Write the final grid to binary file
	char filename[]= "temperature.bin";
	write_grid_binary(src, nx, ny, filename);
	std::cout << "Final grid written to temperature.bin\n"; 
	
	// 8. clean up
	destroy_grids();
	return 0;
}
