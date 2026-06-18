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

	double jacobi_step(void* src, void* dst, int nx, int ny);
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

	// 2. Grid Dimentions
	const int nx = 10;
	const int ny = 10;

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
	
	std::cout << "\nInitial source grid:  \n";
	print_grid(src, nx, ny);

	// 5. Perform one Jacobi step
	double max_change = jacobi_step(src, dst, nx, ny);
	std::cout << "\nAfter one jacobi step, max change = " << max_change << "\n";
	std::cout << "Destination grid: \n";
	print_grid(dst, nx, ny);

	// 6. clean up
	destroy_grids();
	return 0;
}
