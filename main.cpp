#include <iostream>
#include <cstring>

// Declare the Fortran functions as extern "C"
extern "C" {
	void get_version(char version[], int bufsize);
	int compute_square(int n);
	void* allocate_grid(int nx, int ny);
	void fill_initial(void* grid, int nx, int ny);
	void apply_boundary_conditions(void* grid, int nx, int ny, double top, double bottom, double left, double right);
	double get_element(void* grid, int i, int j, int nx, int ny);
	void destroy_grid();
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
	// version
	constexpr int bufsize = 32;
	char version[bufsize] = {};
	get_version(version, bufsize);
	std::cout << "Fortran version: " << version <<std::endl;

	// square test
	int val = 7;
	int sq = compute_square(val);
	std::cout<< val << " squared = " << sq << std::endl;

	//grid test
	const int nx = 10;
	const int ny = 10;

	void* grid = allocate_grid( nx, ny );
	fill_initial(grid, nx, ny);

	std::cout << "\nInitial grid (Gaussian): \n";
	print_grid(grid, nx, ny);

	// Apply boundaries: hot top, cold elsewhere
	apply_boundary_conditions(grid, nx, ny, 
					1.0, 	// top
					0.0,	// bottom
					0.0,	// left
					0.0);	// right
	
	std::cout << "\nGrid after boundary conditions: \n";
	print_grid(grid, nx, ny);
	
	destroy_grid();
	return 0;
}
