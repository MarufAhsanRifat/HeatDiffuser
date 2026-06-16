#include <iostream>
#include <cstring>

// Declare the Fortran functions as extern "C"
extern "C" {
	void get_version(char version[], int bufsize);
	int compute_square(int n);
	void* allocate_grid(int nx, int ny);
	void fill_initial(void* grid, int nx, int ny);
	double get_element(void* grid, int i, int j, int nx, int ny);
	void destroy_grid();
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
	const int nx = 10, ny = 10;
	void* grid = allocate_grid( nx, ny );
	fill_initial(grid, nx, ny);

	// print a few values
	std::cout << "Grid snippet:" << std::endl;
	for(int j = 1; j <= 5 ; ++j) {
		for( int i = 1; i <= 5; ++i){
			std::cout << get_element(grid, i, j, nx, ny) << "\t";
		}
		std::cout << std::endl;
	}
	
	destroy_grid();
	return 0;
}
