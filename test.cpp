//============================================================================
// Name        : test.cpp
// Author      : Hao Xiong
// Date        : 15/02/26
// Description : Test usage of class CircOverlap
//============================================================================
#include "circ_overlap.h"

int main(int argc, char* argv[]) {
	// Test three circles.
	POINT circle_center[]  = {{.5, .5}, {.3, .3}, {.5, .3}};
	double circle_radius[] = {.2, .2, .2};
	// Create class CircOverlap from the given data above.
	CircOverlap circ_overlap(circle_center, circle_radius, sizeof(circle_radius)/sizeof(double));

	// circ_overlap.set_grid_res(0.001); // Set grid resolution.
	circ_overlap.Calc();
	circ_overlap.PrintArea();
	circ_overlap.WriteGrid();

	return 0;
}
