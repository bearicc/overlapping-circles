//============================================================================
// Name        : circ_overlap.cpp
// Author      : Hao Xiong
// Date        : 15/02/26
//               First written in MATLAB on March 8, 2014
// Description : Given the center coordinates and radii of several circles,
//               calculate the overlapped area, print the result to screen,
//               write the corresponding data to file for visualization.
//
// class CircOverlap(POINT* circle_center, double* circle_radius, int circle_num)
//     circle_center: coordinates of the circles. 1D array of POINT, and struct
//                    POINT contains two double members x and y;
//     circle_radius: radii of the circles. 1D array of double;
//     circle_num   : number of the circles.
//
// Algorithm description:
//     (1) Prepare data structures:
//         grid_   : 2D matrix grid_.matrix is the data for visualization.
//         map_    : Map different overlapping case to a binary code.
//         overlap_: store the mapped code of all different overlapping cases.
//     (2) Divide the grid to cells, init grid_.matrix to -1.
//     (3) For each cell, if it is not -1, continue because it has been visited;
//         A cell can be in one circle while not in another. This overlapping case
//         can be mapped to a binary code "map_result". Its n-th bit value indicates
//         whether this cell in the n-th circle.
//         3.1 For index "map_result" of "map_.table", increase its count by 1, meaning
//             a cell of this type found. This can be used to calculate the overlapped area;
//         3.2 Set the value in grid to the number of overlapped circles, meaning it has
//             been visited and can be used for visualization.
//         3.3 If this overlapped type is found for the first time, append "map_result" to
//             overlap[n], which is a list of overlapping cases with n circles overlapped.
//     (4) For the cell processed in (3), find a proper circle around it, mark their value
//         to the same value of this cell. By marking them as visited we can improve efficiency
//         if there are many circles.
//     (5) Loop over all unvisited cells.
//
//============================================================================
#include "circ_overlap.h"

// Constructor, set the default grid resolution to 1e-3, which can be changed later
//     using set_grid_res().
CircOverlap::CircOverlap(POINT* circle_center, double* circle_radius, int circle_num)
{
    int i = 0, map_table_size = int(pow(2, circle_num));

    lim_[0] = 0, lim_[1] = 1, lim_[2] = 0, lim_[3] = 1;
    circle_center_.resize(circle_num);
    circle_center_.assign(circle_center, circle_center+circle_num);
    circle_radius_.resize(circle_num);
    circle_radius_.assign(circle_radius, circle_radius+circle_num);
    grid_.res      = 1e-3;
    grid_.size     = 0;
    map_.weight.resize(circle_num);
    map_.key.resize(circle_num);
    map_.table.assign(map_table_size, 0);
    for(i = 0; i < circle_num; ++i)
    {
        map_.weight[i] = pow(2, circle_num-i-1); // Bit code abcd ---> a*2^3 + b*2^2 + c*2^1 + d*2^0
    }
    overlap_.resize(circle_num);
    circle_num_    = circle_num;
}

void CircOverlap::Calc()
{
    double x = 0, y = 0, min_distance = .0;
    int map_result = 0, i = 0, j = 0;

    // To simplify, we assume the grid is square.
    grid_.size = int((lim_[1] - lim_[0]) / grid_.res);
    // Adjust the resolution value to make sure (grid size)*(grid resolution) = (grid range).
    grid_.res  = (lim_[1] - lim_[0]) / grid_.size;
    grid_.matrix.resize(grid_.size);
    for(i = 0; i < grid_.size; ++i)
    	grid_.matrix[i].assign(grid_.size, -1);

    y = lim_[2] + grid_.res/2;
    x = lim_[0] + grid_.res/2;
    printf("Begin to calculate (res = %g) ...\n", grid_.res);
    for(j = 0; j < grid_.size; ++j)
    {
        y = lim_[2] + grid_.res*(j+0.5);
        for(i = 0; i < grid_.size; ++i)
        {
            if(grid_.matrix[i][j] < 0) // This cell has not been visited.
            {
            	x = lim_[0] + grid_.res*(i+0.5);
            	min_distance = this->Visit(x, y);
                map_result  = std::inner_product(map_.key.begin(), map_.key.end(), map_.weight.begin(), 0);
            	// this->Update(i, j, map_result, x, y, min_distance);
                // Note: Use function Update() if there are many circles.
				grid_.matrix[i][j] = std::accumulate(map_.key.begin(), map_.key.end(), 0);
				if(!map_.table[map_result])
					overlap_[grid_.matrix[i][j]-1].push_back(map_result);
				++(map_.table[map_result]);
            }
        } //end for i
    } // end for j
}

void CircOverlap::PrintCircle(int map_value)
{
    int i = 0;

    std::fill(map_.key.begin(), map_.key.end(), 0);
    for(i = 1; i <= circle_num_; ++i)
    {
        if(map_value & 1)
        {
        	printf("\tx = %f, y = %f, r = %f\n", \
                    circle_center_[circle_num_-i].x, \
                    circle_center_[circle_num_-i].y, \
                    circle_radius_[circle_num_-i]);
        }
        map_value >>= 1;
    }
}

void CircOverlap::PrintArea()
{
    int count = 0, i = 0, j = 0;

    for(i = circle_num_-1; i > 0; --i)
    {
        for(j = 0; j < overlap_[i].size(); ++j)
        {
        	++count;
            printf("(%d) N = %d\tarea = %f\n", count, i+1, map_.table[overlap_[i][j]]*grid_.res*grid_.res);
            this->PrintCircle(overlap_[i][j]);
        }
    }
}

void CircOverlap::WriteGrid()
{
    int i = 0, j = 0;
    ofstream grid_file;

    grid_file.open("circ_overlap.txt", ios::out);
    for(i = 0; i < grid_.size; ++i)
    {
        for(j = 0; j < grid_.size; ++j)
            grid_file << grid_.matrix[j][i] << " ";
        grid_file << endl;
    }
    grid_file.close();
}

double CircOverlap::Visit(double x, double y)
{
    int i = 0;
    double d = .0, min_distance = 2*(lim_[1]-lim_[0]); // Make sure min_distance is larger than the grid range.

    for(i = 0; i < circle_num_; ++i)
    {
    	d = sqrt(pow(circle_center_[i].x-x, 2) + pow(circle_center_[i].y-y, 2)) - circle_radius_[i];
        if(d > 0)
            map_.key[i] = 0;
        else
        {
            map_.key[i] = 1;
            d = -d;
        }
        if(d < min_distance) min_distance = d;
    }

    return min_distance;
}

void CircOverlap::Update(int i, int j, double map_result, double x, double y, double min_distance)
{
	int ii = 0, jj = 0, step = int(min_distance/grid_.res);
	double xx = .0, yy = .0;

	for(ii = max(0, i-step); ii <= min(i+step, grid_.size-1); ++ii)
		for(jj = max(0, j-step); jj <= min(j+step, grid_.size-1); ++jj)
		{
			xx = lim_[0] + grid_.res*(i+0.5);
			yy = lim_[2] + grid_.res*(j+0.5);
			if((xx-x)*(xx-x)+(yy-y)*(yy-y) < min_distance*min_distance)
			{
				grid_.matrix[ii][jj] = std::accumulate(map_.key.begin(), map_.key.end(), 0);
				if(!map_.table[map_result])
					overlap_[grid_.matrix[ii][jj]-1].push_back(map_result);
				++(map_.table[map_result]);
			}
		}
}
