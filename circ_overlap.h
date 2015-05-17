#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstdio>

using namespace std;

#ifndef NULL
#define NULL 0
#endif

typedef struct {
	double x;
	double y;
} POINT;

typedef struct
{
    double res;
    int    size;
    vector< vector<int> >  matrix;
} GRID;

typedef struct
{
    vector<int> weight;
    vector<int> key;
    vector<int> table;
} MAP;

class CircOverlap
{
public:
    vector<POINT>  circle_center_; // Coordinate of circle center
    vector<double> circle_radius_; // Radius of circle
    double  lim_[4];        // Grid range: left, right, bottom, top
    int     circle_num_;    // Number of circles

private:
    GRID grid_;             // Grid struct: resolution, size, grid matrix
    // Treat the overlap type as a bit code, its decimal value as position in the hash table
    MAP map_;

    // The index of outer container represents the number of overlapped circles;
    // The inner container is a list of hashed code of overlapped circles.
    vector<vector<int> > overlap_;

public:
    // Construct from the center coordinates and radii of given number of circles.
    CircOverlap(POINT* circle_center, double* circle_radius, int circle_num);
    // Use this to control the calculation precision.
    void set_grid_res(double grid_res) { grid_.res = grid_res; }
    // Set the bit code of an overlap type and store it in class member hash_.key.
    int* isPointInCircle(POINT point);
    double Visit(double x, double y);
    void Update(int i, int j, double hash_result, double x, double y, double min_distance);
    void PrintCircle(int hash_value);
    void PrintArea();
    void WriteGrid();
    void Calc();
};
