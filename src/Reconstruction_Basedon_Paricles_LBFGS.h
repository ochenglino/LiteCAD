#pragma once
#include <omp.h>
#include "Non_Linear_Optimization.h"
#include "Reconstruction_Parameters.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include "Distance_OBB.h"
#include <Eigen/Sparse>
#include <Eigen/LU>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
using namespace std;

// Derived class for reconstruction based on particles using L-BFGS optimization
class Reconstruction_Basedon_Paricles_LBFGS : public Non_Linear_Optimization
{
protected:

	double m_areaOfWholeSurface;

	double m_factor_for_particle_system;

	int m_nSamples;

	double m_lambdaForBalance;

	double m_sigma;

	double m_lambda_PS;

	double m_lambda_NA;

public:
	Reconstruction_Basedon_Paricles_LBFGS(int nSamples, double factor_for_particle_system,
		double lambdaForBalance, double lambda_PS, double lambda_NA);

protected:

	virtual vector<Point> GetInitialSamplesOnSurface(int nSamples) = 0;

	virtual double GetSurfaceArea() const = 0;

};

