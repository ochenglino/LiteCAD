#pragma once
#include "lbfgs.h"
#include <vector>
#include <iostream>
#include <minmax.h>
#include "fstream"
#include "string"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include "Distance_OBB.h"
#include <algorithm>
#include <filesystem> 

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Triangle_3 Triangle;
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
using namespace std;

// Base class for non-linear optimization using L-BFGS algorithm
class Non_Linear_Optimization
{
protected:
	vector<lbfgsfloatval_t> m_variables;
	int maxIteration;
	clock_t  inital_time;
	lbfgsfloatval_t eps;
	lbfgsfloatval_t ena;
	lbfgsfloatval_t min_eps = 0;
	lbfgsfloatval_t min_ena = 0;
	std::vector<double> prev_variables;

	std::string input_filename;
	std::string output_folder;

	bool early_stopped;

protected:
	std::string getFilenameWithoutExtension(const std::string& filepath) {
		std::filesystem::path p(filepath);
		return p.stem().string();
	}

	virtual int progress(
		const lbfgsfloatval_t* x,
		const lbfgsfloatval_t* g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
	)
	{
		if (k == 1) {
			min_eps = eps;
			min_ena = ena;
			prev_variables = m_variables;
		}
		else {
			if (eps < min_eps) min_eps = eps;
			if (ena < min_ena) min_ena = ena;

			if (eps > min_eps * 1.05) {
				double running_time = ((double)(clock() - inital_time) / CLOCKS_PER_SEC);
				std::string base_name = getFilenameWithoutExtension(input_filename);
				char time_str[50];
				sprintf_s(time_str, "%.3f", running_time);

				// Complete output path with earlystop
				std::string iter_filename = base_name + "_iter" + std::to_string(k - 1) +
					"_time" + std::string(time_str) + "s_nsample" + std::to_string(n / 3) + "_earlystop.xyz";
				std::string iter_path = output_folder + "/" + iter_filename;

				OutResult(iter_path);

				early_stopped = true;

				return 1;
			}
			prev_variables = m_variables;
		}

		/*cerr << " Iteration #" << k << ": ";
		cerr << " Energy = " << fx << ", ";
		cerr << " E_PS = " << eps << ", ";
		cerr << " E_NA = " << ena << ", ";
		cerr << " step = " << step << ", ";
		cerr << " Error = " << gnorm / max(1, xnorm) << endl;
		cerr << "-----------------------------------------\n";*/

		maxIteration = k;

		return 0;
	}

	static int _progress(
		void* instance,
		const lbfgsfloatval_t* x,
		const lbfgsfloatval_t* g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
	)
	{
		return reinterpret_cast<Non_Linear_Optimization*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
	}

	virtual lbfgsfloatval_t evaluate(
		lbfgsfloatval_t* x,
		lbfgsfloatval_t* g,
		const int n,
		const lbfgsfloatval_t step
	) = 0;

	static lbfgsfloatval_t _evaluate(
		void* instance,
		lbfgsfloatval_t* x,
		lbfgsfloatval_t* g,
		const int n,
		const lbfgsfloatval_t step
	)
	{
		return reinterpret_cast<Non_Linear_Optimization*>(instance)->evaluate(x, g, n, step);
	}

public:
	void setInputFilename(const std::string& filename) {
		input_filename = filename;
	}

	void setOutputFolder(const std::string& folder) {
		output_folder = folder;
	}

	// Run the optimization process
	virtual void run(int nSamples, double factor_for_particle_system, double lambdaForBalance)
	{
		lbfgsfloatval_t fx = 0;
		eps = 0, ena = 0;
		early_stopped = false;

		clock_t start = clock();
		inital_time = (clock());
		// Call the L-BFGS optimization algorithm
		int ret = lbfgs(m_variables.size(), &m_variables[0], &fx, _evaluate, _progress, this, NULL);


		double running_time = ((double)(clock() - start) / CLOCKS_PER_SEC);
		cerr << "Total running time = " << running_time << "s" << endl;

		if (early_stopped) {
			return;
		}

		std::string base_name = getFilenameWithoutExtension(input_filename);
		char time_str[50];
		sprintf_s(time_str, "%.3f", running_time);

		// Complete output path with earlystop
		std::string final_filename = base_name + "_iter" + std::to_string(maxIteration) +
			"_time" + std::string(time_str) + "s_nsample" + std::to_string(nSamples) + "_final.xyz";
		std::string final_path = output_folder + "/" + final_filename;

		OutResult(final_path);


	}

	//The result is output as a .xyz point cloud file
	void OutResult(std::string& outFilename)
	{
		std::filesystem::path outPath(outFilename);
		if (outPath.has_parent_path())
		{
			std::filesystem::create_directories(outPath.parent_path());
		}
		ofstream outputModel(outFilename);
		for (vector<double>::iterator it = m_variables.begin(); it != m_variables.end(); it += 3) {
			outputModel << *it << " " << *(it + 1) << " " << *(it + 2) << endl;

		}
		std::cerr << "\n========================================" << std::endl;
		std::cerr << "[Saved] result path: " << outFilename << std::endl;
		std::cerr << "========================================" << std::endl;
		outputModel.close();
	}
};
