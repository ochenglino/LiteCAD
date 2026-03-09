#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include "Inputed_Triangle_Mesh_Reconstruction.h"
using namespace std;

void main(int argc, char** argv)
{
	// default parameters
	int nsample = 2000; // number of samples for reconstruction
	double lambda_PS = 1;
	double lambda_NA = 1;
	double factor_for_balance = 0.95; // decay factor of particle system energy
	double factor_for_particle_system = 0.30; 

	// input, output path
	string modelName = "../../data/wheel.obj";
	string outputFolder = "../../data";


	// cli arguments
	if (argc >= 2) {
		modelName = string(argv[1]);
		cout << "[INFO] model name: " << modelName << endl;
	}

	if (argc >= 3) {
		outputFolder = string(argv[2]);
		cout << "[INFO] output folder: " << outputFolder << endl;
	}

	if (argc >= 4) {
		nsample = atoi(argv[3]);
		cout << "[INFO] nsample: " << nsample << endl;
	}

	if (argc >= 5) {
		factor_for_balance = atof(argv[4]);
		cout << "[INFO] factor for balance: " << factor_for_balance << endl;
	}

	if (argc >= 6) {
		factor_for_particle_system = atof(argv[5]);
		cout << "[INFO] factor for particle system: " << factor_for_particle_system << endl;
	}

	std::cout << "\n[INFO] Starting simplification..." << std::endl;
	try {
		// load model
		if (!std::filesystem::exists(modelName)) {
			throw std::runtime_error("Model is not existed.");
			return;
		}
		CBaseModel inputModel(modelName);
		inputModel.LoadModel();

		// create reconstruction object
		Inputed_Triangle_Mesh_Reconstruction inputed_Triangle_Mesh_Reconstruction(
			inputModel, nsample, factor_for_particle_system,
			factor_for_balance, lambda_PS, lambda_NA
		);

		inputed_Triangle_Mesh_Reconstruction.setInputFilename(modelName);
		inputed_Triangle_Mesh_Reconstruction.setOutputFolder(outputFolder);

		// optimize and output result
		inputed_Triangle_Mesh_Reconstruction.run(
			nsample, factor_for_particle_system,
			factor_for_balance
		);

	}
	catch (const exception& e) {
		cerr << "\n[ERROR] Failed: " << e.what() << endl;
		return;
	}

	if (argc < 2) {
		system("pause");
	}
}
