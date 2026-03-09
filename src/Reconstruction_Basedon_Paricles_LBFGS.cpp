#include "Reconstruction_Basedon_Paricles_LBFGS.h"
#include <algorithm>
using namespace std;

Reconstruction_Basedon_Paricles_LBFGS::Reconstruction_Basedon_Paricles_LBFGS(int nSamples,
	double factor_for_particle_system,
	double lambdaForBalance, double lambda_PS, double lambda_NA)
	:m_factor_for_particle_system(factor_for_particle_system),
	m_nSamples(nSamples),
	m_lambdaForBalance(lambdaForBalance),
	m_lambda_PS(lambda_PS),
	m_lambda_NA(lambda_NA)
{
}




