#include "Inputed_Triangle_Mesh_Reconstruction.h"

#pragma comment(lib, "ANN.lib")

Inputed_Triangle_Mesh_Reconstruction::Inputed_Triangle_Mesh_Reconstruction(CBaseModel& inputModel, int nSamples,
	double factor_for_particle_system, double lambdaForBalance, double lambda_PS, double lambda_NA)
	:Reconstruction_Basedon_Paricles_LBFGS(nSamples, factor_for_particle_system,
		lambdaForBalance, lambda_PS, lambda_NA), m_obbtree(inputModel), m_inputModel(inputModel)
{
	vector<Point> samples = GetInitialSamplesOnSurface(nSamples);

	for (auto& sample_point : samples)
	{
		m_sampleModel.SetVertex(sample_point.x(), sample_point.y(), sample_point.z());
	}

	for (int i = 0; i < samples.size(); ++i)
	{
		m_variables.emplace_back(samples[i].x());
		m_variables.emplace_back(samples[i].y());
		m_variables.emplace_back(samples[i].z());
	}
	m_areaOfWholeSurface = GetSurfaceArea();
	m_sigma = m_factor_for_particle_system * (sqrt(m_areaOfWholeSurface / m_nSamples));
	m_lambda_NA = lambda_NA;
	m_lambda_PS = lambda_PS;

	//build KdTree for ground truth points
	int gtPointNum = m_inputModel.GetNumOfVerts();
	m_gtPointsForKdTree.reserve(gtPointNum);
	m_gtDataPts = annAllocPts(gtPointNum, 3);
	for (int i = 0; i < gtPointNum; i++)
	{
		const auto& pt = m_inputModel.Vert(i);
		m_gtPointsForKdTree.emplace_back(Point(pt.x, pt.y, pt.z));
		m_gtDataPts[i][0] = pt.x;
		m_gtDataPts[i][1] = pt.y;
		m_gtDataPts[i][2] = pt.z;
	}
	m_gtKdTree = new ANNkd_tree(m_gtDataPts, gtPointNum, 3);

	m_workspace.init(m_nSamples);
}

Inputed_Triangle_Mesh_Reconstruction::~Inputed_Triangle_Mesh_Reconstruction()
{
	if (m_gtKdTree)  { delete m_gtKdTree;              m_gtKdTree  = nullptr; }
	if (m_gtDataPts) { annDeallocPts(m_gtDataPts);     m_gtDataPts = nullptr; }
}

//Initialize the sample points by randomly sampling on the surface
vector<Point> Inputed_Triangle_Mesh_Reconstruction::GetInitialSamplesOnSurface(int nSamples)
{
	const int numFaces = m_inputModel.GetNumOfFaces();
	vector<double> areaList;
	areaList.reserve(numFaces);
	double cumArea = 0;
	for (int i = 0; i < numFaces; ++i)
	{
		cumArea += m_inputModel.GetAreaOfTriangle(i);
		areaList.emplace_back(cumArea);
	}

	vector<Point> samplePoints;
	samplePoints.reserve(nSamples);
	for (int i = 0; i < nSamples; ++i)
	{
		double areaQuery = rand() / (double)RAND_MAX * areaList.back();
		int faceID = lower_bound(areaList.begin(), areaList.end(), areaQuery) - areaList.begin();
		if (faceID >= numFaces)
			faceID -= 1;
		double alpha = rand() / (double)RAND_MAX + 1e-6;
		double beta = rand() / (double)RAND_MAX + 1e-6;
		double gama = rand() / (double)RAND_MAX + 1e-6;
		double sum = alpha + beta + gama;
		alpha /= sum;
		beta /= sum;
		gama /= sum;
		CPoint3D sample = m_inputModel.Vert(m_inputModel.Face(faceID)[0]) * alpha
			+ m_inputModel.Vert(m_inputModel.Face(faceID)[1]) * beta
			+ m_inputModel.Vert(m_inputModel.Face(faceID)[2]) * gama;
		samplePoints.push_back(Point(sample.x, sample.y, sample.z));
	}
	return samplePoints;
}

//Triangle point sampling, three points at the midpoints of the edges, three points inside the triangle
array<CPoint3D, 6> Inputed_Triangle_Mesh_Reconstruction::GetSamplesOnTriangle(int faceIndex)
{
	array<CPoint3D, 6> samplesOnTriangle;
	auto& verts = m_inputModel.m_Verts;
	CBaseModel::CFace face = m_inputModel.Face(faceIndex);
	samplesOnTriangle[0] = verts[face[0]] * 0.5 + verts[face[1]] * 0.5;
	samplesOnTriangle[1] = verts[face[0]] * 0.5 + verts[face[2]] * 0.5;
	samplesOnTriangle[2] = verts[face[1]] * 0.5 + verts[face[2]] * 0.5;
	samplesOnTriangle[3] = verts[face[0]] * (1 / 6.0) + verts[face[1]] * (1 / 6.0) + verts[face[2]] * (2 / 3.0);
	samplesOnTriangle[4] = verts[face[0]] * (1 / 6.0) + verts[face[2]] * (1 / 6.0) + verts[face[1]] * (2 / 3.0);
	samplesOnTriangle[5] = verts[face[1]] * (1 / 6.0) + verts[face[2]] * (1 / 6.0) + verts[face[0]] * (2 / 3.0);

	return samplesOnTriangle;
}

double Inputed_Triangle_Mesh_Reconstruction::GetSurfaceArea() const
{
	double totalArea = 0;
	for (int i = 0; i < m_inputModel.GetNumOfFaces(); ++i)
	{
		CPoint3D vec1 = m_inputModel.Vert(m_inputModel.Face(i)[1]) - m_inputModel.Vert(m_inputModel.Face(i)[0]);
		CPoint3D vec2 = m_inputModel.Vert(m_inputModel.Face(i)[2]) - m_inputModel.Vert(m_inputModel.Face(i)[1]);
		double area = (vec1 * vec2).Len() * 0.5;
		totalArea += area;
	}
	return totalArea;
}

lbfgsfloatval_t Inputed_Triangle_Mesh_Reconstruction::evaluate(
	lbfgsfloatval_t* x,
	lbfgsfloatval_t* g,
	const int n,
	const lbfgsfloatval_t step
)
{
	m_iter_count++;

	m_workspace.reset(m_nSamples);

	ParticleSoA& particles = m_workspace.particles;
	//push points to the mesh 
	#pragma omp parallel  for schedule(static)
	for (int i = 0; i < m_nSamples; ++i)
	{
		CPoint3D closedPoint = m_obbtree.Query(CPoint3D(x[i * 3], x[i * 3 + 1], x[i * 3 + 2])).closestPnt;

		particles.xs[i] = closedPoint.x;
		particles.ys[i] = closedPoint.y;
		particles.zs[i] = closedPoint.z;

		x[i * 3] = closedPoint.x;
		x[i * 3 + 1] = closedPoint.y;
		x[i * 3 + 2] = closedPoint.z;


	}
	// create KdTree
	vector<Point> pointsForKdTree = m_workspace.pointsForKdTree;
	for (int i = 0; i < m_nSamples; ++i)
	{
		pointsForKdTree[i] = Point(particles.xs[i], particles.ys[i], particles.zs[i]);
	}

	// Fill the pre-allocated buffer and rebuild only the tree structure
	ANNpointArray annBuf = m_workspace.annDataPts;
	for (int i = 0; i < m_nSamples; ++i)
	{
		annBuf[i][0] = pointsForKdTree[i].x();
		annBuf[i][1] = pointsForKdTree[i].y();
		annBuf[i][2] = pointsForKdTree[i].z();
	}
	m_workspace.kdTree = new ANNkd_tree(annBuf, m_nSamples, 3);

	auto& closetPointsSet = m_workspace.closetPointsSet;
	auto& gtClosetTrianglesSet = m_workspace.gtClosetTrianglesSet;

	// find closet points on the surface and closet triangles on the ground truth for each sample point
	for (int i = 0; i < m_nSamples; ++i)
	{
		closetPointsSet[i] = GetclosetNeighorPoints(m_workspace.kdTree, pointsForKdTree[i], m_sigma, pointsForKdTree);
		vector<std::pair<Point, int>> gtClosetPoints_i = GetclosetNeighorPoints(m_gtKdTree, pointsForKdTree[i], 0.45 * m_sigma, m_gtPointsForKdTree);

		set<int>& triangles = gtClosetTrianglesSet[i];
		for (auto& closetPoint : gtClosetPoints_i)
		{
			for (int face : m_inputModel.GetFacesFromVertex(closetPoint.second))
			{
				triangles.insert(face);
			}
		}
	}

	auto& grad_x = m_workspace.grad_x;
	auto& grad_y = m_workspace.grad_y;
	auto& grad_z = m_workspace.grad_z;

	lbfgsfloatval_t fx = 0.0;
	lbfgsfloatval_t E_PS = 0.0;
	lbfgsfloatval_t E_NA = 0.0;

	const double inv_4sigma2 = 1.0 / (4 * m_sigma * m_sigma);
	const double inv_2sigma2 = 1.0 / (2 * m_sigma * m_sigma);

	#pragma omp parallel  for reduction(+:E_PS, E_NA) schedule(dynamic)
	for (int i = 0; i < m_nSamples; ++i)
	{
		const double xi = particles.xs[i];
		const double yi = particles.ys[i];
		const double zi = particles.zs[i];

		const int idx = i * 3;

		vector <std::pair<Point, int>> closetPoints = closetPointsSet[i];
		set<int> gtClosetTriangles = gtClosetTrianglesSet[i];
		double gi0 = 0, gi1 = 0, gi2 = 0;

		// Compute the particle system energy and its gradient
		for (int j = 0; j < closetPoints.size(); ++j)
		{
			const double cx = closetPoints[j].first.x();
			const double cy = closetPoints[j].first.y();
			const double cz = closetPoints[j].first.z();
			double dist_x = xi - cx;
			double dist_y = yi - cy;
			double dist_z = zi - cz;
			double  dist2 = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
			double  interParticlesEnergy = 1 / exp((dist2 * inv_4sigma2));
			E_PS += m_lambda_PS * interParticlesEnergy;

			double medialValue = interParticlesEnergy * inv_2sigma2;
			gi0 += m_lambda_PS * medialValue * (cx - xi);
			gi1 += m_lambda_PS * medialValue * (cy - yi);
			gi2 += m_lambda_PS * medialValue * (cz - zi);

		}

		// Compute the normal consistency energy and its gradient
		for (int triangle : gtClosetTriangles)
		{
			const CPoint3D  normal = m_inputModel.GetFaceNormal(triangle);
			const array<CPoint3D, 6> pts = GetSamplesOnTriangle(triangle);
			const double area = m_inputModel.GetAreaOfTriangle(triangle);

			for (int k = 0; k < 6; ++k)
			{
				const double weight = (k < 3) ? (1.0 / 30.0) : (9.0 / 30.0);
				const double coeff_na = m_lambda_NA * weight * area;

				const double dist_x = xi - pts[k].x;
				const double dist_y = yi - pts[k].y;
				const double dist_z = zi - pts[k].z;
				const double innerEnergy = normal.x * dist_x + normal.y * dist_y + normal.z * dist_z;

				E_NA += coeff_na * innerEnergy * innerEnergy;

				const double grad_coeff = 2.0 * coeff_na * innerEnergy;
				gi0 += grad_coeff * normal.x;
				gi1 += grad_coeff * normal.y;
				gi2 += grad_coeff * normal.z;
			}
		}

		grad_x[i] = gi0;
		grad_y[i] = gi1;
		grad_z[i] = gi2;
	}

	#pragma omp parallel  for  schedule(static)
	for (int idx = 0; idx < m_nSamples; ++idx)
	{
		g[idx * 3] = grad_x[idx];
		g[idx * 3 + 1] = grad_y[idx];
		g[idx * 3 + 2] = grad_z[idx];
	}

	eps = E_PS / m_lambda_PS, ena = E_NA / m_lambda_NA;

	//Scale the normal consistency energy to the same level as the particle system energy at the first iteration, and then only update the particle system energy by multiplying it with lambdaForBalance in the following iterations.
	if (m_iter_count == 1)
	{
		m_lambda_NA = eps / (ena * 300);
		E_NA = E_NA * m_lambda_NA;
		fx = E_NA + E_PS;
	}
	if (m_iter_count >= 2)
	{
		fx = E_NA + E_PS;
		m_lambda_PS = m_lambda_PS * m_lambdaForBalance;
	}

	return fx;
}

vector<std::pair<Point, int>> Inputed_Triangle_Mesh_Reconstruction::GetclosetNeighorPoints(ANNkd_tree* kdTree, Point queryPoint,
	double sigma, const vector<Point>& pointsForKdTree)const
{
	ANNpoint queryPt = annAllocPt(3);					// allocate query point
	vector <std::pair<Point, int>> queriedPoint;
	queryPt[0] = queryPoint.x();
	queryPt[1] = queryPoint.y();
	queryPt[2] = queryPoint.z();

	const double sqRadius = 25.0 * sigma * sigma;

	ANNidxArray	tmpIdx = new ANNidx[1];					// near neighbor indices
	ANNdistArray tmpDist = new ANNdist[1];					// near neighbor distances
	int pointNumInK = kdTree->annkFRSearch(						// search
		queryPt,						// query point
		sqRadius,                             //squared radius
		0,								//  number of near neighbors to return
		tmpIdx,							// nearest neighbor array (modified)
		tmpDist,							// dist to near neighbors (modified)
		0);                           // error bound

	delete[] tmpIdx;
	delete[] tmpDist;

	if (pointNumInK > 0)
	{
		ANNidxArray nnIdx = new ANNidx[pointNumInK];						// allocate near neigh indices
		ANNdistArray dists = new ANNdist[pointNumInK];						// allocate near neighbor dists
		kdTree->annkFRSearch(						// search
			queryPt,						// query point
			sqRadius,                             //squared radius
			pointNumInK,					//  number of near neighbors to return
			nnIdx,							// nearest neighbor array (modified)
			dists,							// dist to near neighbors (modified)
			0);                           // error bound

		for (int i = 0; i < pointNumInK; ++i)
		{
			queriedPoint.emplace_back(pointsForKdTree[nnIdx[i]], nnIdx[i]);
		}

		delete[] nnIdx;							// clean things up
		delete[] dists;
	}

	annDeallocPt(queryPt);

	return queriedPoint;
}