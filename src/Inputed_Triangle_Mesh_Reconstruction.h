#pragma once
#include "Reconstruction_Basedon_Paricles_LBFGS.h"
#include "ANN.h"
#include "Point3D.h"
#include <algorithm>
#include <vector>
#include <set>
#include <memory>

//Particle array storage structure
struct alignas(64) ParticleSoA
{
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;

    void resize(int n)
    {
        xs.resize(n, 0.0);
        ys.resize(n, 0.0);
        zs.resize(n, 0.0);
    }

    int size() const { return (int)xs.size(); }

    void fromRawArray(const double* x, int nSamples)
    {
        for (int i = 0; i < nSamples; ++i)
        {
            xs[i] = x[i * 3];
            ys[i] = x[i * 3 + 1];
            zs[i] = x[i * 3 + 2];
        }
    }

    void toRawArray(double* x, int nSamples) const
    {
        for (int i = 0; i < nSamples; ++i)
        {
            x[i * 3] = xs[i];
            x[i * 3 + 1] = ys[i];
            x[i * 3 + 2] = zs[i];
        }
    }
};

//The iterative reusable object pool of evaluate()
struct EvaluateWorkspace
{
    ParticleSoA particles;

    std::vector<Point> pointsForKdTree;

    std::vector<std::vector<std::pair<Point, int>>> closetPointsSet;

    std::vector<std::set<int>> gtClosetTrianglesSet;
    std::vector<int> projectedTriangleIds;

    std::vector<double> grad_x;
    std::vector<double> grad_y;
    std::vector<double> grad_z;

    ANNpointArray annDataPts = nullptr;   // reused across iterations
    ANNkd_tree*   kdTree     = nullptr;   // rebuilt each iteration but buffer is reused

    // Initialize the workspace with the number of samples, allocating memory for the data structures needed to implement the target function.
    void init(int nSamples)
    {
        particles.resize(nSamples);
        pointsForKdTree.resize(nSamples);
        closetPointsSet.resize(nSamples);
        gtClosetTrianglesSet.resize(nSamples);
        projectedTriangleIds.resize(nSamples, -1);
        grad_x.resize(nSamples, 0.0);
        grad_y.resize(nSamples, 0.0);
        grad_z.resize(nSamples, 0.0);

        annDataPts = annAllocPts(nSamples, 3);
    }

    void release()
    {
        if (kdTree)     { delete kdTree;             kdTree     = nullptr; }
        if (annDataPts) { annDeallocPts(annDataPts); annDataPts = nullptr; }
    }

    ~EvaluateWorkspace() { release(); }

    void reset(int nSamples)
    {
        if (kdTree) { delete kdTree; kdTree = nullptr; } 
        for (int i = 0; i < nSamples; ++i)
        {
            closetPointsSet[i].clear();
            gtClosetTrianglesSet[i].clear();
            projectedTriangleIds[i] = -1;
            grad_x[i] = 0.0;
            grad_y[i] = 0.0;
            grad_z[i] = 0.0;
        }
    }
};


class Inputed_Triangle_Mesh_Reconstruction : public Reconstruction_Basedon_Paricles_LBFGS
{
private:
    Distance_OBB m_obbtree;                 //The OBB tree, used to find the closest point on the surface for each sample point.
    CBaseModel& m_inputModel = CBaseModel::CBaseModel();
    CBaseModel& m_sampleModel = CBaseModel::CBaseModel();

    ANNkd_tree*   m_gtKdTree  = nullptr;               //The KdTree built on the ground truth points, used to find the closest points on the ground truth for each sample point
    ANNpointArray m_gtDataPts = nullptr;               //Owned data buffer for the gt KdTree
    vector<Point> m_gtPointsForKdTree;

    EvaluateWorkspace m_workspace;

    int m_iter_count = 0;

public:
    Inputed_Triangle_Mesh_Reconstruction(CBaseModel& inputModel, int nSamples,
        double factor_for_particle_system, double lambdaForBalance, double lambda_PS, double lambda_NA);
    ~Inputed_Triangle_Mesh_Reconstruction();
    vector<Point> GetInitialSamplesOnSurface(int nSamples);

protected:
    lbfgsfloatval_t evaluate(
        lbfgsfloatval_t* x,
        lbfgsfloatval_t* g,
        const int n,
        const lbfgsfloatval_t step
    );

    // Get 6 sample points on the triangle specified, which are used to compute the normal consistency energy and its gradient.
    array<CPoint3D, 6> GetSamplesOnTriangle(int faceIndex);

    double GetSurfaceArea() const;

    vector<std::pair<Point, int>> GetclosetNeighorPoints(ANNkd_tree* kdTree, Point queryPoint, double sigma, const vector<Point>& pointsForKdTree) const;
};
