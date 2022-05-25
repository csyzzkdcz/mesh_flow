#include "../include/MeshFlowProcess.h"

Eigen::MatrixXd meshflow::MeshFlowProcess::isoSurfaceErosion(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, const int numIter, const double energyTol, const double dhat, double dt, int quadOrd)
{
	Eigen::MatrixXd erodedIsoPos = isoPos;

	for (int i = 0; i < numIter; i++)
	{
		// step 1: compute the closest points on the reference surfaces
		Eigen::MatrixXd quadPos;

	}
}