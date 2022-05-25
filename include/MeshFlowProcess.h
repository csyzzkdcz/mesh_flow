#pragma once
#include <Eigen/Dense>

namespace meshflow
{
	class MeshFlowProcess
	{
	public:
		MeshFlowProcess() {}
		MeshFlowProcess(const Eigen::MatrixXd& refPos, const Eigen::MatrixXi& refFaces) : _refPos(refPos), _refFaces(refFaces) {}

		Eigen::MatrixXd isoSurfaceErosion(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, const int numIter, const double energyTol, const double dhat, double dt = 1e-3, int quadOrd = 2);
		/*  try to shrink the isosuface, and make sure it close to the reference mesh
		* @paramIn[isoPos]:		the initial vertices postion of the iso surface 
		* @paramIn[isoFaces]:	the initial faces of the iso surface
		* @paramIn[numIter]:	the num of iterations for termination
		* @paramIn[energyTol]:	the energy tolenrance of the termination
		* @paramIn[dhat]:		the IPC tolerance, which is used to control the log barrier sharpness, in our case, the minimum distance between edge-edge. point-edge, point-face. (See IPC paper for details)
		* @paramIn[dt]:			the flowing time step size. 
		* @paramIn[quadOrd]:	the order of quadarature points using to approximate the integration, see "Nasted Cage" paper for details
		* 
		* @return:				the erroded iso-surface, which should satisfy the following properties: (1) it is as close as possible to the reference surface, (2) there is not self-intersection.  
		*/
	private:
		Eigen::Vector3d timeDerivativeOfDisEnergyPerVertex(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, int vertId, int quadOrd = 2); // similar to the eqn (5) of the Nested Cage paper (with sign = 1).
		Eigen::MatrixXd timeDerivativeofDisEnergy(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& is0Faces, int quadId = 2);	

		// ToDo: IPC part.
	private:
		Eigen::MatrixXd _refPos;
		Eigen::MatrixXi _refFaces;
	};
}
