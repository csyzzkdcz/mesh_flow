#include <iostream>
#include <igl/signed_distance.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/normalize_row_lengths.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include "../include/flow.h"
#include "../include/MeshFlowProcess.h"

bool meshflow::MeshFlowProcess::isoSurfaceFlow(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, std::vector<Eigen::MatrixXd>& isoFlow, const int numIter, const double energyTol, const double dhat, double dt, int quadOrd)
{

	using namespace Eigen;
	using namespace std;
	using namespace igl;
	isoFlow.clear();

	MatrixXd V = isoPos;
	isoFlow.push_back(V);


	// while there are intersections or some negative winding number, keep flowing
	SparseMatrix<double> M;
	massmatrix(V, isoFaces, MASSMATRIX_TYPE_BARYCENTRIC, M);
	SparseMatrix<double> M_inv;
	invert_diag(M, M_inv);


	// calculate diameter of the meshes to scale step size
	double diam = diameter(isoPos, _refPos);
	double delta_t = diam * dt;

	VectorXd area_0;
	doublearea(V, isoFaces, area_0);
	area_0 = 0.5 * area_0;
	Eigen::SparseMatrix<double> A_qv;
	gradQ_to_gradV(V, isoFaces, area_0, quadOrd, A_qv);

	int step = 1;
	double energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, NULL);
	std::cout << "initial energy: " << energy << std::endl;

	while (energy > energyTol) {
		cout << "Flow step " << step << ": " << energy << endl;
		Eigen::MatrixXd grad;
		energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, &grad);
		
		V = V - delta_t * M_inv * grad;

		step = step + 1;
		isoFlow.push_back(V);
		// Quit after user-specified steps of the flow and return false
		if (step > numIter) return false;
	}
	return true;
}