#include <iostream>
#include <igl/point_mesh_squared_distance.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/normalize_row_lengths.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/doublearea.h>
#include <igl/edges.h>

#include <ipc/ipc.hpp>
#include <ipc/collision_mesh.hpp>

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
    
    Eigen::MatrixXi E;
    igl::edges(isoFaces, E);
    
    ipc::CollisionMesh mesh(V, E, isoFaces);

    ipc::Constraints constraint_set;
    ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
    
    
    double ipcEnergy = ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
    
    std::cout << "initial ipc energy: " << ipcEnergy << std::endl;

	while (energy > energyTol)
    {
		Eigen::MatrixXd grad;
		energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, &grad);
        
        ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
        ipcEnergy = ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
        Eigen::VectorXd grad_b = ipc::compute_barrier_potential_gradient(
                    mesh, V, constraint_set, dhat);
        
        cout << "Flow step " << step << ": \nflow energy: " << energy << ", ipc energy: " << ipcEnergy << endl;
        cout << "flow grad: " << grad.norm() << ", ipc grad: " << grad_b.norm() << std::endl;
        if(isinf(grad_b.norm()))
            break;
        
        Eigen::MatrixXd unflattend_grad_b = V;
        for(int i = 0; i < V.rows(); i++)
        {
            unflattend_grad_b.row(i) = grad_b.segment<3>(3 * i);
        }
        grad += unflattend_grad_b;
		
		V = V - delta_t * M_inv * grad;

		step = step + 1;
		isoFlow.push_back(V);
		// Quit after user-specified steps of the flow and return false
		if (step > numIter) return false;
	}
	return true;
}
