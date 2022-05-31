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

#include <finitediff.hpp>
#include <ipc/ipc.hpp>
#include <ipc/collision_mesh.hpp>

#include "../include/flow.h"
#include "../include/MeshFlowProcess.h"

bool meshflow::MeshFlowProcess::isoSurfaceFlow(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, std::vector<Eigen::MatrixXd>& isoFlow, const int numIter, const double energyTol, const double dhat, double ipcCoeff, double dt, int quadOrd)
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
	SparseMatrix<double> A_qv;
	gradQ_to_gradV(V, isoFaces, area_0, quadOrd, A_qv);

	int step = 1;
	double energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, NULL);
	cout << "initial energy: " << energy << endl;
    
    MatrixXi E;
    igl::edges(isoFaces, E);
    
    ipc::CollisionMesh mesh(V, E, isoFaces);

    ipc::Constraints constraint_set;
    ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
    
    
    double ipcEnergy = ipcCoeff * ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
    
    cout << "initial ipc energy: " << ipcEnergy << endl;

	while (energy > energyTol)
    {
		MatrixXd grad;
		energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, &grad);
        ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
        ipcEnergy = ipcCoeff * ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
        VectorXd grad_b = ipcCoeff * ipc::compute_barrier_potential_gradient(
                    mesh, V, constraint_set, dhat);
        
        cout << "Flow step " << step << ": \nflow energy: " << energy << ", ipc energy: " << ipcEnergy << endl;
        cout << "flow grad: " << grad.norm() << ", ipc grad: " << grad_b.norm() << endl;
//        SparseMatrix<double> hess_b = ipc::compute_barrier_potential_hessian(mesh, V, constraint_set, dhat);
        
        if(isinf(grad_b.norm()))
            break;
        
        MatrixXd unflattend_grad_b = fd::unflatten(grad_b, V.cols());
        grad += unflattend_grad_b;
        MatrixXd V1 = V - M_inv * grad;
        double maxStep = ipc::compute_collision_free_stepsize(mesh, V, V1);
        
        cout << "max step size: " << maxStep << endl;
    
		
		V = V - maxStep * M_inv * grad;

		step = step + 1;
		isoFlow.push_back(V);
		// Quit after user-specified steps of the flow and return false
		if (step > numIter) return false;
	}
	return true;
}

bool meshflow::MeshFlowProcess::isoSurfaceFlowFullReturn(const Eigen::MatrixXd& isoPos, const Eigen::MatrixXi& isoFaces, std::vector<Eigen::MatrixXd>& isoFlow, const int numIter, const double energyTol, const double dhat, double ipcCoeff, double dt, int quadOrd, std::vector<Eigen::MatrixXd> *flowGrad, std::vector<Eigen::MatrixXd> *ipcGrad)
{

    using namespace Eigen;
    using namespace std;
    using namespace igl;
    isoFlow.clear();
    
    if(flowGrad)
        flowGrad->clear();
    if(ipcGrad)
        ipcGrad->clear();

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
    SparseMatrix<double> A_qv;
    gradQ_to_gradV(V, isoFaces, area_0, quadOrd, A_qv);

    int step = 1;
    double energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, NULL);
    cout << "initial energy: " << energy << endl;
    
    MatrixXi E;
    igl::edges(isoFaces, E);
    ipc::CollisionMesh mesh(V, E, isoFaces);
    ipc::Constraints constraint_set;
    ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
    
    double ipcEnergy = ipcCoeff * ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
    
    cout << "initial ipc energy: " << ipcEnergy << endl;
   

    while (energy > energyTol)
    {
        MatrixXd grad;
        energy = grad_energy(V, isoFaces, _refPos, _refFaces, A_qv, &grad);
        ipc::construct_constraint_set(mesh, V, dhat, constraint_set);
        ipcEnergy = ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
        VectorXd grad_b = ipcCoeff * ipc::compute_barrier_potential_gradient(
                    mesh, V, constraint_set, dhat);
        
        cout << "Flow step " << step << ": \nflow energy: " << energy << ", ipc energy: " << ipcEnergy << endl;
        cout << "flow grad: " << grad.norm() << ", ipc grad: " << grad_b.norm() << endl;
        cout << "check grad size: " << grad.rows() << " " << grad.cols() << " " << grad_b.size() << endl;
        
//        SparseMatrix<double> hess_b = ipc::compute_barrier_potential_hessian(mesh, V, constraint_set, dhat);
        
        if(isinf(grad_b.norm()))
            break;
        
        MatrixXd unflattend_grad_b = fd::unflatten(grad_b, V.cols());
        
        if(flowGrad)
        {
            flowGrad->push_back(grad);
        }
        if(ipcGrad)
        {
            ipcGrad->push_back(unflattend_grad_b);
        }
        
        grad += unflattend_grad_b;
        MatrixXd V1 = V - M_inv * grad;
        double maxStep = ipc::compute_collision_free_stepsize(mesh, V, V1);
        
        cout << "max step size: " << maxStep << endl;
    
        
        V = V - maxStep * M_inv * grad;

        step = step + 1;
        isoFlow.push_back(V);
        // Quit after user-specified steps of the flow and return false
        if (step > numIter) return false;
    }
    return true;
}
