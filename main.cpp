#include "polyscope/polyscope.h"

#include <igl/PI.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/exact_geodesic.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/lscm.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>

#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <unordered_set>
#include <utility>

#include "include/MeshFlowProcess.h"

// The mesh, Eigen representation
Eigen::MatrixXd refV;
Eigen::MatrixXi refF;

Eigen::MatrixXd isoV;
Eigen::MatrixXi isoF;

std::vector<Eigen::MatrixXd> isoFlow;

meshflow::MeshFlowProcess meshFlow;
double energyTol = 1e-8;
int numIter = 100;
double dt = 1e-3;
double dhat = 1e-2;
int quadOrd = 2;
int curFrame = 0;

bool isShowIsoSurface = true;
bool isShowRefSurface = true;

void updateView(int frameId)
{
	polyscope::registerSurfaceMesh("reference mesh", refV, refF);
	polyscope::getSurfaceMesh("reference mesh")->setEnabled(isShowRefSurface);

	polyscope::registerSurfaceMesh("iso mesh", isoFlow[frameId], isoF);
	polyscope::getSurfaceMesh("iso mesh")->setEnabled(isShowIsoSurface);

}

void callback() 
{

	static int numPoints = 2000;
	static float param = 3.14;

	ImGui::PushItemWidth(100);

	if (ImGui::CollapsingHeader("Flow Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ImGui::InputDouble("energy tol", &energyTol))
		{
			if (energyTol <= 0)
				energyTol = 1e-8;
		}
		if (ImGui::InputDouble("IPC dhat", &dhat))
		{
			if (dhat < 0)
				dhat = 1e-2;
		}
		if (ImGui::InputDouble("step size", &dt))
		{
			if (dt < 0)
				dt = 1e-3;
		}
		if (ImGui::InputInt("num of iter", &numIter))
		{
			if (numIter < 0)
				numIter = 100;
		}

	}
	if (ImGui::CollapsingHeader("visualization Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ImGui::Checkbox("Iso Surface", &isShowIsoSurface))
		{
			updateView(curFrame);
		}
		if (ImGui::Checkbox("Ref Surface", &isShowRefSurface))
		{
			updateView(curFrame);
		}
		if (ImGui::SliderInt("current frame", &curFrame, 0, numIter - 1))
		{
			updateView(curFrame);
		}

		if (ImGui::DragInt("current frame", &curFrame, 1, 0, numIter - 1))
		{
			curFrame = curFrame % numIter;
			updateView(curFrame);
		}
	}


	// flow
	if (ImGui::Button("Iso Flow"))
	{
		meshFlow = meshflow::MeshFlowProcess(refV, refF);
		meshFlow.isoSurfaceFlow(isoV, isoF, isoFlow, numIter, energyTol, dhat, dt, quadOrd);
		if (isoFlow.size() < numIter)	// this is just for the visualization convenience
		{
			Eigen::MatrixXd finalIso = isoFlow[isoFlow.size() - 1];
			for (int i = isoFlow.size(); i < numIter; i++)
			{
				isoFlow.push_back(finalIso);
			}
		}
		updateView(curFrame);
	}

	ImGui::PopItemWidth();
}

int main(int argc, char** argv) {
	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();
	if (argc < 3)
	{
		std::cerr << "no enough inputs!" << std::endl;
		std::cerr << "Correct usage: ./mesh-flow_bin [reference mesh path] [iso mesh path]" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string filename = argv[1];
	std::cout << "loading: " << filename << std::endl;

	std::string filePath = filename;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows

	// Read the mesh
	if (!igl::readOBJ(filePath, refV, refF))
	{
		std::cerr << "missing the reference mesh file!" << std::endl;
		exit(EXIT_FAILURE);
	}

	filename = argv[2];
	std::cout << "loading: " << filename << std::endl;
	filePath = filename;
	std::replace(filePath.begin(), filePath.end(), '\\', '/'); // handle the backslash issue for windows
	// Read the mesh
	if (!igl::readOBJ(filePath, isoV, isoF))
	{
		std::cerr << "missing the iso mesh file!" << std::endl;
		exit(EXIT_FAILURE);
	}

	curFrame = 0;
	isoFlow.resize(numIter, isoV);

	// Register the mesh with Polyscope
	updateView(curFrame);

	// Add the callback
	polyscope::state::userCallback = callback;

	// Show the gui
	polyscope::show();

	return 0;
}
