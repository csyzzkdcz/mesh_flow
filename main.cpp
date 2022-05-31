#include "polyscope/polyscope.h"

#include <igl/exact_geodesic.h>
#include <igl/readOBJ.h>

#include "polyscope/messages.h"
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
std::vector<Eigen::MatrixXd> flowGrads;
std::vector<Eigen::MatrixXd> ipcGrads;
std::vector<Eigen::VectorXd> flowGradNorms;
std::vector<Eigen::VectorXd> ipcGradNorms;

meshflow::MeshFlowProcess meshFlow;
double energyTol = 1e-8;
int numIter = 20;
double dt = 1e-3;
double dhat = 1e-2;
double ipcCoeff = 1e3;
int quadOrd = 2;
int curFrame = 0;

bool isShowIsoSurface = true;
bool isShowRefSurface = false;

double globalMinFlow = 0;
double globalMaxFlow = 0;

double globalMinIPC = 0;
double globalMaxIPC = 0;

double vecRatio = 0.01;
bool isRun = false;

void updateView(int frameId)
{
	polyscope::registerSurfaceMesh("reference mesh", refV, refF);
	polyscope::getSurfaceMesh("reference mesh")->setEnabled(isShowRefSurface);

	auto isoMesh = polyscope::registerSurfaceMesh("iso mesh", isoFlow[frameId], isoF);
    isoMesh->setEnabled(isShowIsoSurface);
    
    if(isRun)
    {
        std::cout << "ipc force norm: " << std::endl;
        auto ipcGradNormVis = isoMesh->addVertexScalarQuantity("ipc force norm", ipcGradNorms[frameId]);
        ipcGradNormVis->setMapRange({globalMinIPC, globalMaxIPC});
        ipcGradNormVis->setColorMap("coolwarm");
        
        std::cout << "flow force norm: " << std::endl;
        auto flowGradNormVis = isoMesh->addVertexScalarQuantity("flow force norm", flowGradNorms[frameId]);
        flowGradNormVis->setMapRange({globalMinFlow, globalMaxFlow});
        flowGradNormVis->setColorMap("coolwarm");
        
        std::cout << "ipc force: " << std::endl;
        isoMesh->addVertexVectorQuantity("ipc force", -vecRatio * ipcGrads[frameId], polyscope::VectorType::AMBIENT);
        std::cout << "flow force: " << std::endl;
        isoMesh->addVertexVectorQuantity("flow force", -vecRatio * flowGradNorms[frameId], polyscope::VectorType::AMBIENT);
    }
    
    
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
		if (ImGui::InputDouble("IPC stiffness", &ipcCoeff))
		{
			if (ipcCoeff < 0)
				ipcCoeff = 1e3;
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
        
        if(ImGui::InputDouble("vec ratio", &vecRatio))
        {
            if(vecRatio < 0)
                vecRatio = 0;
            updateView(curFrame);
        }
	}


	// flow
	if (ImGui::Button("Iso Flow"))
	{
		meshFlow = meshflow::MeshFlowProcess(refV, refF);
		meshFlow.isoSurfaceFlowFullReturn(isoV, isoF, isoFlow, numIter, energyTol, dhat, ipcCoeff, dt, quadOrd, &flowGrads, &ipcGrads);
		if (isoFlow.size() < numIter)	// this is just for the visualization convenience
		{
			Eigen::MatrixXd finalIso = isoFlow[isoFlow.size() - 1];
            Eigen::MatrixXd finalFlowGrad = flowGrads[isoFlow.size() - 1];
            Eigen::MatrixXd finalIPCGrad = ipcGrads[isoFlow.size() - 1];
            
			for (int i = isoFlow.size(); i < numIter; i++)
			{
				isoFlow.push_back(finalIso);
                flowGrads.push_back(finalFlowGrad);
                ipcGrads.push_back(finalIPCGrad);
			}
		}
        ipcGradNorms.resize(numIter);
        flowGradNorms.resize(numIter);
        
     
        
        for(int i = 0; i < numIter; i++)
        {
            flowGradNorms[i].resize(flowGrads[i].rows());
            ipcGradNorms[i].resize(ipcGrads[i].rows());
            for(int j = 0; j < flowGrads[i].rows(); j++)
            {
                flowGradNorms[i][j] = flowGrads[i].row(j).norm();
                globalMinFlow = std::min(globalMinFlow, flowGradNorms[i][j]);
                globalMaxFlow = std::max(globalMaxFlow, flowGradNorms[i][j]);
                
                ipcGradNorms[i][j] = ipcGrads[i].row(j).norm();
                globalMinIPC = std::min(globalMinIPC, ipcGradNorms[i][j]);
                globalMaxIPC = std::max(globalMaxIPC, ipcGradNorms[i][j]);
            }
        }
        isRun = true;
		updateView(curFrame);
	}

    if (ImGui::Button("output images", ImVec2(-1, 0)))
    {
        for (int i = 0; i < isoFlow.size(); i++)
        {
            updateView(i);
            //polyscope::options::screenshotExtension = ".jpg";
            std::string name = "output_" + std::to_string(i) + ".jpg";
            polyscope::screenshot(name);
        }
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
    
    isRun = false;

	// Register the mesh with Polyscope
	updateView(curFrame);

	// Add the callback
	polyscope::state::userCallback = callback;

	// Show the gui
	polyscope::show();

	return 0;
}
