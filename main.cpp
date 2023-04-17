#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/arap.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

// Body matrices
Eigen::MatrixXd V_body;
Eigen::MatrixXi F_body;

// Garments matrices
Eigen::MatrixXd V_garment_front; // vertices of the mesh
Eigen::MatrixXi F_garment_front; // faces of the mesh

Eigen::MatrixXd V_garment_back; // vertices of the mesh
Eigen::MatrixXi F_garment_back; // faces of the mesh

// Deformed garments
Eigen::MatrixXd V_garment_front_deformed;
Eigen::MatrixXd V_garment_back_deformed;


igl::ARAPData arap_data;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){}



// ------------ main program ----------------
int main(int argc, char *argv[]) {

    // // Read the body mesh from the OBJ file
	//igl::readOBJ("../data/FinalBaseMesh.obj",V_body, F_body);
	igl::readOFF("../data/FinalBaseMesh.off",V_body, F_body);
	//igl::readPLY("../models/man_15k.ply",V_body, F_body);
	// // read the 2D pattern mesh from file
  	igl::read_triangle_mesh("../data/garment.obj", V_garment_front, F_garment_front);
	igl::read_triangle_mesh("../data/garment.obj", V_garment_back, F_garment_back);
	
	Vector3d t(0.0, 0.0, -4.0);
	V_garment_back.rowwise() += t.transpose();


	igl::opengl::glfw::Viewer viewer;

	viewer.data_list.resize(3);
	viewer.data_list.at(0).set_mesh(V_body, F_body);

    viewer.data_list.at(1).set_mesh(V_garment_front, F_garment_front);
	viewer.data_list.at(2).set_mesh(V_garment_back, F_garment_back);
	

  viewer.launch();
}
