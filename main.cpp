#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/arap.h>
#include <igl/writeDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/is_edge_manifold.h>


#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace Eigen;
using namespace std;

typedef 
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
const RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
const RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);


Eigen::MatrixXd V,U_garment_back;
Eigen::MatrixXd U_sym_garment_back;


Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
igl::ARAPData arap_data;


// Body matrices
Eigen::MatrixXd V_body;
Eigen::MatrixXi F_body;

// Garments matrices
Eigen::MatrixXd V_garment_front; // vertices of the mesh
Eigen::MatrixXi F_garment_front; // faces of the mesh

Eigen::MatrixXd V_garment_back; // vertices of the mesh
Eigen::MatrixXi F_garment_back; // faces of the mesh


MatrixXd V_sym_garment_back;
MatrixXi F_sym_garment_back;

// Deformed garments
Eigen::MatrixXd V_garment_front_deformed;
Eigen::MatrixXd V_garment_back_deformed;

bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{

    MatrixXd bc(b.size(),V_garment_back.cols());
    for(int i = 0;i<b.size();i++)
    {
      bc.row(i) = V_garment_back.row(b(i));

      switch(S(b(i)))
      {
        case 0:
        {
		  if(bc(i, 0) > 0) bc(i,0) = 0;
          break;
        }
        case 1:
        {
          bc(i,2) += 0.5*sin(0.5*anim_t*2.*igl::PI);
          break;
        }
        case 2:
        {
          bc(i,2) += 0.8*cos(igl::PI+0.35*anim_t*2.*igl::PI);
          break;
        }
        default:
          break;
      }
    }
    igl::arap_solve(bc,arap_data,U_garment_back);
    viewer.data_list.at(2).set_vertices(U_garment_back);
    viewer.data_list.at(2).compute_normals();

	U_sym_garment_back = U_garment_back;
	for(int i =0; i < V_sym_garment_back.rows(); i ++)
	{
		U_sym_garment_back.row(i)(0) = - U_sym_garment_back.row(i)(0);
	}
	viewer.data_list.at(3).set_vertices(U_sym_garment_back);
    viewer.data_list.at(3).compute_normals();

	if(viewer.core().is_animating)
	{
	anim_t += anim_t_dir;
	}
	return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      return true;
  }
  return false;
}



int main(int argc, char *argv[])
{
	// read the 3D mesh from input file
	igl::readOFF("../data/FinalBaseMesh.off",V_body, F_body);

	// read the 2D pattern mesh from file
	igl::readOFF("../data/plane.off", V_garment_back, F_garment_back);
	igl::readOFF("../data/plane.off", V_garment_front, F_garment_front);
	//igl::read_triangle_mesh("../data/half_garment.off", V_garment_front, F_garment_front);
	//igl::read_triangle_mesh("../data/half_garment.off", V_garment_back, F_garment_back);
	
	cout << "The input body mesh has " << V_body.rows() << " vertices and " << F_body.rows() << " faces." << endl;
	cout << "The garment pattern has " << V_garment_back.rows() << " vertices and " << F_garment_back.rows() << " faces." << endl;

	igl::readDMAT("../data/plane.dmat", S);

	Vector3d t_front(0.0, 10.0, -4.0);
	V_garment_front.rowwise() += t_front.transpose();


	Vector3d t_back(0.0, 10.0, 2.0);
	V_garment_back.rowwise() += t_back.transpose();

	
	// vertices in selection
	igl::colon<int>(0,V_garment_back.rows()-1,b);
	b.conservativeResize(stable_partition( b.data(), b.data()+b.size(), 
	[](int i)->bool{return S(i)>=0;})-b.data());

	cout << "The number of constrained vectors is" << b.rows() << endl;

	// Centroid
	mid = 0.5*(V_garment_back.colwise().maxCoeff() + V_garment_back.colwise().minCoeff());

	cout << "The centroid of the garment is " << mid << endl;

	// Precomputation
	arap_data.max_iter = 100;
   	igl::arap_precomputation(V_garment_back,F_garment_back,V_garment_back.cols(),b,arap_data);
	igl::arap_precomputation(V_garment_front,F_garment_front,V_garment_front.cols(),b,arap_data);

	cout << "ARAP precomputation done" << endl;

	// Set color based on selection
	MatrixXd C(F_garment_back.rows(),3);

	for(int f = 0;f<F_garment_back.rows();f++)
	{
		if( S(F_garment_back(f,0))>=0 && S(F_garment_back(f,1))>=0 && S(F_garment_back(f,2))>=0)
		{
		C.row(f) = sea_green;
		}else
		{
		C.row(f) = gold;
		}
	}

	// Set the mirror plane axis
	float x_offset = V_garment_back.col(0).maxCoeff();
	Vector3d t_offset(-x_offset, 0.0, 0.0);
	V_garment_back.rowwise() += t_offset.transpose();
	V_garment_front.rowwise() += t_offset.transpose();
	

	MatrixXd V_sym_garment_front;
	MatrixXi F_sym_garment_front;
	V_sym_garment_front = V_garment_front;
	F_sym_garment_front = F_garment_front;
	for(int i =0; i < V_sym_garment_front.rows(); i ++)
	{
		V_sym_garment_front.row(i)(0) = - V_sym_garment_front.row(i)(0);
	}
	for(int i =0; i < F_sym_garment_front.rows(); i ++)
	{
		F_garment_front.row(i)(0) = F_sym_garment_front.row(i)(2);
		F_garment_front.row(i)(2) = F_sym_garment_front.row(i)(0);
	}


	V_sym_garment_back = V_garment_back;
	F_sym_garment_back = F_garment_back;
	for(int i =0; i < V_sym_garment_back.rows(); i ++)
	{
		V_sym_garment_back.row(i)(0) = - V_sym_garment_back.row(i)(0);
	}
	for(int i =0; i < F_sym_garment_back.rows(); i ++)
	{
		F_sym_garment_back.row(i)(0) = F_garment_back.row(i)(2);
		F_sym_garment_back.row(i)(2) = F_garment_back.row(i)(0);
	}

	U_garment_back = V_garment_back;
	U_sym_garment_back = V_sym_garment_back;

	cout << "Symetry pattern generation done" << endl;

	// Plot the mesh with pseudocolors
	igl::opengl::glfw::Viewer viewer;

  	viewer.data_list.resize(5);

	viewer.data_list.at(0).set_mesh(V_body, F_body);
	viewer.data_list.at(0).set_colors(purple);

    viewer.data_list.at(1).set_mesh(V_garment_front, F_garment_front);
	viewer.data_list.at(1).set_colors(C);
	
	viewer.data_list.at(2).set_mesh(U_garment_back, F_garment_back);
	viewer.data_list.at(2).set_colors(C);

	viewer.data_list.at(3).set_mesh(U_sym_garment_back, F_sym_garment_back);
	viewer.data_list.at(3).set_colors(C);

	viewer.data_list.at(4).set_mesh(V_sym_garment_front, F_sym_garment_front);
	viewer.data_list.at(4).set_colors(C);


	viewer.callback_pre_draw = &pre_draw;

	viewer.callback_key_down = &key_down;
	viewer.core().is_animating = false;
	viewer.core().animation_max_fps = 30.;
	cout<<
		"Press [space] to toggle animation"<<endl;
	viewer.launch();


	// Plot the mesh with pseudocolors
	igl::opengl::glfw::Viewer viewer_2D;

	viewer_2D.data_list.resize(4);

	viewer_2D.data_list.at(0).set_mesh(V_garment_back, F_garment_back);
	viewer_2D.data_list.at(0).set_colors(C);

	viewer_2D.data_list.at(1).set_mesh(V_sym_garment_back, F_sym_garment_back);
	viewer_2D.data_list.at(1).set_colors(C);

	Vector3d t_bottom(0.0, -5.0, 6.0);
	V_garment_front.rowwise() += t_bottom.transpose();
	V_sym_garment_front.rowwise() += t_bottom.transpose();
	
	viewer_2D.data_list.at(2).set_mesh(V_garment_front, F_garment_front);
	viewer_2D.data_list.at(2).set_colors(C);

	viewer_2D.data_list.at(3).set_mesh(V_sym_garment_front, F_sym_garment_front);
	viewer_2D.data_list.at(3).set_colors(C);

	viewer_2D.launch();
}