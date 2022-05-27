// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/get_seconds.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/pseudonormal_test.h>
#include <igl/AABB.h>
#include <igl/parallel_for.h>

#include "../include/igl_signed_distance_wrapper.h"
void signed_distance(
        const Eigen::MatrixXd & P,
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const double lower_bound,
        const double upper_bound,
        Eigen::VectorXd & S,
        Eigen::VectorXi & I,
        Eigen::MatrixXd & C,
        Eigen::MatrixXd & N)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;

  const int dim = V.cols();

  assert((V.cols() == 3||V.cols() == 2) && "V should have 3d or 2d positions");
  assert((P.cols() == 3||P.cols() == 2) && "P should have 3d or 2d positions");
  assert(V.cols() == P.cols() && "V should have same dimension as P");

  typedef Eigen::Matrix<double,1,3> RowVector3S;

  // Prepare distance computation
  AABB<Eigen::MatrixXd,3> tree3;
  AABB<Eigen::MatrixXd,2> tree2;
  switch(dim)
  {
    default:
    case 3:
      tree3.init(V,F);
      break;
    case 2:
      tree2.init(V,F);
      break;
  }

  // Need to be Dynamic columns to work with both 2d and 3d
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> FN,VN,EN;
  Eigen::Matrix<int,Eigen::Dynamic,2> E;
  Eigen::Matrix<int,Eigen::Dynamic,1> EMAP;

  switch(dim)
      {
        default:
        case 3:
          // "Signed Distance Computation Using the Angle Weighted Pseudonormal"
          // [Bærentzen & Aanæs 2005]
          per_face_normals(V,F,FN);
          per_vertex_normals(V,F,PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
          per_edge_normals(
            V,F,PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
          break;
        case 2:
          FN.resize(F.rows(),2);
          VN = Eigen::MatrixXd::Zero(V.rows(),2);
          for(int e = 0;e<F.rows();e++)
          {
            // rotate edge vector
            FN(e,0) =  (V(F(e,1),1)-V(F(e,0),1));
            FN(e,1) = -(V(F(e,1),0)-V(F(e,0),0));
            FN.row(e).normalize();
            // add to vertex normal
            VN.row(F(e,1)) += FN.row(e);
            VN.row(F(e,0)) += FN.row(e);
          }
          // normalize to average
          VN.rowwise().normalize();
          break;
      }
    N.resize(P.rows(),dim);
      
  // convert to bounds on (unsiged) squared distances
  typedef double Scalar;
  const Scalar max_abs = std::max(std::abs(lower_bound),std::abs(upper_bound));
  const Scalar up_sqr_d = std::pow(max_abs,2.0);
  const Scalar low_sqr_d = 
    std::pow(std::max(max_abs-(upper_bound-lower_bound),(Scalar)0.0),2.0);

  S.resize(P.rows(),1);
  I.resize(P.rows(),1);
  C.resize(P.rows(),dim);

  parallel_for(P.rows(),[&](const int p)
  //for(int p = 0;p<P.rows();p++)
  {
    RowVector3S q3;
    Eigen::Matrix<double,1,2>  q2;
    switch(P.cols())
    {
      default:
      case 3:
        q3.head(P.row(p).size()) = P.row(p);
        break;
      case 2:
        q2 = P.row(p).head(2);
        break;
    }
    double s=1,sqrd=0;
    Eigen::Matrix<double,1,Eigen::Dynamic>  c;
    Eigen::Matrix<double,1,3> c3;
    Eigen::Matrix<double,1,2>  c2;
    int i=-1;
    // in all cases compute squared unsiged distances
    sqrd = dim==3?
      tree3.squared_distance(V,F,q3,low_sqr_d,up_sqr_d,i,c3):
      tree2.squared_distance(V,F,q2,low_sqr_d,up_sqr_d,i,c2);
    if(sqrd >= up_sqr_d || sqrd < low_sqr_d)
    {
      // Out of bounds gets a nan (nans on grids can be flood filled later using
      // igl::flood_fill)
      S(p) = std::numeric_limits<double>::quiet_NaN();
      I(p) = F.rows()+1;
      C.row(p).setConstant(0);
    }else
    {
      // Determine sign
      RowVector3S n3;
          Eigen::Matrix<double,1,2>  n2;
          dim==3 ?
            pseudonormal_test(V,F,FN,VN,EN,EMAP,q3,i,c3,s,n3):
            // This should use (V,F,FN), not (V,E,EN) since E is auxiliary for
            // 3D case, not the input "F"acets.
            pseudonormal_test(V,F,FN,VN,q2,i,c2,s,n2);
          Eigen::Matrix<double,1,Eigen::Dynamic>  n;
          (dim==3 ? n = n3.template cast<double>() : n = n2.template cast<double>());
          N.row(p) = n.template cast<double>();
      I(p) = i;
      S(p) = s*sqrt(sqrd);
      C.row(p) = (dim==3 ? c=c3 : c=c2).template cast<double>();
    }
  }
  ,10000);
}

 void signed_distance(
         const Eigen::MatrixXd & P,
         const Eigen::MatrixXd & V,
         const Eigen::MatrixXi & F,
         Eigen::VectorXd & S,
         Eigen::VectorXi & I,
         Eigen::MatrixXd & C,
         Eigen::MatrixXd & N)
{
  double lower = std::numeric_limits<double>::min();
  double upper = std::numeric_limits<double>::max();
  return signed_distance(P,V,F,lower,upper,S,I,C,N);
}








