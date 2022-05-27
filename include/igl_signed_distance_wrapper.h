// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once
#include <Eigen/Core>
#include <vector>

  // Computes signed distance to a mesh
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   F  #F by ss list of triangle indices, ss should be 3 unless sign_type ==
  //     SIGNED_DISTANCE_TYPE_UNSIGNED
  //   sign_type  method for computing distance _sign_ S
  //   lower_bound  lower bound of distances needed {std::numeric_limits::min}
  //   upper_bound  lower bound of distances needed {std::numeric_limits::max}
  // Outputs:
  //   S  #P list of smallest signed distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //   N  #P by 3 list of closest normals (only set if
  //     sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
  //
  // Known bugs: This only computes distances to triangles. So unreferenced
  // vertices and degenerate triangles are ignored.
 void signed_distance(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const double lower_bound,
    const double upper_bound,
    Eigen::VectorXd & S,
    Eigen::VectorXi & I,
    Eigen::MatrixXd & C,
    Eigen::MatrixXd & N);
  // Computes signed distance to a mesh, with default bounds
  //
  // Inputs:
  //   P  #P by 3 list of query point positions
  //   V  #V by 3 list of vertex positions
  //   F  #F by ss list of triangle indices, ss should be 3 unless sign_type ==
  //     SIGNED_DISTANCE_TYPE_UNSIGNED
  //   sign_type  method for computing distance _sign_ S
  //   lower_bound  lower bound of distances needed {std::numeric_limits::min}
  //   upper_bound  lower bound of distances needed {std::numeric_limits::max}
  // Outputs:
  //   S  #P list of smallest signed distances
  //   I  #P list of facet indices corresponding to smallest distances
  //   C  #P by 3 list of closest points
  //   N  #P by 3 list of closest normals (only set if
  //     sign_type=SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
void signed_distance(
          const Eigen::MatrixXd & P,
          const Eigen::MatrixXd & V,
          const Eigen::MatrixXi & F,
          Eigen::VectorXd & S,
          Eigen::VectorXi & I,
          Eigen::MatrixXd & C,
          Eigen::MatrixXd & N);


