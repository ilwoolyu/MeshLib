# MeshLib

## Description
Mesh Library (MeshLib) - A collection of fundamental surface processing C++ libraries. The library is very lightweight, which can enable the users to handle surface processing easily. This collection is used for implementation of several surface-based applications (see the publication list below).

## Overview - Libraries
* AABB: fast closest triangle access to general mesh
* AABB_Sphere: fast closest triangle access to sphere mesh (this one is at least x2 faster than AABB above)
* Geom: a collection of general processing classes
  * Vector (single and double precision)
  * Euclidean and spherical coordinate system: conversion, projection, SO(3), etc.
  * Functions: associate Legendre polynomial, normal_pdf, normal_cdf, erf.
  * Basic array operations: mean, std, min/max, correlation, etc.
  * 3X3 matrix solver: eig, trace, det.
* Slicer: fast contour extraction given a plane equation
* Spherical Harmonics: real spherical harmonics basis functions
* SurfaceUtil: a collection of surface processing classes
  * Smoothing: NN smoothing
  * Curvature: tensor-based principal curvature approximation
  * Sphere: conformal/area-preserving sphere mapping
* Geodesic: fast marching solver for static Hamilton-Jacobi PDEs/Eikonal equations (ported from Gabriel Peyré's implementation)

## History
* 2009-2011: initial version - base mesh class (Department of Computer Science at KAIST)
* 2011-2016: extended version - surface utilities (Department of Computer Science at the University of North Carolina at Chapel Hill)
* 2017-Present: maintenance - bug fix, optimization, etc. (Electrical Engineering and Computer Science at Vanderbilt University)

## Publications
* Selected publications using this library
  * Lyu, I., Styner, M., Landman, B., Hierarchical Spherical Deformation for Shape Correspondence, <i>Medical Image Computing and Computer Assisted Intervention (MICCAI) 2018</i>, To appear
  * Lyu, I., Kim, S., Girault, J., Gilmore, J., Styner, M., <a href="https://doi.org/10.1016/j.media.2018.06.009">A Cortical Shape-Adaptive Approach to Local Gyrification Index</a>, <i>Medical Image Analysis</i>, 48, 244-258, 2018
  * Lyu, I., Kim, S., Woodward, N., Styner, M., Landman, B., <a href="http://dx.doi.org/10.1109/TMI.2017.2787589">TRACE: A Topological Graph Representation for Automatic Sulcal Curve Extraction</a>, <i>IEEE Transactions on Medical Imaging</i>, 37(7), 1653-1663, 2018
  * Lyu, I., Kim, S., Bullins, J., Gilmore, J., Styner, M., <a href="http://dx.doi.org/10.1007/978-3-319-66182-7_4">Novel Local Shape-Adaptive Gyrification Index with Application to Brain Development</a>, <i>Medical Image Computing and Computer Assisted Intervention (MICCAI) 2017</i>, LNCS10433, 31-39, 2017
  * Lyu, I., Kim, S., Seong, J., Yoo, S., Evans, A., Shi, Y., Sanchez, M., Niethammer, M., Styner, M., <a href="http://dx.doi.org/10.1007/978-3-642-38868-2_31">Group-wise Cortical Correspondence via Sulcal Curve-constrained Entropy Minimization</a>, <i>Information Processing in Medical Imaging (IPMI) 2013</i>, LNCS7917, 364-375, 2013
