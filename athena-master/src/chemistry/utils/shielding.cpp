//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shielding.cpp
//! \brief implementation of functions in class Shielding

// this class header
#include "shielding.hpp"

// Athena headers
#include "../../utils/interp.hpp"

//CO column density for DB table
const Real Shielding::logNCOvDB_[len_NCO_DB_] = {0, 13, 14, 15, 16, 17, 18, 19};
// H2 column densities for DB table
const Real Shielding::logNH2vDB_[len_NH2_DB_] = {0, 19, 20, 21, 22, 23};
// Tabulated shielding factors
const Real Shielding::ThetavDB_[len_NH2_DB_][len_NCO_DB_] = {
  {1.0, 9.681e-1, 7.764e-1, 3.631e-1, 7.013e-2, 1.295e-2, 1.738e-3, 9.985e-5},
  {8.215e-1, 7.916e-1, 6.160e-1, 2.749e-1, 5.351e-2, 1.065e-2, 1.519e-3, 8.818e-5},
  {7.160e-1, 6.900e-1, 5.360e-1, 2.359e-1, 4.416e-2, 8.769e-3, 1.254e-3, 7.558e-5},
  {3.500e-1, 3.415e-1, 2.863e-1, 1.360e-1, 2.500e-2, 4.983e-3, 7.151e-4, 3.796e-5},
  {4.973e-2, 4.877e-2, 4.296e-2, 2.110e-2, 4.958e-3, 9.245e-4, 1.745e-4, 8.377e-6},
  {1.310e-4, 1.293e-4, 1.160e-4, 6.346e-5, 1.822e-5, 6.842e-6, 3.622e-6, 3.572e-7}
};

//-------------------Visser 2009 Table 5---------------------------
const Real Shielding::logNCOV09_[len_NCO_V09_] = {
  0.000,10.000,10.200,10.400,10.600,10.800,11.000,11.200,11.400,11.600,11.800,
  12.000,12.200,12.400,12.600,12.800,13.000,13.200,13.400,13.600,13.800,14.000,
  14.200,14.400,14.600,14.800,15.000,15.200,15.400,15.600,15.800,16.000,16.200,
  16.400,16.600,16.800,17.000,17.200,17.400,17.600,17.800,18.000,18.200,18.400,
  18.600,18.800,19.000};
const Real Shielding::logNH2V09_[len_NH2_V09_] = {
  0.000,15.000,15.200,15.400,15.600,15.800,16.000,16.200,16.400,16.600,16.800,17.000,
  17.200,17.400,17.600,17.800,18.000,18.200,18.400,18.600,18.800,19.000,19.200,19.400,
  19.600,19.800,20.000,20.200,20.400,20.600,20.800,21.000,21.200,21.400,21.600,21.800,
  22.000,22.200,22.400,22.600,22.800,23.000};
const Real Shielding::ThetaV09_[len_NH2_V09_][len_NCO_V09_] = {
  {1.000e+00,9.997e-01,9.995e-01,9.992e-01,9.988e-01,9.981e-01,9.970e-01,9.953e-01,
   9.926e-01,9.883e-01,9.817e-01,9.716e-01,9.563e-01,9.338e-01,9.021e-01,8.599e-01,
   8.080e-01,7.498e-01,6.900e-01,6.323e-01,5.777e-01,5.250e-01,4.720e-01,4.177e-01,
   3.614e-01,3.028e-01,2.434e-01,1.871e-01,1.387e-01,1.012e-01,7.401e-02,5.467e-02,
   4.075e-02,3.063e-02,2.323e-02,1.775e-02,1.362e-02,1.044e-02,7.963e-03,6.037e-03,
   4.541e-03,3.378e-03,2.470e-03,1.759e-03,1.210e-03,8.046e-04,5.240e-04},
  {8.985e-01,8.983e-01,8.981e-01,8.978e-01,8.974e-01,8.967e-01,8.956e-01,8.939e-01,
   8.913e-01,8.871e-01,8.807e-01,8.707e-01,8.558e-01,8.338e-01,8.030e-01,7.621e-01,
   7.122e-01,6.569e-01,6.014e-01,5.494e-01,5.021e-01,4.578e-01,4.137e-01,3.681e-01,
   3.201e-01,2.694e-01,2.177e-01,1.683e-01,1.256e-01,9.224e-02,6.786e-02,5.034e-02,
   3.767e-02,2.844e-02,2.167e-02,1.664e-02,1.282e-02,9.868e-03,7.556e-03,5.743e-03,
   4.325e-03,3.223e-03,2.363e-03,1.690e-03,1.169e-03,7.815e-04,5.112e-04},
  {8.966e-01,8.963e-01,8.962e-01,8.959e-01,8.955e-01,8.948e-01,8.937e-01,8.920e-01,
   8.894e-01,8.852e-01,8.788e-01,8.688e-01,8.539e-01,8.319e-01,8.011e-01,7.602e-01,
   7.103e-01,6.551e-01,5.996e-01,5.476e-01,5.004e-01,4.562e-01,4.122e-01,3.667e-01,
   3.190e-01,2.685e-01,2.171e-01,1.679e-01,1.254e-01,9.214e-02,6.781e-02,5.031e-02,
   3.765e-02,2.842e-02,2.166e-02,1.663e-02,1.282e-02,9.865e-03,7.554e-03,5.741e-03,
   4.323e-03,3.222e-03,2.362e-03,1.689e-03,1.169e-03,7.811e-04,5.110e-04},
  {8.949e-01,8.946e-01,8.944e-01,8.941e-01,8.937e-01,8.930e-01,8.920e-01,8.903e-01,
   8.876e-01,8.834e-01,8.770e-01,8.671e-01,8.521e-01,8.302e-01,7.993e-01,7.585e-01,
   7.086e-01,6.533e-01,5.979e-01,5.460e-01,4.988e-01,4.546e-01,4.107e-01,3.655e-01,
   3.179e-01,2.677e-01,2.165e-01,1.676e-01,1.252e-01,9.204e-02,6.776e-02,5.028e-02,
   3.763e-02,2.841e-02,2.165e-02,1.662e-02,1.281e-02,9.861e-03,7.551e-03,5.739e-03,
   4.322e-03,3.220e-03,2.361e-03,1.689e-03,1.168e-03,7.808e-04,5.108e-04},
  {8.932e-01,8.929e-01,8.927e-01,8.924e-01,8.920e-01,8.913e-01,8.903e-01,8.886e-01,
   8.859e-01,8.818e-01,8.753e-01,8.654e-01,8.504e-01,8.285e-01,7.976e-01,7.568e-01,
   7.069e-01,6.517e-01,5.962e-01,5.444e-01,4.972e-01,4.531e-01,4.094e-01,3.642e-01,
   3.169e-01,2.669e-01,2.159e-01,1.672e-01,1.250e-01,9.193e-02,6.770e-02,5.025e-02,
   3.761e-02,2.839e-02,2.164e-02,1.661e-02,1.281e-02,9.858e-03,7.549e-03,5.737e-03,
   4.320e-03,3.219e-03,2.360e-03,1.688e-03,1.168e-03,7.805e-04,5.106e-04},
  {8.915e-01,8.912e-01,8.911e-01,8.908e-01,8.904e-01,8.897e-01,8.886e-01,8.869e-01,
   8.843e-01,8.801e-01,8.737e-01,8.637e-01,8.488e-01,8.269e-01,7.960e-01,7.551e-01,
   7.053e-01,6.501e-01,5.947e-01,5.428e-01,4.957e-01,4.517e-01,4.080e-01,3.630e-01,
   3.159e-01,2.661e-01,2.154e-01,1.669e-01,1.248e-01,9.182e-02,6.764e-02,5.022e-02,
   3.759e-02,2.838e-02,2.162e-02,1.661e-02,1.280e-02,9.854e-03,7.546e-03,5.735e-03,
   4.319e-03,3.218e-03,2.359e-03,1.687e-03,1.167e-03,7.802e-04,5.104e-04},
  {8.899e-01,8.896e-01,8.895e-01,8.892e-01,8.888e-01,8.881e-01,8.870e-01,8.853e-01,
   8.827e-01,8.785e-01,8.721e-01,8.621e-01,8.472e-01,8.253e-01,7.944e-01,7.536e-01,
   7.037e-01,6.485e-01,5.931e-01,5.413e-01,4.942e-01,4.503e-01,4.067e-01,3.618e-01,
   3.148e-01,2.653e-01,2.148e-01,1.665e-01,1.246e-01,9.170e-02,6.758e-02,5.018e-02,
   3.757e-02,2.837e-02,2.161e-02,1.660e-02,1.280e-02,9.851e-03,7.544e-03,5.733e-03,
   4.317e-03,3.216e-03,2.358e-03,1.686e-03,1.167e-03,7.799e-04,5.103e-04},
  {8.855e-01,8.852e-01,8.850e-01,8.848e-01,8.843e-01,8.837e-01,8.826e-01,8.809e-01,
   8.782e-01,8.741e-01,8.676e-01,8.577e-01,8.428e-01,8.209e-01,7.900e-01,7.492e-01,
   6.993e-01,6.442e-01,5.888e-01,5.371e-01,4.901e-01,4.463e-01,4.028e-01,3.582e-01,
   3.114e-01,2.622e-01,2.120e-01,1.642e-01,1.227e-01,9.024e-02,6.640e-02,4.917e-02,
   3.667e-02,2.759e-02,2.096e-02,1.608e-02,1.239e-02,9.538e-03,7.308e-03,5.558e-03,
   4.189e-03,3.123e-03,2.290e-03,1.638e-03,1.133e-03,7.572e-04,4.958e-04},
  {8.834e-01,8.831e-01,8.829e-01,8.826e-01,8.822e-01,8.815e-01,8.805e-01,8.788e-01,
   8.761e-01,8.720e-01,8.655e-01,8.556e-01,8.406e-01,8.187e-01,7.879e-01,7.471e-01,
   6.972e-01,6.421e-01,5.867e-01,5.350e-01,4.881e-01,4.443e-01,4.010e-01,3.565e-01,
   3.099e-01,2.609e-01,2.111e-01,1.635e-01,1.223e-01,9.003e-02,6.629e-02,4.911e-02,
   3.664e-02,2.757e-02,2.095e-02,1.607e-02,1.238e-02,9.533e-03,7.305e-03,5.555e-03,
   4.187e-03,3.121e-03,2.289e-03,1.637e-03,1.132e-03,7.567e-04,4.955e-04},
  {8.814e-01,8.811e-01,8.809e-01,8.807e-01,8.802e-01,8.796e-01,8.785e-01,8.768e-01,
   8.741e-01,8.700e-01,8.635e-01,8.536e-01,8.387e-01,8.168e-01,7.859e-01,7.451e-01,
   6.953e-01,6.401e-01,5.848e-01,5.331e-01,4.862e-01,4.425e-01,3.993e-01,3.549e-01,
   3.085e-01,2.597e-01,2.101e-01,1.629e-01,1.219e-01,8.978e-02,6.616e-02,4.905e-02,
   3.661e-02,2.755e-02,2.094e-02,1.606e-02,1.237e-02,9.529e-03,7.302e-03,5.553e-03,
   4.185e-03,3.120e-03,2.288e-03,1.636e-03,1.131e-03,7.562e-04,4.952e-04},
  {8.792e-01,8.789e-01,8.788e-01,8.785e-01,8.781e-01,8.774e-01,8.763e-01,8.746e-01,
   8.720e-01,8.678e-01,8.614e-01,8.515e-01,8.365e-01,8.146e-01,7.838e-01,7.429e-01,
   6.931e-01,6.380e-01,5.827e-01,5.310e-01,4.842e-01,4.405e-01,3.974e-01,3.531e-01,
   3.068e-01,2.583e-01,2.090e-01,1.620e-01,1.214e-01,8.948e-02,6.601e-02,4.897e-02,
   3.657e-02,2.753e-02,2.092e-02,1.605e-02,1.237e-02,9.523e-03,7.297e-03,5.549e-03,
   4.182e-03,3.117e-03,2.286e-03,1.634e-03,1.130e-03,7.557e-04,4.948e-04},
  {8.766e-01,8.764e-01,8.762e-01,8.759e-01,8.755e-01,8.748e-01,8.737e-01,8.720e-01,
   8.694e-01,8.652e-01,8.588e-01,8.489e-01,8.339e-01,8.120e-01,7.812e-01,7.404e-01,
   6.906e-01,6.355e-01,5.802e-01,5.285e-01,4.817e-01,4.381e-01,3.951e-01,3.509e-01,
   3.049e-01,2.566e-01,2.076e-01,1.611e-01,1.207e-01,8.911e-02,6.581e-02,4.887e-02,
   3.652e-02,2.749e-02,2.090e-02,1.603e-02,1.235e-02,9.515e-03,7.291e-03,5.545e-03,
   4.178e-03,3.114e-03,2.284e-03,1.633e-03,1.129e-03,7.549e-04,4.943e-04},
  {8.735e-01,8.732e-01,8.730e-01,8.728e-01,8.723e-01,8.716e-01,8.706e-01,8.689e-01,
   8.662e-01,8.621e-01,8.556e-01,8.457e-01,8.308e-01,8.089e-01,7.781e-01,7.372e-01,
   6.874e-01,6.324e-01,5.771e-01,5.255e-01,4.787e-01,4.352e-01,3.923e-01,3.483e-01,
   3.025e-01,2.546e-01,2.060e-01,1.598e-01,1.199e-01,8.861e-02,6.554e-02,4.874e-02,
   3.644e-02,2.745e-02,2.087e-02,1.601e-02,1.234e-02,9.504e-03,7.284e-03,5.539e-03,
   4.173e-03,3.111e-03,2.281e-03,1.630e-03,1.128e-03,7.540e-04,4.937e-04},
  {8.697e-01,8.694e-01,8.692e-01,8.689e-01,8.685e-01,8.678e-01,8.668e-01,8.651e-01,
   8.624e-01,8.583e-01,8.518e-01,8.419e-01,8.270e-01,8.051e-01,7.743e-01,7.335e-01,
   6.837e-01,6.286e-01,5.734e-01,5.218e-01,4.751e-01,4.317e-01,3.889e-01,3.451e-01,
   2.996e-01,2.520e-01,2.039e-01,1.582e-01,1.188e-01,8.795e-02,6.517e-02,4.854e-02,
   3.634e-02,2.738e-02,2.083e-02,1.598e-02,1.232e-02,9.490e-03,7.273e-03,5.530e-03,
   4.167e-03,3.105e-03,2.277e-03,1.628e-03,1.126e-03,7.527e-04,4.929e-04},
  {8.652e-01,8.649e-01,8.647e-01,8.644e-01,8.640e-01,8.633e-01,8.623e-01,8.606e-01,
   8.579e-01,8.538e-01,8.473e-01,8.374e-01,8.225e-01,8.006e-01,7.698e-01,7.290e-01,
   6.793e-01,6.242e-01,5.690e-01,5.175e-01,4.709e-01,4.276e-01,3.849e-01,3.414e-01,
   2.962e-01,2.490e-01,2.014e-01,1.563e-01,1.175e-01,8.712e-02,6.469e-02,4.828e-02,
   3.619e-02,2.729e-02,2.077e-02,1.594e-02,1.229e-02,9.471e-03,7.259e-03,5.519e-03,
   4.158e-03,3.098e-03,2.271e-03,1.624e-03,1.123e-03,7.509e-04,4.919e-04},
  {8.600e-01,8.597e-01,8.595e-01,8.593e-01,8.588e-01,8.582e-01,8.571e-01,8.554e-01,
   8.528e-01,8.486e-01,8.422e-01,8.323e-01,8.173e-01,7.955e-01,7.647e-01,7.239e-01,
   6.742e-01,6.192e-01,5.640e-01,5.126e-01,4.660e-01,4.229e-01,3.804e-01,3.371e-01,
   2.923e-01,2.456e-01,1.985e-01,1.541e-01,1.159e-01,8.608e-02,6.406e-02,4.791e-02,
   3.598e-02,2.717e-02,2.069e-02,1.589e-02,1.226e-02,9.445e-03,7.239e-03,5.504e-03,
   4.146e-03,3.089e-03,2.264e-03,1.618e-03,1.119e-03,7.486e-04,4.904e-04},
  {8.543e-01,8.540e-01,8.539e-01,8.536e-01,8.532e-01,8.525e-01,8.514e-01,8.497e-01,
   8.471e-01,8.429e-01,8.365e-01,8.266e-01,8.117e-01,7.898e-01,7.591e-01,7.183e-01,
   6.686e-01,6.137e-01,5.586e-01,5.072e-01,4.608e-01,4.178e-01,3.755e-01,3.325e-01,
   2.880e-01,2.418e-01,1.953e-01,1.516e-01,1.140e-01,8.481e-02,6.326e-02,4.743e-02,
   3.569e-02,2.699e-02,2.058e-02,1.582e-02,1.221e-02,9.411e-03,7.213e-03,5.483e-03,
   4.130e-03,3.076e-03,2.254e-03,1.611e-03,1.115e-03,7.456e-04,4.885e-04},
  {8.482e-01,8.479e-01,8.477e-01,8.475e-01,8.470e-01,8.464e-01,8.453e-01,8.436e-01,
   8.410e-01,8.368e-01,8.304e-01,8.205e-01,8.056e-01,7.838e-01,7.531e-01,7.124e-01,
   6.627e-01,6.079e-01,5.529e-01,5.016e-01,4.553e-01,4.125e-01,3.704e-01,3.277e-01,
   2.836e-01,2.379e-01,1.920e-01,1.489e-01,1.120e-01,8.334e-02,6.227e-02,4.680e-02,
   3.530e-02,2.675e-02,2.043e-02,1.572e-02,1.214e-02,9.363e-03,7.177e-03,5.456e-03,
   4.108e-03,3.059e-03,2.241e-03,1.602e-03,1.108e-03,7.415e-04,4.859e-04},
  {8.416e-01,8.414e-01,8.412e-01,8.409e-01,8.405e-01,8.398e-01,8.387e-01,8.371e-01,
   8.344e-01,8.303e-01,8.238e-01,8.140e-01,7.991e-01,7.773e-01,7.466e-01,7.060e-01,
   6.565e-01,6.017e-01,5.468e-01,4.957e-01,4.496e-01,4.069e-01,3.651e-01,3.228e-01,
   2.792e-01,2.340e-01,1.886e-01,1.461e-01,1.098e-01,8.169e-02,6.110e-02,4.602e-02,
   3.479e-02,2.643e-02,2.023e-02,1.559e-02,1.205e-02,9.298e-03,7.129e-03,5.419e-03,
   4.079e-03,3.037e-03,2.225e-03,1.590e-03,1.100e-03,7.364e-04,4.827e-04},
  {8.345e-01,8.342e-01,8.340e-01,8.337e-01,8.333e-01,8.326e-01,8.316e-01,8.299e-01,
   8.273e-01,8.231e-01,8.167e-01,8.069e-01,7.920e-01,7.703e-01,7.397e-01,6.992e-01,
   6.498e-01,5.952e-01,5.405e-01,4.895e-01,4.436e-01,4.012e-01,3.597e-01,3.178e-01,
   2.747e-01,2.300e-01,1.853e-01,1.433e-01,1.075e-01,7.995e-02,5.981e-02,4.510e-02,
   3.417e-02,2.602e-02,1.995e-02,1.541e-02,1.193e-02,9.209e-03,7.063e-03,5.369e-03,
   4.041e-03,3.008e-03,2.204e-03,1.575e-03,1.090e-03,7.299e-04,4.785e-04},
  {8.265e-01,8.262e-01,8.260e-01,8.258e-01,8.254e-01,8.247e-01,8.236e-01,8.220e-01,
   8.193e-01,8.152e-01,8.088e-01,7.990e-01,7.842e-01,7.626e-01,7.321e-01,6.918e-01,
   6.425e-01,5.881e-01,5.337e-01,4.830e-01,4.373e-01,3.952e-01,3.542e-01,3.127e-01,
   2.701e-01,2.261e-01,1.820e-01,1.406e-01,1.053e-01,7.817e-02,5.843e-02,4.407e-02,
   3.343e-02,2.550e-02,1.960e-02,1.516e-02,1.176e-02,9.086e-03,6.973e-03,5.302e-03,
   3.991e-03,2.971e-03,2.176e-03,1.556e-03,1.077e-03,7.215e-04,4.731e-04},
  {8.176e-01,8.173e-01,8.171e-01,8.169e-01,8.164e-01,8.158e-01,8.147e-01,8.130e-01,
   8.104e-01,8.063e-01,8.000e-01,7.903e-01,7.756e-01,7.540e-01,7.237e-01,6.836e-01,
   6.347e-01,5.806e-01,5.265e-01,4.761e-01,4.308e-01,3.891e-01,3.485e-01,3.076e-01,
   2.656e-01,2.222e-01,1.787e-01,1.378e-01,1.031e-01,7.637e-02,5.701e-02,4.297e-02,
   3.260e-02,2.488e-02,1.915e-02,1.484e-02,1.152e-02,8.919e-03,6.851e-03,5.212e-03,
   3.925e-03,2.922e-03,2.141e-03,1.531e-03,1.061e-03,7.108e-04,4.662e-04},
  {8.073e-01,8.070e-01,8.069e-01,8.066e-01,8.062e-01,8.055e-01,8.045e-01,8.028e-01,
   8.002e-01,7.962e-01,7.899e-01,7.802e-01,7.657e-01,7.443e-01,7.143e-01,6.746e-01,
   6.261e-01,5.725e-01,5.189e-01,4.689e-01,4.240e-01,3.828e-01,3.427e-01,3.023e-01,
   2.609e-01,2.182e-01,1.753e-01,1.350e-01,1.008e-01,7.452e-02,5.552e-02,4.179e-02,
   3.167e-02,2.417e-02,1.861e-02,1.443e-02,1.122e-02,8.697e-03,6.689e-03,5.093e-03,
   3.838e-03,2.859e-03,2.096e-03,1.500e-03,1.040e-03,6.969e-04,4.573e-04},
  {7.949e-01,7.946e-01,7.944e-01,7.942e-01,7.938e-01,7.931e-01,7.921e-01,7.904e-01,
   7.879e-01,7.839e-01,7.777e-01,7.682e-01,7.538e-01,7.328e-01,7.032e-01,6.640e-01,
   6.162e-01,5.633e-01,5.104e-01,4.611e-01,4.168e-01,3.760e-01,3.364e-01,2.966e-01,
   2.559e-01,2.138e-01,1.716e-01,1.320e-01,9.827e-02,7.248e-02,5.388e-02,4.048e-02,
   3.064e-02,2.335e-02,1.796e-02,1.393e-02,1.084e-02,8.411e-03,6.477e-03,4.938e-03,
   3.725e-03,2.778e-03,2.038e-03,1.460e-03,1.013e-03,6.792e-04,4.459e-04},
  {7.784e-01,7.782e-01,7.780e-01,7.778e-01,7.774e-01,7.767e-01,7.757e-01,7.741e-01,
   7.716e-01,7.677e-01,7.617e-01,7.524e-01,7.383e-01,7.177e-01,6.888e-01,6.504e-01,
   6.037e-01,5.519e-01,5.000e-01,4.516e-01,4.081e-01,3.680e-01,3.290e-01,2.898e-01,
   2.498e-01,2.085e-01,1.671e-01,1.283e-01,9.528e-02,7.009e-02,5.199e-02,3.899e-02,
   2.946e-02,2.242e-02,1.722e-02,1.334e-02,1.038e-02,8.056e-03,6.209e-03,4.740e-03,
   3.581e-03,2.675e-03,1.965e-03,1.409e-03,9.785e-04,6.566e-04,4.313e-04},
  {7.553e-01,7.550e-01,7.549e-01,7.546e-01,7.543e-01,7.536e-01,7.527e-01,7.511e-01,
   7.487e-01,7.450e-01,7.391e-01,7.301e-01,7.166e-01,6.967e-01,6.688e-01,6.317e-01,
   5.865e-01,5.364e-01,4.861e-01,4.391e-01,3.966e-01,3.575e-01,3.194e-01,2.811e-01,
   2.420e-01,2.017e-01,1.614e-01,1.236e-01,9.158e-02,6.719e-02,4.972e-02,3.723e-02,
   2.808e-02,2.134e-02,1.637e-02,1.267e-02,9.843e-03,7.635e-03,5.885e-03,4.496e-03,
   3.402e-03,2.546e-03,1.874e-03,1.346e-03,9.352e-04,6.280e-04,4.129e-04},
  {7.223e-01,7.220e-01,7.219e-01,7.216e-01,7.213e-01,7.207e-01,7.198e-01,7.183e-01,
   7.160e-01,7.125e-01,7.069e-01,6.984e-01,6.856e-01,6.668e-01,6.404e-01,6.053e-01,
   5.624e-01,5.149e-01,4.670e-01,4.221e-01,3.812e-01,3.434e-01,3.066e-01,2.695e-01,
   2.317e-01,1.929e-01,1.540e-01,1.177e-01,8.696e-02,6.364e-02,4.701e-02,3.515e-02,
   2.649e-02,2.010e-02,1.540e-02,1.190e-02,9.231e-03,7.150e-03,5.506e-03,4.207e-03,
   3.186e-03,2.388e-03,1.761e-03,1.266e-03,8.811e-04,5.923e-04,3.899e-04},
  {6.758e-01,6.756e-01,6.754e-01,6.752e-01,6.749e-01,6.743e-01,6.735e-01,6.722e-01,
   6.701e-01,6.668e-01,6.618e-01,6.540e-01,6.423e-01,6.250e-01,6.008e-01,5.686e-01,
   5.292e-01,4.854e-01,4.410e-01,3.991e-01,3.607e-01,3.249e-01,2.898e-01,2.546e-01,
   2.186e-01,1.816e-01,1.448e-01,1.104e-01,8.136e-02,5.941e-02,4.382e-02,3.274e-02,
   2.465e-02,1.869e-02,1.430e-02,1.103e-02,8.540e-03,6.601e-03,5.074e-03,3.871e-03,
   2.931e-03,2.198e-03,1.623e-03,1.169e-03,8.145e-04,5.484e-04,3.619e-04},
  {6.127e-01,6.125e-01,6.124e-01,6.122e-01,6.119e-01,6.114e-01,6.107e-01,6.095e-01,
   6.077e-01,6.049e-01,6.005e-01,5.937e-01,5.835e-01,5.685e-01,5.473e-01,5.192e-01,
   4.847e-01,4.461e-01,4.067e-01,3.691e-01,3.342e-01,3.012e-01,2.686e-01,2.358e-01,
   2.022e-01,1.678e-01,1.335e-01,1.016e-01,7.474e-02,5.447e-02,4.013e-02,2.996e-02,
   2.255e-02,1.708e-02,1.304e-02,1.004e-02,7.759e-03,5.983e-03,4.587e-03,3.490e-03,
   2.637e-03,1.975e-03,1.459e-03,1.052e-03,7.343e-04,4.956e-04,3.281e-04},
  {5.310e-01,5.309e-01,5.308e-01,5.306e-01,5.304e-01,5.300e-01,5.295e-01,5.285e-01,
   5.271e-01,5.248e-01,5.212e-01,5.158e-01,5.076e-01,4.955e-01,4.784e-01,4.557e-01,
   4.276e-01,3.959e-01,3.632e-01,3.313e-01,3.010e-01,2.718e-01,2.426e-01,2.129e-01,
   1.825e-01,1.513e-01,1.202e-01,9.136e-02,6.707e-02,4.880e-02,3.591e-02,2.678e-02,
   2.014e-02,1.523e-02,1.161e-02,8.915e-03,6.873e-03,5.286e-03,4.039e-03,3.062e-03,
   2.305e-03,1.721e-03,1.270e-03,9.161e-04,6.405e-04,4.336e-04,2.884e-04},
  {4.328e-01,4.327e-01,4.326e-01,4.325e-01,4.323e-01,4.321e-01,4.317e-01,4.310e-01,
   4.300e-01,4.283e-01,4.258e-01,4.220e-01,4.161e-01,4.075e-01,3.952e-01,3.788e-01,
   3.584e-01,3.350e-01,3.103e-01,2.854e-01,2.609e-01,2.364e-01,2.115e-01,1.857e-01,
   1.592e-01,1.320e-01,1.048e-01,7.956e-02,5.830e-02,4.233e-02,3.109e-02,2.315e-02,
   1.737e-02,1.310e-02,9.965e-03,7.636e-03,5.872e-03,4.504e-03,3.430e-03,2.589e-03,
   1.940e-03,1.443e-03,1.061e-03,7.650e-04,5.354e-04,3.635e-04,2.432e-04},
  {3.260e-01,3.260e-01,3.259e-01,3.258e-01,3.258e-01,3.256e-01,3.253e-01,3.250e-01,
   3.243e-01,3.234e-01,3.219e-01,3.196e-01,3.161e-01,3.110e-01,3.036e-01,2.937e-01,
   2.810e-01,2.661e-01,2.497e-01,2.324e-01,2.143e-01,1.953e-01,1.753e-01,1.543e-01,
   1.325e-01,1.098e-01,8.726e-02,6.618e-02,4.841e-02,3.506e-02,2.569e-02,1.907e-02,
   1.426e-02,1.073e-02,8.134e-03,6.216e-03,4.768e-03,3.647e-03,2.768e-03,2.081e-03,
   1.552e-03,1.150e-03,8.429e-04,6.062e-04,4.238e-04,2.882e-04,1.941e-04},
  {2.241e-01,2.241e-01,2.241e-01,2.241e-01,2.240e-01,2.240e-01,2.238e-01,2.237e-01,
   2.234e-01,2.229e-01,2.223e-01,2.212e-01,2.196e-01,2.172e-01,2.137e-01,2.089e-01,
   2.025e-01,1.947e-01,1.854e-01,1.748e-01,1.629e-01,1.496e-01,1.349e-01,1.192e-01,
   1.025e-01,8.509e-02,6.767e-02,5.133e-02,3.752e-02,2.712e-02,1.982e-02,1.466e-02,
   1.092e-02,8.184e-03,6.184e-03,4.713e-03,3.606e-03,2.751e-03,2.083e-03,1.561e-03,
   1.161e-03,8.580e-04,6.273e-04,4.495e-04,3.130e-04,2.126e-04,1.440e-04},
  {1.394e-01,1.394e-01,1.394e-01,1.393e-01,1.393e-01,1.393e-01,1.393e-01,1.392e-01,
   1.391e-01,1.390e-01,1.387e-01,1.384e-01,1.378e-01,1.369e-01,1.356e-01,1.338e-01,
   1.312e-01,1.277e-01,1.232e-01,1.176e-01,1.106e-01,1.024e-01,9.286e-02,8.230e-02,
   7.094e-02,5.903e-02,4.704e-02,3.576e-02,2.616e-02,1.892e-02,1.381e-02,1.020e-02,
   7.568e-03,5.648e-03,4.253e-03,3.233e-03,2.470e-03,1.882e-03,1.423e-03,1.067e-03,
   7.942e-04,5.872e-04,4.290e-04,3.062e-04,2.119e-04,1.432e-04,9.745e-05},
  {7.604e-02,7.604e-02,7.603e-02,7.603e-02,7.603e-02,7.602e-02,7.601e-02,7.599e-02,
   7.596e-02,7.591e-02,7.584e-02,7.573e-02,7.555e-02,7.528e-02,7.486e-02,7.422e-02,
   7.327e-02,7.189e-02,6.995e-02,6.729e-02,6.380e-02,5.941e-02,5.415e-02,4.814e-02,
   4.160e-02,3.472e-02,2.779e-02,2.123e-02,1.561e-02,1.134e-02,8.295e-03,6.124e-03,
   4.536e-03,3.374e-03,2.534e-03,1.924e-03,1.470e-03,1.122e-03,8.510e-04,6.409e-04,
   4.796e-04,3.565e-04,2.612e-04,1.862e-04,1.281e-04,8.610e-05,5.885e-05},
  {3.382e-02,3.382e-02,3.382e-02,3.381e-02,3.381e-02,3.381e-02,3.381e-02,3.380e-02,
   3.379e-02,3.378e-02,3.376e-02,3.372e-02,3.366e-02,3.357e-02,3.343e-02,3.322e-02,
   3.289e-02,3.239e-02,3.165e-02,3.060e-02,2.915e-02,2.725e-02,2.490e-02,2.219e-02,
   1.922e-02,1.612e-02,1.299e-02,1.003e-02,7.466e-03,5.480e-03,4.038e-03,2.990e-03,
   2.213e-03,1.643e-03,1.232e-03,9.352e-04,7.164e-04,5.497e-04,4.205e-04,3.201e-04,
   2.425e-04,1.823e-04,1.348e-04,9.649e-05,6.645e-05,4.483e-05,3.108e-05},
  {1.108e-02,1.108e-02,1.108e-02,1.108e-02,1.108e-02,1.108e-02,1.107e-02,1.107e-02,
   1.107e-02,1.107e-02,1.106e-02,1.105e-02,1.103e-02,1.101e-02,1.097e-02,1.090e-02,
   1.081e-02,1.066e-02,1.043e-02,1.011e-02,9.646e-03,9.033e-03,8.269e-03,7.385e-03,
   6.427e-03,5.434e-03,4.441e-03,3.492e-03,2.653e-03,1.982e-03,1.479e-03,1.102e-03,
   8.166e-04,6.052e-04,4.528e-04,3.439e-04,2.644e-04,2.048e-04,1.589e-04,1.233e-04,
   9.539e-05,7.329e-05,5.529e-05,4.038e-05,2.845e-05,1.982e-05,1.437e-05},
  {2.364e-03,2.364e-03,2.364e-03,2.364e-03,2.364e-03,2.364e-03,2.363e-03,2.363e-03,
   2.363e-03,2.362e-03,2.360e-03,2.358e-03,2.355e-03,2.350e-03,2.342e-03,2.330e-03,
   2.310e-03,2.281e-03,2.236e-03,2.171e-03,2.078e-03,1.955e-03,1.802e-03,1.625e-03,
   1.434e-03,1.236e-03,1.034e-03,8.355e-04,6.521e-04,4.978e-04,3.767e-04,2.830e-04,
   2.105e-04,1.558e-04,1.162e-04,8.813e-05,6.808e-05,5.344e-05,4.250e-05,3.411e-05,
   2.751e-05,2.211e-05,1.752e-05,1.356e-05,1.029e-05,7.862e-06,6.317e-06},
  {2.983e-04,2.982e-04,2.982e-04,2.982e-04,2.982e-04,2.982e-04,2.982e-04,2.982e-04,
   2.981e-04,2.980e-04,2.979e-04,2.977e-04,2.974e-04,2.970e-04,2.962e-04,2.951e-04,
   2.933e-04,2.905e-04,2.863e-04,2.801e-04,2.713e-04,2.595e-04,2.446e-04,2.268e-04,
   2.066e-04,1.841e-04,1.592e-04,1.325e-04,1.059e-04,8.236e-05,6.326e-05,4.825e-05,
   3.638e-05,2.723e-05,2.054e-05,1.583e-05,1.253e-05,1.023e-05,8.603e-06,7.410e-06,
   6.479e-06,5.689e-06,4.977e-06,4.328e-06,3.767e-06,3.331e-06,3.040e-06},
  {2.297e-05,2.297e-05,2.297e-05,2.296e-05,2.296e-05,2.296e-05,2.296e-05,2.296e-05,
   2.296e-05,2.296e-05,2.295e-05,2.295e-05,2.294e-05,2.292e-05,2.290e-05,2.286e-05,
   2.280e-05,2.270e-05,2.255e-05,2.233e-05,2.200e-05,2.154e-05,2.089e-05,2.002e-05,
   1.886e-05,1.736e-05,1.548e-05,1.326e-05,1.093e-05,8.801e-06,7.086e-06,5.751e-06,
   4.678e-06,3.829e-06,3.200e-06,2.755e-06,2.446e-06,2.236e-06,2.094e-06,1.996e-06,
   1.922e-06,1.857e-06,1.796e-06,1.737e-06,1.684e-06,1.640e-06,1.608e-06},
  {1.546e-06,1.546e-06,1.546e-06,1.546e-06,1.546e-06,1.546e-06,1.545e-06,1.545e-06,
   1.545e-06,1.545e-06,1.545e-06,1.545e-06,1.545e-06,1.545e-06,1.544e-06,1.544e-06,
   1.542e-06,1.541e-06,1.538e-06,1.534e-06,1.527e-06,1.516e-06,1.500e-06,1.477e-06,
   1.442e-06,1.393e-06,1.329e-06,1.251e-06,1.168e-06,1.092e-06,1.032e-06,9.867e-07,
   9.497e-07,9.197e-07,8.973e-07,8.813e-07,8.701e-07,8.624e-07,8.574e-07,8.539e-07,
   8.513e-07,8.490e-07,8.467e-07,8.444e-07,8.421e-07,8.397e-07,8.374e-07},
  {3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,
   3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,3.938e-07,
   3.938e-07,3.937e-07,3.937e-07,3.937e-07,3.937e-07,3.936e-07,3.935e-07,3.933e-07,
   3.931e-07,3.928e-07,3.923e-07,3.918e-07,3.913e-07,3.908e-07,3.904e-07,3.901e-07,
   3.898e-07,3.896e-07,3.894e-07,3.893e-07,3.893e-07,3.892e-07,3.891e-07,3.891e-07,
   3.890e-07,3.890e-07,3.889e-07,3.887e-07,3.885e-07,3.881e-07,3.875e-07}
};

//----------------------------------------------------------------------------------------
//! constructor for class Shielding
Shielding::Shielding() {}

//----------------------------------------------------------------------------------------
//! \fn Real Shielding::fShield_CO_vDB(const Real NCO, const Real NH2)
//! \brief CO self-sheilding and sheilding by H2,
//!         from van Dishoeck & Black's shielding function.
//!
//! NCO, NH2: column densith of CO and H2 in cm^-2.
Real Shielding::fShield_CO_vDB(const Real NCO, const Real NH2) {
  const Real logNCO = std::log10(NCO);
  const Real logNH2 = std::log10(NH2);
  int iCO0, iCO1;
  int iH20, iH21;
  //interpretation value of first, second line, and final value
  Real fl1, fl2, fl;
  //find which two points are we interpolate/extroplating
  iCO0 = Interpolation::LinearInterpIndex(len_NCO_DB_, logNCOvDB_, logNCO);
  iCO1 = iCO0+1;
  iH20 = Interpolation::LinearInterpIndex(len_NH2_DB_, logNH2vDB_, logNH2);
  iH21 = iH20+1;
  //linear interpretation of the rows
  fl1 = Interpolation::LinearInterp(logNCOvDB_[iCO0], logNCOvDB_[iCO1],
                                    std::log(ThetavDB_[iH20][iCO0]),
                                    std::log(ThetavDB_[iH20][iCO1]), logNCO);
  fl2 = Interpolation::LinearInterp(logNCOvDB_[iCO0], logNCOvDB_[iCO1],
                                    std::log(ThetavDB_[iH21][iCO0]),
                                    std::log(ThetavDB_[iH21][iCO1]), logNCO);
  //linear interpretation of the column
  fl = Interpolation::LinearInterp(logNH2vDB_[iH20], logNH2vDB_[iH21],
                    fl1, fl2, logNH2);
  return std::exp(fl);
}

//----------------------------------------------------------------------------------------
//! \fn Real Shielding::fShield_CO_V09(const Real NCO, const Real NH2)
//! \brief CO self-sheilding and sheilding by H2, from Visser+2009, Table 5
//!
//! NCO, NH2: column densith of CO and H2 in cm^-2.
Real Shielding::fShield_CO_V09(const Real NCO, const Real NH2) {
  const Real N_small_ = 1.0e10;
  Real logNCO, logNH2;
  //restrain values on the table
  if (NCO < N_small_ && NH2 < N_small_) {
    return 1.0;
  }
  if (NCO < N_small_) {
    logNCO = 10.;
    logNH2 = std::log10(NH2);
  } else if (NH2 < N_small_) {
    logNH2 = 10.;
    logNCO = std::log10(NCO);
  } else {
    logNCO = std::log10(NCO);
    logNH2 = std::log10(NH2);
  }

  int iCO0, iCO1;
  int iH20, iH21;
  //interpretation value of first, second line, and final value
  Real fl1, fl2, fl;
  //find which two points are we interpolate/extroplating
  iCO0 = Interpolation::LinearInterpIndex(len_NCO_V09_, logNCOV09_, logNCO);
  iCO1 = iCO0+1;
  iH20 = Interpolation::LinearInterpIndex(len_NH2_V09_, logNH2V09_, logNH2);
  iH21 = iH20+1;
  //linear interpretation of the rows
  fl1 = Interpolation::LinearInterp(logNCOV09_[iCO0], logNCOV09_[iCO1],
                                    std::log(ThetaV09_[iH20][iCO0]),
                                    std::log(ThetaV09_[iH20][iCO1]), logNCO);
  fl2 = Interpolation::LinearInterp(logNCOV09_[iCO0], logNCOV09_[iCO1],
                                    std::log(ThetaV09_[iH21][iCO0]),
                                    std::log(ThetaV09_[iH21][iCO1]), logNCO);
  //linear interpretation of the column
  fl = Interpolation::LinearInterp(logNH2V09_[iH20], logNH2V09_[iH21],
                    fl1, fl2, logNH2);
  return std::exp(fl);
}

//----------------------------------------------------------------------------------------
//! \fn Real Shielding::fShield_H2(const Real NH2, const Real bH2)
//! \brief H2 self shielding from Draine+Bertoldi1996
Real Shielding::fShield_H2(const Real NH2, const Real bH2) {
  const Real N_small_ = 10.;
  if (NH2 < N_small_) {
    return 1.;
  }
  const Real b5 = bH2 / 1.0e5;
  const Real x = NH2 / 5.0e14;
  Real p1, p2, term; //first and second term
  p1 = 0.965 / ( (1.+x/b5) * (1.+x/b5) );
  term = std::sqrt(1. + x);
  p2 = 0.035/term * std::exp(-8.5e-4*term);
  return p1 + p2;
}

//----------------------------------------------------------------------------------------
//! \fn Real Shielding::fShield_C(const Real NC, const Real NH2)
//! \brief CI self shielding.
Real Shielding::fShield_C(const Real NC, const Real NH2) {
  //Tielens+Hollenbach1985 (A6) and (A7). Also see Wolfire's email on April 18,
  //2016.
  Real NCi, NH2i;
  if (NC < 0) {
    NCi = 0;
  } else {
    NCi = NC;
  }
  if (NH2 < 0) {
    NH2i = 0;
  } else {
    NH2i = NH2;
  }
  const Real AH2 = 1.17e-8;
  const Real tau_H2 = 1.2e-14 * 2 * NH2i;
  const Real y = AH2 * tau_H2;
  const Real ry = std::exp(-y) / (1. + y);
  const Real rc = std::exp(-1.6e-17*NCi);
  return rc * ry;
}

//----------------------------------------------------------------------------------------
//! \fn Real Shielding::fShield_CO_C(const Real NC)
//! \brief CI shielding of CO
Real Shielding::fShield_CO_C(const Real NC) {
  if (NC < 0) {
    return 0.;
  } else {
    return std::exp(-1.6e-17*NC);
  }
}
