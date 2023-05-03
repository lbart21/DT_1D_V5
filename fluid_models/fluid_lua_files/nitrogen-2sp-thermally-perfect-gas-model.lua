-- Auto-generated by prep-gas on: 04-Apr-2023 13:03:48

model = 'CompositeGas'
species = {'N2', 'N', }

physical_model = 'thermally-perfect-gas'
energyModes = {'equilibrium'}
db = {}
db['N2'] = {}
db['N2'].atomicConstituents = { N=2, }
db['N2'].charge = 0
db['N2'].M = 2.80134000e-02
db['N2'].sigma = 3.62100000
db['N2'].epsilon = 97.53000000
db['N2'].Lewis = 1.15200000
db['N2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
    2.210371497e+04,
   -3.818461820e+02,
    6.082738360e+00,
   -8.530914410e-03,
    1.384646189e-05,
   -9.625793620e-09,
    2.519705809e-12,
    7.108460860e+02,
   -1.076003744e+01,
  },
  segment1 = {
    5.877124060e+05,
   -2.239249073e+03,
    6.066949220e+00,
   -6.139685500e-04,
    1.491806679e-07,
   -1.923105485e-11,
    1.061954386e-15,
    1.283210415e+04,
   -1.586640027e+01,
  },
  segment2 = {
    8.310139160e+08,
   -6.420733540e+05,
    2.020264635e+02,
   -3.065092046e-02,
    2.486903333e-06,
   -9.705954110e-11,
    1.437538881e-15,
    4.938707040e+06,
   -1.672099740e+03,
  },
}
db['N2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['N2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}
db['N'] = {}
db['N'].atomicConstituents = { N=1, }
db['N'].charge = 0
db['N'].M = 1.40067000e-02
db['N'].sigma = 3.29800000
db['N'].epsilon = 71.40000000
db['N'].Lewis = 1.15200000
db['N'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 4, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, 50000.00, },
  T_blend_ranges = { 400.0, 1000.0, 1000.0, },
  segment0 = {
    0.000000000e+00,
    0.000000000e+00,
    2.500000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    5.610610630e+04,
    4.194251390e+00,
  },
  segment1 = {
   -2.270732770e+05,
    8.140529440e+02,
    1.327051370e+00,
    8.627217310e-04,
   -3.357470890e-07,
    6.290106870e-11,
   -3.906745870e-15,
    5.109431410e+04,
    1.228237130e+01,
  },
  segment2 = {
   -2.047389940e+09,
    1.458428470e+06,
   -4.188338240e+02,
    6.259944070e-02,
   -4.965428220e-06,
    1.982470520e-10,
   -3.054701940e-15,
   -1.127277300e+07,
    3.584874170e+03,
  },
  segment3 = {
    5.742919020e+11,
   -1.290392940e+08,
    1.153814670e+04,
   -5.250785680e-01,
    1.292190900e-05,
   -1.639742310e-10,
    8.418785850e-16,
    1.152618360e+09,
   -1.116492320e+05,
  },
}
db['N'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3724737e-01,
      B =  4.3997150e+02,
      C = -1.7450753e+05,
      D =  1.0365689e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.9986588e-01,
      B =  1.4112801e+03,
      C = -1.8200478e+06,
      D = -5.5811716e-01,
   },
}
db['N'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3771661e-01,
      B =  4.4243270e+02,
      C = -1.7578446e+05,
      D =  8.9942915e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  9.0001710e-01,
      B =  1.4141175e+03,
      C = -1.8262403e+06,
      D =  2.4048513e-01,
   },
}