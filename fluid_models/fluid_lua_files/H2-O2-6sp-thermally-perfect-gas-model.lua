-- Auto-generated by prep-gas on: 04-Apr-2023 13:02:26

model = 'CompositeGas'
species = {'H2', 'O2', 'OH', 'H2O', 'H', 'O', }

physical_model = 'thermally-perfect-gas'
energyModes = {'equilibrium'}
db = {}
db['H2'] = {}
db['H2'].atomicConstituents = { H=2, }
db['H2'].charge = 0
db['H2'].M = 2.01588000e-03
db['H2'].sigma = 2.92000000
db['H2'].epsilon = 38.00000000
db['H2'].Lewis = 0.31700000
db['H2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
    4.078323210e+04,
   -8.009186040e+02,
    8.214702010e+00,
   -1.269714457e-02,
    1.753605076e-05,
   -1.202860270e-08,
    3.368093490e-12,
    2.682484665e+03,
   -3.043788844e+01,
  },
  segment1 = {
    5.608128010e+05,
   -8.371504740e+02,
    2.975364532e+00,
    1.252249124e-03,
   -3.740716190e-07,
    5.936625200e-11,
   -3.606994100e-15,
    5.339824410e+03,
   -2.202774769e+00,
  },
  segment2 = {
    4.966884120e+08,
   -3.147547149e+05,
    7.984121880e+01,
   -8.414789210e-03,
    4.753248350e-07,
   -1.371873492e-11,
    1.605461756e-16,
    2.488433516e+06,
   -6.695728110e+02,
  },
}
db['H2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.4553182e-01,
      B =  4.3555109e+01,
      C = -3.2579340e+03,
      D =  1.3556243e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.6730605e-01,
      B =  6.7931897e+02,
      C = -2.1025179e+05,
      D = -1.8251697e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  1.0126129e+00,
      B =  1.4973739e+03,
      C = -1.4428484e+06,
      D = -2.3254928e+00,
   },
}
db['H2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  1.0059461e+00,
      B =  2.7951262e+02,
      C = -2.9792018e+04,
      D =  1.1996252e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  1.0582450e+00,
      B =  2.4875372e+02,
      C =  1.1736907e+04,
      D =  8.2758695e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -2.2364420e-01,
      B = -6.9650442e+03,
      C = -7.7771313e+04,
      D =  1.3189369e+01,
   },
}
db['O2'] = {}
db['O2'].atomicConstituents = { O=2, }
db['O2'].charge = 0
db['O2'].M = 3.19988000e-02
db['O2'].sigma = 3.45800000
db['O2'].epsilon = 107.40000000
db['O2'].Lewis = 1.08600000
db['O2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
   -3.425563420e+04,
    4.847000970e+02,
    1.119010961e+00,
    4.293889240e-03,
   -6.836300520e-07,
   -2.023372700e-09,
    1.039040018e-12,
   -3.391454870e+03,
    1.849699470e+01,
  },
  segment1 = {
   -1.037939022e+06,
    2.344830282e+03,
    1.819732036e+00,
    1.267847582e-03,
   -2.188067988e-07,
    2.053719572e-11,
   -8.193467050e-16,
   -1.689010929e+04,
    1.738716506e+01,
  },
  segment2 = {
    4.975294300e+08,
   -2.866106874e+05,
    6.690352250e+01,
   -6.169959020e-03,
    3.016396027e-07,
   -7.421416600e-12,
    7.278175770e-17,
    2.293554027e+06,
   -5.530621610e+02,
  },
}
db['O2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01,
   },
}
db['O2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01,
   },
}
db['OH'] = {}
db['OH'].atomicConstituents = { H=1, O=1, }
db['OH'].charge = 0
db['OH'].M = 1.70073400e-02
db['OH'].sigma = 2.75000000
db['OH'].epsilon = 80.00000000
db['OH'].Lewis = 0.73600000
db['OH'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
   -1.998858990e+03,
    9.300136160e+01,
    3.050854229e+00,
    1.529529288e-03,
   -3.157890998e-06,
    3.315446180e-09,
   -1.138762683e-12,
    2.991214235e+03,
    4.674110790e+00,
  },
  segment1 = {
    1.017393379e+06,
   -2.509957276e+03,
    5.116547860e+00,
    1.305299930e-04,
   -8.284322260e-08,
    2.006475941e-11,
   -1.556993656e-15,
    2.019640206e+04,
   -1.101282337e+01,
  },
  segment2 = {
    2.847234193e+08,
   -1.859532612e+05,
    5.008240900e+01,
   -5.142374980e-03,
    2.875536589e-07,
   -8.228817960e-12,
    9.567229020e-17,
    1.468393908e+06,
   -4.023555580e+02,
  },
}
db['OH'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  5.9711536e-01,
      B = -4.6100678e+02,
      C =  3.7606286e+04,
      D =  2.4041761e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4287721e-01,
      B = -1.8173747e+02,
      C = -8.8543767e+04,
      D =  1.9636057e+00,
   },
}
db['OH'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  6.8627561e-01,
      B = -7.4033274e+02,
      C =  2.7559033e+04,
      D =  2.8308741e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.7918112e-01,
      B = -9.3769908e+03,
      C =  7.0509952e+06,
      D =  1.4203688e+01,
   },
}
db['H2O'] = {}
db['H2O'].atomicConstituents = { H=2, O=1, }
db['H2O'].charge = 0
db['H2O'].M = 1.80152800e-02
db['H2O'].sigma = 2.60500000
db['H2O'].epsilon = 572.40000000
db['H2O'].Lewis = 0.85400000
db['H2O'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 2, 
  T_break_points = { 200.00, 1000.00, 6000.00, },
  T_blend_ranges = { 400.0, },
  segment0 = {
   -3.947960830e+04,
    5.755731020e+02,
    9.317826530e-01,
    7.222712860e-03,
   -7.342557370e-06,
    4.955043490e-09,
   -1.336933246e-12,
   -3.303974310e+04,
    1.724205775e+01,
  },
  segment1 = {
    1.034972096e+06,
   -2.412698562e+03,
    4.646110780e+00,
    2.291998307e-03,
   -6.836830480e-07,
    9.426468930e-11,
   -4.822380530e-15,
   -1.384286509e+04,
   -7.978148510e+00,
  },
}
db['H2O'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  5.0019557e-01,
      B = -6.9712796e+02,
      C =  8.8163892e+04,
      D =  3.0836508e+00,
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  5.8988538e-01,
      B = -5.3769814e+02,
      C =  5.4263513e+04,
      D =  2.3386375e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  6.4330087e-01,
      B = -9.5668913e+01,
      C = -3.7742283e+05,
      D =  1.8125190e+00,
   },
}
db['H2O'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 373.2,
      T_upper = 1073.2,
      A =  1.0966389e+00,
      B = -5.5513429e+02,
      C =  1.0623408e+05,
      D = -2.4664550e-01,
   },
   segment1 = {
      T_lower = 1073.2,
      T_upper = 5000.0,
      A =  3.9367933e-01,
      B = -2.2524226e+03,
      C =  6.1217458e+05,
      D =  5.8011317e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -4.1858737e-01,
      B = -1.4096649e+04,
      C =  1.9179190e+07,
      D =  1.4345613e+01,
   },
}
db['H'] = {}
db['H'].atomicConstituents = { H=1, }
db['H'].charge = 0
db['H'].M = 1.00794000e-03
db['H'].sigma = 2.05000000
db['H'].epsilon = 145.00000000
db['H'].Lewis = 0.18900000
db['H'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
    0.000000000e+00,
    0.000000000e+00,
    2.500000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    0.000000000e+00,
    2.547370801e+04,
   -4.466828530e-01,
  },
  segment1 = {
    6.078774250e+01,
   -1.819354417e-01,
    2.500211817e+00,
   -1.226512864e-07,
    3.732876330e-11,
   -5.687744560e-15,
    3.410210197e-19,
    2.547486398e+04,
   -4.481917770e-01,
  },
  segment2 = {
    2.173757694e+08,
   -1.312035403e+05,
    3.399174200e+01,
   -3.813999680e-03,
    2.432854837e-07,
   -7.694275540e-12,
    9.644105630e-17,
    1.067638086e+06,
   -2.742301051e+02,
  },
}
db['H'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4226149e-01,
      B = -4.0132865e+02,
      C =  1.8554165e+05,
      D =  4.6741844e-02,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7486623e-01,
      B = -2.5022902e+03,
      C =  7.0955048e+06,
      D = -9.3888455e-01,
   },
}
db['H'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.4166119e-01,
      B = -4.0487203e+02,
      C =  1.8775642e+05,
      D =  3.4843121e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7447639e-01,
      B = -2.5089452e+03,
      C =  7.1081294e+06,
      D =  2.4970991e+00,
   },
}
db['O'] = {}
db['O'].atomicConstituents = { O=1, }
db['O'].charge = 0
db['O'].M = 1.59994000e-02
db['O'].sigma = 2.75000000
db['O'].epsilon = 80.00000000
db['O'].Lewis = 0.71200000
db['O'].thermoCoeffs = {
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
    2.927275780e+04,
    5.204407390e+00,
  },
  segment1 = {
   -4.594930000e+05,
    1.360042990e+03,
    9.307798890e-01,
    9.081572960e-04,
   -2.817466540e-07,
    4.494516750e-11,
   -2.685644900e-15,
    2.061677860e+04,
    1.639289840e+01,
  },
  segment2 = {
   -1.361804880e+09,
    9.176836840e+05,
   -2.468849420e+02,
    3.485201650e-02,
   -2.617134300e-06,
    9.915441590e-11,
   -1.441920000e-15,
   -7.153289770e+06,
    2.140150020e+03,
  },
  segment3 = {
    4.381213090e+11,
   -8.983106230e+07,
    7.336265460e+03,
   -3.053691850e-01,
    6.903180540e-06,
   -8.092578220e-11,
    3.863224590e-16,
    8.105768910e+08,
   -7.164011440e+04,
  },
}
db['O'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7269241e-01,
      B =  8.3842977e+01,
      C = -5.8502098e+04,
      D =  8.5100827e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7669586e-01,
      B =  1.0158420e+03,
      C = -1.0884566e+06,
      D = -1.8001077e-01,
   },
}
db['O'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7271664e-01,
      B =  8.3989100e+01,
      C = -5.8580966e+04,
      D =  1.5179900e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7676666e-01,
      B =  1.0170744e+03,
      C = -1.0906690e+06,
      D =  4.8644232e-01,
   },
}
