species = {[0]='H2', [1]='O2', [2]='OH', [3]='H2O', [4]='H', [5]='O', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='alpha-qss', eps1= 1.000000e-03, eps2= 5.000000e-04, delta= 1.000000e-10, maxIters=10},
  tightTempCoupling = true,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "H + O2 <=> OH + O",
  type = "elementary",
  frc = {model='Arrhenius', A=1.250000000000e+08, n=0.000000, C=8.202441517830e+03, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 1, 4,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 5,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[2] = {
  equation = "O + H2 <=> OH + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.960000000000e+07, n=0.000000, C=4.931529256119e+03, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[3] = {
  equation = "H2 + OH <=> H2O + H",
  type = "elementary",
  frc = {model='Arrhenius', A=2.100000000000e+07, n=0.000000, C=2.566408082266e+03, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 2,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 3, 4,},
  prodCoeffs = { 1.000000e+00, 1.000000e+00,},
}

reaction[4] = {
  equation = "O + H2O <=> OH + OH",
  type = "elementary",
  frc = {model='Arrhenius', A=5.750000000000e+07, n=0.000000, C=9.057910878586e+03, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 5,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 2,},
  prodCoeffs = { 2.000000e+00,},
}

