species = {[0]='N2', [1]='N', }
config = {
  tempLimits = {lower=300.000000, upper=50000.000000},
  odeStep = {method='rkf', errTol=1.000000e-03},
  tightTempCoupling = false,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "N2 + N2 <=> N + N + N2",
  type = "elementary",
  frc = {model='Arrhenius', A=7.000000000000e+15, n=-1.600000, C=1.132000000000e+05, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 2.000000e+00,},
  prodIdx = { 0, 1,},
  prodCoeffs = { 1.000000e+00, 2.000000e+00,},
}

reaction[2] = {
  equation = "N2 + N <=> N + N + N",
  type = "elementary",
  frc = {model='Arrhenius', A=3.000000000000e+16, n=-1.600000, C=1.132000000000e+05, rctIndex=-1},
  brc = {model='fromEqConst'},
  ec = {},
  reacIdx = { 0, 1,},
  reacCoeffs = { 1.000000e+00, 1.000000e+00,},
  prodIdx = { 1,},
  prodCoeffs = { 3.000000e+00,},
}

