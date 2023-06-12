-- H2 + O2 6sp chain reaction and dissociation/recombination mechanism 
-- Author: Luke Bartholomew
-- Date: 5/06/23
-- Place: UQ, Saint Lucia
-- Reference: Bittker-Scullin reaction mechanism model
-- url: https://onlinelibrary.wiley.com/doi/epdf/10.1002/kin.20026
--
-- History:
--   

C_factor = 8.3145 / 4.184

Config{
    odeStep = {method='alpha-qss', eps1=0.001}, --default eps1=0.001
    tightTempCoupling = true
}

Reaction{
    'H + O2 <=> OH + O',
    fr = {'Arrhenius', A=1.25e14,  n=0.0, C=16300.0 / C_factor}
}

Reaction{
    'O + H2 <=> OH + H',
    fr = {'Arrhenius', A=2.96e13,  n=0.0, C=9800.0 / C_factor}
}

Reaction{
    'H2 + OH <=> H2O + H',
    fr = {'Arrhenius', A=2.1e13,  n=0.0, C=5100.0 / C_factor}
}

Reaction{
    'O + H2O <=> OH + OH',
    fr = {'Arrhenius', A=5.75e13,  n=0.0, C=18000.0 / C_factor}
}

Reaction{
    'H2 + O2 <=> OH + OH',
    fr = {'Arrhenius', A=1.0e13, n=0.0, C=43000.0 / C_factor}
}

Reaction{
    'H + H + M <=> H2 + M',
    fr = {'Arrhenius', A=1.0e18, n=-1.0, C=0.0},
    efficiencies = {O2=2.0, H2=5.0, H2O=15.0}
}

Reaction{
    'O + O + M <=> O2 + M',
    fr = {'Arrhenius', A=1.38e18, n=-1.0, C=340.0/C_factor}
}

Reaction{
    'H + OH + M <=> H2O + M',
    fr = {'Arrhenius', A=7.5e23, n=-2.6, C=0.0,
    efficiencies={H2=4.0, H2O=20.0, O2=1.6}}
}