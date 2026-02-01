within IDEAS.Fluid.Geothermal.Borefields.Validation;
model ConvectionCalculatorDynamic
  // Change Medium to the medium package you want to use (IDEAS.Media.Water or TemperatureDependentDensity)
  replaceable package Medium = IDEAS.Media.Specialized.Water.TemperatureDependentDensity;

  // Inputs you must connect to live signals:
  input Real T_fluid "Bulk fluid temperature (K) — connect a temperature sensor output";
  input Real m_flow "instantaneous mass flow (kg/s)";
  input Real m_flow_nominal "nominal mass flow (kg/s)";

  // geometry / bore info (set as parameters or connect if dynamic)
  parameter Real hSeg = 1 "segment height (m)"; // set to your actual value per segment
  parameter Real rTub = 0.01 "tube radius (m)"; // your rTub
  parameter Real eTub = 0.002 "tube wall thickness (m)"; // your eTub

  // outputs (time varying)
  output Real RFluPip;
  output Real Nu;
  output Real h "convective heat transfer coeff (W/m2.K)";
  output Real Re;
  output Real NuTurb;

protected
  Medium.ThermodynamicState medState;
  Real muMed;
  Real kMed;
  Real cpMed;
  parameter Real p_nominal = 101330 "Use nominal pressure if actual is not available";

equation
  // Build a medium state from current T and nominal p
  medState = Medium.setState_pTX(p = p_nominal, T = T_fluid);

  // Get dynamic properties from medium at current state
  muMed = Medium.dynamicViscosity(medState); // dynamic viscosity [Pa.s]
  kMed  = Medium.thermalConductivity(medState); // thermal conductivity [W/m.K]
  cpMed = Medium.specificHeatCapacityCp(medState); // specific heat [J/kg.K]

  // Call the convection function (dynamic evaluation)
  (RFluPip, Nu, h, Re, NuTurb) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
      hSeg = hSeg,
      rTub = rTub,
      eTub = eTub,
      kMed = kMed,
      muMed = muMed,
      cpMed = cpMed,
      m_flow = m_flow,
      m_flow_nominal = m_flow_nominal);
end ConvectionCalculatorDynamic;
