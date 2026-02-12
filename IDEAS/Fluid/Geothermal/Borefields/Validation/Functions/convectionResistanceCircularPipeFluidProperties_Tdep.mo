within IDEAS.Fluid.Geothermal.Borefields.Validation.Functions;
function convectionResistanceCircularPipeFluidProperties_Tdep
  "Convection resistance, with fluid properties computed from current T and p"
  extends Modelica.Icons.Function;

  // Allow selecting any compatible medium (default to your temperature-dependent water)
  replaceable package Medium = IDEAS.Media.Specialized.Water.TemperatureDependentDensity
    "Medium used to evaluate thermophysical properties";
  // Geometry of the borehole
  input Modelica.Units.SI.Height hSeg "Height of the element";
  input Modelica.Units.SI.Radius rTub "Tube radius";
  input Modelica.Units.SI.Length eTub "Tube thickness";

  // Current state (make p and T inputs so function is purely algebraic)
  input Modelica.Units.SI.Temperature T "Fluid temperature (K) at which to evaluate properties";
  input Modelica.Units.SI.Pressure    p = Medium.p_default "Pressure (Pa)";

  // Mass flow inputs (same as before)
  input Modelica.Units.SI.MassFlowRate m_flow "Mass flow rate";
  input Modelica.Units.SI.MassFlowRate m_flow_nominal "Nominal mass flow rate";

  // Outputs
  output Modelica.Units.SI.ThermalResistance RFluPip "Convection resistance";
  output Real Nu "Nusselt";
  output Modelica.Units.SI.CoefficientOfHeatTransfer h "Convective heat transfer coeff";
  output Real Re "Reynolds number";
  output Real NuTurb "Nusselt at Re=2400";

protected
  Modelica.Units.SI.Radius rTub_in = rTub - eTub "Pipe inner radius";
  Real k_coef "aux coefficient used in Re computation";
  // compute properties from current T and p
  Modelica.Units.SI.SpecificHeatCapacity cpMed;
  Modelica.Units.SI.ThermalConductivity       kMed;
  Modelica.Units.SI.DynamicViscosity          muMed;
  Modelica.Units.SI.MassFlowRate m_flow_abs;
algorithm
  // Evaluate medium properties at (p,T)
  cpMed := Medium.specificHeatCapacityCp(Medium.setState_pTX(p, T, Medium.X_default));
  kMed  := Medium.thermalConductivity(     Medium.setState_pTX(p, T, Medium.X_default));
  muMed := Medium.dynamicViscosity(        Medium.setState_pTX(p, T, Medium.X_default));

  // Convection resistance and Reynolds number (same logic as original)
  k_coef := 2/(muMed*Modelica.Constants.pi*rTub_in);
  m_flow_abs := IDEAS.Utilities.Math.Functions.spliceFunction(
                   m_flow, -m_flow, m_flow, m_flow_nominal/30);

  Re := m_flow_abs * k_coef;

  if Re >= 2400 then
    Nu := 0.023*(cpMed*muMed/kMed)^(0.35) *
          IDEAS.Utilities.Math.Functions.regNonZeroPower(x=Re, n=0.8, delta=0.01*m_flow_nominal*k_coef);
  else
    NuTurb := 0.023*(cpMed*muMed/kMed)^(0.35)*(2400)^(0.8);
    Nu := IDEAS.Utilities.Math.Functions.spliceFunction(NuTurb, 3.66, Re-(2300+2400)/2, ((2300+2400)/2)-2300);
  end if;

  h := Nu * kMed / (2*rTub_in);
  RFluPip := 1/(2*Modelica.Constants.pi*rTub_in*hSeg*h);

end convectionResistanceCircularPipeFluidProperties_Tdep;
