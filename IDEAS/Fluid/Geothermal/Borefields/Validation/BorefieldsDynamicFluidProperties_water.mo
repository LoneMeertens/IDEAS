within IDEAS.Fluid.Geothermal.Borefields.Validation;
model BorefieldsDynamicFluidProperties_water
  "Validation model of borefields with different media specifications operating simultaneously"
  extends Modelica.Icons.Example;

  package Medium = IDEAS.Media.Water(reference_T=273.15);
  package Medium1 = IDEAS.Media.Specialized.Water.TemperatureDependentDensity;
  parameter Modelica.Units.SI.Time tLoaAgg=300
    "Time resolution of load aggregation";

  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal/2;

  parameter Modelica.Units.SI.Temperature TGro=283.15 "Ground temperature";
  parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Example borFieUTubDat(
    filDat=IDEAS.Fluid.Geothermal.Borefields.Data.Filling.Bentonite(
    steadyState=true),
    conDat=IDEAS.Fluid.Geothermal.Borefields.Data.Configuration.Example(
    borCon=IDEAS.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.SingleUTube))
    annotation (Placement(transformation(extent={{-90,76},{-70,96}})));

  Real RFluPip_value, Nu_value, h_value, Re_value, NuTurb_value;
  Real RFluPip_value1, Nu_value1, h_value1, Re_value1, NuTurb_value1;

    // --- runtime medium-state & properties for Medium1 (TemperatureDependentDensity) ---
  Medium1.ThermodynamicState medState1;
  Real muMed1_dyn "dynamic viscosity (Pa.s) computed from Medium1 at runtime";
  Real kMed1_dyn  "thermal conductivity (W/m.K) computed from Medium1 at runtime";
  Real cpMed1_dyn "specific heat capacity cp (J/kg.K) computed from Medium1 at runtime";

  IDEAS.Fluid.Geothermal.Borefields.OneUTube borFieUTub(
    redeclare package Medium = Medium,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro)
    "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,50},{10,70}})));
  IDEAS.Fluid.Sources.MassFlowSource_T sou(
    redeclare package Medium = Medium,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,50},{-72,
            70}},      rotation=0)));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubIn(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,50},{-40,70}})));
  IDEAS.Fluid.Sources.Boundary_pT sin(
    redeclare package Medium = Medium,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,50},{70,70}},
                   rotation=0)));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubOut(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,50},{60,70}})));

   Modelica.Blocks.Sources.RealExpression Nu(y=Nu_value) "Nusselt number"
    annotation (Placement(transformation(extent={{128,56},{148,76}})));
  Modelica.Blocks.Sources.RealExpression RFluPip(y=RFluPip_value)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,40},{124,60}})));
  Modelica.Blocks.Sources.RealExpression h(y=h_value)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,40},{148,60}})));
  Modelica.Blocks.Sources.RealExpression Re(y=Re_value) "Reynolds number"
    annotation (Placement(transformation(extent={{104,56},{124,76}})));
  Modelica.Blocks.Sources.RealExpression NuTurb(y=NuTurb_value)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,56},{172,76}})));
   Modelica.Blocks.Sources.RealExpression Nu1(y=Nu_value1)
                                                         "Nusselt number"
    annotation (Placement(transformation(extent={{128,2},{148,22}})));
  Modelica.Blocks.Sources.RealExpression RFluPip1(y=RFluPip_value1)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,-14},{124,6}})));
  Modelica.Blocks.Sources.RealExpression h1(y=h_value1)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,-14},{148,6}})));
  Modelica.Blocks.Sources.RealExpression Re1(y=Re_value1)
                                                        "Reynolds number"
    annotation (Placement(transformation(extent={{104,2},{124,22}})));
  Modelica.Blocks.Sources.RealExpression NuTurb1(y=NuTurb_value1)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,2},{172,22}})));
  OneUTube                                   borFieUTub1(
    redeclare package Medium = Medium1,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro)
    "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,6},{10,26}})));
  Sources.MassFlowSource_T             sou1(
    redeclare package Medium = Medium1,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,6},{-72,
            26}},      rotation=0)));
  Sensors.TemperatureTwoPort             TUTubIn1(
    redeclare package Medium = Medium1,
    m_flow_nominal=m_flow_nominal,
    tau=0)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,6},{-40,26}})));
  Sources.Boundary_pT             sin1(
    redeclare package Medium = Medium1,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,6},{70,26}},
                   rotation=0)));
  Sensors.TemperatureTwoPort             TUTubOut1(
    redeclare package Medium = Medium1,
    m_flow_nominal=m_flow_nominal,
    tau=0)
    "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,6},{60,26}})));
equation
  (RFluPip_value, Nu_value, h_value, Re_value, NuTurb_value) = IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borFieUTub.borHol.intHex[1].hSeg,
    rTub=borFieUTub.borFieDat.conDat.rTub,
    eTub=borFieUTub.borFieDat.conDat.eTub,
    kMed=borFieUTub.borHol.intHex[1].kMed,
    muMed=borFieUTub.borHol.intHex[1].muMed,
    cpMed=borFieUTub.borHol.intHex[1].cpMed,
    m_flow=borFieUTub.borHol.intHex[1].m1_flow,
    m_flow_nominal=borFieUTub.borHol.intHex[1].m1_flow_nominal);

  // --- compute Medium1 properties from current outlet temperature (K) ---
  medState1 = Medium1.setState_pTX(p=101330, T = TUTubOut1.T, X = Medium1.X_default);

  muMed1_dyn = Medium1.dynamicViscosity(medState1);          // Pa.s
  kMed1_dyn  = Medium1.thermalConductivity(medState1);       // W/m.K
  cpMed1_dyn = Medium1.specificHeatCapacityCp(medState1);    // J/kg.K

  // --- call convection function with runtime properties and live mass flow ---
  (RFluPip_value1, Nu_value1, h_value1, Re_value1, NuTurb_value1) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
      hSeg = borFieUTub1.borHol.intHex[1].hSeg,
      rTub = borFieUTub1.borFieDat.conDat.rTub,
      eTub = borFieUTub1.borFieDat.conDat.eTub,
      kMed = kMed1_dyn,
      muMed = muMed1_dyn,
      cpMed = cpMed1_dyn,
      m_flow = borFieUTub1.borHol.intHex[1].m1_flow,
      m_flow_nominal = m_flow_nominal);


  connect(TUTubIn.port_b, borFieUTub.port_a)
    annotation (Line(points={{-40,60},{-10,60}},   color={0,127,255}));
  connect(borFieUTub.port_b, TUTubOut.port_a)
    annotation (Line(points={{10,60},{40,60}},            color={0,127,255}));
  connect(TUTubOut.port_b, sin.ports[1])
    annotation (Line(points={{60,60},{70,60}},            color={0,127,255}));

  connect(TUTubIn1.port_b, borFieUTub1.port_a)
    annotation (Line(points={{-40,16},{-10,16}}, color={0,127,255}));
  connect(borFieUTub1.port_b, TUTubOut1.port_a)
    annotation (Line(points={{10,16},{40,16}}, color={0,127,255}));
  connect(TUTubOut1.port_b, sin1.ports[1])
    annotation (Line(points={{60,16},{70,16}}, color={0,127,255}));
  connect(sou.ports[1], TUTubIn.port_a)
    annotation (Line(points={{-72,60},{-60,60}}, color={0,127,255}));
  connect(TUTubIn1.port_a, sou1.ports[1])
    annotation (Line(points={{-60,16},{-72,16}}, color={0,127,255}));
  annotation (__Dymola_Commands(file="modelica://IDEAS/Resources/Scripts/Dymola/Fluid/Geothermal/Borefields/Examples/Borefields.mos"
        "Simulate and plot"),
  Documentation(info="<html>
<p>
This example shows three different borefields, each with a different configuration
(single U-tube, double U-tube in parallel, and double U-tube in series) and compares
the thermal behaviour of the circulating fluid in each case.
</p>
</html>",
revisions="<html>
<ul>
<li>
May 17, 2024, by Michael Wetter:<br/>
Updated model due to removal of parameter <code>dynFil</code>.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1885\">IBPSA, #1885</a>.
</li>
<li>
April 8, 2021, by Michael Wetter:<br/>
Added missing <code>parameter</code> keyword.<br/>
For <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1464\">IBPSA, issue 1464</a>.
</li>
<li>
June 2018, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</html>"),
    experiment(
      StopTime=36000,Tolerance=1e-6),
    Diagram(coordinateSystem(extent={{-100,-100},{180,100}}), graphics={
          Rectangle(extent={{98,80},{176,42}}, lineColor={28,108,200}),
          Rectangle(extent={{98,26},{176,-12}},lineColor={28,108,200})}),
    Icon(coordinateSystem(extent={{-100,-100},{180,100}})));
end BorefieldsDynamicFluidProperties_water;
