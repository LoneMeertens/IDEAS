within IDEAS.Fluid.Geothermal.Borefields.Validation;
model BorefieldsDynamicFluidProperties_water_injection
  "Validation model of borefields with different media specifications operating simultaneously"
  extends Modelica.Icons.Example;


  package Medium10 = IDEAS.Media.Specialized.Water.ConstantProperties_pT(T_nominal=273.15+10);
  package Medium15 = IDEAS.Media.Specialized.Water.ConstantProperties_pT(T_nominal=273.15+15);
  package Medium20 = IDEAS.Media.Specialized.Water.ConstantProperties_pT(T_nominal=273.15+20);
  package Medium30 = IDEAS.Media.Specialized.Water.ConstantProperties_pT(T_nominal=273.15+30);
  package Mediumdyn = IDEAS.Media.Specialized.Water.TemperatureDependentDensity;
  parameter Modelica.Units.SI.Time tLoaAgg=300
    "Time resolution of load aggregation";

  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal/4;

  parameter Modelica.Units.SI.Temperature TGro=283.15 "Ground temperature";
  parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Example borFieUTubDat(
    filDat=IDEAS.Fluid.Geothermal.Borefields.Data.Filling.Bentonite(
    steadyState=true),
    conDat=IDEAS.Fluid.Geothermal.Borefields.Data.Configuration.Example(
    borCon=IDEAS.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.SingleUTube))
    annotation (Placement(transformation(extent={{-90,76},{-70,96}})));

  Real RFluPip_value30, Nu_value30, h_value30, Re_value30, NuTurb_value30;
  Real RFluPip_value10, Nu_value10, h_value10, Re_value10, NuTurb_value10;
  Real RFluPip_value15, Nu_value15, h_value15, Re_value15, NuTurb_value15;
  Real RFluPip_value20, Nu_value20, h_value20, Re_value20, NuTurb_value20;
  Real RFluPip_value_dyn10, Nu_value_dyn10, h_value_dyn10, Re_value_dyn10, NuTurb_value_dyn10;
  Real RFluPip_value_dyn30, Nu_value_dyn30, h_value_dyn30, Re_value_dyn30, NuTurb_value_dyn30;

   // --- runtime medium-state & properties for Medium1 (TemperatureDependentDensity) ---
  Medium10.ThermodynamicState medState10;
  Real muMed10_dyn "dynamic viscosity (Pa.s) computed from Medium1 at runtime";
  Real kMed10_dyn  "thermal conductivity (W/m.K) computed from Medium1 at runtime";
  Real cpMed10_dyn "specific heat capacity cp (J/kg.K) computed from Medium1 at runtime";

   // --- runtime medium-state & properties for Medium1 (TemperatureDependentDensity) ---
  Medium30.ThermodynamicState medState30;
  Real muMed30_dyn "dynamic viscosity (Pa.s) computed from Medium1 at runtime";
  Real kMed30_dyn  "thermal conductivity (W/m.K) computed from Medium1 at runtime";
  Real cpMed30_dyn "specific heat capacity cp (J/kg.K) computed from Medium1 at runtime";

  IDEAS.Fluid.Geothermal.Borefields.OneUTube borFieUTub30(
    redeclare package Medium = Medium30,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
  IDEAS.Fluid.Sources.MassFlowSource_T sou30(
    redeclare package Medium = Medium30,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,-70},{-72,
            -50}}, rotation=0)));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubIn30(
    redeclare package Medium = Medium30,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,-70},{-40,-50}})));
  IDEAS.Fluid.Sources.Boundary_pT sin30(
    redeclare package Medium = Medium30,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,-70},{70,-50}},
          rotation=0)));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubOut30(
    redeclare package Medium = Medium30,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,-70},{60,-50}})));

   Modelica.Blocks.Sources.RealExpression Nu30(y=Nu_value30) "Nusselt number"
    annotation (Placement(transformation(extent={{128,-64},{148,-44}})));
  Modelica.Blocks.Sources.RealExpression RFluPip30(y=RFluPip_value30)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,-80},{124,-60}})));
  Modelica.Blocks.Sources.RealExpression h30(y=h_value30)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,-80},{148,-60}})));
  Modelica.Blocks.Sources.RealExpression Re30(y=Re_value30) "Reynolds number"
    annotation (Placement(transformation(extent={{104,-64},{124,-44}})));
  Modelica.Blocks.Sources.RealExpression NuTurb30(y=NuTurb_value30)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,-64},{172,-44}})));
   Modelica.Blocks.Sources.RealExpression Nu10(y=Nu_value10)
                                                         "Nusselt number"
    annotation (Placement(transformation(extent={{128,64},{148,84}})));
  Modelica.Blocks.Sources.RealExpression RFluPip10(y=RFluPip_value10)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,48},{124,68}})));
  Modelica.Blocks.Sources.RealExpression h10(y=h_value10)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,48},{148,68}})));
  Modelica.Blocks.Sources.RealExpression Re10(y=Re_value10)
                                                        "Reynolds number"
    annotation (Placement(transformation(extent={{104,64},{124,84}})));
  Modelica.Blocks.Sources.RealExpression NuTurb10(y=NuTurb_value10)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,64},{172,84}})));
  OneUTube borFieUTub10(
    redeclare package Medium = Medium10,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,50},{10,70}})));
  Sources.MassFlowSource_T sou10(
    redeclare package Medium = Medium10,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,50},{-72,
            70}}, rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn10(
    redeclare package Medium = Medium10,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,50},{-40,70}})));
  Sources.Boundary_pT sin10(
    redeclare package Medium = Medium10,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,50},{70,70}},
          rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut10(
    redeclare package Medium = Medium10,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,50},{60,70}})));
   Modelica.Blocks.Sources.RealExpression Nu15(y=Nu_value15)
                                                         "Nusselt number"
    annotation (Placement(transformation(extent={{128,22},{148,42}})));
  Modelica.Blocks.Sources.RealExpression RFluPip15(y=RFluPip_value15)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,6},{124,26}})));
  Modelica.Blocks.Sources.RealExpression h15(y=h_value15)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,6},{148,26}})));
  Modelica.Blocks.Sources.RealExpression Re15(y=Re_value15)
                                                        "Reynolds number"
    annotation (Placement(transformation(extent={{104,22},{124,42}})));
  Modelica.Blocks.Sources.RealExpression NuTurb15(y=NuTurb_value15)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,22},{172,42}})));
  OneUTube borFieUTub15(
    redeclare package Medium = Medium15,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,10},{10,30}})));
  Sources.MassFlowSource_T sou15(
    redeclare package Medium = Medium15,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,10},{-72,
            30}},  rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn15(
    redeclare package Medium = Medium15,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,10},{-40,30}})));
  Sources.Boundary_pT sin15(
    redeclare package Medium = Medium15,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,10},{70,30}},
          rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut15(
    redeclare package Medium = Medium15,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,10},{60,30}})));
   Modelica.Blocks.Sources.RealExpression Nu20(y=Nu_value20)
                                                         "Nusselt number"
    annotation (Placement(transformation(extent={{128,-20},{148,0}})));
  Modelica.Blocks.Sources.RealExpression RFluPip20(y=RFluPip_value20)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{104,-36},{124,-16}})));
  Modelica.Blocks.Sources.RealExpression h20(y=h_value20)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{128,-36},{148,-16}})));
  Modelica.Blocks.Sources.RealExpression Re20(y=Re_value20)
                                                        "Reynolds number"
    annotation (Placement(transformation(extent={{104,-20},{124,0}})));
  Modelica.Blocks.Sources.RealExpression NuTurb20(y=NuTurb_value20)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{152,-20},{172,0}})));
  OneUTube borFieUTub20(
    redeclare package Medium = Medium20,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
  Sources.MassFlowSource_T sou20(
    redeclare package Medium = Medium20,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-92,-30},{-72,
            -10}}, rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn20(
    redeclare package Medium = Medium20,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
  Sources.Boundary_pT sin20(
    redeclare package Medium = Medium20,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{90,-30},{70,-10}},
          rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut20(
    redeclare package Medium = Medium20,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{40,-30},{60,-10}})));
equation
  (RFluPip_value30, Nu_value30, h_value30, Re_value30, NuTurb_value30) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borFieUTub30.borHol.intHex[1].hSeg,
    rTub=borFieUTub30.borFieDat.conDat.rTub,
    eTub=borFieUTub30.borFieDat.conDat.eTub,
    kMed=borFieUTub30.borHol.intHex[1].kMed,
    muMed=borFieUTub30.borHol.intHex[1].muMed,
    cpMed=borFieUTub30.borHol.intHex[1].cpMed,
    m_flow=borFieUTub30.borHol.intHex[1].m1_flow,
    m_flow_nominal=borFieUTub30.borHol.intHex[1].m1_flow_nominal);

  (RFluPip_value10, Nu_value10, h_value10, Re_value10, NuTurb_value10) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borFieUTub10.borHol.intHex[1].hSeg,
    rTub=borFieUTub10.borFieDat.conDat.rTub,
    eTub=borFieUTub10.borFieDat.conDat.eTub,
    kMed=borFieUTub10.borHol.intHex[1].kMed,
    muMed=borFieUTub10.borHol.intHex[1].muMed,
    cpMed=borFieUTub10.borHol.intHex[1].cpMed,
    m_flow=borFieUTub10.borHol.intHex[1].m1_flow,
    m_flow_nominal=borFieUTub10.borHol.intHex[1].m1_flow_nominal);

  (RFluPip_value15, Nu_value15, h_value15, Re_value15, NuTurb_value15) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borFieUTub15.borHol.intHex[1].hSeg,
    rTub=borFieUTub15.borFieDat.conDat.rTub,
    eTub=borFieUTub15.borFieDat.conDat.eTub,
    kMed=borFieUTub15.borHol.intHex[1].kMed,
    muMed=borFieUTub15.borHol.intHex[1].muMed,
    cpMed=borFieUTub15.borHol.intHex[1].cpMed,
    m_flow=borFieUTub15.borHol.intHex[1].m1_flow,
    m_flow_nominal=borFieUTub15.borHol.intHex[1].m1_flow_nominal);

   (RFluPip_value20, Nu_value20, h_value20, Re_value20, NuTurb_value20) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borFieUTub20.borHol.intHex[1].hSeg,
    rTub=borFieUTub20.borFieDat.conDat.rTub,
    eTub=borFieUTub20.borFieDat.conDat.eTub,
    kMed=borFieUTub20.borHol.intHex[1].kMed,
    muMed=borFieUTub20.borHol.intHex[1].muMed,
    cpMed=borFieUTub20.borHol.intHex[1].cpMed,
    m_flow=borFieUTub20.borHol.intHex[1].m1_flow,
    m_flow_nominal=borFieUTub20.borHol.intHex[1].m1_flow_nominal);

    // --- compute dynamic properties from current outlet temperature (K) ---
  medState10 = Mediumdyn.setState_pTX(p=101330, T = TUTubOut10.T, X = Mediumdyn.X_default);

  muMed10_dyn = Mediumdyn.dynamicViscosity(medState10);          // Pa.s
  kMed10_dyn  = Mediumdyn.thermalConductivity(medState10);       // W/m.K
  cpMed10_dyn = Mediumdyn.specificHeatCapacityCp(medState10);    // J/kg.K

  medState30 = Mediumdyn.setState_pTX(p=101330, T = TUTubOut30.T, X = Mediumdyn.X_default);

  muMed30_dyn = Mediumdyn.dynamicViscosity(medState30);          // Pa.s
  kMed30_dyn  = Mediumdyn.thermalConductivity(medState30);       // W/m.K
  cpMed30_dyn = Mediumdyn.specificHeatCapacityCp(medState30);    // J/kg.K

  // --- call convection function with runtime properties and live mass flow ---
  (RFluPip_value_dyn10, Nu_value_dyn10, h_value_dyn10, Re_value_dyn10, NuTurb_value_dyn10) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
      hSeg = borFieUTub10.borHol.intHex[1].hSeg,
      rTub = borFieUTub10.borFieDat.conDat.rTub,
      eTub = borFieUTub10.borFieDat.conDat.eTub,
      kMed = kMed10_dyn,
      muMed = muMed10_dyn,
      cpMed = cpMed10_dyn,
      m_flow = borFieUTub10.borHol.intHex[1].m1_flow,
      m_flow_nominal = m_flow_nominal);

   (RFluPip_value_dyn30, Nu_value_dyn30, h_value_dyn30, Re_value_dyn30, NuTurb_value_dyn30) =
    IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
      hSeg = borFieUTub30.borHol.intHex[1].hSeg,
      rTub = borFieUTub30.borFieDat.conDat.rTub,
      eTub = borFieUTub30.borFieDat.conDat.eTub,
      kMed = kMed30_dyn,
      muMed = muMed30_dyn,
      cpMed = cpMed30_dyn,
      m_flow = borFieUTub30.borHol.intHex[1].m1_flow,
      m_flow_nominal = m_flow_nominal);

  connect(TUTubIn30.port_b, borFieUTub30.port_a)
    annotation (Line(points={{-40,-60},{-10,-60}}, color={0,127,255}));
  connect(borFieUTub30.port_b, TUTubOut30.port_a)
    annotation (Line(points={{10,-60},{40,-60}}, color={0,127,255}));
  connect(TUTubOut30.port_b, sin30.ports[1])
    annotation (Line(points={{60,-60},{70,-60}}, color={0,127,255}));

  connect(TUTubIn10.port_b, borFieUTub10.port_a)
    annotation (Line(points={{-40,60},{-10,60}}, color={0,127,255}));
  connect(borFieUTub10.port_b, TUTubOut10.port_a)
    annotation (Line(points={{10,60},{40,60}}, color={0,127,255}));
  connect(TUTubOut10.port_b, sin10.ports[1])
    annotation (Line(points={{60,60},{70,60}}, color={0,127,255}));
  connect(sou30.ports[1], TUTubIn30.port_a)
    annotation (Line(points={{-72,-60},{-60,-60}}, color={0,127,255}));
  connect(TUTubIn10.port_a, sou10.ports[1])
    annotation (Line(points={{-60,60},{-72,60}}, color={0,127,255}));
  connect(TUTubIn15.port_b, borFieUTub15.port_a)
    annotation (Line(points={{-40,20},{-10,20}},   color={0,127,255}));
  connect(borFieUTub15.port_b, TUTubOut15.port_a)
    annotation (Line(points={{10,20},{40,20}},   color={0,127,255}));
  connect(TUTubOut15.port_b, sin15.ports[1])
    annotation (Line(points={{60,20},{70,20}},   color={0,127,255}));
  connect(TUTubIn15.port_a, sou15.ports[1])
    annotation (Line(points={{-60,20},{-72,20}},   color={0,127,255}));
  connect(TUTubIn20.port_b, borFieUTub20.port_a)
    annotation (Line(points={{-40,-20},{-10,-20}}, color={0,127,255}));
  connect(borFieUTub20.port_b, TUTubOut20.port_a)
    annotation (Line(points={{10,-20},{40,-20}}, color={0,127,255}));
  connect(TUTubOut20.port_b, sin20.ports[1])
    annotation (Line(points={{60,-20},{70,-20}}, color={0,127,255}));
  connect(TUTubIn20.port_a, sou20.ports[1])
    annotation (Line(points={{-60,-20},{-72,-20}}, color={0,127,255}));
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
          Rectangle(extent={{98,-40},{176,-78}},
                                               lineColor={28,108,200}),
          Rectangle(extent={{98,88},{176,50}}, lineColor={28,108,200}),
          Rectangle(extent={{98,46},{176,8}},  lineColor={28,108,200}),
          Rectangle(extent={{98,4},{176,-34}}, lineColor={28,108,200})}),
    Icon(coordinateSystem(extent={{-100,-100},{180,100}})));
end BorefieldsDynamicFluidProperties_water_injection;
