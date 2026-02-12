within IDEAS.Fluid.Geothermal.Borefields.Validation;
model BorefieldsDynamicFluidProperties_water_injection_var
  "Validation model of borefields with different media specifications operating simultaneously"
  extends Modelica.Icons.Example;

  package MediumTdep = IDEAS.Media.Specialized.Water.TemperatureDependentDensity;
  package Medium = IDEAS.Media.Water;
  parameter Modelica.Units.SI.Time tLoaAgg=300
    "Time resolution of load aggregation";

  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal/4;
  parameter Real freq = 1/3600 "sinusoid frequency in Hz (1/3600 = 1 cycle per hour)";
  parameter Real phase = 0 "phase shift (rad)";

  parameter Modelica.Units.SI.Temperature TGro=283.15 "Ground temperature";
  parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Example borFieUTubDat(
    filDat=IDEAS.Fluid.Geothermal.Borefields.Data.Filling.Bentonite(
    steadyState=true),
    conDat=IDEAS.Fluid.Geothermal.Borefields.Data.Configuration.Example(
    borCon=IDEAS.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.SingleUTube))
    annotation (Placement(transformation(extent={{-88,72},{-68,92}})));


  IDEAS.Fluid.Geothermal.Borefields.Validation.OneUTube borFieUTub1(
    redeclare package Medium = MediumTdep,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{30,50},{50,70}})));
  Sources.MassFlowSource_T sou1(
    redeclare package Medium = MediumTdep,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-52,50},{-32,
            70}}, rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn1(
    redeclare package Medium = MediumTdep,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-20,50},{0,70}})));
  Sources.Boundary_pT sin1(
    redeclare package Medium = MediumTdep,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{130,50},{110,
            70}},
          rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut1(
    redeclare package Medium = MediumTdep,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{80,50},{100,70}})));
  IDEAS.Fluid.Geothermal.Borefields.OneUTube borFieUTub4(
    redeclare package Medium = Medium,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{34,-86},{54,-66}})));
  Sources.MassFlowSource_T sou4(
    redeclare package Medium = Medium,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-48,-86},{
            -28,-66}},
                   rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn4(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-16,-86},{4,-66}})));
  Sources.Boundary_pT sin4(
    redeclare package Medium = Medium,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{134,-86},{
            114,-66}},
                   rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut4(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{84,-86},{104,-66}})));
  OneUTube borFieUTub2(
    redeclare package Medium = Medium,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{34,4},{54,24}})));
  Sources.MassFlowSource_T sou2(
    redeclare package Medium = Medium,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-48,4},{
            -28,24}},
                   rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn2(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-16,4},{4,24}})));
  Sources.Boundary_pT sin2(
    redeclare package Medium = Medium,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{134,4},{114,
            24}},  rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut2(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{84,4},{104,24}})));
  IDEAS.Fluid.Geothermal.Borefields.OneUTube borFieUTub3(
    redeclare package Medium = Medium,
    borFieDat=borFieUTubDat,
    tLoaAgg=tLoaAgg,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
    annotation (Placement(transformation(extent={{34,-40},{54,-20}})));
  Sources.MassFlowSource_T sou3(
    redeclare package Medium = Medium,
    nPorts=1,
    use_T_in=false,
    m_flow=m_flow_nominal,
    T=303.15) "Source" annotation (Placement(transformation(extent={{-48,-40},{
            -28,-20}},
                   rotation=0)));
  Sensors.TemperatureTwoPort TUTubIn3(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
  Sources.Boundary_pT sin3(
    redeclare package Medium = Medium,
    use_p_in=false,
    use_T_in=false,
    nPorts=1,
    p=101330,
    T=283.15) "Sink" annotation (Placement(transformation(extent={{134,-40},{
            114,-20}},
                   rotation=0)));
  Sensors.TemperatureTwoPort TUTubOut3(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    tau=0) "Inlet temperature of the borefield with UTube configuration"
    annotation (Placement(transformation(extent={{84,-40},{104,-20}})));
equation

  connect(TUTubIn1.port_b,borFieUTub1. port_a)
    annotation (Line(points={{0,60},{30,60}}, color={0,127,255}));
  connect(borFieUTub1.port_b,TUTubOut1. port_a)
    annotation (Line(points={{50,60},{80,60}}, color={0,127,255}));
  connect(TUTubOut1.port_b,sin1. ports[1])
    annotation (Line(points={{100,60},{110,60}}, color={0,127,255}));
  connect(TUTubIn1.port_a,sou1. ports[1])
    annotation (Line(points={{-20,60},{-32,60}}, color={0,127,255}));
  connect(TUTubIn4.port_b,borFieUTub4. port_a)
    annotation (Line(points={{4,-76},{34,-76}}, color={0,127,255}));
  connect(borFieUTub4.port_b,TUTubOut4. port_a)
    annotation (Line(points={{54,-76},{84,-76}}, color={0,127,255}));
  connect(TUTubOut4.port_b,sin4. ports[1])
    annotation (Line(points={{104,-76},{114,-76}}, color={0,127,255}));
  connect(TUTubIn4.port_a,sou4. ports[1])
    annotation (Line(points={{-16,-76},{-28,-76}}, color={0,127,255}));
  connect(TUTubIn2.port_b,borFieUTub2. port_a)
    annotation (Line(points={{4,14},{34,14}},   color={0,127,255}));
  connect(borFieUTub2.port_b,TUTubOut2. port_a)
    annotation (Line(points={{54,14},{84,14}},   color={0,127,255}));
  connect(TUTubOut2.port_b,sin2. ports[1])
    annotation (Line(points={{104,14},{114,14}},   color={0,127,255}));
  connect(TUTubIn2.port_a,sou2. ports[1])
    annotation (Line(points={{-16,14},{-28,14}},   color={0,127,255}));
  connect(TUTubIn3.port_b,borFieUTub3. port_a)
    annotation (Line(points={{4,-30},{34,-30}}, color={0,127,255}));
  connect(borFieUTub3.port_b,TUTubOut3. port_a)
    annotation (Line(points={{54,-30},{84,-30}}, color={0,127,255}));
  connect(TUTubOut3.port_b,sin3. ports[1])
    annotation (Line(points={{104,-30},{114,-30}}, color={0,127,255}));
  connect(TUTubIn3.port_a,sou3. ports[1])
    annotation (Line(points={{-16,-30},{-28,-30}}, color={0,127,255}));
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
        Rectangle(extent={{-60,76},{140,42}}, lineColor={28,108,200}),
        Rectangle(extent={{-60,-60},{140,-94}}, lineColor={28,108,200}),
        Text(
          extent={{-58,84},{168,78}},
          textColor={28,108,200},
          textString=
              "Temperature-dependent fluid properties (MediumTdep = IDEAS.Media.Specialized.Water.TemperatureDependentDensity)"),
        Text(
          extent={{-54,-54},{48,-54}},
          textColor={28,108,200},
          textString="Fixed fluid properties (Medium = IDEAS.Media.Water)"),
        Rectangle(extent={{-60,30},{140,-4}},  lineColor={28,108,200}),
        Text(
          extent={{-56,36},{156,38}},
          textColor={28,108,200},
          textString=
              "Fixed fluid properties for Medium, temperature varying heat transfer fluid properties (Medium = IDEAS.Media.Water)"),
        Rectangle(extent={{-60,-14},{140,-48}}, lineColor={28,108,200}),
        Text(
          extent={{-54,-8},{168,-10}},
          textColor={28,108,200},
          textString=
              "Temperature-dependent fluid properties (MediumTdep = IDEAS.Media.Specialized.Water.TemperatureDependentDensity)")}),
    Icon(coordinateSystem(extent={{-100,-100},{180,100}})));
end BorefieldsDynamicFluidProperties_water_injection_var;
