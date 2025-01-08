within IDEAS.Examples.Tutorial.DetailedHouse;
model DetailedHouse6
  "Extension of example 5 that adds a heating system"
  extends DetailedHouse5;
  package MediumWater = IDEAS.Media.Water "Water Medium";

  Fluid.HeatPumps.ScrollWaterToWater heaPum(
    m2_flow_nominal=pumpPrim.m_flow_nominal,
    enable_variable_speed=false,
    m1_flow_nominal=pumpSec.m_flow_nominal,
    redeclare package Medium1 = MediumWater,
    redeclare package Medium2 = MediumWater,
    datHeaPum=
        Fluid.HeatPumps.Data.ScrollWaterToWater.Heating.Viessmann_BW301A21_28kW_5_94COP_R410A(),
    scaling_factor=0.025,
    dp1_nominal=10000,
    dp2_nominal=10000,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
    "Heat pump model, rescaled for low thermal powers"
                                            annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={190,10})));

  Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(
    redeclare package Medium = MediumWater,
    Q_flow_nominal=500,
    T_a_nominal=318.15,
    T_b_nominal=308.15,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
                        "Radiator for zone 1" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={50,-10})));
  Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad1(
    redeclare package Medium = MediumWater,
    Q_flow_nominal=500,
    T_a_nominal=318.15,
    T_b_nominal=308.15,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
                        "Radiator for zone 2" annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={90,-10})));
  Fluid.Actuators.Valves.TwoWayTRV val(
    m_flow_nominal=rad.m_flow_nominal,
    dpValve_nominal=20000,
    redeclare package Medium = MediumWater) "Thermostatic valve for first zone"
                           annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={50,30})));
  Fluid.Movers.FlowControlled_dp pumpSec(
    dp_nominal=20000,
    inputType=IDEAS.Fluid.Types.InputType.Constant,
    m_flow_nominal=rad.m_flow_nominal + rad1.m_flow_nominal,
    redeclare package Medium = MediumWater,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
    "Circulation pump at secondary side"
    annotation (Placement(transformation(extent={{120,50},{100,70}})));
  Fluid.Sources.Boundary_pT bou(
    nPorts=2,
    redeclare package Medium = MediumWater,
    T=283.15) "Cold water source for heat pump"
              annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={270,10})));
  Fluid.Actuators.Valves.TwoWayTRV val1(
    dpValve_nominal=20000,
    m_flow_nominal=rad1.m_flow_nominal,
    redeclare package Medium = MediumWater)
    "Thermostatic valve for second zone"
                                        annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={90,30})));
  Fluid.Movers.FlowControlled_dp pumpPrim(
    inputType=IDEAS.Fluid.Types.InputType.Constant,
    dp_nominal=10000,
    m_flow_nominal=rad.m_flow_nominal + rad1.m_flow_nominal,
    redeclare package Medium = MediumWater,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial)
    "Circulation pump at primary side"
    annotation (Placement(transformation(extent={{240,50},{220,70}})));
  Modelica.Blocks.Sources.IntegerConstant HPOn(k=1) "Heat pump is on"
    annotation (Placement(transformation(extent={{160,-62},{180,-42}})));
  Modelica.Blocks.Continuous.Integrator EEl(k=1/3600000)
    "Electrical energy meter with conversion to kWh"
    annotation (Placement(transformation(extent={{280,40},{300,60}})));
  Fluid.Sensors.TemperatureTwoPort senTemSup(redeclare package Medium =
        MediumWater, m_flow_nominal=pumpSec.m_flow_nominal)
    "Supply water temperature sensor"
    annotation (Placement(transformation(extent={{146,70},{126,50}})));
  Fluid.Sources.Boundary_pT bou1(
    nPorts=1,
    redeclare package Medium = MediumWater,
    T=283.15) "Cold water source for heat pump"
              annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={110,90})));
  Fluid.Storage.Stratified tan(
    redeclare package Medium = MediumWater,
    m_flow_nominal=pumpSec.m_flow_nominal,
    VTan=0.1,
    hTan=0.5,
    dIns=0.1) "Buffer tank for avoiding excessive heat pump on/off switches"
    annotation (Placement(transformation(extent={{178,50},{158,70}})));
equation
  connect(val.port_b, rad.port_a)
    annotation (Line(points={{50,20},{50,0}}, color={0,127,255}));
  connect(val1.port_b, rad1.port_a)
    annotation (Line(points={{90,20},{90,0}}, color={0,127,255}));
  connect(val.port_a, pumpSec.port_b)
    annotation (Line(points={{50,40},{50,60},{100,60}}, color={0,127,255}));
  connect(val1.port_a, pumpSec.port_b)
    annotation (Line(points={{90,40},{90,60},{100,60}}, color={0,127,255}));
  connect(heaPum.port_a2, pumpPrim.port_b)
    annotation (Line(points={{196,20},{196,60},{220,60}}, color={0,127,255}));
  connect(heaPum.port_b2, bou.ports[1]) annotation (Line(points={{196,0},{196,
          -30},{248,-30},{248,11},{260,11}},      color={0,127,255}));
  connect(pumpPrim.port_a, bou.ports[2]) annotation (Line(points={{240,60},{248,
          60},{248,9},{260,9}},   color={0,127,255}));
  connect(rad.heatPortCon, rectangularZoneTemplate.gainCon) annotation (Line(
        points={{42.8,-8},{20,-8},{20,27},{10,27}}, color={191,0,0}));
  connect(rad.heatPortRad, rectangularZoneTemplate.gainRad) annotation (Line(
        points={{42.8,-12},{16,-12},{16,24},{10,24}}, color={191,0,0}));
  connect(rad1.heatPortCon, rectangularZoneTemplate1.gainCon) annotation (Line(
        points={{82.8,-8},{66,-8},{66,-33},{10,-33}}, color={191,0,0}));
  connect(rad1.heatPortRad, rectangularZoneTemplate1.gainRad) annotation (Line(
        points={{82.8,-12},{70,-12},{70,-36},{10,-36}}, color={191,0,0}));
  connect(HPOn.y, heaPum.stage) annotation (Line(points={{181,-52},{187,-52},{
          187,-2}},
                color={255,127,0}));
  connect(heaPum.P, EEl.u)
    annotation (Line(points={{190,21},{190,50},{278,50}}, color={0,0,127}));
  connect(rectangularZoneTemplate.TSensor, val.T) annotation (Line(points={{11,
          32},{26,32},{26,30},{39.4,30}}, color={0,0,127}));
  connect(rectangularZoneTemplate1.TSensor, val1.T) annotation (Line(points={{
          11,-28},{32,-28},{32,12},{79.4,12},{79.4,30}}, color={0,0,127}));
  connect(bou1.ports[1], pumpSec.port_b)
    annotation (Line(points={{110,80},{100,80},{100,60}}, color={0,127,255}));
  connect(senTemSup.port_b, pumpSec.port_a)
    annotation (Line(points={{126,60},{120,60}}, color={0,127,255}));
  connect(senTemSup.port_a, tan.port_b)
    annotation (Line(points={{146,60},{158,60},{158,50},{168,50}},
                                                 color={0,127,255}));
  connect(tan.port_a, heaPum.port_b1)
    annotation (Line(points={{168,70},{184,70},{184,20}}, color={0,127,255}));
  connect(rad1.port_b, heaPum.port_a1) annotation (Line(points={{90,-20},{90,
          -30},{184,-30},{184,0}}, color={0,127,255}));
  connect(rad.port_b, heaPum.port_a1) annotation (Line(points={{50,-20},{50,-30},
          {184,-30},{184,0}}, color={0,127,255}));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{280,100}},
          initialScale=0.1), graphics={Text(
          extent={{138,98},{224,90}},
          lineColor={28,108,200},
          textString="This sets the absolute pressure only"), Line(points={{126,
              86},{134,92}}, color={28,108,200})}), Icon(coordinateSystem(
          extent={{-100,-100},{100,100}}, initialScale=0.1)),
    __Dymola_Commands(file=
          "Resources/Scripts/Dymola/Examples/Tutorial/DetailedHouse/DetailedHouse6.mos"
        "Simulate and plot"),
    experiment(
      StartTime=10000000,
      StopTime=11000000,
      __Dymola_NumberOfIntervals=5000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Lsodar"),
    Documentation(revisions="<html>
<ul>
<li>
September 18, 2019 by Filip Jorissen:<br/>
First implementation for the IDEAS crash course.
</li>
</ul>
</html>",
info="<html>
<p>
This model extends <a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse5\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse5</a> by adding a HVAC system.
The system consists of a water-water heat pump, radiators, a storage tank, circulation
pumps and a low temperature heat source for the heat pump. The model includes constant 
control set points for the heat pump and pumps. A system is incorporated to measure 
the electrical energy consumption of the heat pump.
</p>
<h4>Required models</h4>
<ul>
<li>
<a href=\"modelica://IDEAS.Fluid.HeatPumps.ScrollWaterToWater\">
IDEAS.Fluid.HeatPumps.ScrollWaterToWater</a>
</li>
<li>
<a href=\"modelica://IDEAS.Fluid.HeatExchangers.Radiators.RadiatorEN442_2\">
IDEAS.Fluid.HeatExchangers.Radiators.RadiatorEN442_2</a>
</li>
<li>
<a href=\"modelica://IDEAS.Fluid.Actuators.Valves.TwoWayTRV\">
IDEAS.Fluid.Actuators.Valves.TwoWayTRV</a>
</li>
<li>
<a href=\"modelica://IDEAS.Fluid.Movers.FlowControlled_dp\">
IDEAS.Fluid.Movers.FlowControlled_dp</a>
</li>
<li>
<a href=\"modelica://IDEAS.Fluid.Sources.Boundary_pT\">
IDEAS.Fluid.Sources.Boundary_pT</a>
</li>
<li>
<a href=\"modelica://IDEAS.Fluid.Storage.Stratified\">
IDEAS.Fluid.Storage.Stratified</a>
</li>
<li>
<a href=\"modelica://Modelica.Blocks.Sources.IntegerConstant\">
Modelica.Blocks.Sources.IntegerConstant</a>
</li>
<li>
<a href=\"modelica://Modelica.Blocks.Continuous.Integrator\">
Modelica.Blocks.Continuous.Integrator</a>
</li>
</ul>
<h4>Connection instructions</h4>
<p>
The model components of <a href=\"modelica://IDEAS.Buildings\"> IDEAS.Buildings</a> have been
integrated where possible. E.g. a zone can always have an occupancy schedule, or lighting models. However,
for HVAC, many configurations exist in reality. Some zones have a radiator, which physically resides within the
zone, while others may use floor heating, which is a property of the floor. Some zones may even have multiple
radiators with different sizes, etc. Due to this multitude of options, HVAC components have not been integrated
into the IDEAS.Buildings components (yet). Instead, they can be added manually by the user. This requires
some knowledge by the user about how to interconnect the different components, and even about the relevant
physical principles that are at play.
</p>
<p>
One notable example is that in each fluid loop the <i>absolute</i> pressure of that loop has to be
defined somewhere: pumps and valves only provide information about <i>differential</i> pressures. For this purpose
use <a href=\"modelica://IDEAS.Fluid.Sources.Boundary_pT\"> IDEAS.Fluid.Sources.Boundary_pT</a> and connect it to the loop, 
which will set the absolute pressure at the connection point.
</p>
<p>
A reference implementation for this example is shown in the figure below. 
</p>
<p align=\"center\">
<img alt=\" The schematic of Example 6.\"
src=\"modelica://IDEAS/Resources/Images/Examples/Tutorial/DetailedHouse/Schematic6.png\" width=\"700\"/>
</p>
<p>
<h4>Reference result</h4>
<p>
The following figures plot plot the zone temperatures <code>rectangularZoneTemplate.TSensor<\code> and <code>rectangularZoneTemplate1<\code>.
TSensor, the radiator heat flow rates <code>rad.Q_flow<\code> and <code>rad1.Q_flow<\code>, the heat pump condenser temperature
<code>heaPum.con.T</code> and the heat pump heat flow rate <code>heaPum.QCon_flow</code>.
</p>
<p align=\"center\">
<img alt=\" The schematic of Example 6.\"
src=\"modelica://IDEAS/Resources/Images/Examples/Tutorial/DetailedHouse/Example6.png\" width=\"700\"/>
</p>
<p>
<p align=\"center\">
<img alt=\" The schematic of Example 6.\"
src=\"modelica://IDEAS/Resources/Images/Examples/Tutorial/DetailedHouse/Example6_bis.png\" width=\"700\"/>
</p>
<p>
This example illustrates the importance of control, which is currently not modelled. All pumps and the heat
pump are assumed to be active continuously, which is detrimental for the system performance. The COP
(heaPum.com.COP) is only about 2.9.
</p>
</html>"));
end DetailedHouse6;
