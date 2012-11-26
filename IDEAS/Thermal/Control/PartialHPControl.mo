within IDEAS.Thermal.Control;
partial model PartialHPControl "Basis of the heat Pump control algorithm"

  /* 
  This partial class contains the temperature control algorithm. It has to be extended
  in order to be complete controller.  
  

  */

  input Modelica.SIunits.Temperature TTankTop
    "Top (or near top) tank temperature";
  input Modelica.SIunits.Temperature TTankBot
    "Bottom (or near bottom) tank temperature";

  Modelica.SIunits.Temperature TBotSet(start=283.15)
    "Bottom temperature setpoint";
  Modelica.SIunits.Temperature TTopSet(start=283.15) "Top temperature setpoint";
  Modelica.SIunits.Temperature TBotEmpty(start=283.15)
    "Temperature in bottom corresponding to SOC = 0";

  Real SOC(start = 0);

  //output SI.Temperature THPSet(start = 283.15) "Heat pump set temperature";
  //output Real onOff(start=0) "onoff signal as Real";

  parameter Modelica.SIunits.TemperatureDifference dTSupRetNom=10
    "Nominal difference between supply and return water temperatures";
  parameter Modelica.SIunits.TemperatureDifference dTSafetyTop=3
    "Safety margin on top temperature setpoint" annotation(Evaluate=false);
  parameter Modelica.SIunits.TemperatureDifference dTSafetyBot=
      dTSafetyTop "Safety margin on bottom temperature setpoint";
  parameter Modelica.SIunits.TemperatureDifference dTHPTankSet(min=0)=2
    "Difference between tank setpoint and heat pump setpoint";

  parameter Boolean DHW = false "if true, the system has to foresee DHW";
  parameter Modelica.SIunits.Temperature TDHWSet=0
    "Setpoint temperature for the DHW outlet";
  parameter Modelica.SIunits.Temperature TColdWaterNom=273.15 + 10
    "Nominal cold water temperature";
  parameter Modelica.SIunits.Time timeFilter=43200
    "Time constant for filter on ambient temperature";

  HeatingCurve heatingCurve(
    timeFilter = timeFilter,
    dTOutHeaBal = 0,
    TSup_nominal = 273.15+45,
    TRet_nominal = 273.15+35,
    TRoo_nominal = 273.15+21,
    TOut_nominal = 273.15-8,
    redeclare IDEAS.BaseClasses.Math.MovingAverage filter(period=timeFilter))
    annotation (Placement(transformation(extent={{-54,44},{-34,64}})));
  outer IDEAS.SimInfoManager         sim
    annotation (Placement(transformation(extent={{24,50},{44,70}})));
  Modelica.Blocks.Interfaces.RealOutput THPSet(start = 283.15)
    "Heat pump set temperature"
    annotation (Placement(transformation(extent={{90,-10},{110,10}}),
        iconTransformation(extent={{90,-10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput onOff(start=0) "onoff signal as Real"
    annotation (Placement(transformation(extent={{90,-50},{110,-30}}),
        iconTransformation(extent={{90,-50},{110,-30}})));
  Modelica.Blocks.Interfaces.RealOutput THeaCur "Heating curve setpoint"
    annotation (Placement(transformation(extent={{92,-90},{112,-70}}),
        iconTransformation(extent={{92,-90},{112,-70}})));
initial equation
  //der(onOff) = 0;

equation
  heatingCurve.TOut = sim.Te;
  TBotEmpty = if DHW then TColdWaterNom else TTopSet - dTSupRetNom;
  //tankSOC is intentionally computed based only on 2 temperature sensors for practical reasons.  It is computed
  // with regard to TTopSet and TBotSet and a reference temperature (TBotEmpty)
  SOC=0.5*(TTankBot-TBotEmpty)/(TBotSet+dTSafetyBot-TBotEmpty)+0.5*(TTankTop-(TTopSet+dTSafetyTop))/(dTSupRetNom-dTSafetyTop);

  annotation(Icon(graphics={
        Line(
          points={{-100,80},{-100,-80}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{-100,0},{-40,0},{60,60}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{0,24},{0,48},{10,48},{10,70},{0,70},{0,88}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{40,0},{40,-40},{100,-40}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{60,0},{100,0}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{20,0},{20,-80},{100,-80}},
          color={0,0,127},
          smooth=Smooth.None)}),
      Diagram(graphics));
end PartialHPControl;
