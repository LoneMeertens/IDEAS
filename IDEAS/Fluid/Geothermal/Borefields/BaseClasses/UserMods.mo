within IDEAS.Fluid.Geothermal.Borefields.BaseClasses;
package UserMods

model GroundTemperatureResponse
  "Model calculating discrete load aggregation with optional .mat input"

  parameter Modelica.Units.SI.Time tLoaAgg(final min=Modelica.Constants.eps)=3600
    "Time resolution of load aggregation";
  parameter Integer nCel(min=1)=5 "Number of cells per aggregation level";
  parameter Integer nSeg=12 "Number of segments per borehole for g-function calculation";
  parameter Integer nClu=5 "Number of clusters for g-function calculation";
  parameter Boolean forceGFunCalc = false
    "Set to true to force the thermal response to be calculated at the start instead of checking whether it has been pre-computed";
  parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Template borFieDat
    "Record containing all the parameters of the borefield model"
    annotation (choicesAllMatching=true, Placement(transformation(extent={{-80,-80},{-60,-60}})));

  // New parameters to read an existing .mat file instead of computing the g-function
  parameter Boolean temResMatInput = false
    "If true, read matrix 'TStep' from file specified by temResMatPath instead of computing it";
  parameter String temResMatPath = ""
    "Path to .mat file containing variable 'TStep' with size nTimTot x 2. Can be absolute, relative or a resource URI (tool-dependent).";

  Modelica.Blocks.Interfaces.RealOutput delTBor(unit="K")
    "Temperature difference current borehole wall temperature minus initial borehole wall temperature"
    annotation (Placement(transformation(extent={{100,-14},{126,12}}),
        iconTransformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealInput QBor_flow(unit="W")
    "Heat flow from all boreholes combined (positive if heat from fluid into soil)"
    annotation (Placement(transformation(extent={{-120,-10},{-100,10}}),
        iconTransformation(extent={{-120,-10},{-100,10}})));

protected
  constant Integer nTimSho = 26 "Number of time steps in short time region";
  constant Integer nTimLon = 50 "Number of time steps in long time region";
  constant Real ttsMax = exp(5) "Maximum non-dimensional time for g-function calculation";
  constant Integer nTimTot = nTimSho+nTimLon "Total length of g-function vector";
  constant Real lvlBas = 2 "Base for exponential cell growth between levels";

  parameter String SHAgfun=
    IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.ThermalResponseFactors.shaGFunction(
      nBor=borFieDat.conDat.nBor,
      cooBor=borFieDat.conDat.cooBor,
      hBor=borFieDat.conDat.hBor,
      dBor=borFieDat.conDat.dBor,
      rBor=borFieDat.conDat.rBor,
      aSoi=borFieDat.soiDat.aSoi,
      nSeg=nSeg,
      nClu=nClu,
      nTimSho=nTimSho,
      nTimLon=nTimLon,
      ttsMax=ttsMax) "String with encrypted g-function arguments";
  parameter Modelica.Units.SI.Time timFin=(borFieDat.conDat.hBor^2/(9*borFieDat.soiDat.aSoi)) *ttsMax
    "Final time for g-function calculation";
  parameter Integer i(min=1)=
    IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.countAggregationCells(
      lvlBas=lvlBas,
      nCel=nCel,
      timFin=timFin,
      tLoaAgg=tLoaAgg)
      "Number of aggregation cells";
  final parameter Real[nTimTot,2] timSer(each fixed=false)
    "g-function input from matrix, with the second column as temperature Tstep";
  final parameter Modelica.Units.SI.Time t_start(fixed=false)
    "Simulation start time";
  final parameter Modelica.Units.SI.Time[i] nu(each fixed=false)
    "Time vector for load aggregation";
  final parameter Real[i] kappa(each fixed=false)
    "Weight factor for each aggregation cell";
  final parameter Real[i] rCel(each fixed=false) "Cell widths";

  discrete Modelica.Units.SI.HeatFlowRate[i] QAgg_flow
    "Vector of aggregated loads";
  discrete Modelica.Units.SI.HeatFlowRate[i] QAggShi_flow
    "Shifted vector of aggregated loads";
  discrete Integer curCel "Current occupied cell";

  discrete Modelica.Units.SI.TemperatureDifference delTBor0
    "Previous time step's temperature difference current borehole wall temperature minus initial borehole temperature";
  discrete Real derDelTBor0(unit="K/s")
    "Derivative of wall temperature change from previous time steps";
  final parameter Real dTStepdt(fixed=false)
    "Time derivative of g/(2*pi*H*Nb*ks) within most recent cell";

  Modelica.Units.SI.Heat U "Accumulated heat flow from all boreholes";
  discrete Modelica.Units.SI.Heat U_old
    "Accumulated heat flow from all boreholes at last aggregation step";

initial equation
  QAgg_flow = zeros(i);
  curCel = 1;
  delTBor = 0;
  QAggShi_flow = QAgg_flow;
  delTBor0 = 0;
  U = 0;
  U_old = 0;
  derDelTBor0 = 0;

  (nu,rCel) = IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationCellTimes(
    i=i,
    lvlBas=lvlBas,
    nCel=nCel,
    tLoaAgg=tLoaAgg,
    timFin=timFin);

  t_start = time;

  kappa = IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.aggregationWeightingFactors(
    i=i,
    nTimTot=nTimTot,
    TStep=timSer,
    nu=nu);

  dTStepdt = kappa[1]/tLoaAgg;

  // --- NEW: read timSer from a .mat file instead of calculating ---
  timSer =  Modelica.Utilities.Streams.readRealMatrix(
    fileName=temResMatPath,
    matrixName="TStep",
    nrow=nTimTot,
    ncol=2);


equation
  der(delTBor) = dTStepdt*QBor_flow + derDelTBor0;
  der(U) = QBor_flow;

  when sample(t_start, tLoaAgg) then
    // Assign average load since last aggregation step to the first cell of the
    // aggregation vector
    U_old = U;

    // Store (U - pre(U_old))/tLoaAgg in QAgg_flow[1], and pre(QAggShi_flow) in the other elements
    QAgg_flow = cat(1, {(U - pre(U_old))/tLoaAgg}, pre(QAggShi_flow[2:end]));
    // Shift loads in aggregation cells
    (curCel,QAggShi_flow) = IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.shiftAggregationCells(
      i=i,
      QAgg_flow=QAgg_flow,
      rCel=rCel,
      nu=nu,
      curTim=(time - t_start));

    // Determine the temperature change at the next aggregation step (assuming
    // no loads until then)
    delTBor0 = IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.LoadAggregation.temporalSuperposition(
      i=i,
      QAgg_flow=QAggShi_flow,
      kappa=kappa,
      curCel=curCel);

    derDelTBor0 = (delTBor0-delTBor)/tLoaAgg;
  end when;

end GroundTemperatureResponse;


// -----------------------------------------------------------------------------
// Modified PartialBorefield
// - Adds temResMatInput and temResMatPath parameters and passes them to the
//   GroundTemperatureResponse instance (groTemRes). This makes the option
//   available at the "top model" level (e.g. OneUTube extends this and will
//   inherit the parameter).
// -----------------------------------------------------------------------------

partial model PartialBorefield
  "Borefield model using single U-tube borehole heat exchanger configuration.Calculates the average fluid temperature T_fts of the borefield for a given (time dependent) load Q_flow"

  extends IDEAS.Fluid.Interfaces.PartialTwoPortInterface(
    final m_flow_nominal=borFieDat.conDat.mBorFie_flow_nominal);

  extends IDEAS.Fluid.Interfaces.TwoPortFlowResistanceParameters(
    final dp_nominal=borFieDat.conDat.dp_nominal,
    final computeFlowResistance=(borFieDat.conDat.dp_nominal > Modelica.Constants.eps));

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component"
      annotation (choices(
        choice(redeclare package Medium = IDEAS.Media.Water "Water"),
        choice(redeclare package Medium =
            IDEAS.Media.Antifreeze.PropyleneGlycolWater (
              property_T=293.15,
              X_a=0.40)
              "Propylene glycol water, 40% mass fraction")));

  constant Real mSenFac(min=1)=1
    "Factor for scaling the sensible thermal mass of the volume";

  // Assumptions
  parameter Modelica.Fluid.Types.Dynamics energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
    "Type of energy balance: dynamic (3 initialization options) or steady state"
    annotation(Evaluate=true, Dialog(tab = "Dynamics", group="Conservation equations"));

  // Initialization
  parameter Medium.AbsolutePressure p_start = Medium.p_default
    "Start value of pressure"
    annotation(Dialog(tab = "Initialization"));

  // Simulation parameters
  parameter Modelica.Units.SI.Time tLoaAgg=300
    "Time resolution of load aggregation";
  parameter Integer nCel(min=1)=5 "Number of cells per aggregation level";
  parameter Integer nSeg(min=1)=10
    "Number of segments to use in vertical discretization of the boreholes";
  parameter Boolean forceGFunCalc = false
    "Set to true to force the thermal response to be calculated at the start instead of checking whether this has been pre-computed"
    annotation (Dialog(tab="Advanced", group="g-function"));
  parameter Integer nSegGFun(min=1)=12 "Number of segments to use in the calculation of the g-function"
    annotation (Dialog(tab="Advanced", group="g-function"));
  parameter Integer nClu(min=1)=5
    "Number of borehole clusters to use in the calculation of the g-function"
    annotation (Dialog(tab="Advanced", group="g-function"));

  // New top-level options to choose an external .mat file
  parameter String temResMatPath = "" "Path to .mat file (contains variable TStep)";

  // General parameters of borefield
  parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Template borFieDat "Borefield data"
    annotation (choicesAllMatching=true,Placement(transformation(extent={{-80,-80},{-60,-60}})));

  // Temperature gradient in undisturbed soil
  parameter Modelica.Units.SI.Temperature TExt0_start=283.15
    "Initial far field temperature"
    annotation (Dialog(tab="Initialization", group="Soil"));
  parameter Modelica.Units.SI.Temperature TExt_start[nSeg]={if z[i] >= z0 then
      TExt0_start + (z[i] - z0)*dT_dz else TExt0_start for i in 1:nSeg}
    "Temperature of the undisturbed ground"
    annotation (Dialog(tab="Initialization", group="Soil"));

  parameter Modelica.Units.SI.Temperature TGro_start[nSeg]=TExt_start
    "Start value of grout temperature"
    annotation (Dialog(tab="Initialization", group="Filling material"));

  parameter Modelica.Units.SI.Temperature TFlu_start[nSeg]=TGro_start
    "Start value of fluid temperature" annotation (Dialog(tab="Initialization"));

  parameter Modelica.Units.SI.Height z0=10
    "Depth below which the temperature gradient starts"
    annotation (Dialog(tab="Initialization", group="Temperature profile"));
  parameter Real dT_dz(final unit="K/m", min=0) = 0.01
    "Vertical temperature gradient of the undisturbed soil for h below z0"
    annotation (Dialog(tab="Initialization", group="Temperature profile"));

  Modelica.Blocks.Interfaces.RealOutput TBorAve(final quantity="ThermodynamicTemperature",
                                                final unit="K",
                                                displayUnit = "degC",
                                                start=TExt0_start)
    "Average borehole wall temperature in the borefield"
    annotation (Placement(transformation(extent={{100,34},{120,54}})));

  // Instantiate modified GroundTemperatureResponse and forward temResMat parameters
  UserMods.GroundTemperatureResponse groTemRes(
    final tLoaAgg=tLoaAgg,
    final nCel=nCel,
    final nSeg=nSegGFun,
    final nClu=nClu,
    final borFieDat=borFieDat,
    final forceGFunCalc=forceGFunCalc,
    final temResMatPath=temResMatPath)
    "Ground temperature response (optionally read from .mat)"
    annotation (Placement(transformation(extent={{20,70},{40,90}})));

  replaceable IDEAS.Fluid.Geothermal.Borefields.BaseClasses.Boreholes.BaseClasses.PartialBorehole borHol constrainedby
    IDEAS.Fluid.Geothermal.Borefields.BaseClasses.Boreholes.BaseClasses.PartialBorehole(
    redeclare final package Medium = Medium,
    final borFieDat=borFieDat,
    final nSeg=nSeg,
    final m_flow_nominal=m_flow_nominal/borFieDat.conDat.nBor,
    final dp_nominal=dp_nominal,
    final allowFlowReversal=allowFlowReversal,
    final m_flow_small=m_flow_small,
    final show_T=show_T,
    final computeFlowResistance=computeFlowResistance,
    final from_dp=from_dp,
    final linearizeFlowResistance=linearizeFlowResistance,
    final deltaM=deltaM,
    final energyDynamics=energyDynamics,
    final p_start=p_start,
    final mSenFac=mSenFac,
    final TFlu_start=TFlu_start,
    final TGro_start=TGro_start) "Borehole"
    annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

protected
  parameter Modelica.Units.SI.Height z[nSeg]={borFieDat.conDat.hBor/nSeg*(i -
      0.5) for i in 1:nSeg}
    "Distance from the surface to the considered segment";

  IDEAS.Fluid.BaseClasses.MassFlowRateMultiplier masFloDiv(
    redeclare final package Medium = Medium,
    final allowFlowReversal=allowFlowReversal,
    final use_input=false,
    final k=1/borFieDat.conDat.nBor)
                                   "Division of flow rate"
    annotation (Placement(transformation(extent={{-80,-50},{-60,-30}})));

  IDEAS.Fluid.BaseClasses.MassFlowRateMultiplier masFloMul(
    redeclare final package Medium = Medium,
    final allowFlowReversal=allowFlowReversal,
    final use_input=false,
    final k=borFieDat.conDat.nBor) "Mass flow multiplier"
    annotation (Placement(transformation(extent={{60,-50},{80,-30}})));

  Modelica.Blocks.Math.Gain gaiQ_flow(k=borFieDat.conDat.nBor)
    "Gain to multiply the heat extracted by one borehole by the number of boreholes"
    annotation (Placement(transformation(extent={{-20,70},{0,90}})));
  IDEAS.Utilities.Math.Average AveTBor(nin=nSeg)
    "Average temperature of all the borehole segments"
    annotation (Placement(transformation(extent={{50,34},{70,54}})));

  Modelica.Blocks.Sources.Constant TSoiUnd[nSeg](
    k = TExt_start,
    each y(unit="K",
           displayUnit="degC"))
    "Undisturbed soil temperature"
    annotation (Placement(transformation(extent={{-40,14},{-20,34}})));

  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor QBorHol[nSeg]
    "Heat flow rate of all segments of the borehole"
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=90,
        origin={0,-10})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature TemBorWal[nSeg]
    "Borewall temperature"
    annotation (Placement(transformation(extent={{50,6},{70,26}})));

  Modelica.Blocks.Math.Add TSoiDis[nSeg](each final k1=1, each final k2=1)
    "Addition of undisturbed soil temperature and change of soil temperature"
    annotation (Placement(transformation(extent={{10,20},{30,40}})));

  Modelica.Blocks.Math.Sum QTotSeg_flow(
    final nin=nSeg,
    final k = ones(nSeg))
    "Total heat flow rate for all segments of this borehole"
    annotation (Placement(transformation(extent={{-60,70},{-40,90}})));

  Modelica.Blocks.Routing.Replicator repDelTBor(final nout=nSeg)
    "Signal replicator for temperature difference of the borehole"
    annotation (Placement(transformation(extent={{60,70},{80,90}})));

equation
  connect(masFloMul.port_b, port_b)
    annotation (Line(points={{80,-40},{90,-40},{90,0},{100,0}},
                                                     color={0,127,255}));
  connect(masFloDiv.port_a, port_a)
    annotation (Line(points={{-80,-40},{-90,-40},{-90,0},{-100,0}},
                                                color={0,127,255}));
  connect(masFloDiv.port_b, borHol.port_a)
    annotation (Line(points={{-60,-40},{-10,-40}},     color={0,127,255}));
  connect(borHol.port_b, masFloMul.port_a)
    annotation (Line(points={{10,-40},{60,-40}},    color={0,127,255}));
  connect(QBorHol.port_a, borHol.port_wall)
    annotation (Line(points={{-4.44089e-16,-20},{0,-20},{0,-30}},
                                                        color={191,0,0}));
  connect(QBorHol.Q_flow, QTotSeg_flow.u)
    annotation (Line(points={{-11,-10},{-86,-10},{-86,80},{-62,80}},
                                                          color={0,0,127}));
  connect(groTemRes.delTBor, repDelTBor.u)
    annotation (Line(points={{41,80},{58,80}}, color={0,0,127}));
  connect(TSoiDis.u1, repDelTBor.y) annotation (Line(points={{8,36},{0,36},{0,
          60},{90,60},{90,80},{81,80}},
                        color={0,0,127}));
  connect(TSoiDis.u2, TSoiUnd.y) annotation (Line(points={{8,24},{-19,24}},
                         color={0,0,127}));
  connect(QTotSeg_flow.y, gaiQ_flow.u)
    annotation (Line(points={{-39,80},{-22,80}}, color={0,0,127}));
  connect(gaiQ_flow.y, groTemRes.QBor_flow)
    annotation (Line(points={{1,80},{19,80}}, color={0,0,127}));
  connect(TSoiDis.y, TemBorWal.T)
    annotation (Line(points={{31,30},{36,30},{36,16},{48,16}},
                                               color={0,0,127}));
  connect(QBorHol.port_b, TemBorWal.port) annotation (Line(points={{4.44089e-16,
          0},{0,0},{0,4},{80,4},{80,16},{70,16}},   color={191,0,0}));
  connect(TSoiDis.y, AveTBor.u) annotation (Line(points={{31,30},{36,30},{36,44},
          {48,44}}, color={0,0,127}));
  connect(AveTBor.y, TBorAve)
    annotation (Line(points={{71,44},{110,44}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={
        Rectangle(
          extent={{-100,60},{100,-66}},
          lineColor={0,0,0},
          fillColor={234,210,210},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-88,-6},{-32,-62}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-82,-12},{-38,-56}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-88,54},{-32,-2}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-82,48},{-38,4}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-26,54},{30,-2}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-20,48},{24,4}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-28,-6},{28,-62}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{-22,-12},{22,-56}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{36,56},{92,0}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{42,50},{86,6}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{38,-4},{94,-60}},
          lineColor={0,0,0},
          fillColor={223,188,190},
          fillPattern=FillPattern.Forward),
        Ellipse(
          extent={{44,-10},{88,-54}},
          lineColor={0,0,0},
          fillColor={0,0,255},
          fillPattern=FillPattern.Forward)}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}})),Documentation(info="<html>
<p>
This model simulates a borefield containing one or multiple boreholes
using the parameters in the <code>borFieDat</code> record.
</p>
<p>
Heat transfer to the soil is modeled using only one borehole heat exchanger
(To be added in an extended model). The
fluid mass flow rate into the borehole is divided to reflect the per-borehole
fluid mass flow rate. The borehole model calculates the dynamics within the
borehole itself using an axial discretization and a resistance-capacitance
network for the internal thermal resistances between the individual pipes and
between each pipe and the borehole wall.
</p>
<p>
The thermal interaction between the borehole wall and the surrounding soil
is modeled using
<a href=\"modelica://IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse\">
IDEAS.Fluid.Geothermal.Borefields.BaseClasses.HeatTransfer.GroundTemperatureResponse</a>,
which uses a cell-shifting load aggregation technique to calculate the borehole wall
temperature after calculating and/or read (from a previous calculation) the borefield's thermal response factor.
</p>
</html>", revisions="<html>
<ul>
<li>
April 9, 2021, by Michael Wetter:<br/>
Corrected placement of <code>each</code> keyword.<br/>
See <a href=\"https://github.com/lbl-srg/modelica-IDEAS/pull/2440\">IDEAS, PR #2440</a>.
</li>
<li>
August 25, 2020, by Filip Jorissen:<br/>
Switched port connections for <code>masFloDiv</code>.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/41\">#41</a>.
</li>
<li>
March 24, 2020, by Damien Picard:<br/>
Propagated flowReversal into <code>masFloDiv</code> and <code>masFloMul</code>.
</li>
<li>
June 7, 2019, by Massimo Cimmino:<br/>
Converted instances that are not of interest to user to be <code>protected</code>.
</li>
<li>
June 4, 2019, by Massimo Cimmino:<br/>
Added an output for the average borehole wall temperature.
See
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1107\">#1107</a>.
</li>
<li>
April 11, 2019, by Filip Jorissen:<br/>
Added <code>choicesAllMatching</code> for <code>borFieDat</code>.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1117\">#1117</a>.
</li>
<li>
January 18, 2019, by Jianjun Hu:<br/>
Limited the media choice to water and glycolWater.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1050\">#1050</a>.
</li>
<li>
July 2018, by Alex Laferri&egrave;re:<br/>
Changed into a partial model and changed documentation to reflect the new approach
used by the borefield models.
</li>
<li>
July 2014, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</html>"));
end PartialBorefield;

  model OneUTube "Borefield model containing single U-tube boreholes"
    extends IDEAS.Fluid.Geothermal.Borefields.BaseClasses.UserMods.PartialBorefield(
      redeclare IDEAS.Fluid.Geothermal.Borefields.BaseClasses.Boreholes.OneUTube borHol);

    annotation (
    defaultComponentName="borFie",
    Documentation(info="<html>
<p>
This model simulates a borefield containing one or many single U-tube boreholes
using the parameters in the <code>borFieDat</code> record.
</p>
<p>
Heat transfer to the soil is modeled using only one borehole heat exchanger. The
fluid mass flow rate into the borehole is divided to reflect the per-borehole
fluid mass flow rate. The borehole model calculates the dynamics within the
borehole itself using an axial discretization and a resistance-capacitance
network for the internal thermal resistances between the individual pipes and
between each pipe and the borehole wall.
</p>
</html>",   revisions="<html>
<ul>
<li>
July 2018, by Alex Laferri&egrave;re:<br/>
Extended partial model and changed documentation to reflect the new approach
used by the borefield models.
</li>
<li>
July 2014, by Damien Picard:<br/>
First implementation.
</li>
</ul>
</html>"));
  end OneUTube;

  package Examples "Example models for IDEAS.Fluid.Geothermal.Borefields"
  extends Modelica.Icons.ExamplesPackage;

    model Borefields
      "Example model with several borefield configurations operating simultaneously"
      extends Modelica.Icons.Example;

      package Medium = IDEAS.Media.Water;

      parameter Modelica.Units.SI.Time tLoaAgg=300
        "Time resolution of load aggregation";

      parameter Modelica.Units.SI.Temperature TGro=283.15 "Ground temperature";
      parameter IDEAS.Fluid.Geothermal.Borefields.Data.Borefield.Example borFieUTubDat(
        filDat=IDEAS.Fluid.Geothermal.Borefields.Data.Filling.Bentonite(
        steadyState=true),
        conDat=IDEAS.Fluid.Geothermal.Borefields.Data.Configuration.Template(
          borCon = Types.BoreholeConfiguration.SingleUTube,
          cooBor={{0,0},{6,0},{0,6},{6,6},{12,0},{12,6},{18,0},{18,6}},
          mBor_flow_nominal=0.3,
          dp_nominal=5e4,
          hBor=150.0,
          rBor=0.075,
          dBor=4.0,
          rTub=0.02,
          kTub=0.5,
          eTub=0.002,
          xC=0.05))
        annotation (Placement(transformation(extent={{70,-92},{90,-72}})));

      IDEAS.Fluid.Geothermal.Borefields.BaseClasses.UserMods.OneUTube borFieUTubCircularImport(
        redeclare package Medium = Medium,
        temResMatPath=Modelica.Utilities.Files.loadResource("modelica://IDEAS/Resources/Data/Fluid/Geothermal/Borefields/HeatTransfer/temperatureResponseMatrix/customTStepInclined2.mat"),
        borFieDat=borFieUTubDat,
        tLoaAgg=tLoaAgg,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        TExt0_start=TGro)
        "Inclined Borefield with a U-tube boreholes in circular configuration"
        annotation (Placement(transformation(extent={{-10,68},{10,88}})));
      IDEAS.Fluid.Sources.MassFlowSource_T sou(
        redeclare package Medium = Medium,
        nPorts=1,
        use_T_in=false,
        m_flow=borFieUTubDat.conDat.mBorFie_flow_nominal,
        T=303.15) "Source" annotation (Placement(transformation(extent={{-92,68},{-72,
                88}},      rotation=0)));
      IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubIn(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{-60,68},{-40,88}})));
      IDEAS.Fluid.Sources.Boundary_pT sin(
        redeclare package Medium = Medium,
        use_p_in=false,
        use_T_in=false,
        nPorts=1,
        p=101330,
        T=283.15) "Sink" annotation (Placement(transformation(extent={{90,68},{70,88}},
                       rotation=0)));
      IDEAS.Fluid.Sensors.TemperatureTwoPort TUTubOut(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{40,68},{60,88}})));

      OneUTube borFieUTubOptimalImport(
        redeclare package Medium = Medium,
        temResMatPath=Modelica.Utilities.Files.loadResource("modelica://IDEAS/Resources/Data/Fluid/Geothermal/Borefields/HeatTransfer/temperatureResponseMatrix/customTStepInclined1.mat"),
        borFieDat=borFieUTubDat,
        tLoaAgg=tLoaAgg,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        TExt0_start=TGro)
        "Inclined borefield with a U-tube boreholes in 'Optimal' configuration"
        annotation (Placement(transformation(extent={{-10,34},{10,54}})));
      Sources.MassFlowSource_T             sou1(
        redeclare package Medium = Medium,
        nPorts=1,
        use_T_in=false,
        m_flow=borFieUTubDat.conDat.mBorFie_flow_nominal,
        T=303.15) "Source" annotation (Placement(transformation(extent={{-92,34},{-72,
                54}},      rotation=0)));
      Sensors.TemperatureTwoPort             TUTubIn1(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{-60,34},{-40,54}})));
      Sources.Boundary_pT             sin1(
        redeclare package Medium = Medium,
        use_p_in=false,
        use_T_in=false,
        nPorts=1,
        p=101330,
        T=283.15) "Sink" annotation (Placement(transformation(extent={{90,34},{70,54}},
                       rotation=0)));
      Sensors.TemperatureTwoPort             TUTubOut1(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{40,34},{60,54}})));
      OneUTube borFieUTubRectImport(
        redeclare package Medium = Medium,
        temResMatPath=Modelica.Utilities.Files.loadResource("modelica://IDEAS/Resources/Data/Fluid/Geothermal/Borefields/HeatTransfer/temperatureResponseMatrix/customTStepRect.mat"),
        borFieDat=borFieUTubDat,
        tLoaAgg=tLoaAgg,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        TExt0_start=TGro)
        "Inclined borefield with a U-tube boreholes in 'Optimal' configuration"
        annotation (Placement(transformation(extent={{-10,0},{10,20}})));
      Sources.MassFlowSource_T             sou2(
        redeclare package Medium = Medium,
        nPorts=1,
        use_T_in=false,
        m_flow=borFieUTubDat.conDat.mBorFie_flow_nominal,
        T=303.15) "Source" annotation (Placement(transformation(extent={{-92,0},{-72,
                20}},      rotation=0)));
      Sensors.TemperatureTwoPort             TUTubIn2(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{-60,0},{-40,20}})));
      Sources.Boundary_pT             sin2(
        redeclare package Medium = Medium,
        use_p_in=false,
        use_T_in=false,
        nPorts=1,
        p=101330,
        T=283.15) "Sink" annotation (Placement(transformation(extent={{90,0},{70,20}},
                       rotation=0)));
      Sensors.TemperatureTwoPort             TUTubOut2(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{40,0},{60,20}})));
      IDEAS.Fluid.Geothermal.Borefields.OneUTube borFieUTubRect(
        redeclare package Medium = Medium,
        borFieDat=borFieUTubDat,
        tLoaAgg=tLoaAgg,
        energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
        TExt0_start=TGro) "Borefield with a U-tube borehole configuration"
        annotation (Placement(transformation(extent={{-10,-34},{10,-14}})));
      Sources.MassFlowSource_T             sou3(
        redeclare package Medium = Medium,
        nPorts=1,
        use_T_in=false,
        m_flow=borFieUTubDat.conDat.mBorFie_flow_nominal,
        T=303.15) "Source" annotation (Placement(transformation(extent={{-92,-34},{-72,
                -14}},     rotation=0)));
      Sensors.TemperatureTwoPort             TUTubIn3(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{-60,-34},{-40,-14}})));
      Sources.Boundary_pT             sin3(
        redeclare package Medium = Medium,
        use_p_in=false,
        use_T_in=false,
        nPorts=1,
        p=101330,
        T=283.15) "Sink" annotation (Placement(transformation(extent={{90,-34},{70,-14}},
                       rotation=0)));
      Sensors.TemperatureTwoPort             TUTubOut3(
        redeclare package Medium = Medium,
        m_flow_nominal=borFieUTubDat.conDat.mBorFie_flow_nominal,
        tau=0)
        "Inlet temperature of the borefield with UTube configuration"
        annotation (Placement(transformation(extent={{40,-34},{60,-14}})));
    equation
      connect(sou.ports[1], TUTubIn.port_a)
        annotation (Line(points={{-72,78},{-60,78}},   color={0,127,255}));
      connect(TUTubIn.port_b, borFieUTubCircularImport.port_a)
        annotation (Line(points={{-40,78},{-10,78}}, color={0,127,255}));
      connect(borFieUTubCircularImport.port_b, TUTubOut.port_a)
        annotation (Line(points={{10,78},{40,78}}, color={0,127,255}));
      connect(TUTubOut.port_b, sin.ports[1])
        annotation (Line(points={{60,78},{70,78}},            color={0,127,255}));
      connect(sou1.ports[1], TUTubIn1.port_a)
        annotation (Line(points={{-72,44},{-60,44}}, color={0,127,255}));
      connect(TUTubIn1.port_b, borFieUTubOptimalImport.port_a)
        annotation (Line(points={{-40,44},{-10,44}}, color={0,127,255}));
      connect(borFieUTubOptimalImport.port_b, TUTubOut1.port_a)
        annotation (Line(points={{10,44},{40,44}}, color={0,127,255}));
      connect(TUTubOut1.port_b, sin1.ports[1])
        annotation (Line(points={{60,44},{70,44}}, color={0,127,255}));
      connect(sou2.ports[1],TUTubIn2. port_a)
        annotation (Line(points={{-72,10},{-60,10}}, color={0,127,255}));
      connect(TUTubIn2.port_b, borFieUTubRectImport.port_a)
        annotation (Line(points={{-40,10},{-10,10}}, color={0,127,255}));
      connect(borFieUTubRectImport.port_b, TUTubOut2.port_a)
        annotation (Line(points={{10,10},{40,10}}, color={0,127,255}));
      connect(TUTubOut2.port_b,sin2. ports[1])
        annotation (Line(points={{60,10},{70,10}}, color={0,127,255}));
      connect(sou3.ports[1], TUTubIn3.port_a)
        annotation (Line(points={{-72,-24},{-60,-24}}, color={0,127,255}));
      connect(TUTubIn3.port_b, borFieUTubRect.port_a)
        annotation (Line(points={{-40,-24},{-10,-24}}, color={0,127,255}));
      connect(borFieUTubRect.port_b, TUTubOut3.port_a)
        annotation (Line(points={{10,-24},{40,-24}}, color={0,127,255}));
      connect(TUTubOut3.port_b, sin3.ports[1])
        annotation (Line(points={{60,-24},{70,-24}}, color={0,127,255}));
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
          StopTime=1800000000,
          Tolerance=1e-06,
          __Dymola_Algorithm="Dassl"));
    end Borefields;
  annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains examples for the use of functions that can be found in
<a href=\"modelica://IDEAS.Fluid.Geothermal.Borefields\">
IDEAS.Fluid.Geothermal.Borefields</a>.
</p>
</html>"));
  end Examples;
end UserMods;
