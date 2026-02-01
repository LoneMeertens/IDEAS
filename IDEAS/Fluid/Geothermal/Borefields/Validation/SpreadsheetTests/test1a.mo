within IDEAS.Fluid.Geothermal.Borefields.Validation.SpreadsheetTests;
model test1a
  import Buildings;
  extends Modelica.Icons.Example;
  package Medium =
      BorSizingMod.Validation.Media.PropyleneGlycolWater_test1 (
       property_T=273.15,X_a=0);
  parameter Integer nSeg = 10;
  parameter Modelica.Units.SI.Temperature T_startAll = 273.15 + 17.5;
  parameter Modelica.Units.SI.Temperature TExt0_start=T_startAll;
  parameter Modelica.Units.SI.Length z0=4;
  parameter Real dT_dz(final unit="K/m", min=0) = 0;
  parameter Modelica.Units.SI.Height z[nSeg]={borefield.borFieDat.conDat.hBor/nSeg*(i -
      0.5) for i in 1:nSeg};
  Real RFluPip_value, Nu_value, h_value, Re_value, NuTurb_value;

  IDEAS.Fluid.Sensors.TemperatureTwoPort TBorFieIn(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    tau=0) "Inlet temperature of the borefield"
    annotation (Placement(transformation(extent={{20,38},{40,18}})));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TBorFieOut(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    tau=0)
    "Outlet temperature of the borefield"
    annotation (Placement(transformation(extent={{80,38},{100,18}})));
  IDEAS.Fluid.HeatExchangers.HeaterCooler_u hea(
    redeclare package Medium = Medium,
    dp_nominal=10000,
    show_T=true,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    T_start=T_startAll,
    Q_flow_nominal=1,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    m_flow(start=borefield.borFieDat.conDat.mBorFie_flow_nominal))
                    "Heater"
    annotation (Placement(transformation(extent={{-40,18},{-20,38}})));

  IDEAS.Fluid.Sources.Boundary_pT bou(
    redeclare package Medium = Medium,
    T=T_startAll,
    nPorts=1)
    annotation (Placement(transformation(extent={{72,50},{92,70}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
    tableOnFile=true,
    tableName="data",
    offset={0},
    columns={2},
    fileName=Modelica.Utilities.Files.loadResource("modelica://BorSizingMod/Resources/Data/test1a.txt"))
    annotation (Placement(transformation(extent={{-74,50},{-54,70}})));

  IDEAS.Fluid.Sensors.TemperatureTwoPort TheaIn(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    tau=0) "Inlet temperature of the borefield"
    annotation (Placement(transformation(extent={{-72,18},{-52,38}})));
  Buildings.Fluid.Movers.FlowControlled_m_flow mov(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    addPowerToMedium=false,
    use_riseTime=false,
    energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    nominalValuesDefineDefaultPressureCurve=true,
    inputType=Buildings.Fluid.Types.InputType.Constant,
    dp_nominal=60E3)
    annotation (Placement(transformation(extent={{-12,18},{8,38}})));
  Buildings.Fluid.Geothermal.Borefields.OneUTube borefield(
    redeclare package Medium = Medium,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    tLoaAgg=300,
    forceGFunCalc=false,
    borFieDat(
      filDat(
        kFil=1.4,
        dFil=1925,
        cFil=2026),
      soiDat(
        kSoi=1.8,
        cSoi=1150,
        dSoi=1801),
      conDat(
        borCon=Buildings.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.SingleUTube,
        use_Rb=true,
        Rb=0.13,
        mBor_flow_nominal=0.44,
        hBor=60,
        rBor=0.075,
        dBor=4,
        nBor=1,
        cooBor=[0,0],
        rTub=0.0167,
        kTub=0.43,
        eTub=0.003,
        xC=0.075/2,
        dp_nominal=5e4)),
    TExt0_start=TExt0_start,
    TExt_start={TExt0_start for i in 1:nSeg},
    z0=z0,
    dT_dz=dT_dz)
    annotation (Placement(transformation(extent={{48,18},{68,38}})));

   Modelica.Blocks.Sources.RealExpression Nu(y=Nu_value) "Nusselt number"
    annotation (Placement(transformation(extent={{132,60},{152,80}})));
  Modelica.Blocks.Sources.RealExpression RFluPip(y=RFluPip_value)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{132,40},{152,60}})));
  Modelica.Blocks.Sources.RealExpression h(y=h_value)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{164,42},{184,62}})));
  Modelica.Blocks.Sources.RealExpression Re(y=Re_value) "Reynolds number"
    annotation (Placement(transformation(extent={{192,42},{212,62}})));
  Modelica.Blocks.Sources.RealExpression NuTurb(y=NuTurb_value)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{164,60},{184,80}})));
algorithm
  (RFluPip_value, Nu_value, h_value, Re_value, NuTurb_value) := IDEAS.Fluid.Geothermal.Borefields.Validation.Functions.convectionResistanceCircularPipe(
    hSeg=borefield.borHol.intHex[1].hSeg,
    rTub=borefield.borFieDat.conDat.rTub,
    eTub=borefield.borFieDat.conDat.eTub,
    kMed=borefield.borHol.intHex[1].kMed,
    muMed=borefield.borHol.intHex[1].muMed,
    cpMed=borefield.borHol.intHex[1].cpMed,
    m_flow=borefield.borHol.intHex[1].m1_flow,
    m_flow_nominal=borefield.borHol.intHex[1].m1_flow_nominal);
equation
  connect(bou.ports[1],TBorFieOut. port_b) annotation (Line(points={{92,60},{
          110,60},{110,28},{100,28}},
                                   color={0,127,255}));
  connect(combiTimeTable.y[1],hea. u) annotation (Line(points={{-53,60},{-48,60},
          {-48,34},{-42,34}},  color={0,0,127}));
  connect(hea.port_a, TheaIn.port_b)
    annotation (Line(points={{-40,28},{-52,28}},   color={0,127,255}));
  connect(TheaIn.port_a, TBorFieOut.port_b) annotation (Line(points={{-72,28},{
          -90,28},{-90,88},{110,88},{110,28},{100,28}},      color={0,127,255}));
  connect(hea.port_b, mov.port_a)
    annotation (Line(points={{-20,28},{-12,28}},   color={0,127,255}));
  connect(TBorFieIn.port_a, mov.port_b)
    annotation (Line(points={{20,28},{8,28}},    color={0,127,255}));
  connect(borefield.port_a, TBorFieIn.port_b)
    annotation (Line(points={{48,28},{40,28}},   color={0,127,255}));
  connect(TBorFieOut.port_a, borefield.port_b)
    annotation (Line(points={{80,28},{68,28}},   color={0,127,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{220,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            220,100}})),
    experiment(
      StopTime=864000,
      __Dymola_fixedstepsize=5,
      __Dymola_Algorithm="Dassl"));
end test1a;
