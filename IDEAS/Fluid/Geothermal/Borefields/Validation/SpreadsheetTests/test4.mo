within IDEAS.Fluid.Geothermal.Borefields.Validation.SpreadsheetTests;
model test4
  import Buildings;
  extends Modelica.Icons.Example;
  package Medium =
      BorSizingMod.Validation.Media.PropyleneGlycolWater_test2 (
       property_T=273.15,X_a=0);
  parameter Integer nSeg = 10;
  parameter Modelica.Units.SI.Temperature T_startAll = 273.15 + 15;
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
    annotation (Placement(transformation(extent={{10,-10},{30,-30}})));
  IDEAS.Fluid.Sensors.TemperatureTwoPort TBorFieOut(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    tau=0)
    "Outlet temperature of the borefield"
    annotation (Placement(transformation(extent={{70,-10},{90,-30}})));
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
    annotation (Placement(transformation(extent={{-50,-30},{-30,-10}})));

  Buildings.Fluid.Geothermal.Borefields.OneUTube borefield(
    redeclare package Medium = Medium,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    tLoaAgg=300,
    forceGFunCalc=false,
    nClu=4,
    borFieDat(
      filDat(
        kFil=0.69,
        dFil=1925,
        cFil=2026),
      soiDat(
        kSoi=1.9,
        cSoi=1145,
        dSoi=1792),
      conDat(
        borCon=Buildings.Fluid.Geothermal.Borefields.Types.BoreholeConfiguration.SingleUTube,
        use_Rb=true,
        Rb=0.2,
        mBor_flow_nominal=0.074*139.731/25,
        hBor=100,
        rBor=0.075,
        dBor=4,
        cooBor=[0,0; 8,0; 16,0; 24,0; 32,0; 0,8; 8,8; 16,8; 24,8; 32,8; 0,16; 8,
            16; 16,16; 24,16; 32,16; 0,24; 8,24; 16,24; 24,24; 32,24; 0,32; 8,32;
            16,32; 24,32; 32,32],
        rTub=0.0167,
        kTub=0.4,
        eTub=0.0037,
        xC=0.083/2,
        dp_nominal=5e4)),
    TExt0_start=TExt0_start,
    TExt_start={TExt0_start for i in 1:nSeg},
    z0=z0,
    dT_dz=dT_dz)
    annotation (Placement(transformation(extent={{40,-30},{60,-10}})));

  IDEAS.Fluid.Sources.Boundary_pT bou(
    redeclare package Medium = Medium,
    T=T_startAll,
    nPorts=1)
    annotation (Placement(transformation(extent={{62,2},{82,22}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
    tableOnFile=true,
    tableName="data",
    offset={0},
    columns={2},
    fileName=Modelica.Utilities.Files.loadResource("modelica://BorSizingMod/Resources/Data/test4.txt"))
    annotation (Placement(transformation(extent={{-84,2},{-64,22}})));

  IDEAS.Fluid.Sensors.TemperatureTwoPort TheaIn(
    redeclare package Medium = Medium,
    T_start=T_startAll,
    m_flow_nominal=borefield.borFieDat.conDat.mBorFie_flow_nominal,
    tau=0) "Inlet temperature of the borefield"
    annotation (Placement(transformation(extent={{-82,-30},{-62,-10}})));
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
    annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
   Modelica.Blocks.Sources.RealExpression Nu(y=Nu_value) "Nusselt number"
    annotation (Placement(transformation(extent={{-68,-72},{-48,-52}})));
  Modelica.Blocks.Sources.RealExpression RFluPip(y=RFluPip_value)
    "Convection resistance (or conduction in fluid if no mass flow)"
    annotation (Placement(transformation(extent={{-68,-94},{-48,-74}})));
  Modelica.Blocks.Sources.RealExpression h(y=h_value)
    "Convective heat transfer coefficient of the fluid"
    annotation (Placement(transformation(extent={{-32,-94},{-12,-74}})));
  Modelica.Blocks.Sources.RealExpression Re(y=Re_value) "Reynolds number"
    annotation (Placement(transformation(extent={{6,-94},{26,-74}})));
  Modelica.Blocks.Sources.RealExpression NuTurb(y=NuTurb_value)
    "Nusselt at Re=2400"
    annotation (Placement(transformation(extent={{-32,-72},{-12,-52}})));
algorithm
  (RFluPip_value, Nu_value, h_value, Re_value, NuTurb_value) := BorSizingMod.Functions.convectionResistanceCircularPipe(
    hSeg=borefield.borHol.intHex[1].hSeg,
    rTub=borefield.borFieDat.conDat.rTub,
    eTub=borefield.borFieDat.conDat.eTub,
    kMed=borefield.borHol.intHex[1].kMed,
    muMed=borefield.borHol.intHex[1].muMed,
    cpMed=borefield.borHol.intHex[1].cpMed,
    m_flow=borefield.borHol.intHex[1].m1_flow,
    m_flow_nominal=borefield.borHol.intHex[1].m1_flow_nominal);
equation
  connect(TBorFieIn.port_b, borefield.port_a)
    annotation (Line(points={{30,-20},{40,-20}}, color={0,127,255}));
  connect(TBorFieOut.port_a, borefield.port_b)
    annotation (Line(points={{70,-20},{60,-20}}, color={0,127,255}));
  connect(bou.ports[1],TBorFieOut. port_b) annotation (Line(points={{82,12},{
          100,12},{100,-20},{90,-20}},
                                   color={0,127,255}));
  connect(combiTimeTable.y[1],hea. u) annotation (Line(points={{-63,12},{-58,12},
          {-58,-14},{-52,-14}},color={0,0,127}));
  connect(hea.port_a, TheaIn.port_b)
    annotation (Line(points={{-50,-20},{-62,-20}}, color={0,127,255}));
  connect(TheaIn.port_a, TBorFieOut.port_b) annotation (Line(points={{-82,-20},
          {-100,-20},{-100,40},{100,40},{100,-20},{90,-20}}, color={0,127,255}));
  connect(mov.port_a, hea.port_b)
    annotation (Line(points={{-20,-20},{-30,-20}}, color={0,127,255}));
  connect(mov.port_b, TBorFieIn.port_a)
    annotation (Line(points={{0,-20},{10,-20}}, color={0,127,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{120,
            100}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            120,100}})),
    experiment(
      StopTime=864000,
      Interval=1.00000224,
      __Dymola_fixedstepsize=30,
      __Dymola_Algorithm="Dassl"));
end test4;
