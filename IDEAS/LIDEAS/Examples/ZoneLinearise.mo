within IDEAS.LIDEAS.Examples;
model ZoneLinearise "Test for linearisation"
  extends LIDEAS.Components.LinearisationInterface(sim(nWindow=1));
  parameter Integer nZones=2 "Number of zone";
  parameter Integer nEmb=1 "Number of embbed systems";
  package Medium = IDEAS.Media.Air;
  extends Modelica.Icons.Example;
  LIDEAS.Components.LinZone zone(
    redeclare package Medium = Medium,
    allowFlowReversal=true,
    V=20,
    nSurf=5) "First zone"
    annotation (Placement(transformation(extent={{20,-20},{40,0}})));
  IDEAS.Buildings.Components.BoundaryWall commonWall(
    redeclare parameter IDEAS.Buildings.Validation.Data.Constructions.HeavyWall
      constructionType,
    A=10,
    azi=0,
    inc=1.5707963267949,
    use_T_in=true) "Common wall model"
    annotation (Placement(transformation(extent={{-54,-2},{-44,18}})));
  IDEAS.Buildings.Components.InternalWall floor(
    A=10,
    azi=0,
    inc=IDEAS.Types.Tilt.Floor,
    redeclare IDEAS.Buildings.Data.Constructions.TABS constructionType)
    "TABS Floor between the two zones." annotation (Placement(transformation(
        extent={{-5,-10},{5,10}},
        rotation=90,
        origin={9,-38})));

  LIDEAS.Components.LinWindow window(
    A=1,
    redeclare parameter IDEAS.Buildings.Data.Glazing.Ins2 glazing,
    redeclare IDEAS.Buildings.Components.Shading.Screen shaType,
    redeclare IDEAS.Buildings.Data.Frames.Pvc fraType,
    inc=IDEAS.Types.Tilt.Wall,
    azi=IDEAS.Types.Azimuth.S,
    indexWindow=1) "Window model"
    annotation (Placement(transformation(extent={{-54,-82},{-44,-62}})));
  IDEAS.Buildings.Components.SlabOnGround slabOnGround(
    redeclare parameter IDEAS.Buildings.Validation.Data.Constructions.LightWall
      constructionType,
    A=20,
    PWall=3,
    inc=0,
    azi=0) "Floor model"
    annotation (Placement(transformation(extent={{-54,20},{-44,40}})));
  IDEAS.Buildings.Components.OuterWall outerWall(
    azi=0,
    redeclare parameter IDEAS.Buildings.Validation.Data.Constructions.HeavyWall
      constructionType,
    A=10,
    inc=1.5707963267949) "Outer wall model"
    annotation (Placement(transformation(extent={{-54,-58},{-44,-38}})));
  LIDEAS.Components.LinZone zone1(
    nSurf=2,
    redeclare package Medium = Medium,
    allowFlowReversal=true,
    V=20) "Second zone"
    annotation (Placement(transformation(extent={{20,-70},{40,-50}})));
  IDEAS.Buildings.Components.Shading.ShadingControl shadingControl
    "Shading control model"
    annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
  IDEAS.Buildings.Components.BoundaryWall commonWall1(
    redeclare parameter IDEAS.Buildings.Validation.Data.Constructions.HeavyWall
      constructionType,
    A=10,
    azi=0,
    inc=1.5707963267949,
    use_T_in=false,
    use_Q_in=true) "Common wall model"
    annotation (Placement(transformation(extent={{-54,-26},{-44,-6}})));

  Modelica.Blocks.Sources.Constant const(k=273.15 + 20) "Temperature at one side of the common wall"
    annotation (Placement(transformation(extent={{-100,0},{-80,20}})));
  Modelica.Blocks.Sources.Constant const1(k=10) "Constant heat flow at one side of the common wall 1"
    annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
equation
  connect(commonWall.propsBus_a, zone.propsBus[1]) annotation (Line(
      points={{-44.8333,10},{-12,10},{-12,-4.4},{20,-4.4}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(floor.propsBus_a, zone.propsBus[2]) annotation (Line(
      points={{7,-33.8333},{6,-33.8333},{6,-5.2},{20,-5.2}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(window.propsBus_a, zone1.propsBus[2]) annotation (Line(
      points={{-44.8333,-70},{6,-70},{6,-57},{20,-57}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(outerWall.propsBus_a, zone.propsBus[4]) annotation (Line(
      points={{-44.8333,-46},{-12,-46},{-12,-6.8},{20,-6.8}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(slabOnGround.propsBus_a, zone.propsBus[3]) annotation (Line(
      points={{-44.8333,32},{-12,32},{-12,-6},{20,-6}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(floor.propsBus_b, zone1.propsBus[1]) annotation (Line(
      points={{7,-42.1667},{6.5,-42.1667},{6.5,-55},{20,-55}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shadingControl.y, window.Ctrl) annotation (Line(
      points={{-60,-84},{-58,-84},{-58,-86},{-52.3333,-86},{-52.3333,-82}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(commonWall1.propsBus_a, zone.propsBus[5]) annotation (Line(
      points={{-44.8333,-14},{-12,-14},{-12,-7.6},{20,-7.6}},
      color={255,204,51},
      thickness=0.5));
  connect(const.y, commonWall.T)
    annotation (Line(points={{-79,10},{-57.8333,10}},
                                                   color={0,0,127}));
  connect(const1.y, commonWall1.Q_flow) annotation (Line(points={{-79,-20},{
          -57.8333,-20},{-57.8333,-18}},
                                   color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),           __Dymola_Commands(file=
          "Scripts/linearize_ZoneLinearise.mos" "Linearise"),
    Documentation(revisions="<html>
<ul>
<li>April 10, 2018 by Damien Picard: <br>Add documentation.</li>
<li>July 18, 2016 by Filip Jorissen:<br>Cleaned up code and implementation. </li>
<li>By Filip Jorissen:<br>First implementation. </li>
</ul>
</html>", info="<html>
<p>This model represents a very simple two-zone building. The model can be linearised by executing the command <i>Linearise</i> 
in the <i>Commands</i> menu of Dymola. This executes the script at LIDEAS/LIDEAS/Resources/Scripts/linearize.mos, 
which contains the following code:</p>
<p>
<br>OutputCPUtime:=false;    </br>
<br>  </br>
<br>// linearise model and save results in \"re\" </br>
<br>re=Modelica_LinearSystems2.ModelAnalysis.Linearize(\"IDEAS.LIDEAS.Examples.ZoneLinearise\"); </br>
<br>  </br>
<br>// extract state space matrices from 're' and save them to a mat file in the current working directory</br>
<br>writeMatrix(fileName=\"ssm_ZoneLinearise.mat\",matrixName=\"A\",matrix=re.A);        </br>
<br>writeMatrix(fileName=\"ssm_ZoneLinearise.mat\",matrixName=\"B\",matrix=re.B, append=true);     </br>
<br>writeMatrix(fileName=\"ssm_ZoneLinearise.mat\",matrixName=\"C\",matrix=re.C, append=true);     </br>
<br>writeMatrix(fileName=\"ssm_ZoneLinearise.mat\",matrixName=\"D\",matrix=re.D, append=true);     </br>
<br>  </br>
<br>// save the variable names of the inputs, outputs and states in the current working directory      </br>
<br>Modelica.Utilities.Files.remove(\"uNames_ZoneLinearise.txt\");</br>
<br>for i in 1:size(re.uNames,1) loop Modelica.Utilities.Streams.print(re.uNames[i], \"uNames_ZoneLinearise.txt\"); end for;     </br>
<br>Modelica.Utilities.Files.remove(\"yNames_ZoneLinearise.txt\");</br>
<br>for i in 1:size(re.yNames,1) loop Modelica.Utilities.Streams.print(re.yNames[i], \"yNames_ZoneLinearise.txt\"); end for;     </br>
<br>Modelica.Utilities.Files.remove(\"xNames_ZoneLinearise.txt\");</br>
<br>for i in 1:size(re.xNames,1) loop Modelica.Utilities.Streams.print(re.xNames[i], \"xNames_ZoneLinearise.txt\"); end for;      </br>
<br>  </br>
<br>OutputCPUtime:=true;     </br>
</p> 
<p>The first and last line of the code avoid that the CPU time is also linearised. 
The second line performs the actual linearisation and stores the resulting state space model in variable 
<i>re</i>. The next four lines are used to save the four state space matrices contained by <i>re</i> 
to a file <i>ssm.mat</i>. The remaining four lines of code also export the variable names of the 
state space inputs (u), states (x) and outputs (y) to three .txt files, which can then be used externally.</p>
</html> "));
end ZoneLinearise;
