within IDEAS.Utilities.IO.Files;
model WeeklySchedule "Weekly schedule model"
  extends Modelica.Blocks.Icons.DiscreteBlock;
  parameter Integer[:] columns
    "Columns that should be loaded from file";
  parameter String fileName
    "Filename";
  parameter Modelica.SIunits.Time t_offset = 0
    "Timestamp that corresponds to Monday midnight";

  Modelica.Blocks.Interfaces.RealOutput[n_columns] y = {getCalendarValue(cal, iCol-1, time) for iCol in columns} "Outputs"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  IDEAS.Utilities.IO.Files.BaseClasses.WeeklyScheduleObject cal=
      IDEAS.Utilities.IO.Files.BaseClasses.WeeklyScheduleObject(fileName, t_offset)
    "Schedule object";
protected
  parameter Integer n_columns = size(columns,1) "Number of columns";
  function getCalendarValue
    "Returns the interpolated (zero order hold) value"
    extends Modelica.Icons.Function;
    input IDEAS.Utilities.IO.Files.BaseClasses.WeeklyScheduleObject ID "Pointer to file writer object";
    input Integer iCol;
    input Real timeIn;
    output Real y;
    external "C" y=getScheduleValue(ID, iCol, timeIn)
    annotation(Include=" #include <WeeklySchedule.c>",
    IncludeDirectory="modelica://IDEAS/Resources/C-Sources");
  end getCalendarValue;

  annotation (experiment(
      StartTime=-10000,
      StopTime=1000000,
      Interval=100),                                     Documentation(
        revisions="<html>
<ul>
<li>
March 9 2022, by Filip Jorissen:<br/>
First implementation.
</li>
</ul>
</html>", info="<html>
<p>
This model interprets a specific file and performs a cyclic (weekly) extrapolation on the source data.
See <a href=\"modelica://IDEAS/Resources/Data/schedule.txt\">IDEAS/Resources/Data/schedule.txt</a> 
for an example of the supported file format.
</p>
</html>"));
end WeeklySchedule;
