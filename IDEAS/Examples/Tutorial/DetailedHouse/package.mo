within IDEAS.Examples.Tutorial;
package DetailedHouse "Package with example for how to built a simple house, using the detailed building envelope component
  models, HVAC and other components"
  extends Modelica.Icons.ExamplesPackage;

  annotation (Documentation(info="<html>
<p>
This package contains examples with step-by-step instructions for how to build a system model
for a simple house with detailed building envelope component models, HVAC and other components.
It serves as a demonstration case of how the <code>IDEAS</code> library can be used.
</p>
<p>
The goal of this exercise is to become familiar with Modelica and the IDEAS library.
Since the IDEAS library components are typically used by combining several components graphically,
the use of equations falls outside of the scope of this exercise.
</p>
<p>
For this exercise you will create a model of a simple house,
consisting of a heating system, detailed building envelope component models, and a ventilation model.
The exercise starts from a template file that should not produce any errors.
This file will be extended in several steps, adding complexity.
In between each step the user should be able to simulate the model,
i.e., no errors should be produced and simulation results may be compared.
</p>
<p>
The model has been created in the following stages:
</p>
<ol>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse1\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse1</a>
represents a basic building model that includes a single zone, four walls, a window, a floor, and a ceiling
and serves as a starting model to implement the entire <code>DetailedHouse</code> model.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse2\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse2</a>
implements solar shading extensions by adding a screen to the window model.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse2\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse3</a>
adds continuous occupancy of 1 person and LED lighting for the zone.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse4\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse4</a>
adds an occupancy schedule that assigns two occupants during office hours and 
zero otherwise. This schedule will be created using extensions of the IDEAS 
templates and replaceable models.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse5\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse5</a>
This example teaches the use of the <code>RectangularZoneTemplate</code> .
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse6\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse6</a>
adds a HVAC system. The system consists of a water-water heat pump, radiators, a storage tank, circulation
pumps and a low temperature heat source for the heat pump.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse7\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse7</a>
adds a controller that disables the heat pump when the supply water
temperature exceeds 45◦C.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse8\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse8</a>
adds workflow automation, generating a CSV file at a userdefined location that contains data 
for each of the inputs of the block .
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse9\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse9</a>
adds a CO2-controlled ventilaation system. The occupancy model from <a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse4\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse4</a> is added to one zone and a 
fixed occupancy of 1 person to the other zone. The ventilation system
consists of two fans, two supply and two return air VAVs (Variable Air Volume), a heat recovery unit and an
outdoor air source. The control consists of PI controllers with a set point of 1000 ppm.
</li>
<li>
<a href=\"modelica://IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse10\">
IDEAS.Examples.Tutorial.DetailedHouse.DetailedHouse10</a>
removes small time constants from the problem to decrease computation time.
</li>
</ol>
<p>
For each stage, firstly the model part is qualitatively explained.
Next, the names of the required Modelica models (from the Modelica Standard Library and/or IDEAS library) are listed.
Finally, we provide high-level instructions of how to set up the model.
If these instructions are not clear immediately, have a look at the model documentation and at the type of connectors the model has,
try out some things, make an educated guess, etc.
Finally, we provide reference results that allow you to check if your implementation is correct.
Depending on the parameter values that you choose, results may differ.
</p>
<p>
The graphical representation of the final model is given below.
</p>
<p align=\"center\">
<img alt=\"Graphical representation of the final simple house model.\"
src=\"modelica://IDEAS/Resources/Images/Examples/Tutorial/DetailedHouse/DetailedHouse10.png\"/>
</p>
</html>"));
end DetailedHouse;
