# AstroResearch

This program is designed to input a .txt file of parameters and create a triple system and a star. The program is used to simulate accretion of the star by the system.
For just a binary system, use `systemrk.jl`. For more particles, use `Main.jl`. For a nested binary system, use `NestedBinary.jl`.

Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.

## systemrk.jl

  With `systemrk.jl`, the `Plot` function will plot the orbits of two bodies interacting with each other gravitationally. You can also plot the system's energy and/or angular momentum versus time to ensure that the simulation is precise enough, as well as the timestep used by the numerical integration during the course of the simulation.
  
  ```
  Plot("Input.txt","blue")
  ```
  will yield the plot of the orbits. The color can be changed to any color recognized by matplotlib, and will be the color of one of the orbits.
  
  ```
  julia> Plot("Input.txt","E")
  
  julia> Plot("Input.txt","L")
  
  julia> Plot("Input.txt","EL")
  ```
  will yield the energy vs. time plot, the angular momentum vs. time plot, or both on the same axes, respectively.
  
  ```
  Plot("Input.txt","time")
  ```
  will yield the timestep vs. iteration plot.
  
  Additionally, `Plot` can save a .txt file that can reproduce a graph later from the simulation (saving you from having to run it again). To do this, enter `1` as a third parameter. The .txt file produced will have a name consisting of the initial conditions. `ExternalPlotter.jl` can recreate the plots from the .txt file. This is especially useful to use alongside with `ExternalPlot` if you want to see all the graphs, as `ExternalPlot` can create any of the above mentioned graphs with the .txt file created this way.
  
  ```
  julia> Plot("Input.txt","blue",1)
  
  julia> ExternalPlot("Output.txt","time")
  
  julia> ExternalPlot("Output.txt","EL")
  ```
  will produce both the orbital graph from `Plot`, then the timestep graph and the E+L graph from the two `ExternalPlot`s.


### The .txt file must be in the following format for systemrk.jl:
  
  There are six elements in one line, each seperated by a comma:
  1. Mass of the primary body
  2. Mass of the secondary body
  3. Semi-major axis of the system
  4. Eccentricity
  5. Time that the simulation runs for
     - (if t=0, then the simulation will run for 100 periods)
  6. Timestep parameter, which determines the length of the timestep the program uses
  
  Example (black hole-star system): 
  ```
  1,8,4.20984,0.2,0,0.001
  ```

## Main.jl

  `Plot` has the same functionality here as `systemrk.jl`, but can include three massive bodies and a massless test particle. The syntax is exactly the same. Notice that 3D plotting may occur here, since two bodies' orbits can always be displayed in a plane, while three bodies may not be. `Plot` detects this automatically, and will plot the orbits on a plane if possible.
  
  `Main` also has the option to manually chose your .txt save file for `ExternalPlot`.
  
  ```
  Plot("Input.txt","blue","CoolOrbits.txt")
  ```
  would plot the orbits, as well as return the .txt data file `CoolOrbits.txt`.
  
  ```
  ExternalPlot("CoolOrbits.txt","blue")
  ```
  would recreate this.

### The .txt file must be in the following format for Main.jl:

  The elements are in three lines, each element being seperated by a comma:
  1. List of masses of the bodies
  2. List of initial conditions: X1, V1, X2, V2, X3, V3, X4, V4
     - X's are the x, y, and z components of a body's initial position
     - V's are the x, y, and z components of a body's initial velocity
     - X4 and V4 are only present if a test particle is present
  3. Total time the simulation runs for, and the timestep parameter (see above)

  Example:
  ```
  1,0.000003003,0.00000003694
  0,0,0,0,0,0,215.032,0,0,0,3.7,0,215.585,0,0,0,3.826,0,0,1,0,0,0,50 
  365,0.00001
  #Sun-Earth-Moon system with a test particle with a polar orbit around the Sun
  ```

## NestedBinary.jl

  `NestedBinary` is essentially `Main`, except it's designed with specific three-body setup: two bodies form an inner binary, while the third body orbits from a certain distance away. 
  
  `Plot` also has the option `"none"`, where nothing is plotted.
  
  `NestedBinary` also comes with the `AutomaticTester` function, which  simulates a variety of initial conditions that are typical of observed triple systems. It then saves the data in an Excel spreadsheet, as well as compressed .txt files for later graphing. The .txt files are in the format `Name_n.txt`, where `n` is the Excel row number where that simulation was recorded.
  
  ```
  Plot("Input.txt","none","CoolOrbits.txt")
  ```
  would do the same thing as the above `Main` example.
  
  ```
  AutomaticTester("AutoRun")
  ```
  would produce a spreadsheet of simulations, whose .txt files would have the name `AutoRun_n.txt`.

### The .txt file must be in the following format for NestedBinary.jl:

  The elements are in three lines, each element being seperated by a comma:
  1. The first line is the list of masses of the bodies.
  2. The second line is the list of initial conditions: a₁, e₁, a₂, e₂, i, and Θ.
     - a₁ and e₁ are the semi-major axis and eccentricity of the inner binary's orbit.
     - a₂ and e₂ are the semi-major axis and eccentricity of the outer binary's orbit. 
     - i is the inclination of the third body's orbit. It's measured from the plane formed by the inner binary.
     - Θ is the angle between the line formed by the two inner bodies and the line formed from the inner binary's center of mass to the projection of the third body onto the above mentioned plane. It is also the longitude of the ascending node minus 270 degrees.
  3. The third line consists of the time the simulation should run for and the timestep parameter (see above).

  Example:
  ```
  8,8,1
  2,0,20,0.99,0,0
  10,0.0001
  ```
  
## ExternalPlotter.jl
 
  This contains the function `ExternalPlot`, which can recreate any plot created by `Plot` in `systemrk`, `Main`, or `NestedBinary`. The syntax is the same as `Plot`.
  
  ```
  ExternalPlot("Output.txt","blue")
  ```
  would return the orbits from the simulation that created `Output.txt`.

## Calculations

  The equations we used, along with our derivations that we used in this program, are explained in `Equations.pdf`. The LaTeX file for `Equations.pdf` is also avaliable in `main.tex`.
  
 ## Credits
 This code was produced for Prof. Ann Esin's summer research in astrophysics.
 This code was written by Brennen Quigley (@bquigley528) and Benjamin Khoury (@bbk001).
 Thank you to @timholy for their `ProgressMeter` package and @felipenoris for their `XLSX` package.
