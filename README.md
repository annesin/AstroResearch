# AstroResearch

This program is designed to input a .txt file of parameters and create a triple system and a star. The program is used to simulate accretion of the star by the system.
For just a binary system, use systemrk.jl. For more particles, use Main.jl. For a nested binary system, use NestedBinary.jl.

Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.

## systemrk.jl

### The .txt file must be in the following format for the systemrk.jl:
  
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

