# AstroResearch

This program is designed to input a .txt file of parameters and create a triple system and a star. The program is used to simulate accretion of the star by the system.
For just a binary system, use systemrk.jl. For more particles, use Main.jl.

The .txt file must be in the following format for the systemrk.jl:
  
  There are six elements in one line, each seperated by a comma.
  
  Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.
  
  The first parameter is the mass of the primary body.
  
  The second parameter is the mass of the secondary body.
  
  The third parameter is the semi-major axis of the system.
  
  The fourth parameter is the eccentricity.
  
  The fifth parameter is the time that the simulation runs for. (if t=0, then the simulation will run for 100 periods)

  The sixth parameter is the timestep parameter, which determines the length of the timestep the program uses
  
  Example (black hole-star system): 
  1,8,4.20984,0.2,0,0.001
  (See systemrk_input.txt)

The .txt file must be in the following format for the Main.jl:

  The elements are in three lines, each element being seperated by a comma.

  Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.

  The first line is the list of masses of the bodies.

  The second line is the list of initial conditions: X1, V1, X2, V2, X3, V3, X4, V4
    (X are the x, y, and z components of the body's initial position, while V are the x, y, and z components of the body's initial velocity)
    (X4 and V4 are only present if a test particle is present)

  The third line consists of the total time the simulation runs for and the timestep parameter (see above)

  Example:
  1,0.000003003,0.00000003694
  0,0,0,0,0,0,215.032,0,0,0,3.7,0,215.585,0,0,0,3.826,0,0,1,0,0,0,50 
  365,0.00001
  Sun-Earth-Moon system with a test particle with a polar orbit around the Sun
