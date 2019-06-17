# AstroResearch

This program is designed to input a .txt file of parameters and create a triple system and a star. The program is used to simulate accretion of the star by the system.
Currently, the program only can simulate a binary system without a test particle.

The .txt file must be in the following format:
  
  There are five elements, each seperated by a comma.
  
  Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.
  
  The first parameter is the mass of the primary body.
  
  The second parameter is the mass of the secondary body.
  
  The third parameter is the semi-major axis of the system.
  
  The fourth parameter is the eccentricity.
  
  The fifth parameter is the time that the simulation runs for. (if t=0, then the simulation will run for 100 periods)
  
  Example (black hole-star system): 
  1,8,4.20984,0.2,0
