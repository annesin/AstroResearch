# AstroResearch

This program is designed to input a .txt file of parameters and create a triple system and a star. The program is used to simulate accretion of the star by the system.
Currently, the program only can simulate a binary system without a test particle.

The .txt file must be in the following format:

  There is only one line.
  
  There are twelve elements (not in array, so no brackets) each seperated by a comma.
  
  Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.
  
  The first parameter is the mass of the primary body.
  
  The second parameter is the mass of the secondary body.
  
  The third parameter is the initial x-coordinate of the primary body.
  
  The fourth parameter is the initial y-coordinate of the primary body.
  
  The fifth parameter is the initial x-component of the velocity of the primary body.
  
  The sixth parameter is the initial y-component of the velcoity of the secondary body.
  
  The seventh parameter is the initial x-coordinate of the secondary body.
  
  The eighth parameter is the initial y-coordinate of the secondary body.
  
  The ninth parameter is the initial x-component of the velocity of the secondary body.
  
  The tenth parameter is the initial y-component of the velocity of the secondary body.
  
  The eleventh parameter is the time that the simulation runs for.
  
  The twelfth parameter is the time-step value.
  
  Example (using Sun-Jupiter system): 
  1,0.0009543,-1.068,0,0,-0.001548,1119.12,0,0,1.62318,4332.71,1
