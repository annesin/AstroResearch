# AstroResearch

This branch consists of code that is currently needing to be completed. As of now, this file only consists of one program: `AutomaticTesterIslands.jl`

Masses are in units of solar mass, distances are in the units of solar radii, and time is in the units of days.

## AutomaticTesterIslands.jl

Over the course of our examinations of what separations of the inner and outer binaries yield stable systems, we discovered that there exists several systems that remained stable that had separations smaller than one would expect for stability. 

  With `AutomaticTesterIslands.jl`, the program should identitify any of these islands of stablility. This program would keep track of exactly what parameters yielded a stable system where one was not expected. These parameters would then be stored in `Islands.xlsx` for later inspection.

  
 ## Credits
 This code was produced for Prof. Ann Esin's summer research in astrophysics.
 This code was written by Brennen Quigley (@bquigley528) and Benjamin Khoury (@bbk001).
 Thank you to @timholy for their `ProgressMeter` package and @felipenoris for their `XLSX` package.
