import LinearAlgebra.cross
#= Uncomment this to install matplotlib =#
using Pkg
Pkg.add("PyPlot")
using PyPlot 
#import matplotlib.pyplot
#plt = matplotlib.pyplot Not sure why the program works without actually defining plt, but it does
G = 2945.49

function fileInput(file)
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
       fileText = parse.(Float64,split(readlines(file)[1],","))
       m₁ = fileText[1] 
       m₂ = fileText[2]
       x₁0 = fileText[3]
       y₁0 = fileText[4]
       v₁0 = fileText[5]
       w₁0 = fileText[6]
       x₂0 = fileText[7]
       y₂0 = fileText[8]
       v₂0 = fileText[9]
       w₂0 = fileText[10]
       t = fileText[11]
       h = fileText[12]
       return m₁, m₂, x₁0, y₁0, v₁0, w₁0, x₂0, y₂0, v₂0, w₂0, t, h
       end
    
function SystemRK(file)

                            m₁, m₂, x₁0, y₁0, v₁0, w₁0, x₂0, y₂0, v₂0, w₂0, t, h = fileInput(file)

	 	      	      #= This first section of code serves to set the initial conditions and prepare the variables for the upcoming loop=#

                                 x₁ = x₁0
                                 y₁ = y₁0
                                 v₁ = v₁0
                                 w₁ = w₁0
                                 x₂ = x₂0
                                 y₂ = y₂0
                                 v₂ = v₂0
                                 w₂ = w₂0
				     t0 = 0
  
                            Llist = [] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
                            Elist = [] #same, but for energy
                            X1 = [] #keeps track of the first body's x coordinate
                            X2 = [] #similar for these
                            Y1 = []
                            Y2 = []

				for i = 1:(t/h) #=The number of iterations depends on how many steps it takes to reach the desired time t=#
                                	x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, L, E = RungeKutta(m₁, m₂, x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, h, t0) #=numerically integrates the system=#
                                   push!(Llist,L[3]) #we only add the z component to the list, since by right hand rule, the rotational momentum will always be in the z-direction
                                   push!(Elist,E)
                                   push!(X1,x₁)
                                   push!(X2,x₂)
                                   push!(Y1,y₁)
                                   push!(Y2,y₂)
				end
                                   

                                return x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, Llist, Elist, X1, X2, Y1, Y2 #returns the desired values
end

function dEdt(x₁, x₂, y₁, y₂, m₁, m₂) #DE for v1
       g = G * m₂ *(x₂-x₁)
       r = ((x₂-x₁)^2 + (y₂-y₁)^2)^0.5
       return g/(r^3)
       end

 function dFdt(x₁, x₂, y₁, y₂, m₁, m₂) #DE for v2
       g = G * m₁ *(x₁-x₂)
       r = ((x₂-x₁)^2 + (y₂-y₁)^2)^0.5
       return g/(r^3)
       end

function dGdt(x₁, x₂, y₁, y₂, m₁, m₂) #DE for w1
       g = G * m₂ *(y₂-y₁)
       r = ((x₂-x₁)^2 + (y₂-y₁)^2)^0.5
       return g/(r^3)
       end

function dHdt(x₁, x₂, y₁, y₂, m₁, m₂) #DE for w2
       g = G * m₁ *(y₁-y₂)
       r = ((x₂-x₁)^2 + (y₂-y₁)^2)^0.5
       return g/(r^3)
       end

function RungeKutta(m₁, m₂, x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, h, t0) #numeric integrator

                                 A1 = h * v₁ #change in x1
                                 B1 = h * v₂ #change in x2
                                 C1 = h * w₁ #change in y1
                                 D1 = h * w₂ #change in y2
                                 E1 = h * dEdt(x₁, x₂, y₁, y₂, m₁, m₂) #change in v1
                                 F1 = h * dFdt(x₁, x₂, y₁, y₂, m₁, m₂) #change in v2
                                 G1 = h * dGdt(x₁, x₂, y₁, y₂, m₁, m₂) #change in w1
                                 H1 = h * dHdt(x₁, x₂, y₁, y₂, m₁, m₂) #change in w2
                                 A2 = h * (v₁ + 0.5*E1) #the next lines follow Runge Kutta formulism
                                 B2 = h * (v₂ + 0.5*F1)
                                 C2 = h * (w₁ + 0.5*G1)
                                 D2 = h * (w₂ + 0.5*H1)
                                 E2 = h * dEdt(x₁ + 0.5*A1, x₂ + 0.5*B1, y₁ + 0.5*C1, y₂ + 0.5*D1, m₁, m₂)
                                 F2 = h * dFdt(x₁ + 0.5*A1, x₂ + 0.5*B1, y₁ + 0.5*C1, y₂ + 0.5*D1, m₁, m₂)
                                 G2 = h * dGdt(x₁ + 0.5*A1, x₂ + 0.5*B1, y₁ + 0.5*C1, y₂ + 0.5*D1, m₁, m₂)
                                 H2 = h * dHdt(x₁ + 0.5*A1, x₂ + 0.5*B1, y₁ + 0.5*C1, y₂ + 0.5*D1, m₁, m₂)
                                 A3 = h * (v₁ + 0.5*E2)
                                 B3 = h * (v₂ + 0.5*F2)
                                 C3 = h * (w₁ + 0.5*G2)
                                 D3 = h * (w₂ + 0.5*H2)
                                 E3 = h * dEdt(x₁ + 0.5*A2, x₂ + 0.5*B2, y₁ + 0.5*C2, y₂ + 0.5*D2, m₁, m₂)
                                 F3 = h * dFdt(x₁ + 0.5*A2, x₂ + 0.5*B2, y₁ + 0.5*C2, y₂ + 0.5*D2, m₁, m₂)
                                 G3 = h * dGdt(x₁ + 0.5*A2, x₂ + 0.5*B2, y₁ + 0.5*C2, y₂ + 0.5*D2, m₁, m₂)
                                 H3 = h * dHdt(x₁ + 0.5*A2, x₂ + 0.5*B2, y₁ + 0.5*C2, y₂ + 0.5*D2, m₁, m₂)
                                 A4 = h * (v₁ + E3)
                                 B4 = h * (v₂ + F3)
                                 C4 = h * (w₁ + G3)
                                 D4 = h * (w₂ + H3)
                                 E4 = h * dEdt(x₁ + A3, x₂ + B3, y₁ + C3, y₂ + D3, m₁, m₂)
                                 F4 = h * dFdt(x₁ + A3, x₂ + B3, y₁ + C3, y₂ + D3, m₁, m₂)
                                 G4 = h * dGdt(x₁ + A3, x₂ + B3, y₁ + C3, y₂ + D3, m₁, m₂)
                                 H4 = h * dHdt(x₁ + A3, x₂ + B3, y₁ + C3, y₂ + D3, m₁, m₂)

                                 x₁ = x₁ + (1.0/6.0)*(A1 + 2*A2 + 2*A3 + A4) #=weighted sums of the differences in variable values=#
                                 x₂ = x₂ + (1.0/6.0)*(B1 + 2*B2 + 2*B3 + B4)
                                 y₁ = y₁ + (1.0/6.0)*(C1 + 2*C2 + 2*C3 + C4)
                                 y₂ = y₂ + (1.0/6.0)*(D1 + 2*D2 + 2*D3 + D4)
                                 v₁ = v₁ + (1.0/6.0)*(E1 + 2*E2 + 2*E3 + E4)
                                 v₂ = v₂ + (1.0/6.0)*(F1 + 2*F2 + 2*F3 + F4)
                                 w₁ = w₁ + (1.0/6.0)*(G1 + 2*G2 + 2*G3 + G4)
                                 w₂ = w₂ + (1.0/6.0)*(H1 + 2*H2 + 2*H3 + H4)

                            r₁ = [x₁,y₁,0]
                            r₂ = [x₂,y₂,0]
                            r = sqrt((x₁-x₂)^2+(y₁-y₂)^2) #distance between two bodies

                            V₁ = [v₁,w₁,0] #vector form of velocity
                            V₂ = [v₂,w₂,0]

                            L₁ = m₁ * (cross(r₁,V₁))
                            L₂ = m₂ * (cross(r₂,V₂))
                            L = L₁+L₂
                            K₁ = .5*m₁*(v₁^2+w₁^2) #kinetic energies
                            K₂ = .5*m₂*(v₂^2+w₂^2)
                            P = (-G*m₁*m₂)/r #gravitational potential
                            E = K₁ + K₂ + P #total energy

                            t0 = t0 + h #time step
                                 
                            return x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, L, E
                               
       end

function Plot(file, color) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
       x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, Llist, Elist, X1, X2, Y1, Y2 = SystemRK(file)
       if color == "L"
              plt.plot(Llist) 
       elseif color == "E"
              plt.plot(Elist)
       else
              plt.plot(X1,Y1,linestyle="solid",color="yellow")
              plt.plot(X2,Y2,linestyle="solid",color=color)
       end
end