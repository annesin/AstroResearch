import LinearAlgebra.cross
#= Uncomment this to install matplotlib =#
using Pkg
Pkg.add("PyPlot")
using PyPlot 
#import matplotlib.pyplot
#plt = matplotlib.pyplot Not sure why the program works without actually defining plt, but it does=#
G = 2945.49
hParam = 0.001 #this helps decide what timestep to use when integrating values

function fileInput(file) #change initial conditions to m1, m2, semi-major axis, e, 
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
       fileText = parse.(Float64,split(readlines(file)[1],",")) #change to 2
       m₁ = fileText[1] 
       m₂ = fileText[2]
       a = fileText[3]
       e = fileText[4]
       t = fileText[5]
       M = m₁+m₂ #total mass
       a₁ = (a*m₂)/(m₁+m₂)
       a₂ = (a*m₁)/(m₁+m₂)
       x₁0 = -(a₁+a₁*e)
       y₁0 = 0
       v₁0 = 0
       x₂0 = (a₂+a₂*e)
       y₂0 = 0
       v₂0 = 0
       r = x₂0-x₁0
       v = sqrt(G*M*((2/r)-(1/a))) #reduced mass velocity
       w₁0 = (v*m₂)/M
       w₂0 = -(v*m₁)/M
       if t == 0.0
              t = 100*sqrt((a^3*4*pi^2)/(G*M)) #if no time is entered, then it just goes for 100 periods
       end
       return m₁, m₂, x₁0, y₁0, v₁0, w₁0, x₂0, y₂0, v₂0, w₂0, t, a
       end
    
function SystemRK(file)
       m₁, m₂, x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, a = fileInput(file) #these are initial conditions read from the .txt file
  
       Llist = [] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
       Elist = [] #same, but for energy
       X1 = [] #keeps track of the first body's x coordinate
       X2 = [] #similar for these
       Y1 = []
       Y2 = []

       t0 = 0

       r = sqrt((x₁-x₂)^2+(y₁-y₂)^2) #distance between two bodies
       v = sqrt(G*(m₁+m₂)*((2/r)-(1/a))) #velocity of reduced mass
       h = hParam*(r/v) #this calculates the initial timestep

       while t0<t #we now iterate t0 until we get to t. The number of steps we need to get there is unknown, since our timestep varies
              x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂ = RungeKutta(m₁, m₂, x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, h) #=numerically integrates the system=#

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

              r = sqrt((x₁-x₂)^2+(y₁-y₂)^2) #distance between two bodies
              v = sqrt(G*(m₁+m₂)*(abs((2/r)-(1/a)))) #velocity of reduced mass, needed because if r is about 2a, then because of numerical error, the system may crash because of an imaginary result

              h = hParam*(r/v) #this calculates the initial timestep

              t0 += h

              push!(Llist,L[3]) #we only add the z component to the list, since by right hand rule, the rotational momentum will always be in the z-direction
              push!(Elist,E)
              push!(X1,x₁)
              push!(X2,x₂)
              push!(Y1,y₁) #move L ane E calculation up here
              push!(Y2,y₂)
       end
              
              println("$t0 days later...")
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

function RungeKutta(m₁, m₂, x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, h) #numeric integrator

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
       return x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂
                               
end

function Plot(file, color) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
       x₁, y₁, v₁, w₁, x₂, y₂, v₂, w₂, t, Llist, Elist, X1, X2, Y1, Y2 = SystemRK(file)
       if color == "L"
              L0 = Llist[1]
              Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
              plt.plot(Llist) 
       elseif color == "E"
              E0 = Elist[1]
              Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
              plt.plot(Elist) #find out what the scale things are, actually change to deltaE/E0
       elseif color == "EL"
              L0 = Llist[1]
              Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
              plt.plot(Llist,linestyle="solid",color="green") 
              E0 = Elist[1]
              Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
              plt.plot(Elist,linestyle="solid",color="red") #find out what the scale things are, actually change to deltaE/E0
       else
              plt.plot(X1,Y1,linestyle="solid",color="red")
              plt.plot(X2,Y2,linestyle="solid",color=color)
              plt.axis("equal") #makes axes equal, especially helpful if orbits are highly elliptical
              E0 = Elist[1]
              L0 = Llist[1]
              Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E'
              Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
              println("The angular momentum varied by $(minimum(Llist)) to $(maximum(Llist)) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).")
       end
end