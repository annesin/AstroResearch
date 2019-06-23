using Pkg
Pkg.add("PyPlot")
using PyPlot 

function ExternalPlot(file, color) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
    file = "h≈(r÷v) data files/$file"
    Llist = parse.(Float64,split(readlines(file)[1],","))
    Elist = parse.(Float64,split(readlines(file)[2],",")) 
    Tlist = parse.(Float64,split(readlines(file)[3],",")) 
    lList = parse.(Float64,split(readlines(file)[4],",")) 
    X1 = parse.(Float64,split(readlines(file)[5],",")) 
    X2 = parse.(Float64,split(readlines(file)[6],",")) 
    Y1 = parse.(Float64,split(readlines(file)[7],",")) 
    Y2 = parse.(Float64,split(readlines(file)[8],",")) 
    numBodies = parse.(Int64,readlines(file)[9])
    if color == "L"
           L0 = Llist[1]
           Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
           plt.plot(Tlist,Llist) 
    elseif color == "E"
           E0 = Elist[1]
           Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
           plt.plot(Tlist,Elist) #find out what the scale things are, actually change to deltaE/E0
    elseif color == "EL"
           L0 = Llist[1]
           Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
           plt.plot(Tlist,Llist,linestyle="solid",color="green") 
           E0 = Elist[1]
           Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
           plt.plot(Tlist,Elist,linestyle="solid",color="red") #find out what the scale things are, actually change to deltaE/E0
           println("The angular momentum varied by $(minimum(Llist)) to $(maximum(Llist)) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).")
    elseif color == "time"
           plt.plot(lList,linestyle="solid",color="green")
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