import LinearAlgebra.norm
using Pkg
Pkg.add("PyPlot")
using PyPlot 

function ExternalPlot(file, color, equal=0) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
       file = "h≈(r÷v) data files/$file"
       numBodies = parse.(Int64,readlines(file)[1])
       Llist = parse.(Float64,split(readlines(file)[2],","))
       Elist = parse.(Float64,split(readlines(file)[3],",")) 
       Tlist = parse.(Float64,split(readlines(file)[4],",")) 
       lList = parse.(Float64,split(readlines(file)[5],",")) 
       X1 = parse.(Float64,split(readlines(file)[6],",")) 
       X2 = parse.(Float64,split(readlines(file)[7],",")) 
       Y1 = parse.(Float64,split(readlines(file)[8],",")) 
       Y2 = parse.(Float64,split(readlines(file)[9],",")) 
       if numBodies > 2
              Z1 = parse.(Float64,split(readlines(file)[10],","))
              Z2 = parse.(Float64,split(readlines(file)[11],","))
              X3 = parse.(Float64,split(readlines(file)[12],","))
              Y3 = parse.(Float64,split(readlines(file)[13],","))
              Z3 = parse.(Float64,split(readlines(file)[14],","))
              if numBodies > 3
                     X4 = parse.(Float64,split(readlines(file)[15],","))
                     Y4 = parse.(Float64,split(readlines(file)[16],","))
                     Z4 = parse.(Float64,split(readlines(file)[17],","))
              else
                     Z4 = [] #for deciding whether or not to plot 3D
              end
       end
       println("The timestep varied from $(minimum(lList)) to $(maximum(lList)).")
	println("The angular momentum varied by $(minimum(Llist)) to $(maximum(Llist)) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).")
       if color == "L"
              plt.plot(Tlist,Llist) 
       elseif color == "E"
              plt.plot(Tlist,Elist) #find out what the scale things are, actually change to deltaE/E0
       elseif color == "EL"
              plt.plot(Tlist,Llist,linestyle="solid",color="green") 
              plt.plot(Tlist,Elist,linestyle="solid",color="red") #find out what the scale things are, actually change to deltaE/E0
       elseif color == "time"
              plt.plot(lList,linestyle="solid",color="green")
       else
              if numBodies == 2
                     plt.plot(X1,Y1,linestyle="solid",color="red")
                     plt.plot(X2,Y2,linestyle="solid",color=color)
                     plt.axis("equal") #makes axes equal, especially helpful if orbits are highly elliptical
              else
                     if maximum(vcat(Z1,Z2,Z3,Z4))<0.0001 #checks to see if it has to graphed in 3D
                            plt.plot(X1,Y1,linestyle="solid",color="red")
                            plt.plot(X2,Y2,linestyle="solid",color=color)
                            plt.plot(X3,Y3,linestyle="solid",color="green")
                            if numBodies == 4
                                   plt.plot(X4,Y4,linestyle="solid",color="orange")
                            end
                            if equal == 0 #this will equalize the axes by default. Otherwise, it'll just plot only what it needs too
                                   plt.axis("equal")
                            end
                     else
                            plot3D(X1,Y1,Z1,linestyle="solid",color="red")
                            plot3D(X2,Y2,Z2,linestyle="solid",color=color)
                            plot3D(X3,Y3,Z3,linestyle="solid",color="green")
                            if numBodies == 4
                                   plot3D(X4,Y4,Z4,linestyle="solid",color="orange")
                            end
                            if equal == 0
                                   maxPoint = maximum(vcat(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4))
                                   minPoint = minimum(vcat(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4))
                                   scatter3D(minPoint,minPoint,minPoint,alpha=0)
                                   scatter3D(minPoint,minPoint,maxPoint,alpha=0)
                                   scatter3D(minPoint,maxPoint,minPoint,alpha=0)
                                   scatter3D(minPoint,maxPoint,maxPoint,alpha=0)
                                   scatter3D(maxPoint,minPoint,minPoint,alpha=0)
                                   scatter3D(maxPoint,minPoint,maxPoint,alpha=0)
                                   scatter3D(maxPoint,maxPoint,minPoint,alpha=0)
                                   scatter3D(maxPoint,maxPoint,maxPoint,alpha=0) #creates cube made out of extrema so it includes everything and so that the axes are all equal. I think this will work.
                            end
                     end
              end
       end
end