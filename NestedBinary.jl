import LinearAlgebra.cross
import LinearAlgebra.norm
using Pkg
Pkg.add("PyPlot")
using PyPlot
Pkg.clone("https://github.com/timholy/ProgressMeter.jl.git")
Pkg.add("ProgressMeter")
using ProgressMeter
Pkg.clone("https://github.com/felipenoris/XLSX.jl.git")
Pkg.add("XLSX")
import XLSX
#!!look at how julia runs with multiple cores
G = 2945.49 #gravitational constant

"Inputs a file and retrieves the necessary information from the file. This includes the masses of the bodies and their initial conditions."
function fileInput(file) #change initial conditions to m1, m2, semi-major axis, e, 
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
	fArray = [] #we need differential functions to calculate new positions and velocities with each timestep. These functions will be stored in this array. What functions are entered depends on how many bodies we are working with.

	mArray = parse.(Float64,split(readlines(file)[1],",")) #this inputs the masses of the bodies
	XArray = parse.(Float64,split(readlines(file)[2],",")) #this inputs the initial positions and velocities of the bodies

	if length(XArray)%6 != 0 #this makes sure that each body has 6 entries: a position and velocity in the x, y, and z direction.
		error("Invalid input: not enough positions and velocities for given number of masses.")
	end

	numBodies = (length(XArray)/3)+1 #this is different than the number of masses, because the test particle counts as a massless body here

	if numBodies == 3 #so, here, we either have three massive bodies or two massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C) #these are the functions we need for that
		if length(mArray) == 2 #test particle with 2 numBodies
			push!(mArray,0) #including the test particle as a third body with zero mass
		end
	elseif numBodies == 4 #Main.jl cannot handle four massive bodies, so presumably, this is three massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C, f1D, f2D, f3D, f4D, f5D, f6D) #these are the functions we need
	else
		error("Check how many bodies you have. Number of position vectors are $numBodies while number of masses are $(length(mArray)).") #if we only have two bodies, then we run systemrk.jl, so not this. 
	end
	t, notPeriods = Pfunction(file,XArray,mArray)
	hParam = parse.(Float64,split(readlines(file)[3],",")[2]) #these should be the elements of the third line of the .txt file
	return fArray, XArray, mArray, t, hParam, numBodies, notPeriods
end

Pfunction(file,XArray,mArray) = try
	[parse.(Float64,split(readlines(file)[3],","))[1],true]
catch
	[sqrt((XArray[3]/(XArray[4]+1))^3*(4*pi^2)/(G*(mArray[1]+mArray[2]+mArray[3])))*parse(Float64,split(readlines(file)[3],",")[1][1:end-1]),parse(Float64,split(readlines(file)[3],",")[1][1:end-1])] #number of periods
end

"Inputs a file (that is a triple system) and numerically calculates the system's energy and angular momentum versus time, as well as the bodies' positions versus time."
function System(file)
	#this is the main function that integrates with RK4 and returns the final positions (as well as arrays with information we can plot)
	f, x, m, t, hParam, numBodies, periods = fileInput(file) #gets info from file

	OriginalX = x

	M1 = m[1]
	M2 = m[2]
	M3 = m[3]
	M = M1+M2+M3

	A1 = x[1]
	A2 = x[3]
	e1 = x[2]
	e2 = x[4]
	i = x[5]
	Θ = x[6]

	X1 = [(-(A1*M2)/(M1+M2))-cosd(Θ)*cosd(i)*A2*(M3/M)] #keeps track of the first body's x coordinate
	X2 = [((M1*A1)/(M1 + M2))-cosd(Θ)*cosd(i)*A2*(M3/M)] #similar for these
	X3 = [cosd(Θ)*A2*cosd(i)*(M1+M2)/M]
	Y1 = [-sind(Θ)*cosd(i)*A2*(M3/M)]
	Y2 = [-sind(Θ)*cosd(i)*A2*(M3/M)]
	Y3 = [sind(Θ)*A2*cosd(i)*(M1+M2)/M]
	Z1 = [-sind(i)*A2*(M3/M)]
	Z2 = [-sind(i)*A2*(M3/M)]
	Z3 = [A2*sind(i)*(M1+M2)/M]
	if numBodies == 4
		X4 = [x[7]]
		Y4 = [x[8]]
		Z4 = [x[9]]
	end

	R₁ = [X1[1],Y1[1],Z1[1]]
	R₂ = [X2[1],Y2[1],Z2[1]]
	R₃ = [X3[1],Y3[1],Z3[1]]
	CM₁₂ = (M1*R₁+M2*R₂)/(M1+M2)
	R₁₂ = R₁-R₂
	R₁₃ = R₁-R₃
	R₂₃ = R₂-R₃
	velocityM3 = sqrt(G*(M1+M2)^2*(1-e2)/(A2*M))
	velocityM1M2 = sqrt((G*(M3^2)*(1-e2))/(A2*M)) #velocity of inner CM
	V₁ = [velocityM1M2*sind(Θ), -sqrt(G*(M2^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cosd(Θ),0.0]
	V₂ = [-velocityM1M2*sind(Θ), sqrt(G*(M1^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cosd(Θ),0.0]
	V₃ = [velocityM3*sind(Θ), velocityM3*cosd(Θ),0.0]
	if numBodies>3
		V₄ = [x[7],x[8],x[9]]
	end
	CMv₁₂ = [-velocityM1M2*sind(Θ), -velocityM1M2*cosd(Θ),0]
	V₁₂ = V₁-V₂
	V₁₃ = V₁-V₃
	V₂₃ = V₂-V₃
	VCM = (M1*V₁+M2*V₂+M3*V₃)/M
	V₁ -= VCM
	V₂ -= VCM
	V₃ -= VCM
	VINCM = (m[2]*V₁+m[2]*V₂)/(m[1]+m[2]) #velocity of inner center of mass

	K = .5*m[1]*norm(V₁)^2+.5*m[2]*norm(V₂)^2+.5*m[3]*norm(V₃)^2 #overall kinetic energy
	U = G*m[1]*m[2]/norm(R₁₂)+G*m[1]*m[3]/norm(R₁₃)+G*m[2]*m[3]/norm(R₂₃) #total gravitational potential energy
	E = K - U #total energy
	L = cross(R₁,m[1]*V₁)+cross(R₂,m[2]*V₂)+cross(R₃,m[3]*V₃) #finds angular momentum of system

	E₁ = .5*m[1]*norm(V₁-VINCM)^2+.5*m[2]*norm(V₂-VINCM)^2 - G*m[1]*m[2]/norm(R₁₂)#Energy of inner binary
	E₂ = .5*(m[1]+m[2])*norm(VINCM)^2+.5*m[3]*norm(V₃)^2 - G*(m[1]+m[2])*m[3]/norm(R₃-CM₁₂)#Energy of outer binary
	L₁ = cross(R₁-CM₁₂,m[1]*(V₁-VINCM))+cross(R₂-CM₁₂,m[2]*(V₂-VINCM)) #momentum of inner binary
	L₂ = cross(CM₁₂,(m[1]+m[2])*VINCM)+cross(R₃,m[3]*V₃) #momentum of outer Binary

	t0 = 0.0

	Llist = [L] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
	Elist = [E] #same, but for energy
	E₁list = [E₁]
	E₂list = [E₂]
	L₁list = [L₁]
	L₂list = [L₂]
	Tlist = [t0] #keeps track of time independent of timestep

	h = hParam*(minimum([norm(R₁₂)/norm(V₁₂),norm(R₁₃)/norm(V₁₃),norm(R₂₃)/norm(V₂₃)])) #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented

	lList = [h] #testing length of timestep

	x=vcat(X1,Y1,Z1,V₁,X2,Y2,Z2,V₂,X3,Y3,Z3,V₃)
	if numBodies>3
		push!(x,X4,Y4,Z4,V₄)
	end
	#until the desired time has been reached, the code runs RK4

	prog = Progress(convert(Int,ceil(t)),0.5)
	counter = 0
	while t0 < t
		#we will add an adaptive timestep later
		x = RK4(f, x, m, h)
		R₁ = [x[1],x[2],x[3]]
		R₂ = [x[7],x[8],x[9]]
		R₃ = [x[13],x[14],x[15]]
		R₁₂ = R₁-R₂
		R₁₃ = R₁-R₃
		R₂₃ = R₂-R₃
		CM₁₂ = (M1*R₁+M2*R₂)/(M1+M2)
		V₁ = [x[4],x[5],x[6]]
		V₂ = [x[10],x[11],x[12]]
		V₃ = [x[16],x[17],x[18]]
		V₁₂ = V₁-V₂
		V₁₃ = V₁-V₃
		V₂₃ = V₂-V₃
		VINCM = (m[2]*V₁+m[2]*V₂)/(m[1]+m[2])

		K = .5*m[1]*norm(V₁)^2+.5*m[2]*norm(V₂)^2+.5*m[3]*norm(V₃)^2 #overall kinetic energy
		U = G*m[1]*m[2]/norm(R₁₂)+G*m[1]*m[3]/norm(R₁₃)+G*m[2]*m[3]/norm(R₂₃) #total gravitational potential energy
		E = K - U #total energy
		L = cross(R₁,m[1]*V₁)+cross(R₂,m[2]*V₂)+cross(R₃,m[3]*V₃) #finds angular momentum of system
		
		E₁ = .5*m[1]*norm(V₁-VINCM)^2+.5*m[2]*norm(V₂-VINCM)^2 - G*m[1]*m[2]/norm(R₁₂)#Energy of inner binary
		E₂ = .5*(m[1]+m[2])*norm(VINCM)^2+.5*m[3]*norm(V₃)^2 - G*(m[1]+m[2])*m[3]/norm(R₃-CM₁₂)#Energy of outer binary
		L₁ = cross(R₁-CM₁₂,m[1]*(V₁-VINCM))+cross(R₂-CM₁₂,m[2]*(V₂-VINCM)) #momentum of inner binary
		L₂ = cross(CM₁₂,(m[1]+m[2])*VINCM)+cross(R₃,m[3]*V₃) #momentum of outer Binary

		t0 = t0 + h #advances time
		
		h = hParam*(minimum([norm(R₁₂)/norm(V₁₂),norm(R₁₃)/norm(V₁₃),norm(R₂₃)/norm(V₂₃)])) #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented
		push!(Elist,E)
		push!(Llist,L)
		push!(E₁list,E₁)
		push!(E₂list,E₂)
		push!(L₁list,L₁)
		push!(L₂list,L₂)
	   	push!(lList,h)
		push!(Tlist,t0)
		push!(X1,x[1])
		push!(X2,x[7])
		push!(X3,x[13])
		push!(Y1,x[2]) 
		push!(Y2,x[8])
		push!(Y3,x[14])
		push!(Z1,x[3])
		push!(Z2,x[9])
		push!(Z3,x[15])
		if numBodies == 4
			push!(X4,x[19])
			push!(Y4,x[20])
			push!(Z4,x[21])
		end
		#if counter%10000==0
			update!(prog,convert(Int64,floor(t0)))
		#end
		counter += 1
	end
	if numBodies == 3
		v4x, v4y, v4z, X4, Y4, Z4 = [0,0,0,[],[],[]] #so, if we only have three bodies, we just return empty values for the "test particle" to make julia happy
	else 
		v4x, v4y, v4z = [x[22],x[23],x[24]]
	end
	return m, x, Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies, hParam, x[4], x[5], x[6], x[10], x[11], x[12], x[16], x[17], x[18], v4x, v4y, v4z, OriginalX, t0, E₁list, E₂list, L₁list, L₂list, periods #need initial velocities to write filename for saving
end

function RK4(f,x,m,h)
#Inputs are initial position arr=ay, mass array and step size h

    d=length(f)

#Setting the runge kutta parameters to zero

	k1=zeros(length(x))
	k2=zeros(length(x))
	k3=zeros(length(x))
	k4=zeros(length(x))

#calculating the runge kutta parameters for every function in f

	for i in 1:d
		k1[i]=h*f[i](x, m)
	end
	for i in 1:d
		k2[i]=h*f[i](x+k1/2, m)
	end
	for i in 1:d
		k3[i]=h*f[i](x+k2/2, m)
	end
	for i in 1:d
		k4[i]=h*f[i](x+k3, m)
	end

#returns the desired step
	return  x+(k1+2*k2+2*k3+k4)/6

end 

"""
Takes a file with the data, then plots what you choose.

Plot(file, object[, fileSave, writeData, equal])

The input file must be a .txt file in the format described in README.md.

The object is what is plotted. If "E" is typed, then the energy will be plotted versus time. If "L" is typed, then angular momentum will be plotted versus time. If "EL" is typed, then both are plotted. "time" plots the timestep of the integration versus iteration. Finally, a color accepted by matplotlib will plot the trajectories of the bodies, one of them having a path with the color specified.

writeData is optional. Unless something other than 0 is its input, it will write the data of the simulation to the NestedBinaryData spreadsheet.

fileSave is optional. However, if a string is entered, for example, "Sample.txt", then a .txt file will be created that will store the system's data. This file can then be plotted using ExternalPlotter.jl without needing to recalculate the system again.

equal is also optional. Plot() equalizes the axes of the trajectories by default. If anything besides 0 is its input, it will not do this.
"""
function Plot(file, color="none", fileSave=0, writeData=0, equal=0) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
	firstTime = time() #measures runtime of program
	m, x, Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies, hParam, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, v4x, v4y, v4z, OriginalX, t0, E₁list, E₂list, L₁list, L₂list, periods = System(file)
	NowTime = time()
	#=stability calculation=#
	if E₁list[end]>0
		stability = 0
		println("This is an unstable system.")
	elseif E₂list[end]>0
		stability = 0.5
		println("This is a partially unstable system.")
	else
		stability = 1
		println("This is a stable system.")
	end
	L0 = Llist[1]
	Llist = map(x -> (x-L0)/norm(L0),Llist) #plotting ΔL, not L
	E0 = Elist[1]
	Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
	E10 = E₁list[1]
	E₁list = map(x -> (x-E10)/E10,E₁list)
	E20 = E₂list[1]
	E₂list = map(x -> (x-E20)/E20,E₂list)
	L10 = L₁list[1]
	L₁list = map(x -> (x-L10)/norm(L10),L₁list) #plotting ΔL, not L
	L20 = L₂list[1]
	L₂list = map(x -> (x-L20)/norm(L20),L₂list) #plotting ΔL, not L
	println("The timestep varied from $(minimum(lList)) to $(maximum(lList)).")
	println("The angular momentum varied by $(minimum(map(x -> norm(x),Llist))) to $(maximum(map(x -> norm(x),Llist))) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).") #magnitude of angular momentum here for simplicity
	println("This ran in $(NowTime-firstTime) seconds.")
	if color == "L"
		   plt.plot(Tlist,map(x -> x[1],Llist),linestyle="solid",color="red")
		   plt.plot(Tlist,map(x -> x[2],Llist),linestyle="solid",color="green")
		   plt.plot(Tlist,map(x -> x[3],Llist),linestyle="solid",color="blue") #plotting change in individual components
	elseif color == "E"
		   plt.plot(Tlist,Elist) #find out what the scale things are, actually change to deltaE/E0
	elseif color == "EL"
		plt.plot(Tlist,map(x -> x[1],Llist),linestyle="solid",color="green")
		plt.plot(Tlist,map(x -> x[2],Llist),linestyle="solid",color="blue") 
		plt.plot(Tlist,map(x -> x[3],Llist),linestyle="solid",color="darkblue") #plotting change in individual components
		plt.plot(Tlist,Elist,linestyle="solid",color="red") #find out what the scale things are, actually change to deltaE/E0
	elseif color == "E1"
		plt.plot(Tlist,E₁list)
	elseif color == "L1"
		plt.plot(Tlist,map(x -> x[1],L₁list),linestyle="solid",color="red")
		plt.plot(Tlist,map(x -> x[2],L₁list),linestyle="solid",color="green") 
		plt.plot(Tlist,map(x -> x[3],L₁list),linestyle="solid",color="blue") #plotting change in individual components
	elseif color == "E2"
		plt.plot(Tlist,E₂list)
	elseif color == "L2"
		plt.plot(Tlist,map(x -> x[1],L₂list),linestyle="solid",color="red")
		plt.plot(Tlist,map(x -> x[2],L₂list),linestyle="solid",color="green") 
		plt.plot(Tlist,map(x -> x[3],L₂list),linestyle="solid",color="blue") #plotting change in individual components
	elseif color == "Es"
		plt.plot(Tlist,E₁list,linestyle="solid",color="green")
		plt.plot(Tlist,E₂list,linestyle="solid",color="red")
	elseif color == "Ls"
		plt.plot(Tlist,map(x -> x[1],L₁list),linestyle="solid",color="brown")
		plt.plot(Tlist,map(x -> x[2],L₁list),linestyle="solid",color="red") 
		plt.plot(Tlist,map(x -> x[3],L₁list),linestyle="solid",color="maroon") #plotting change in individual components
		plt.plot(Tlist,map(x -> x[1],L₂list),linestyle="solid",color="royalblue")
		plt.plot(Tlist,map(x -> x[2],L₂list),linestyle="solid",color="blue") 
		plt.plot(Tlist,map(x -> x[3],L₂list),linestyle="solid",color="darkblue") #plotting change in individual components
	elseif color == "time"
		plt.plot(lList,linestyle="solid",color="green")
	elseif color == "none" #used for automatic testing
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
			ax = plt.axes(projection="3d")
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
	if fileSave != 0
		bigArray = [Llist,Elist,Tlist,lList,X1,X2,Y1,Y2,Z1,Z2,X3,Y3,Z3] #stores the arrays we want to write into a 2-dimensional array
		if numBodies == 4
			push!(bigArray,X4,Y4,Z4)
		end
		push!(bigArray,L₁list,L₂list,[NowTime-firstTime],E₁list,E₂list)
		tracker = [] #this will store a data point for each 1D array
		for i in 1:length(bigArray)
			if occursin("Array{Float64,1}","$(bigArray[i])")
				push!(tracker,18)
			elseif occursin("Any","$(bigArray[i])") #we loop through each stringed version of the array to see if "Any[" is at the beginning of any of them
				push!(tracker,5) #if there is, we store a 5
			else
				push!(tracker,2) #if there isn't, we store a 2
			end
			if floor(NowTime-firstTime) < 1
				bigArray[i]=bigArray[i][1:convert(Int,ceil(NowTime-firstTime)):end]#we only save every nth term, where n is the amount of time it took to run the code
			else
				bigArray[i]=bigArray[i][1:convert(Int,floor(NowTime-firstTime)):end]
			end
		end
		if numBodies == 3
			#=I'm saving these in case we need them later
			open("h≈(r÷v) data files/$(m[1]), $(m[2]), $(m[3]), $(X1[1]), $(Y1[1]), $(Z1[1]), $v1x, $v1y, $v1z, $(X2[1]), $(Y2[1]), $(Z2[1]), $v2x, $v2y, $v2z, $(X3[1]), $(Y3[1]), $(Z3[1]), $v3x, $v3y, $v3z, $hParam.txt","w") do f =#
			open("h≈(r÷v) data files/$fileSave","w") do f
				write(f,"3","\n")
				for i in 1:length(bigArray)
					write(f, "$(bigArray[i])"[tracker[i]:end-1],"\n") #now, we loop through, cutting off either the first 1 or 4 characters of the stringed array, depending on if it had that Any[, and also we cut off the last character, which is ].
				end
				write(f,"$m"[2:end-1],"\n","$OriginalX"[2:end-1],"\n")
				write(f,"N")
			end
		else
			#open("h≈(r÷v) data files/$(round(m[1];digits=5)), $(round(m[2];digits=5)), $(round(m[3];digits=5)), $(round(X1[1];digits=5)), $(round(Y1[1];digits=5)), $(round(Z1[1];digits=5)), $(round(v1x;digits=5)), $(round(v1y;digits=5)) $(round(v1z;digits=5)), $(round(X2[1];digits=5)), $(round(Y2[1];digits=5)), $(round(Z2[1];digits=5)), $(round(v2x;digits=5)), $(round(v2y;digits=5)), $(round(v2z;digits=5)), $(round(X3[1];digits=5)), $(round(Y3[1];digits=5)), $(round(Z3[1];digits=5)), $(round(v3x;digits=5)), $(round(v3y;digits=5)), $(round(v3z;digits=5)), $(round(X4[1];digits=5)), $(round(Y4[1];digits=5)), $(round(Z4[1];digits=5)), $(round(v4x;digits=1)), $(round(v4y;digits=5)), $(round(v4z;digits=5)), $hParam.txt","w") do f #UGH I had to do this because otherwise the file name would've been too long (ノಠ益ಠ)ノ彡┻━┻
			open("h≈(r÷v) data files/$fileSave","w") do f
				write(f,"4","\n")
				for i in 1:length(bigArray)
					write(f, "$(bigArray[i])"[tracker[i]:end-1],"\n") #now, we loop through, cutting off either the first 1 or 4 characters of the stringed array, depending on if it had that Any[, and also we cut off the last character, which is ].
				end
				write(f,"$m","\n","$OriginalX","\n")
				write(f,"N")
			end
		end
	end
	record = true
	rowNumber = 0
	if writeData == 0
		XLSX.openxlsx("NestedBinaryData.xlsm",mode="rw") do xf
			sheet = xf[1]
			i = 1
			while typeof(sheet["A$i"]) != Missing #gets next blank row
				if [sheet["A$i"],sheet["B$i"],sheet["C$i"],sheet["D$i"],sheet["E$i"],sheet["F$i"],sheet["G$i"],sheet["H$i"],sheet["I$i"],sheet["J$i"],sheet["K$i"]]==[m[1],m[2],m[3],OriginalX[1],OriginalX[2],OriginalX[3],OriginalX[4],OriginalX[5],OriginalX[6],t0,hParam]
					println("Not saving to spreadsheet: This data already has an entry at line $i.")
					record = false
				end
				i += 1
			end
			rowNumber = i
			if record
				sheet["A$i"] = m[1]
				sheet["B$i"] = m[2]
				sheet["C$i"] = m[3]
				sheet["D$i"] = OriginalX[1]
				sheet["E$i"] = OriginalX[2]
				sheet["F$i"] = OriginalX[3]
				sheet["G$i"] = OriginalX[4]
				sheet["H$i"] = OriginalX[5]
				sheet["I$i"] = OriginalX[6]
				if periods == true
					sheet["J$i"] = t0
				else
					sheet["J$i"] = "$periods periods"
				end
				sheet["K$i"] = hParam
				sheet["L$i"] = "[$(minimum(Elist)),$(maximum(Elist))]"
				Lx = map(x -> x[1],Llist)
				Ly = map(x -> x[2],Llist)
				Lz = map(x -> x[3],Llist)
				sheet["M$i"] = "[[$(minimum(Lx)),$(minimum(Ly)),$(minimum(Lz))],[$(maximum(Lx)),$(maximum(Ly)),$(maximum(Lz))]]"
				sheet["N$i"] = "[$(minimum(E₁list)),$(maximum(E₁list))]"
				sheet["O$i"] = "[$(minimum(E₂list)),$(maximum(E₂list))]"
				L1x = map(x -> x[1],L₁list)
				L1y = map(x -> x[2],L₁list)
				L1z = map(x -> x[3],L₁list)
				L2x = map(x -> x[1],L₂list)
				L2y = map(x -> x[2],L₂list)
				L2z = map(x -> x[3],L₂list)
				sheet["P$i"] = "[[$(minimum(L1x)),$(minimum(L1y)),$(minimum(L1z))],[$(maximum(L1x)),$(maximum(L1y)),$(maximum(L1z))]]"
				sheet["Q$i"] = "[[$(minimum(L2x)),$(minimum(L2y)),$(minimum(L2z))],[$(maximum(L2x)),$(maximum(L2y)),$(maximum(L2z))]]"
				sheet["R$i"] = NowTime-firstTime
				sheet["S$i"] = "[$(minimum(lList)),$(maximum(lList))]"
				sheet["T$i"] = stability
				sheet["U$i"] = rowNumber
				if minimum(Elist)<-10.0^-3 || maximum(Elist)>10.0^-3
					sheet["V$i"] = 0
				else
					sheet["V$i"] = 1
				end
			end
		end
	end
	return record, rowNumber #used for autotester
end

function AutomaticTester(fileSave,iI=1,iJ=-2,iK=1,iL=2)
	file = "AutomaticTester.txt"
	Masses = ["8,8,1","1.5,1.5,1","8,5,1","8,1.5,1"]
	Eccentricities = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
	for i in iI:4, j in 0:-0.1:iJ, k in iK:11, l in 10:-1:iL 
		hParam = 0.01 #check and change this
		x = open("$file","w")
		write(x,"$(Masses[i])\n")
		write(x,"$(10^j),$(Eccentricities[k]),$l,0,0,0\n") 
		write(x,"100P,$hParam")
		close(x)
		record, rowNumber = Plot(file, "none", "$fileSave"*"_0.txt")
		rm(file)
		mv("h≈(r÷v) data files/$fileSave"*"_0.txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt")
		if record == false
			rm("h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt") #deletes text file if data wasn't recorded in spreadsheet
		end
	end
end

function f0(x::Array{Float64,1}) #DE for Masses
	return 0
end

function f1A(x::Array{Float64,1}, m::Array{Float64,1}) #x1
	return x[4]
end

function f1B(x::Array{Float64,1}, m::Array{Float64,1}) #x2
	return x[10]
end

function f1C(x::Array{Float64,1}, m::Array{Float64,1}) #x3
	return x[16]
end

function f1D(x::Array{Float64,1}, m::Array{Float64,1}) #x4(TestParticle)
	return x[22]
end

function f2A(x::Array{Float64,1}, m::Array{Float64,1}) #y1
	return x[5]
end

function f2B(x::Array{Float64,1}, m::Array{Float64,1}) #y2
	return x[11]
end

function f2C(x::Array{Float64,1}, m::Array{Float64,1}) #y3
	return x[17]
end

function f2D(x::Array{Float64,1}, m::Array{Float64,1}) #y4(TestParticle)
	return x[23]
end	

function f3A(x::Array{Float64,1}, m::Array{Float64,1}) #z1
	return x[6]
end

function f3B(x::Array{Float64,1}, m::Array{Float64,1}) #z2
	return x[12]
end

function f3C(x::Array{Float64,1}, m::Array{Float64,1}) #z3
	return x[18]
end

function f3D(x::Array{Float64,1}, m::Array{Float64,1}) #z4(TestParticle)
	return x[24]
end

function f4A(x::Array{Float64,1}, m::Array{Float64,1}) #v1
	g1 = G * m[2] * (x[7] - x[1])
	g2 = G * m[3] * (x[13] - x[1])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f4B(x::Array{Float64,1}, m::Array{Float64,1}) #v2
	g1 = G * m[1] * (x[1] - x[7])
	g2 = G * m[3] * (x[13] - x[7])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f4C(x::Array{Float64,1}, m::Array{Float64,1}) #v3
	g1 = G * m[1] * (x[1] - x[13])
	g2 = G * m[2] * (x[7] - x[13])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f4D(x::Array{Float64,1}, m::Array{Float64,1}) #v4
	g1 = G * m[1] * (x[1] - x[19])
	g2 = G * m[2] * (x[7] - x[19])
	g3 = G * m[3] * (x[13] - x[19])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end
	
function f5A(x::Array{Float64,1}, m::Array{Float64,1}) #w1
	g1 = G * m[2] * (x[8] - x[2])
	g2 = G * m[3] * (x[14] - x[2])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f5B(x::Array{Float64,1}, m::Array{Float64,1}) #w2
	g1 = G * m[1] * (x[2] - x[8])
	g2 = G * m[3] * (x[14] - x[8])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f5C(x::Array{Float64,1}, m::Array{Float64,1}) #w3
	g1 = G * m[1] * (x[2] - x[14])
	g2 = G * m[2] * (x[8] - x[14])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f5D(x::Array{Float64,1}, m::Array{Float64,1}) #w4
	g1 = G * m[1] * (x[2] - x[20])
	g2 = G * m[2] * (x[8] - x[20])
	g3 = G * m[3] * (x[14] - x[20])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end

function f6A(x::Array{Float64,1}, m::Array{Float64,1}) #u1
	g1 = G * m[2] * (x[9] - x[3])
	g2 = G * m[3] * (x[15] - x[3])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f6B(x::Array{Float64,1}, m::Array{Float64,1}) #u2
	g1 = G * m[1] * (x[3] - x[9])
	g2 = G * m[3] * (x[15] - x[9])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f6C(x::Array{Float64,1}, m::Array{Float64,1}) #u3
	g1 = G * m[1] * (x[3] - x[15])
	g2 = G * m[2] * (x[9] - x[15])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f6D(x::Array{Float64,1}, m::Array{Float64,1}) #w4
	g1 = G * m[1] * (x[3] - x[21])
	g2 = G * m[2] * (x[9] - x[21])
	g3 = G * m[3] * (x[15] - x[21])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end