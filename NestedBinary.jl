import LinearAlgebra.cross
import LinearAlgebra.norm
using Pkg
Pkg.add("PyPlot")
using PyPlot

G = 2945.49 #gravitational constant


function fileInput(file) #change initial conditions to m1, m2, m3, A1, e1, A2, e2  
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
	fArray = [] #we need differential functions to calculate new positions and velocities with each timestep. These functions will be stored in this array. What functions are entered depends on how many bodies we are working with.

	mArray = parse.(Float64,split(readlines(file)[1],",")) #this inputs the masses of the bodies
	XArray = parse.(Float64,split(readlines(file)[2],",")) #this inputs the initial separations and eccentricities of the bodies

	numBodies = length(mArray) #this is different than the number of masses, so input a test particle with m=0

	if numBodies == 3 #so, here, we either have three massive bodies or two massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C) #these are the functions we need for that
		if length(mArray) == 2 #test particle with 2 numBodies
			push!(mArray,0) #including the test particle as a third body with zero mass
		end
	elseif numBodies == 4 #Main.jl cannot handle four massive bodies, so presumably, this is three massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C, f1D, f2D, f3D, f4D, f5D, f6D) #these are the functions we need
	else
		error("Check how many bodies you have. Number of position vectors are $numBodies while number of masses are $(length(mArray)).") #if we only have two bodies, then we run systemrk.jl, so not this. !!this should be implemented
	end
	t = parse.(Float64,split(readlines(file)[3],","))[1]
	hParam = parse.(Float64,split(readlines(file)[3],","))[2] #these should be the elements of the third line of the .txt file
	return fArray, XArray, mArray, t, hParam, numBodies
end

function System(file)
	#this is the main function that integrates with RK4 and returns the final positions (as well as arrays with information we can plot)
	f, x, m, t, hParam, numBodies = fileInput(file) #gets info from file

	#Masses: NOTE:FOR NOW WE ARE ASSUMING M1 IS ABOUT AS MASSIVE AS M2
	M1 = m[1]
	M2 = m[2]
	M3 = m[3]
	M = M1 + M2 + M3

	#Positions
	A1 = x[1] #keeps track of the separation between bodies 1 and  2
	A2 = x[3] #separations between the cm of bodies 1 and 2 and the mass 3
	e1 = x[2] #eccentricity of orbit 1
	e2 = x[4] #eccentricity of orbit 2
	i = x[5] #inclination of third body orbit - angle between its orbital plane and the orbital plane of the inner binary
	Θ = x[6] #offset angle - angle between the line formed by the x-axis and the x-y position of the third body. I think that this is the longitude of the ascending node-270 degrees, but I'm not sure.
	X1 = [(-(A1*M2)/(M1+M2))-cos(Θ)*cos(i)*A2*(M3/M)]
	X2 = [((M1*A1)/(M1 + M2))-cos(Θ)*cos(i)*A2*(M3/M)]#initial position of center mass
	X3 = [cos(Θ)*A2*cos(i)*(M1+M2)/M] #initial position of rightmost mass
	Y1 = [-sin(Θ)*cos(i)*A2*(M3/M)]
	Y2 = [-sin(Θ)*cos(i)*A2*(M3/M)]
	Y3 = [sin(Θ)*A2*cos(i)*(M1+M2)/M]
	Z1 = [-sin(i)*A2*(M3/M)]
	Z2 = [-sin(i)*A2*(M3/M)]
	Z3 = [A2*sin(i)*(M1+M2)/M]
	if numBodies == 4
		X4 = [x[7]]
		Y4 = [x[8]]
		Z4 = [x[9]]
	end

	#separations
	R₁ = [X1[1], Y1[1], Z1[1]] #Distance of body 1 from origin
	R₂ = [X2[1], Y2[1], Z2[1]] #similar
	R₃ = [X3[1], Y3[1], Z3[1]]
	CM₁₂ = (M1*R₁+M2*R₂)/(M1+M2)
	R₁₂ = R₁-R₂ #Distance from body 1 to 2
	R₁₃ = R₁-R₃
	R₂₃ = R₂-R₃
	
	V₁ = [velocityM1M2*sin(Θ), -sqrt(G*(M2^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cos(Θ),0.0]
	V₂ = [-velocityM1M2*sin(Θ), sqrt(G*(M1^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cos(Θ),0.0]
	V₃ = [velocityM3*sin(Θ), velocityM3*cos(Θ),0.0]
	CMv₁₂ = [-velocityM1M2*sin(Θ), -velocityM1M2*cos(Θ),0]
	V₁₂ = V₁-V₂
	V₁₃ = V₁-V₃
	V₂₃ = V₂-V₃
	p = [R₁[1], R₁[2], R₁[3], V₁[1], V₁[2], V₁[3], R₂[1], R₂[2], R₂[3], V₂[1], V₂[2], V₂[3], R₃[1], R₃[2], R₃[3], V₃[1], V₃[2], V₃[3]]

	K = .5*M1*norm(V₁)^2+.5*M2*norm(V₂)^2+.5*M3*norm(V₃)^2 #overall kinetic energy
	U = G*M1*M2/norm(R₁₂)+G*M1*M3/norm(R₁₃)+G*M2*M3/norm(R₂₃) #total gravitational potential energy
	E = K - U #total energy
	E₁ = .5*M1*norm(V₁)^2+.5*M2*norm(V₂)^2 - G*M1*M2/norm(R₁₂) #Energy of the Inner Binary
	E₂ = .5*(M1+M2)*(G*(M1+M2)/A2) + .5*M3*norm(V₃)^2 - G*(M1+M2)*M3/norm(A2) #Energy of the Outer Binary

	L = cross(R₁,M1*V₁)+cross(R₂,M2*V₂)+cross(R₃,M3*V₃) #finds angular momentum of system

	L₁ = cross(R₁,M1*V₁)+cross(R₂,M2*V₂) #Angular Momentum of the Inner binary
	L₂ = cross(CM₁₂,(M1+M2)*CMv₁₂) + cross(R₃,M3*V₃) #Momentum of the Outer Binary

	t0 = 0.0

	Llist = [L] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
	L₁list = [L₁]
	L₂list = [L₂]
	Elist = [E] #same, but for energy
	E₁list = [E₁]
	E₂list = [E₂]
	Tlist = [t0] #keeps track of time independent of timestep
	println([R₁,R₂,R₃,V₁,V₂,V₃,CMv₁₂])
	h = hParam*(minimum([norm(R₁₂)/norm(V₁₂),norm(R₁₃)/norm(V₁₃),norm(R₂₃)/norm(V₂₃)])) #=this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented=#
	
	lList = [h] #testing length of timestep

	#until the desired time has been reached, the code runs RK4

	while t0 < t
		#we will add an adaptive timestep later
		p = RK4(f, p, m, h)

		R₁ = [p[1], p[2], p[3]] #Distance of body 1 from origin
		R₂ = [p[7], p[8], p[9]] #similar
		R₃ = [p[13], p[14], p[15]]
		CM₁₂ = (M1*R₁+M2*R₂)/(M1+M2) #updated inner center of mass
		R₁₂ = R₁-R₂
		R₁₃ = R₁-R₃
		R₂₃ = R₂-R₃
		V₁ = [p[4],p[5],p[6]]
		V₂ = [p[10],p[11],p[12]]
		V₃ = [p[16],p[17],p[18]]
		CMv₁₂ = (M1*V₁+M2*V₂)/(M1+M2)
		
		V₁₂ = V₁-V₂
		V₁₃ = V₁-V₃
		V₂₃ = V₂-V₃

		K = .5*m[1]*norm(V₁)^2+.5*m[2]*norm(V₂)^2+.5*m[3]*norm(V₃)^2 #overall kinetic energy
		U = G*m[1]*m[2]/norm(R₁₂)+G*m[1]*m[3]/norm(R₁₃)+G*m[2]*m[3]/norm(R₂₃) #total gravitational potential energy
		E = K - U #total energy
		E₁ = .5*M1*norm(V₁)^2+.5*M2*norm(V₂)^2 - G*M1*M2/norm(R₁₂) #Energy of the Inner Binary
		E₂ = .5*(M1+M2)*(G*(M1+M2)/A2) + .5*M3*norm(V₃)^2 - G*(M1+M2)*M3/norm(A2)#Energy of the Outer Binary

		L = cross(R₁,m[1]*V₁)+cross(R₂,m[2]*V₂)+cross(R₃,m[3]*V₃) #finds angular momentum of system

		L₁ = cross(R₁,M1*V₁)+cross(R₂,M2*V₂) #Angular Momentum of the Inner binary
		L₂ = cross(CM₁₂,(M1+M2)*CMv₁₂) + cross(R₃,M3*V₃) #Momentum of the Outer Binary

		t0 = t0 + h #advances time

		h = hParam*(minimum([norm(R₁₂)/norm(V₁₂),norm(R₁₃)/norm(V₁₃),norm(R₂₃)/norm(V₂₃)])) #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented

		push!(Elist,E)
		push!(E₁list, E₁)
		push!(E₂list, E₂)
		push!(Llist,L)
		push!(L₁list, L₁)
		push!(L₂list, L₂)
		push!(lList,h)
		push!(Tlist,t0)
		push!(X1,p[1])
		push!(X2,p[7])
		push!(X3,p[13])
		push!(Y1,p[2]) 
		push!(Y2,p[8])
		push!(Y3,p[14])
		push!(Z1,p[3])
		push!(Z2,p[9])
		push!(Z3,p[15])
		if numBodies == 4
			push!(X4,p[19])
			push!(Y4,p[20])
			push!(Z4,p[21])
		end
	end
	if numBodies == 3
		v4x, v4y, v4z, X4, Y4, Z4 = [0,0,0,[],[],[]] #so, if we only have three bodies, we just return empty values for the "test particle" to make julia happy
	else 
		v4x, v4y, v4z = [p[22],p[23],p[24]]
	end
	return m, x, Elist, E₁list, E₂list, Llist, L₁list, L₂list, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies, hParam
end

function RK4(f,x,m,h)
#Inputs are initial position array, mass array and step size h

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

function Plot(file, color, fileSave=0, equal=0) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
	m, x, Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies, hParam = System(file)
	L0 = Llist[1]
	Llist = map(x -> norm(x-L0),Llist) #plotting the magntiude of the difference vector from L0
	E0 = Elist[1]
	Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
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
	if fileSave != 0
		bigArray = [Llist,Elist,Tlist,lList,X1,X2,Y1,Y2,Z1,Z2,X3,Y3,Z3] #stores the arrays we want to write into a 2-dimensional array
		if numBodies == 4
			push!(bigArray,X4,Y4,Z4)
		end
		tracker = [] #this will store a data point for each 1D array
		for i in 1:length(bigArray)
			if occursin("Any","$(bigArray[i])") #we loop through each stringed version of the array to see if "Any[" is at the beginning of any of them
				push!(tracker,5) #if there is, we store a 5
			else
				push!(tracker,2) #if there isn't, we store a 2
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
			end
		else
			#open("h≈(r÷v) data files/$(round(m[1];digits=5)), $(round(m[2];digits=5)), $(round(m[3];digits=5)), $(round(X1[1];digits=5)), $(round(Y1[1];digits=5)), $(round(Z1[1];digits=5)), $(round(v1x;digits=5)), $(round(v1y;digits=5)) $(round(v1z;digits=5)), $(round(X2[1];digits=5)), $(round(Y2[1];digits=5)), $(round(Z2[1];digits=5)), $(round(v2x;digits=5)), $(round(v2y;digits=5)), $(round(v2z;digits=5)), $(round(X3[1];digits=5)), $(round(Y3[1];digits=5)), $(round(Z3[1];digits=5)), $(round(v3x;digits=5)), $(round(v3y;digits=5)), $(round(v3z;digits=5)), $(round(X4[1];digits=5)), $(round(Y4[1];digits=5)), $(round(Z4[1];digits=5)), $(round(v4x;digits=1)), $(round(v4y;digits=5)), $(round(v4z;digits=5)), $hParam.txt","w") do f #UGH I had to do this because otherwise the file name would've been too long (ノಠ益ಠ)ノ彡┻━┻
			open("h≈(r÷v) data files/$fileSave","w") do f
				write(f,"4","\n")
				for i in 1:length(bigArray)
					write(f, "$(bigArray[i])"[tracker[i]:end-1],"\n") #now, we loop through, cutting off either the first 1 or 4 characters of the stringed array, depending on if it had that Any[, and also we cut off the last character, which is ].
				end
				write(f,"4")
			end
		end
	end
end
#Begin functions for the potentials

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
	g1 = 2945.49 * m[2] * (x[7] - x[1])
	g2 = 2945.49 * m[3] * (x[13] - x[1])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f4B(x::Array{Float64,1}, m::Array{Float64,1}) #v2
	g1 = 2945.49 * m[1] * (x[1] - x[7])
	g2 = 2945.49 * m[3] * (x[13] - x[7])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f4C(x::Array{Float64,1}, m::Array{Float64,1}) #v3
	g1 = 2945.49 * m[1] * (x[1] - x[13])
	g2 = 2945.49 * m[2] * (x[7] - x[13])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f4D(x::Array{Float64,1}, m::Array{Float64,1}) #v4
	g1 = 2945.49 * m[1] * (x[1] - x[19])
	g2 = 2945.49 * m[2] * (x[7] - x[19])
	g3 = 2945.49 * m[3] * (x[13] - x[19])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end
	
function f5A(x::Array{Float64,1}, m::Array{Float64,1}) #w1
	g1 = 2945.49 * m[2] * (x[8] - x[2])
	g2 = 2945.49 * m[3] * (x[14] - x[2])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f5B(x::Array{Float64,1}, m::Array{Float64,1}) #w2
	g1 = 2945.49 * m[1] * (x[2] - x[8])
	g2 = 2945.49 * m[3] * (x[14] - x[8])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f5C(x::Array{Float64,1}, m::Array{Float64,1}) #w3
	g1 = 2945.49 * m[1] * (x[2] - x[14])
	g2 = 2945.49 * m[2] * (x[8] - x[14])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f5D(x::Array{Float64,1}, m::Array{Float64,1}) #w4
	g1 = 2945.49 * m[1] * (x[2] - x[20])
	g2 = 2945.49 * m[2] * (x[8] - x[20])
	g3 = 2945.49 * m[3] * (x[14] - x[20])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end

function f6A(x::Array{Float64,1}, m::Array{Float64,1}) #u1
	g1 = 2945.49 * m[2] * (x[9] - x[3])
	g2 = 2945.49 * m[3] * (x[15] - x[3])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3)
end

function f6B(x::Array{Float64,1}, m::Array{Float64,1}) #u2
	g1 = 2945.49 * m[1] * (x[3] - x[9])
	g2 = 2945.49 * m[3] * (x[15] - x[9])
	r1 = ((x[1] - x[7])^2 + (x[2] - x[8])^2 + (x[3] - x[9])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r3^3)
end

function f6C(x::Array{Float64,1}, m::Array{Float64,1}) #u3
	g1 = 2945.49 * m[1] * (x[3] - x[15])
	g2 = 2945.49 * m[2] * (x[9] - x[15])
	r2 = ((x[1] - x[13])^2 + (x[2] - x[14])^2 +(x[3] - x[15])^2)^0.5
	r3 = ((x[7] - x[13])^2 + (x[8] - x[14])^2 +(x[9] - x[15])^2)^0.5
	return g1/(r2^3) + g2/(r3^3)
end

function f6D(x::Array{Float64,1}, m::Array{Float64,1}) #w4
	g1 = 2945.49 * m[1] * (x[3] - x[21])
	g2 = 2945.49 * m[2] * (x[9] - x[21])
	g3 = 2945.49 * m[3] * (x[15] - x[21])
	r1 = ((x[1] - x[19])^2 + (x[2] - x[20])^2 + (x[3] - x[21])^2)^0.5
	r2 = ((x[7] - x[19])^2 + (x[8] - x[20])^2 +(x[9] - x[21])^2)^0.5
	r3 = ((x[19] - x[13])^2 + (x[20] - x[14])^2 +(x[21] - x[15])^2)^0.5
	return g1/(r1^3) + g2/(r2^3) + g3/(r3^3)
end
