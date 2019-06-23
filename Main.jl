import LinearAlgebra.cross
import LinearAlgebra.norm
using Pkg
Pkg.add("PyPlot")
using PyPlot

G = 2945.49 #gravitational constant

function fileInput(file) #change initial conditions to m1, m2, semi-major axis, e, 
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
	fArray = [] #we need differential functions to calculate new positions and velocities with each timestep. These functions will be stored in this array. What functions are entered depends on how many bodies we are working with.

	mArray = parse.(Float64,split(readlines(file)[1],",")) #this inputs the masses of the bodies
	XArray = parse.(Float64,split(readlines(file)[2],",")) #this inputs the initial positions and velocities of the bodies

	if length(XArray)%6 != 0 #this makes sure that each body has 6 entries: a position and velocity in the x, y, and z direction.
		error("Invalid input: not enough positions and velocities for given number of masses.")
	end

	numBodies = length(XArray)/6 #this is different than the number of masses, because the test particle counts as a massless body here

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

	X1 = [x[1]] #keeps track of the first body's x coordinate
	X2 = [x[7]] #similar for these
	X3 = [x[13]]
	Y1 = [x[2]]
	Y2 = [x[8]]
	Y3 = [x[14]]
	Z1 = [x[3]]
	Z2 = [x[9]]
	Z3 = [x[15]]
	if numBodies == 4
		X4 = [x[19]]
		Y4 = [x[20]]
		Z4 = [x[21]]
	end

	R₁ = [x[1],x[2],x[3]]
	R₂ = [x[7],x[8],x[9]]
	R₃ = [x[13],x[14],x[15]]
	R₁₂ = R₁-R₂
	R₁₃ = R₁-R₃
	R₂₃ = R₂-R₃
	V₁ = [x[4],x[5],x[6]]
	V₂ = [x[10],x[11],x[12]]
	V₃ = [x[16],x[17],x[18]]
	V₁₂ = V₁-V₂
	V₁₃ = V₁-V₃
	V₂₃ = V₂-V₃

	K = .5*m[1]*norm(V₁)^2+.5*m[2]*norm(V₂)^2+.5*m[3]*norm(V₃)^2 #overall kinetic energy
	U = G*m[1]*m[2]/norm(R₁₂)+G*m[1]*m[3]/norm(R₁₃)+G*m[2]*m[3]/norm(R₂₃) #total gravitational potential energy
	E = K - U #total energy

	L = cross(R₁,m[1]*V₁)+cross(R₂,m[2]*V₂)+cross(R₃,m[3]*V₃) #finds angular momentum of system

	t0 = 0

	Llist = [L] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
	Elist = [E] #same, but for energy
	Tlist = [t0] #keeps track of time independent of timestep

	h = hParam*(maximum(norm(R₁₂)/norm(V₁₂),norm(R₁₃)/norm(V₁₃),norm(R₂₃)/norm(V₂₃))) #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented
	hMax = h
	hMin = h #these find the minima and maxima of the timestep interval

	lList = [h] #testing length of timestep

	#until the desired time has been reached, the code runs RK4

	while t0 < t
		#we will add an adaptive timestep later
		x = RK4(f, x, m, h)

		R₁ = [x[1],x[2],x[3]]
		R₂ = [x[7],x[8],x[9]]
		R₃ = [x[13],x[14],x[15]]
		R₁₂ = R₁-R₂
		R₁₃ = R₁-R₃
		R₂₃ = R₂-R₃
		V₁ = [x[4],x[5],x[6]]
		V₂ = [x[10],x[11],x[12]]
		V₃ = [x[16],x[17],x[18]]
		V₁₂ = V₁-V₂
		V₁₃ = V₁-V₃
		V₂₃ = V₂-V₃

		K = .5*m[1]*norm(V₁)^2+.5*m[2]*norm(V₂)^2+.5*m[3]*norm(V₃)^2 #overall kinetic energy
		U = G*m[1]*m[2]/norm(R₁₂)+G*m[1]*m[3]/norm(R₁₃)+G*m[2]*m[3]/norm(R₂₃) #total gravitational potential energy
		E = K - U #total energy

		L = cross(R₁,m[1]*V₁)+cross(R₂,m[2]*V₂)+cross(R₃,m[3]*V₃) #finds angular momentum of system

		t0 = t0 + h #advances time
		push!(Elist,E)
		push!(Llist,L)
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
	end
	if numBodies == 3
		X4, Y4, Z4 = [[],[],[]] #so, if we only have three bodies, we just return empty arrays for these to make julia happy
	end
	return x, Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies
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

function Plot(file, color, equal=0) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
	x, Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies = System(file)
	#= these will be implemented once we figure out how to calculate E and L=#
	L0 = Llist[1]
	Llist = map(x -> norm(x-L0),Llist) #plotting the magntiude of the difference vector from L0
	E0 = Elist[1]
	Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
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
				plt.plot(X4,Y4,linestype="solid",color="green")
			end
			if equal == 0 #this will equalize the axes by default. Otherwise, it'll just plot only what it needs too
				plt.axis("equal")
			end
		else
			plot3D(X1,Y1,Z1,linestyle="solid",color="red")
			plot3D(X2,Y2,Z2,linestyle="solid",color=color)
			plot3D(X3,Y3,Z3,linestyle="solid",color="green")
			if numBodies == 4
				plot3D(X4,Y4,linestype="solid",color="green")
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