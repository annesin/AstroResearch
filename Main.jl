using Pkg
Pkg.add("PyPlot")
using PyPlot

G = 2945.49

function fileInput(file) #change initial conditions to m1, m2, semi-major axis, e, 
#= This function inputs a .txt file and extracts data from it to get the inputs needed for SystemRK =#
	fArray = []

	mArray = parse.(Float64,split(readlines(file)[1],",")) #change to 2
	XArray = parse.(Float64,split(readlines(file)[2],","))

	if length(XArray)%6 != 0
		error("Invalid input: not enough positions and velocities for given number of masses.")
	end

	numBodies = length(XArray)/6

	if numBodies == 3
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C)
		if length(mArray) == 2 #test particle with 2 numBodies
			push!(mArray,0) #this is really three bodies, with the third body having 0 mass for the test particle
		end
	elseif numBodies == 4
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C, f1D, f2D, f3D, f4D, f5D, f6D)
	else
		error("Check how many bodies you have. Something's not right. Number of positions are $numBodies while number of masses are $(length(mArray)).")
	end
		t = parse.(Float64,split(readlines(file)[3],","))[1]
	dev = parse.(Float64,split(readlines(file)[3],","))[2]
	return fArray, XArray, mArray, t, dev, numBodies
end

function System(file)
	f, x, m, t, dev, numBodies = fileInput(file)

	Llist = [] #this keeps track of the system's rotational momentum over time, each entry is an L at time t
	Elist = [] #same, but for energy
	Tlist = [] #keeps track of time independent of timestep
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

	h = 0.0001 #this calculates the initial timestep
	hMax = h
	hMin = h #these find the minima and maxima of the timestep interval

	lList = [h] #testing length of timestep

	counter = 0 #this counter is in order to break the h calculation in case it gets stuck

	t0 = 0

	#until the desired time has been reached, the code runs RK4

	while t0 < t
		#we will add an adaptive timestep later
		x = RK4(f, x, m, h)
		t0 = t0 + h
	    push!(lList,h)
		push!(Tlist,t0)
		push!(X1,x[1])
		push!(X2,x[7])
		push!(X3,x[13])
		push!(Y1,x[2]) #move L ane E calculation up here
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
		X4, Y4, Z4 = [[],[],[]]
	end
	return x, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies
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
	x, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies = System(file)
	#=if color == "L"
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
		   plt.plot(#=Tlist,=#Llist,linestyle="solid",color="green") 
		   E0 = Elist[1]
		   Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E
		   plt.plot(#=Tlist,=#Elist,linestyle="solid",color="red") #find out what the scale things are, actually change to deltaE/E0
		   println("The angular momentum varied by $(minimum(Llist)) to $(maximum(Llist)) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).")
	elseif color == "time"
		   plt.plot(lList,linestyle="solid",color="green")
	else=#
		if maximum(vcat(Z1,Z2,Z3,Z4))<0.0001 #checks to see if it has to graphed in 3D
			plt.plot(X1,Y1,linestyle="solid",color="red")
			plt.plot(X2,Y2,linestyle="solid",color=color)
			plt.plot(X3,Y3,linestyle="solid",color="green")
			if numBodies == 4
				plt.plot(X4,Y4,linestype="solid",color="green")
			end
			if equal != 0
				plt.axis("equal")
			end
		else
			Axes3D.plot(X1,Y1,Z1,linestyle="solid",color="red")
			Axes3D.plot(X2,Y2,Z2,linestyle="solid",color=color)
			Axes3D.plot(X3,Y3,Z3,linestyle="solid",color="green")
			if numBodies == 4
				Axes3D.plot(X4,Y4,linestype="solid",color="green")
			end
			if equal != 0
				maxPoint = maximum(vcat(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4))
				minPoint = minimum(vcat(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4))
				Axes3D.plot(minPoint,minPoint,minPoint)
				Axes3D.plot(minPoint,minPoint,maxPoint)
				Axes3D.plot(minPoint,maxPoint,minPoint)
				Axes3D.plot(minPoint,maxPoint,maxPoint)
				Axes3D.plot(maxPoint,minPoint,minPoint)
				Axes3D.plot(maxPoint,minPoint,maxPoint)
				Axes3D.plot(maxPoint,maxPoint,minPoint)
				Axes3D.plot(maxPoint,maxPoint,maxPoint) #creates cube made out of extrema so it includes everything and so that the axes are all equal. I think this will work.
			end
		end
			

		#=E0 = Elist[1]
		L0 = Llist[1]
		Elist = map(x -> (x-E0)/E0,Elist) #plotting ΔE, not E'
		Llist = map(x -> (x-L0)/L0,Llist) #plotting ΔL, not L
		println("The angular momentum varied by $(minimum(Llist)) to $(maximum(Llist)) while the energy varied by $(minimum(Elist)) to $(maximum(Elist)).")=#
	#end
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