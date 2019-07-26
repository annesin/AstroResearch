using Pkg
Pkg.add("PyPlot")
using PyPlot
Pkg.clone("https://github.com/timholy/ProgressMeter.jl.git")
Pkg.add("ProgressMeter")
using ProgressMeter
Pkg.clone("https://github.com/felipenoris/XLSX.jl.git")
Pkg.add("XLSX")
import XLSX

const G = 2945.49 #gravitational constant, making this const greatly reduces memory usage

"Inputs a file and retrieves the necessary information from the file. This includes the masses of the bodies and their initial conditions."
function fileInput(file) #change initial conditions to m1, m2, semi-major axis, e, 
#= This function inputs a .txt file and extracts data from it to get the inputs needed for NestedBinary =#
	fArray = Function[] #we need differential functions to calculate new positions and velocities with each timestep. These functions will be stored in this array. What functions are entered depends on how many bodies we are working with.
	#println(fArray)
	#!!make this array a defined type!! Along with the other arrays
	mArray = parse.(Float64,split(readlines(file)[1],",")) #this inputs the masses of the bodies
	XArray = parse.(Float64,split(readlines(file)[2],",")) #this inputs the initial conditions of the bodies

	if length(XArray) != 6 #this makes sure that each body has 6 entries: a position and velocity in the x, y, and z direction.
		error("Invalid input: wrong number of initial conditions; go back and check the input.")
	end

	numBodies = (length(XArray)/6)+2 #this is different than the number of masses, because the test particle counts as a massless body here

	if numBodies == 3 #so, here, we either have three massive bodies or two massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C) #these are the functions we need for that
		if length(mArray) == 2 #test particle with 2 numBodies
			push!(mArray,0) #including the test particle as a third body with zero mass
		end
	elseif numBodies == 4 #NestedBinary.jl cannot handle four massive bodies, so presumably, this is three massive bodies with a test particle
		push!(fArray, f1A, f2A, f3A, f4A, f5A, f6A, f1B, f2B, f3B, f4B, f5B, f6B, f1C, f2C, f3C, f4C, f5C, f6C, f1D, f2D, f3D, f4D, f5D, f6D) #these are the functions we need
	else
		error("Check how many bodies you have. Number of position vectors are $numBodies while number of masses are $(length(mArray)).") #if we only have two bodies, then we run systemrk.jl, so not this. 
	end
	Parse = split(readlines(file)[3],",")
	if Parse[1][end] == "P"[1]
		t = sqrt((XArray[3]/(XArray[4]+1))^3*(4*pi^2)/(G*(mArray[1]+mArray[2]+mArray[3])))*parse(Float64,Parse[1][1:end-1])
		notPeriods = Parse[1]
	else
		t = parse(Float64,Parse[1])
		notPeriods = true
	end
	hParam = parse.(Float64,split(readlines(file)[3],",")[2]) #these should be the elements of the third line of the .txt file
	return fArray, XArray, mArray, t, hParam, numBodies, notPeriods
end

"Inputs a file (that is a triple system) and numerically calculates the system's energy and angular momentum versus time, as well as the bodies' positions versus time."
function System(file, fileSave, MemorySave=true)
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

	X1 = Float64[(-(A1*M2)/(M1+M2))-cosd(Θ)*cosd(i)*A2*(M3/M)] #keeps track of the first body's x coordinate
	X2 = Float64[((M1*A1)/(M1 + M2))-cosd(Θ)*cosd(i)*A2*(M3/M)] #similar for these
	X3 = Float64[cosd(Θ)*A2*cosd(i)*(M1+M2)/M] #!!these don't need to arrays anymore, just updated
	Y1 = Float64[-sind(Θ)*cosd(i)*A2*(M3/M)]
	Y2 = Float64[-sind(Θ)*cosd(i)*A2*(M3/M)]
	Y3 = Float64[sind(Θ)*A2*cosd(i)*(M1+M2)/M]
	Z1 = Float64[-sind(i)*A2*(M3/M)]
	Z2 = Float64[-sind(i)*A2*(M3/M)]
	Z3 = Float64[A2*sind(i)*(M1+M2)/M]
	if numBodies == 4
		X4 = Float64[x[7]]
		Y4 = Float64[x[8]]
		Z4 = Float64[x[9]]
	end

	R₁X = X1[1]
	R₁Y = Y1[1]
	R₁Z = Z1[1]
	R₂X = X2[1]
	R₂Y = Y2[1]
	R₂Z = Z2[1]
	R₃X = X3[1] 
	R₃Y = Y3[1]
	R₃Z = Z3[1]
	CM₁₂X = (M1*R₁X+M2*R₂X)/(M1+M2)
	CM₁₂Y = (M1*R₁Y+M2*R₂Y)/(M1+M2)
	CM₁₂Z = (M1*R₁Z+M2*R₂Z)/(M1+M2)
	R₁₂X = R₁X-R₂X
	R₁₂Y = R₁Y-R₂Y
	R₁₂Z = R₁Z-R₂Z
	R₁₃X = R₁X-R₃X
	R₁₃Y = R₁Y-R₃Y
	R₁₃Z = R₁Z-R₃Z
	R₂₃X = R₂X-R₃X
	R₂₃Y = R₂Y-R₃Y
	R₂₃Z = R₂Z-R₃Z
	velocityM3 = sqrt(G*(M1+M2)^2*(1-e2)/(A2*M)) 
	velocityM1M2 = sqrt((G*(M3^2)*(1-e2))/(A2*M)) #velocity of inner CM
	V₁X = -velocityM1M2*sind(Θ)
	V₁Y = -sqrt(G*(M2^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cosd(Θ)
	V₁Z = 0.0
	V₂X = -velocityM1M2*sind(Θ)
	V₂Y = sqrt(G*(M1^2)*(1-e1)/(A1*(M2+M1)))-velocityM1M2*cosd(Θ)
	V₂Z = 0.0
	V₃X = velocityM3*sind(Θ)
	V₃Y = velocityM3*cosd(Θ)
	V₃Z = 0.0
	if numBodies>3
		V₄X = x[7]
		V₄Y = x[8]
		V₄Z = x[9]
	end
	V₁₂X = V₁X-V₂X
	V₁₂Y = V₁Y-V₂Y
	V₁₂Z = V₁Z-V₂Z
	V₁₃X = V₁X-V₃X
	V₁₃Y = V₁Y-V₃Y
	V₁₃Z = V₁Z-V₃Z
	V₂₃X = V₂X-V₃X
	V₂₃Y = V₂Y-V₃Y
	V₂₃Z = V₂Z-V₃Z
	VINCMX = (m[2]*V₁X+m[2]*V₂X)/(m[1]+m[2]) #velocity of inner center of mass
	VINCMY = (m[2]*V₁Y+m[2]*V₂Y)/(m[1]+m[2]) #velocity of inner center of mass
	VINCMZ = (m[2]*V₁Z+m[2]*V₂Z)/(m[1]+m[2]) #velocity of inner center of mass

	K = .5*m[1]*sqrt(V₁X^2+V₁Y^2+V₁Z^2)^2+.5*m[2]*sqrt(V₂X^2+V₂Y^2+V₂Z^2)^2+.5*m[3]*sqrt(V₃X^2+V₃Y^2+V₃Z^2)^2 #overall kinetic energy
	U = -(G*m[1]*m[2]/sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)+G*m[1]*m[3]/sqrt(R₁₃X^2+R₁₃Y^2+R₁₃Z^2)+G*m[2]*m[3]/sqrt(R₂₃X^2+R₂₃Y^2+R₂₃Z^2)) #total gravitational potential energy
	E = K + U #total energy 
	LX = m[1]*(R₁Y*V₁Z-R₁Z*V₁Y)+m[2]*(R₂Y*V₂Z-R₂Z*V₂Y)+m[3]*(R₃Y*V₃Z-R₃Z*V₃Y)
	LY = m[1]*(R₁Z*V₁X-R₁X*V₁Z)+m[2]*(R₂Z*V₂X-R₂X*V₂Z)+m[3]*(R₃Z*V₃X-R₃X*V₃Z)
	LZ = m[1]*(R₁X*V₁Y-R₁Y*V₁X)+m[2]*(R₂X*V₂Y-R₂Y*V₂X)+m[3]*(R₃X*V₃Y-R₃Y*V₃X)

	E₁ = .5*m[1]*sqrt((V₁X-VINCMX)^2+(V₁Y-VINCMY)^2+(V₁Z-VINCMZ)^2)^2+.5*m[2]*sqrt((V₂X-VINCMX)^2+(V₂Y-VINCMY)^2+(V₂Z-VINCMZ)^2)^2 - G*m[1]*m[2]/sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)#Energy of inner binary
	E₂ = .5*(m[1]+m[2])*sqrt(VINCMX^2+VINCMY^2+VINCMZ^2)^2+.5*m[3]*sqrt(V₃X^2+V₃Y^2+V₃Z^2)^2 - G*(m[1]+m[2])*m[3]/sqrt((R₃X-CM₁₂X)^2+(R₃X-CM₁₂X)^2+(R₃X-CM₁₂X)^2)#Energy of outer binary
	L₁X = m[1]*((R₁Y-CM₁₂Y)*(V₁Z-VINCMZ)-(R₁Z-CM₁₂Z)*(V₁Y-VINCMY))+m[2]*((R₂Y-CM₁₂Y)*(V₂Z-VINCMZ)-(R₂Z-CM₁₂Z)*(V₂Y-VINCMY))
	L₁Y = m[1]*((R₁Z-CM₁₂Z)*(V₁X-VINCMX)-(R₁X-CM₁₂X)*(V₁Z-VINCMZ))+m[2]*((R₂Z-CM₁₂Z)*(V₂X-VINCMX)-(R₂X-CM₁₂X)*(V₂Z-VINCMZ))
	L₁Z = m[1]*((R₁X-CM₁₂X)*(V₁Y-VINCMY)-(R₁Y-CM₁₂Y)*(V₁X-VINCMX))+m[2]*((R₂X-CM₁₂X)*(V₂Y-VINCMY)-(R₂Y-CM₁₂Y)*(V₂X-VINCMX))
	L₂X = (m[1]+m[2])*(CM₁₂Y*VINCMZ-CM₁₂Z*VINCMY)+m[3]*(R₃Y*V₃Z-R₃Z*V₃Y)
	L₂Y = (m[1]+m[2])*(CM₁₂Z*VINCMX-CM₁₂X*VINCMZ)+m[3]*(R₃Z*V₃X-R₃X*V₃Z)
	L₂Z = (m[1]+m[2])*(CM₁₂X*VINCMY-CM₁₂Y*VINCMX)+m[3]*(R₃X*V₃Y-R₃Y*V₃X)

	t0 = 0.0

	Lmax = sqrt(LX^2+LY^2+LZ^2)
	Lmin = sqrt(LX^2+LY^2+LZ^2)
	L0 = sqrt(LX^2+LY^2+LZ^2)
	Emax = E
	Emin = E
	E0 = E
	E₁max = E₁
	E₁min = E₁
	E₁0 = E₁
	E₂max = E₂
	E₂min = E₂
	E₂0 = E₂
	L₁max = sqrt(L₁X^2+L₁Y^2+L₁Z^2)
	L₁min = sqrt(L₁X^2+L₁Y^2+L₁Z^2)
	L₁0 = sqrt(L₁X^2+L₁Y^2+L₁Z^2)
	L₂max = sqrt(L₂X^2+L₂Y^2+L₂Z^2)
	L₂min = sqrt(L₂X^2+L₂Y^2+L₂Z^2)
	L₂0 = sqrt(L₂X^2+L₂Y^2+L₂Z^2)

	h1 = hParam*sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)/sqrt(V₁₂X^2+V₁₂Y^2+V₁₂Z^2)
	h2 = hParam*sqrt(R₁₃X^2+R₁₃Y^2+R₁₃Z^2)/sqrt(V₁₃X^2+V₁₃Y^2+V₁₃Z^2)
	h3 = hParam*sqrt(R₂₃X^2+R₂₃Y^2+R₂₃Z^2)/sqrt(V₂₃X^2+V₂₃Y^2+V₂₃Z^2)
	if h1 < h2
		if h1 < h3
			h = h1
		else
			h = h3
		end
	else
		if h2 < h3
			h = h2
		else
			h = h3
		end
	end

	lmax = h
	lmin = h

	#calculates inner binary period
	Iperiod = sqrt(4*pi^2*A1/(1+e1)/G*(M1+M2))

	x=[X1[1],Y1[1],Z1[1],V₁X,V₁Y,V₁Z,X2[1],Y2[1],Z2[1],V₂X,V₂Y,V₂Y,X3[1],Y3[1],Z3[1],V₃X,V₃Y,V₃Z]
	if numBodies>3
		push!(x,X4[1],Y4[1],Z4[1],V₄X,V₄Y,V₄Z)
	end
	x1=x
	#until the desired time has been reached, the code runs RK4
	counter=0
	#sees how many timesteps are in one inner binary orbital period
	if MemorySave
		while t0 < Iperiod
			x = RK4(f, x, m, h)
			R₁₂X = x[1]-x[7]
			R₁₂Y = x[2]-x[8]
			R₁₂Z = x[3]-x[9]
			R₁₃X = x[1]-x[13]
			R₁₃Y = x[2]-x[14]
			R₁₃Z = x[3]-x[15]
			R₂₃X = x[7]-x[13]
			R₂₃Y = x[8]-x[14]
			R₂₃Z = x[9]-x[15]
			V₁₂X = x[4]-x[10]
			V₁₂Y = x[5]-x[11]
			V₁₂Z = x[6]-x[12]
			V₁₃X = x[4]-x[16]
			V₁₃Y = x[5]-x[17] 
			V₁₃Z = x[6]-x[18]
			V₂₃X = x[10]-x[16]
			V₂₃Y = x[11]-x[17] 
			V₂₃Z = x[12]-x[18]

			t0 = t0 + h #advances time
	
			h1 = hParam*sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)/sqrt(V₁₂X^2+V₁₂Y^2+V₁₂Z^2)
			h2 = hParam*sqrt(R₁₃X^2+R₁₃Y^2+R₁₃Z^2)/sqrt(V₁₃X^2+V₁₃Y^2+V₁₃Z^2)
			h3 = hParam*sqrt(R₂₃X^2+R₂₃Y^2+R₂₃Z^2)/sqrt(V₂₃X^2+V₂₃Y^2+V₂₃Z^2)
			if h1 < h2
				if h1 < h3
					h = h1
				else
					h = h3
				end
			else
				if h2 < h3
					h = h2
				else
					h = h3
				end
			end #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented
			counter += 1
		end

		stepSave=convert(Int64,round(counter/100))
		println("The period of the inner binary is $(Iperiod), which takes $counter timesteps to simulate, so we'll save every $(stepSave)th timestep.")
		t0=0
		x=x1
		counter = 0
	else
		stepSave = 1
	end 
	prog = Progress(convert(Int,ceil(t)),0.5)
	open("h≈(r÷v) data files/$fileSave"*".txt","w") do datafile
		write(datafile,"$numBodies","\n")
		firstTime = time()
		while t0 < t
			#we will add an adaptive timestep later
			x = RK4(f, x, m, h)
			R₁₂X = x[1]-x[7]
			R₁₂Y = x[2]-x[8]
			R₁₂Z = x[3]-x[9]
			R₁₃X = x[1]-x[13]
			R₁₃Y = x[2]-x[14]
			R₁₃Z = x[3]-x[15]
			R₂₃X = x[7]-x[13]
			R₂₃Y = x[8]-x[14]
			R₂₃Z = x[9]-x[15]
			V₁₂X = x[4]-x[10]
			V₁₂Y = x[5]-x[11]
			V₁₂Z = x[6]-x[12]
			V₁₃X = x[4]-x[16]
			V₁₃Y = x[5]-x[17] 
			V₁₃Z = x[6]-x[18]
			V₂₃X = x[10]-x[16]
			V₂₃Y = x[11]-x[17] 
			V₂₃Z = x[12]-x[18]

			K = .5*m[1]*sqrt(V₁X^2+V₁Y^2+V₁Z^2)^2+.5*m[2]*sqrt(V₂X^2+V₂Y^2+V₂Z^2)^2+.5*m[3]*sqrt(V₃X^2+V₃Y^2+V₃Z^2)^2 #overall kinetic energy
			U = -(G*m[1]*m[2]/sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)+G*m[1]*m[3]/sqrt(R₁₃X^2+R₁₃Y^2+R₁₃Z^2)+G*m[2]*m[3]/sqrt(R₂₃X^2+R₂₃Y^2+R₂₃Z^2)) #total gravitational potential energy
			E = K + U #total energy 
			LX = m[1]*(R₁Y*V₁Z-R₁Z*V₁Y)+m[2]*(R₂Y*V₂Z-R₂Z*V₂Y)+m[3]*(R₃Y*V₃Z-R₃Z*V₃Y)
			LY = m[1]*(R₁Z*V₁X-R₁X*V₁Z)+m[2]*(R₂Z*V₂X-R₂X*V₂Z)+m[3]*(R₃Z*V₃X-R₃X*V₃Z)
			LZ = m[1]*(R₁X*V₁Y-R₁Y*V₁X)+m[2]*(R₂X*V₂Y-R₂Y*V₂X)+m[3]*(R₃X*V₃Y-R₃Y*V₃X)
			t0 = t0 + h #advances time
			if counter%stepSave == 0 || sqrt(LX^2+LY^2+LZ^2) > Lmax || sqrt(LX^2+LY^2+LZ^2) < Lmin || E > Emax || E < Emin || h > lmax || h < lmin || E₁ > E₁max || E₁ < E₁min || E₂ > E₂max || E₂ < E₂min || sqrt(L₁X^2+L₁Y^2+L₁Z^2) > L₁max || sqrt(L₁X^2+L₁Y^2+L₁Z^2) < L₁min || sqrt(L₂X^2+L₂Y^2+L₂Z^2) > L₂max || sqrt(L₂X^2+L₂Y^2+L₂Z^2) < L₂min
				VINCMX = (m[2]*V₁X+m[2]*V₂X)/(m[1]+m[2]) #velocity of inner center of mass
				VINCMY = (m[2]*V₁Y+m[2]*V₂Y)/(m[1]+m[2]) #velocity of inner center of mass
				VINCMZ = (m[2]*V₁Z+m[2]*V₂Z)/(m[1]+m[2]) #velocity of inner center of mass
				CM₁₂X = (M1*R₁X+M2*R₂X)/(M1+M2)
				CM₁₂Y = (M1*R₁Y+M2*R₂Y)/(M1+M2)
				CM₁₂Z = (M1*R₁Z+M2*R₂Z)/(M1+M2)
				E₁ = .5*m[1]*sqrt((V₁X-VINCMX)^2+(V₁Y-VINCMY)^2+(V₁Z-VINCMZ)^2)^2+.5*m[2]*sqrt((V₂X-VINCMX)^2+(V₂Y-VINCMY)^2+(V₂Z-VINCMZ)^2)^2 - G*m[1]*m[2]/sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)#Energy of inner binary
				E₂ = .5*(m[1]+m[2])*sqrt(VINCMX^2+VINCMY^2+VINCMZ^2)^2+.5*m[3]*sqrt(V₃X^2+V₃Y^2+V₃Z^2)^2 - G*(m[1]+m[2])*m[3]/sqrt((R₃X-CM₁₂X)^2+(R₃X-CM₁₂X)^2+(R₃X-CM₁₂X)^2)#Energy of outer binary
				L₁X = m[1]*((R₁Y-CM₁₂Y)*(V₁Z-VINCMZ)-(R₁Z-CM₁₂Z)*(V₁Y-VINCMY))+m[2]*((R₂Y-CM₁₂Y)*(V₂Z-VINCMZ)-(R₂Z-CM₁₂Z)*(V₂Y-VINCMY))
				L₁Y = m[1]*((R₁Z-CM₁₂Z)*(V₁X-VINCMX)-(R₁X-CM₁₂X)*(V₁Z-VINCMZ))+m[2]*((R₂Z-CM₁₂Z)*(V₂X-VINCMX)-(R₂X-CM₁₂X)*(V₂Z-VINCMZ))
				L₁Z = m[1]*((R₁X-CM₁₂X)*(V₁Y-VINCMY)-(R₁Y-CM₁₂Y)*(V₁X-VINCMX))+m[2]*((R₂X-CM₁₂X)*(V₂Y-VINCMY)-(R₂Y-CM₁₂Y)*(V₂X-VINCMX))
				L₂X = (m[1]+m[2])*(CM₁₂Y*VINCMZ-CM₁₂Z*VINCMY)+m[3]*(R₃Y*V₃Z-R₃Z*V₃Y)
				L₂Y = (m[1]+m[2])*(CM₁₂Z*VINCMX-CM₁₂X*VINCMZ)+m[3]*(R₃Z*V₃X-R₃X*V₃Z)
				L₂Z = (m[1]+m[2])*(CM₁₂X*VINCMY-CM₁₂Y*VINCMX)+m[3]*(R₃X*V₃Y-R₃Y*V₃X)
				if numBodies == 3 #!!make X1 arrays and stuff work again
					write(datafile,"$E, $(LX), $(LY), $(LZ), $E₁, $E₂, $(L₁X), $(L₁Y), $(L₁Z), $(L₂X), $(L₂Y), $(L₂Z), $h, $t0, $(X1[1]), $(X2[1]), $(X3[1]), $(Y1[1]), $(Y2[1]), $(Y3[1]), $(Z1[1]), $(Z2[1]), $(Z3[1])","\n")
				else
					write(datafile,"$E, $(LX), $(LY), $(LZ), $E₁, $E₂, $(L₁X), $(L₁Y), $(L₁Z), $(L₂X), $(L₂Y), $(L₂Z), $h, $t0, $(X1[1]), $(X2[1]), $(X3[1]), $(Y1[1]), $(Y2[1]), $(Y3[1]), $(Z1[1]), $(Z2[1]), $(Z3[1]), $(X4[1]), $(Y4[1]), $(Z4[1])","\n")
				end
				if sqrt(LX^2+LY^2+LZ^2) > Lmax
					Lmax = sqrt(LX^2+LY^2+LZ^2)
				elseif sqrt(LX^2+LY^2+LZ^2) < Lmin
					Lmin = sqrt(LX^2+LY^2+LZ^2)
				end
				if E > Emax
					Emax = E
				elseif E < Emin
					Emin = E
				end
				if h > lmax
					lmax = h
				elseif h < lmin
					lmin = h
				end
				if E₁ > E₁max
					E₁max = E₁
				elseif E₁ < E₁min
					E₁min = E₁
				end
				if E₂ > E₂max
					E₂max = E₂
				elseif E₂ < E₂min
					E₂min = E₂
				end
				if sqrt(L₁X^2+L₁Y^2+L₁Z^2) > L₁max
					L₁max = sqrt(L₁X^2+L₁Y^2+L₁Z^2)
				elseif sqrt(L₁X^2+L₁Y^2+L₁Z^2) < L₁min
					L₁min = sqrt(L₁X^2+L₁Y^2+L₁Z^2)
				end
				if sqrt(L₂X^2+L₂Y^2+L₂Z^2) > L₂max
					L₂max = sqrt(L₂X^2+L₂Y^2+L₂Z^2)
				elseif sqrt(L₂X^2+L₂Y^2+L₂Z^2) < L₂min
					L₂min = sqrt(L₂X^2+L₂Y^2+L₂Z^2)
				end
			end
			h1 = hParam*sqrt(R₁₂X^2+R₁₂Y^2+R₁₂Z^2)/sqrt(V₁₂X^2+V₁₂Y^2+V₁₂Z^2)
			h2 = hParam*sqrt(R₁₃X^2+R₁₃Y^2+R₁₃Z^2)/sqrt(V₁₃X^2+V₁₃Y^2+V₁₃Z^2)
			h3 = hParam*sqrt(R₂₃X^2+R₂₃Y^2+R₂₃Z^2)/sqrt(V₂₃X^2+V₂₃Y^2+V₂₃Z^2)
			if h1 < h2
				if h1 < h3
					h = h1
				else
					h = h3
				end
			else
				if h2 < h3
					h = h2
				else
					h = h3
				end
			end #this calculates the initial timestep, later this will tie into the energy of the system, once that's implemented
			if counter%10000 == 0
				update!(prog,convert(Int64,floor(t0)))
			end
			counter += 1
		end
		NowTime = time()
		write(datafile,"$m"[2:end-1],"\n","$OriginalX"[2:end-1],"\n")
		write(datafile,"$counter","\n")
		write(datafile,"$(NowTime-firstTime)","\n")
		write(datafile,"N")
	end
	if E₁>0
		stability = 0
		println("This is an unstable system.")
	elseif E₂>0
		stability = 0.5
		println("This is a partially unstable system.")
	else
		stability = 1
		println("This is a stable system.")
	end
	return hParam, t0, periods, counter, Emin, Emax, Lmin, Lmax, E₁min, E₁max, E₂min, E₂max, L₁min, L₁max, L₂min, L₂max, lmin, lmax, stability, E0, L0, E₁0, E₂0, L₁0, L₂0 #need initial velocities to write filename for saving
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
		k1[i]=h*f[i](x, m,)
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
	y=x+(k1+2*k2+2*k3+k4)/6
	return y

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
function Master(file, fileSave="AutoSave", writeData=0, MemorySave=true) #plotting L, E, or positions over time, type "L" or "E" to plot those and type a color to plot the orbits
	#Elist, Llist, lList, Tlist, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, numBodies, hParam, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, v4x, v4y, v4z, OriginalX, t0, E₁list, E₂list, L₁list, L₂list, periods, timesteps = System(file, fileSave, MemorySave)
	hParam, t0, periods, timesteps, Emin, Emax, Lmin, Lmax, E₁min, E₁max, E₂min, E₂max, L₁min, L₁max, L₂min, L₂max, lmin, lmax, stability, E0, L0, E₁0, E₂0, L₁0, L₂0 = System(file, fileSave, MemorySave)
	datafile = "h≈(r÷v) data files/$fileSave"*".txt"
	m = parse.(Float64,split(readlines(datafile)[end-4],",")) #keep
	OriginalX = parse.(Float64,split(readlines(datafile)[end-3],",")) #keep
	numBodies = parse.(Float64,split(readlines(datafile)[1],","))[1]
	timeTaken = parse(Float64,readlines(datafile)[end-1])
	#=stability calculation=#
	println("\n")
	Lmin = (Lmin-L0)/(L0)
	Lmax = (Lmax-L0)/(L0)
	Emin = (Emin-E0)/(E0)
	Emax = (Emax-E0)/(E0)
	L₁min = (L₁min-L₁0)/(L₁0)
	L₁max = (L₁max-L₁0)/(L₁0)
	L₂min = (L₂min-L₂0)/(L₂0)
	L₂max = (L₂max-L₂0)/(L₂0)
	E₁min = (E₁min-E₁0)/(E₁0)
	E₁max = (E₁max-E₁0)/(E₁0)
	E₂min = (E₂min-E₂0)/(E₂0)
	E₂max = (E₂max-E₂0)/(E₂0)
	println("The timestep varied from $lmin to $lmax.")
	println("The angular momentum varied by $((Lmin)) to $((Lmax)) while the energy varied by $Emin to $Emax.") #magnitude of angular momentum here for simplicity
	println("This ran in $timeTaken seconds.")
	println("This took $timesteps timesteps to simulate.")
	record = true
	rowNumber = 0
	if writeData == 0
		XLSX.openxlsx("NestedBinaryData.xlsm",mode="rw") do xf
			sheet = xf[1]
			i = 1
			while typeof(sheet["A$i"]) != Missing #gets next blank row
				if [sheet["A$i"],sheet["B$i"],sheet["C$i"],sheet["D$i"],sheet["E$i"],sheet["F$i"],sheet["G$i"],sheet["H$i"],sheet["I$i"],sheet["J$i"],sheet["K$i"]] == [m[1],m[2],m[3],OriginalX[1],OriginalX[2],OriginalX[3],OriginalX[4],OriginalX[5],OriginalX[6],t0,hParam] || [sheet["A$i"],sheet["B$i"],sheet["C$i"],sheet["D$i"],sheet["E$i"],sheet["F$i"],sheet["G$i"],sheet["H$i"],sheet["I$i"],sheet["J$i"],sheet["K$i"]] == [m[1],m[2],m[3],OriginalX[1],OriginalX[2],OriginalX[3],OriginalX[4],OriginalX[5],OriginalX[6],"$periods periods",hParam]
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
				sheet["L$i"] = "[$Emin,$Emax]"
				sheet["M$i"] = "[[$(Lmax[1]),$(Lmin[2]),$(Lmin[3])],[$(Lmax[1]),$(Lmax[2]),$(Lmax[3])]]"
				sheet["N$i"] = "[$E₁min,$E₁max]"
				sheet["O$i"] = "[$E₂min,$E₂max]"
				sheet["P$i"] = "[[$(L₁min[1]),$(L₁min[2]),$(L₁min[3])],[$(L₁max[1]),$(L₁max[2]),$(L₁max[3])]]"
				sheet["Q$i"] = "[[$(L₂min[1]),$(L₂min[2]),$(L₂min[3])],[$(L₂max[1]),$(L₂max[2]),$(L₂max[3])]]"
				sheet["R$i"] = timeTaken
				sheet["S$i"] = "[$lmin,$lmax]"
				sheet["T$i"] = stability
				sheet["U$i"] = rowNumber
				if Emin<-10.0^-3 || Emax>10.0^-3
					sheet["V$i"] = 0
				else
					sheet["V$i"] = 1
				end
				sheet["W$i"] = timesteps
			end
		end
	end
	return record, rowNumber #used for autotester
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