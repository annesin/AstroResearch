function System(f::Array{Function,1},x::Array{Float64,1}, m::Array{Float64,1}, t::Float64,  h::Float64)

	 t0 = 0

#until the desired time has been reached, the code runs RK4

	 while t0 < t

	       x = RK4(f, x, m, h)
	       t0 = t0 + h
	       
	 end
	 return x

end
