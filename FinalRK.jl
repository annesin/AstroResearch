function RK4(f::Array{Function,1},x::Array{Float64,1}, m::Array{Float64,1}, h::Float64)
#Inputs are initial position array, mass array and step size h

                  d=length(f)

#Setting the runge kutta parameters to zero

                  k1=zeros(length(x))
                  k2=zeros(length(x))
                  k3=zeros(length(x))
                  k4=zeros(length(x))

#caluclating the runge kutta parameters for every function in f

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
                      
