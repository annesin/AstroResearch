#welcome to I/O hell
#This automatically runs conditions, but using multi-core threading, significantly reducing runtime.
include("NestedBinary.jl") #we need this code, obviously

#=These are channels. You can think of them as inboxes or outboxes, where objects are stored sequentially.
So, if a function takes an object from the channel, it will take the one "on top".
The objects in the jobs channel consist of a four element array, which are the parameters of the simulation.=#
const jobs = Channel(8316)

Masses = ["8,8,1","1.5,1.5,1","8,5,1","8,1.5,1"] #instead of putting these strings in the job, we'll just put those indices in the jobs
Eccentricities = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]

function do_work(fileSave) #this function takes a job from the jobs channel and runs Plot() on it
    for job_id in jobs #this goes through every job in the channel
        hParam = 0.01
        x = open("Test_$job_id"*".txt","w") #this opens a unique file in the format needed for NestedBinary.jl to run
        write(x,"$(Masses[convert(Int64,job_id[1])])\n")
		write(x,"$(10^(job_id[2])),$(Eccentricities[convert(Int64,job_id[3])]),$(-job_id[4]+11),0,0,0\n") 
		write(x,"100P,$hParam")
        close(x)
        println("The parameters here are [$(Masses[convert(Int64,job_id[1])])] masses, $(10^(job_id[2])) and $(-job_id[4]+11) seperations, and $(Eccentricities[convert(Int64,job_id[3])]) eccentricity.") #prints what parameters we're running
        record, rowNumber = Master("Test_$job_id"*".txt", "none", "$fileSave"*"_$job_id.txt") #This, finally, runs the simulation. It also records what row Excel is saving the data to (or if its saving at all)
		rm("Test_$job_id"*".txt") #deletes the input file
		if record
			mv("h≈(r÷v) data files/$fileSave"*"_$job_id.txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true) #if data was recorded in spreadsheet, it renames the .txt file to the corresponding row number
		else
			rm("h≈(r÷v) data files/$fileSave"*"_$job_id"*".txt")#deletes text file if data wasn't recorded in spreadsheet
        end
    end
end

#this function goes through all the permutations in the parameter space and creates a job from each one
function make_jobs(a=1,b=0,c=1,d=1)
    for i in 1:4, j in 0:-0.1:-2, k in 1:11, l in 1:9
        #=this mess of if/elseif statements are needed because we want to have a starting point, not a parameter space.
            For example, if we were to simply have a nested for loop with the below restrictions, it wouldn't loop back around.
            Think about it like a clock. If you say 'start at 6:30', you want the clock to loop back to 7:00, not 7:30.
            So, we check to see if the test we have is "above" our minimum requirements in this clock-like fashion in the manner below.=#
        if i>a
            put!(jobs,[i,j,k,l])
        elseif i==a
            if j<b
                put!(jobs,[i,j,k,l])
            elseif j==b
                if k>c
                    put!(jobs,[i,j,k,l])
                elseif k==c
                    if l >= d
                        put!(jobs,[i,j,k,l])
                    end
                end
            end
        end
    end
end

#the master function
function Speed(fileSave,a=1,b=0,c=1,d=1)
    make_jobs(a,b,c,d) #loops through and makes jobs from the parameter space we give it
    numCores = length(Sys.cpu_info()) #this isn't really the number of cores, but the number of parallel processes avaliable. It varies by machine.
    println("Running with $numCores cores...")
    for i in 1:numCores # start n tasks to process requests in parallel
        @async do_work(fileSave) #In parallel, runs simulations
    end
end
