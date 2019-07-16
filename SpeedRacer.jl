#welcome to I/O hell
include("NestedBinary.jl")

const jobs = Channel(8316)

const results = Channel(8316)

Masses = ["8,8,1","1.5,1.5,1","8,5,1","8,1.5,1"]
Eccentricities = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]

function do_work(fileSave)
    for job_id in jobs
        hParam = 0.01
        #job_id=fetch(jobs)
        println(job_id)
        x = open("Test_$job_id"*".txt","w")
        write(x,"$(Masses[convert(Int64,job_id[1])])\n")
		write(x,"$(10^(job_id[2])),$(Eccentricities[convert(Int64,job_id[3])]),$(-job_id[4]+11),0,0,0\n") 
		write(x,"100P,$hParam")
        close(x)
        println("The parameters here are [$(Masses[convert(Int64,job_id[1])])] masses, $(10^(job_id[2])) and $(-job_id[4]+11) seperations, and $(Eccentricities[convert(Int64,job_id[3])]) eccentricity.")
        record, rowNumber = Plot("Test_$job_id"*".txt", "none", "$fileSave"*"_$job_id.txt")
        put!(results,job_id)
		rm("Test_$job_id"*".txt")
		if record
			mv("h≈(r÷v) data files/$fileSave"*"_$job_id.txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
		else
			rm("h≈(r÷v) data files/$fileSave"*"_$job_id"*".txt")#deletes text file if data wasn't recorded in spreadsheet
        end
    end
end

function make_jobs(a=1,b=0,c=1,d=1)
    for i in 1:4, j in 0:-0.1:-2, k in 1:11, l in 1:9
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

function Speed(fileSave,a=1,b=0,c=1,d=1)
    make_jobs(a,b,c,d)
    counter = sizeof(jobs)
    numCores = length(Sys.cpu_info())
    println("Running with $numCores cores...")
    for i in 1:numCores # start 12 tasks to process requests in parallel
        @async do_work(fileSave)
    end
    while counter > 0 # print out results
        job_id = take!(results)
        println("Done!")
        counter = counter - 1
    end
end
