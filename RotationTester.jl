include("NestedBinaryFinal.jl")
include("AutomaticTester.jl")
ENV["MPLBACKEND"]="qt5agg" #switching the environment in which the packages are imported and run in. If you don't do this, your computer will crash and will not be able to run Julia until it is restarted
using Pkg
Pkg.add("PyPlot")
using PyPlot

#this function takes a certain fixed a1, a2, and m array and runs the simulation with varying amount of 
#rotations, that is, the three bodies will no longer always initially start in a line

function rotationStabilityFinder(m, a, angles, percent, ask=true, t="1000P", hParam=0.01, fileSave="AutoSave")
    if ask
        #we ask the user what they mean by segments
        println("Would you like $angles different ANGLES, or should each angle be spaced $angles DEGREES apart?")
        segmentType=readline()
        if segmentType == "angles" || segmentType == "ANGLES"
            angles = 360/angles
        end
    else
        segmentType = "angles"
    end
    if segmentType == "degrees" || segmentType == "DEGREES" || segmentType == "angles" || segmentType == "ANGLES"
        a1=a[1]
        a2=a[2]
        segments = floor(360/angles)
        for i in 1:segments
            theta = (i-1)*angles
            x = open("Test_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
            write(x,"$m"[2:end-1],"\n")
            write(x,"$a1,0,$a2,0,$theta,0\n")
            write(x,"$t,$hParam, $percent")
            close(x)
            println("The parameters here are masses $m, separations $a1 and $a2, angular offset of $theta degrees, with zero eccentricity or inclination.")
            record, rowNumber, stability = Master("Test_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSave_$([m[1],m[2],m[3],a1,a2,theta])") #run simulation, get stability. Note that we're saving stuff for the heck of it
            rm("Test_$([m[1],m[2],m[3],a1,a2,theta]).txt") #deleting the input .txt file
            if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
                mv("data_files/AutoSave_$([m[1],m[2],m[3],a1,a2,theta]).txt","data_files/$fileSave"*"_$rowNumber"*".txt", force=true)
            else
                rm("data_files/AutoSave_$([m[1],m[2],m[3],a1,a2,theta]).txt")#deletes text file if data wasn't recorded in spreadsheet
            end
        end
    else
        println("Please enter 'angles' or 'degrees'.")
    end
end
        
function autoRotator(m, a1, angles, percent, precision=2, t="1000P", hParam=0.01, fileSave="AutoSave")
    #this function finds the minimum a2 first at theta=0, then sees how stability changes at different angles
    a2 = StabilityFinder(m, a1, 0, percent, precision, t, hParam, fileSave)
    rotationStabilityFinder(m, [a1,a2], angles, percent, false, t, hParam, fileSave)
end

function minAOfTheta(m, a1, angles, percent=5, precision=2, t="1000P", hParam=0.01, fileSave="AutoSave", graph=false, saveGraph=true)
    #finds the stable a2 of all the angles
    angleStep = floor(360/angles)
    thetas = []
    for i in range(0,length=angles,step=angleStep)
        append!(thetas,i)
    end
    a2 = StabilityFinder(m, a1, 0, percent, precision, t, hParam, fileSave)
    a2s = [a2]
    #then loop through and find each
    for i in 2:length(thetas)
        append!(a2s,StabilityFinder(m, a1, thetas[i], percent, precision, t, hParam, fileSave, a2s[i-1]*0.9))
    end
    if graph
        theta = range(0.0, length=convert(Int32,angles), stop=2*pi)
        radii = a2s
        width = (2pi) / angles
        ax = plt.subplot(polar=true)
        bars = ax.bar(theta, radii, width=width)
        plt.show()
    end
    if saveGraph
        open("data_files/Rotation_Graph_$([m[1],m[2],m[3],a1,angles]).txt","w") do io#creating the input file
            write(io,"$(length(a2s))","\n")
            for i in 1:length(a2s)
                write(io, "$(thetas[i]), $(a2s[i])","\n")
            end
        end
    end
end