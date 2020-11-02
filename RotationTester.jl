include("NestedBinaryFinal.jl")

#this function takes a certain fixed a1, a2, and m array and runs the simulation with varying amount of 
#rotations, that is, the three bodies will no longer always initially start in a line

function rotationStabilityFinder(m, a, angles, percent, t="1000P", hParam=0.01, fileSave="AutoSave")
    #we ask the user what they mean by segments
    println("Would you like $angles different ANGLES, or should each angle be spaced $angles DEGREES apart?")
    segmentType=readline()
    if segmentType == "angles" || segmentType == "ANGLES"
        angles = 360/angles
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
        
