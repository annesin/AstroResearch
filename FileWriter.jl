#If you just want to run a specific system as a test, this will write the file for 
#NestedBinaryFinal.jl for you!

include("NestedBinaryFinal_Fast.jl")

function fileWriter(m, a, percent,t = "1000P",hParam = 0.01,fileSave = "AutoSave",MemorySave="hybrid")
    x = open("Test_$([m[1],m[2],m[3],a[1],a[3]]).txt","w") #creating the input file
    write(x,"$m"[2:end-1],"\n")
    write(x,"$(a[1]),$(a[2]),$(a[3]),$(a[4]),$(a[5]),$(a[6])\n")
    write(x,"$t,$hParam, $percent")
    close(x)
    println("The parameters here are masses $m, separations $(a[1]) and $(a[3]), with zero eccentricity, inclination, or angular offset.")
    record, dataFileExists, rowNumber, stability = Master("Test_$([m[1],m[2],m[3],a[1],a[3]]).txt",true,"AutoSave_$([m[1],m[2],m[3],a[1],a[3]])",0,MemorySave) #run simulation, get stability. Note that we're saving stuff for the heck of it
    rm("Test_$([m[1],m[2],m[3],a[1],a[3]]).txt") #deleting the input .txt file
    if MemorySave!="all"
        if record && dataFileExists #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            mv("data_files/AutoSave_$([m[1],m[2],m[3],a[1],a[3]]).txt","data_files/$fileSave"*"_$rowNumber"*".txt", force=true)
        elseif dataFileExists
            rm("data_files/AutoSave_$([m[1],m[2],m[3],a[1],a[3]]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
    end
end
