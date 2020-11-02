include("NestedBinaryFinal.jl") 

function StabilityFinder2(m, a2, percent, t="1000P", hParam=0.01, fileSave="AutoSave")
    #here, m is an stringed array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
    #gotta get t here for the GR condition, even though NestedBinary will calculate it later
    timelist = split(t, "P") #note: split removes the P
    if t[end] == 'P'
        t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
        notPeriods = timelist[1]*timelist[2]
    else
        t = parse(Float64,timelist[1])
        notPeriods = true
    end
    a1= 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in gigayears?
    stability = 0 
    sizeStep = 0 #this determines the amount we'll change a2 by with each test
    counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
    println("test before while 1")
    while sizeStep < 2 #This line indicates the number of significant figures that will be evaluated. The higher, the longer the runtime
        println("test after while 1")
        x = open("Test_$([m[1],m[2],m[3],a1,a2]).txt","w") #creating the input file
        write(x,"$m"[2:end-1],"\n")
        println("m array thing is ", m[2:end-1])
        write(x,"$a1,0,$a2,0,0,0\n")
        write(x,"$t,$hParam, $percent")
        close(x)
        println("The parameters here are masses $m, separations $a1 and $a2, with zero eccentricity, inclination, or angular offset.")
        println("before Master") #gets here
        println(m)
        record, rowNumber, stability = Master("Test_$([m[1],m[2],m[3],a1,a2]).txt",true,"AutoSave_$([m[1],m[2],m[3],a1,a2])") #run simulation, get stability. Note that we're saving stuff for the heck of it
        println("after Master") #but not here
        rm("Test_$([m[1],m[2],a1,a2]).txt") #deleting the input .txt file
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
        if stability == 1 && counter >= 4 #this is the case if we have a stable system after 5 consecutive checks
            sizeStep += 1 #we're going to narrow our search case
            a1 -= 5*10.0^(-(sizeStep-1))-10.0^-sizeStep #we cancel those 5 checks we had before to get back to the most recent unstable system, then increase increase a2
            counter = 0 #we reset the counter too
        else #in this case, we either have an unstable system or a stable system that we just encountered
            a1 += 10.0^-sizeStep  #we increase a1
            if stability == 1 #if we do have a stable system
                counter += 1 #we increase the consecutive number of systems we've found to be stable
            else #if the system's actually not stable
                counter = 0 #we reset the number of consecutive stable system we've found
            end
        end
        println("test end of while 1")
    end
    println("after while 1")
    XLSX.openxlsx("StabilityConditions.xlsx",mode="rw") do xf #then we save in this nifty spreadsheet
        sheet = xf[1]
        i = 1
        record = true
        while typeof(sheet["A$i"]) != Missing #gets next blank row
            if [sheet["A$i"],sheet["B$i"],sheet["C$i"],sheet["D$i"],sheet["F$i"],sheet["G$i"]] == [m[1],m[2],a1,t,hParam]
                println("Not saving to spreadsheet: This data already has an entry at line $i.")
                record = false
            end
            i += 1
        end
        rowNumber = i
        if record
            sheet["A$i"] = m[1]
            sheet["B$i"] = m[2]
            sheet["C$i"] = round(a1 + 10.0^(-(sizeStep-1));digits=2)
            sheet["D$i"] = a2
            sheet["E$i"] = t
            sheet["F$i"] = hParam
        end
    end
    return round(a1 + 10.0^(-(sizeStep-1));digits=2)  #returns smallest a1 that is stable, digits denotes the number of sig figs after the decimal evaluated
end

