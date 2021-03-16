#include("NestedBinary_updated.jl") 
include("NestedBinaryFinal.jl")
#include("NestedBinary_Ben.jl")

function StabilityFinder2(m, a2, theta, percent, t="1000P", hParam=0.01, fileSave="AutoSaveIn")
    #Note: Tester_Inner and Automatic tester have different meanings of percent and precision
    #here, m is an stringed array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
    #gotta get t here for the GR condition, even though NestedBinary will calculate it later
    timelist = split(t, "P") #note: split removes the P
    if t[end] == 'P' #check  that this does not depend on theta
        t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
        notPeriods = timelist[1]*timelist[2]
    else
        t = parse(Float64,timelist[1])
        notPeriods = true
    end
    #gr = 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in gigayears?
    a1 = .5*a2
    stability = 0 
    diff = a1 #the differences between consecutive guesses
    precision = 1 
    divfactor = 1.1 #for 10% decreases, making sizesteps logarithmic
    counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
    i = 0 #counter to reset in while loops
    minunstable = a2
    maxstable = 0
    avg = 0
    islandcheck = false
    #removed finding stability
    while precision > percent || islandcheck == false #islandcheck only evaluated if precision <= percent
    #while a1 > .5
        println("before first Master\n") #call TestIn to not mixup with A_T files
        x = open("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
        write(x,"$m"[2:end-1],"\n")
        write(x,"$a1,0,$a2,0,$theta,0\n")
        write(x,"$t,$hParam, $percent")
        close(x)
        record, rowNumber, stability = Master("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta])") #run simulation, get stability. Note that we're saving stuff for the heck of it
        println("after first Master\n")
        rm("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt") #deleting the input .txt file
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            #use h≈(r÷v) data files on copernicus, data_files locally
            mv("h≈(r÷v) data files/AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
        y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
        write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2,\n")
        close(y) #by opening and closing each time, should record even if pipe is broken
        if stability == 1 && counter > 4 #this is the case if we have a stable system after 5 consecutive checks 
            println("test stable big count")
            #can this just be replaced with minunstable?
            #what about precision? is that calculated?
            #=
            while i <= 4 #reset to the original stable condition
                a1 = a1*divfactor
                i += 1
             #   println("while a1 = $a1")
            end
            =#
            a1 = minunstable
            #i = 0
            println("after while a1 is $a1")
            divfactor = 1+(divfactor-1)/10
            a1 = a1/divfactor # do we want overestimate or underestimate?
            println("                 stab 1 and counter 4")
            println("new a1 is $a1")
            println("before stable count 4 master\n")
            #removed finding stability
            counter = 0 #reset
            islandcheck = true
            println("stable, counter >= 4")
            println("precision is $precision")
        elseif stability == 1 && islandcheck == true #we don't have to keep checking islands after the first time
            maxstable = a1
            divfactor = 1+(divfactor-1)/10
            a1 = a1*divfactor # do we want overestimate or underestimate?
            #maxstable = a1
            println("                 stab 1 and islandcheck")
            println("new a1 is $a1")
            println("before stable islandcheck master\n")
            #removed finding stability
            println("stable, islandcheck = true")
            println("precision is $precision")
        else
            #islandcheck = false
            #a1 = a1/divfactor
            #removed finding stability
            println("small count, before if/else")
            if stability == 1 && counter <= 4 #stable system, but not enough consecutive steps yet
                println("test stable but small count")
                if counter == 0
                    maxstable = a1
                end
                counter+=1
                a1 = a1/divfactor
                println("                    stab 1 counter <4")
                println("new a1 is $a1")
            else #in this case, we have an unstable system
                println("test unstable")
                counter = 0
                minunstable = a1 #a1 will reset to minunstable after five checks, will this cause problem?
                if minunstable < maxstable #found an island of stability
                    maxstable = 0
                end
                a1 = a1/divfactor
                println("                       stab 0")
                println("new a1 is $a1")
            end
            println("precision is $precision")
        end
        #=
        rm("Test_$([m[1],m[2],m[3],a1,a2]).txt") #deleting the input .txt file
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
        =#
        diff = minunstable - maxstable
        println("minunstable = $minunstable, maxstable = $maxstable, diff is $diff")
        avg = (minunstable+maxstable)/2
        precision = diff/avg
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
            sheet["C$i"] = m[3]
            sheet["D$i"] = round(a1;digits=2)
            sheet["E$i"] = a2
            sheet["F$i"] = t
            sheet["G$i"] = hParam
            sheet["H$i"] = theta
        end
    end
    #return round(a1 + 10.0^(-(sizeStep-1));digits=2)  #returns smallest a1 that is stable, digits denotes the number of sig figs after the decimal evaluated
    return round(a1;digits=2)  #changed to just a1
end

#top comment
#= if a1<gr
            stability = 0
            println("a1 is $a1, which is smaller than the GR condition, $gr, so this system is unstable")
            break
        end 
        =#
        #=
        x = open("Test_$([m[1],m[2],m[3],a1,a2]).txt","w") #creating the input file
        write(x,"$m"[2:end-1],"\n")
        #println("m array thing is ", m[2:end-1])
        write(x,"$a1,0,$a2,0,0,0\n")
        write(x,"$t,$hParam, $percent")
        close(x)
        println("The parameters here are masses $m, separations $a1 and $a2, with zero eccentricity, inclination, or angular offset.")
        println("before Master") #gets here
        #println(m)
        record, rowNumber, stability = Master("Test_$([m[1],m[2],m[3],a1,a2]).txt",true,"AutoSave_$([m[1],m[2],m[3],a1,a2])") #run simulation, get stability. Note that we're saving stuff for the heck of it
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
        =#



#= 
       if stability == 1 && counter >= 4 #this is the case if we have a stable system after 5 consecutive checks
            sizeStep += 1 #we're going to narrow our search case
            #changed from - to +
            if sizeStep == maxOrder
                a1 += 4*10.0^(-(sizeStep-1)) #because we are looking for the max a1, we'd rather have an underestimate than overestimate
            else
                a1 += 5*10.0^(-(sizeStep-1))-10.0^-sizeStep #we cancel those 5 checks we had before to get back to the most recent unstable system, then increase increase a2 
                counter = 0 #we reset the counter too
            end
        else #in this case, we either have an unstable system or a stable system that we just encountered
            a1 -= 10.0^-sizeStep  #we increase a1
            if stability == 1 #if we do have a stable system
                counter += 1 #we increase the consecutive number of systems we've found to be stable
            else #if the system's actually not stable
                counter = 0 #we reset the number of consecutive stable system we've found
            end
        end
        =#
