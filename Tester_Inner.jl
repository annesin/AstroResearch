include("NestedBinaryFinal.jl")

function StabilityFinder2(m, a2, theta, percent, t="1000P", hParam=0.01, fileSave="AutoSaveIn")
    #Note: Tester_Inner and Automatic tester have different meanings of percent and precision
    #here, m is an array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
    #percent is such that 5 is 5%

    #have to get t here for the GR condition, even though NestedBinary will calculate it later
       #since we don't calculate gr, we don't need the t calculation
       #but if we need to calculate gr in a later iteration, just uncomment

    # timelist = split(t, "P") #note: split removes the P
    # if t[end] == 'P' #check  that this does not depend on theta
    #     t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
    #     notPeriods = timelist[1]*timelist[2]
    # else
    #     t = parse(Float64,timelist[1])
    #     notPeriods = true
    # end
 
    #gr = 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in Megayears
    
    ### Initialization ###
    a1 = .5*a2 #set initial a1
    stability = 0 #initially unstable
    diff = a1 #the differences between consecutive guesses
    precision = 100
    divfactor = 1.1 #for 10% decreases, making sizesteps logarithmic
    counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
    #i = 0 #counter to reset in while loops
    minunstable = a2 #minimum unstable a1 found so far
    maxstable = 0 #maximum stable a1 found so far
    avg = 0
    islandcheck = false #we have not stepped 5 times yet to confirm we aren't in an island of stability
    #removed finding stability

    ### Finding precise a1 ###
    while precision > percent || islandcheck == false #islandcheck only evaluated if precision <= percent
        println("before first Master\n") #call TestIn to not mixup with A_T files
        
        ### Making file for Master ###
        x = open("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
        write(x,"$m"[2:end-1],"\n")
        write(x,"$a1,0,$a2,0,$theta,0\n")
        write(x,"$t,$hParam, $percent")
        close(x)

        ### Find stability ###
        record, rowNumber, stability = Master("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta])") #run simulation, get stability. Note that we're saving stuff for the heck of it
        println("after first Master\n")
        
        ### Removing file and recording information ###
        rm("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt") #deleting the input .txt file
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            #use data_files instead of h≈(r÷v) data files if this doesn't work
            mv("h≈(r÷v) data files/AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end

        ### Recording current a1, stability, counter, and a2 ###
        y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
        write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
        write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
        close(y) #by opening and closing each time, should record even if pipe is broken
        
        ### Stability Cases ###

        ## Stable, just finished checking islands ##
        if stability == 1 && counter > 4 #this is the case if we have a stable system after 5 consecutive checks 
            println("test stable big count")
            a1 = minunstable #reset a1 to the minimum unstable we've found
            println("after while a1 is $a1")
            divfactor = 1+(divfactor-1)/10 #decrease divfactor to become more precise, eg. from 1.1 to 1.01
            a1 = a1/divfactor #decrease a1 to search for stable case
            println("                 stab 1 and counter 4")
            println("new a1 is $a1")
            println("before stable count 4 master\n")
            #removed finding stability
            counter = 0 #reset
            islandcheck = true #can skip the 5 step search now
            println("stable, counter >= 4")
            println("precision is $precision")
            y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
            write(y,"Stable, just finished island check\n")
            close(y)
        
          ## Stable, has checked islands - want to go back up to unstable now ##
          elseif stability == 1 && islandcheck == true #we don't have to keep checking islands after the first time
            maxstable = a1 #update the most recent maxstable to current a1
            divfactor = 1+(divfactor-1)/10
            #a1 = a1*divfactor
            a1 = minunstable #reset a1 to minunstable
            a1 = a1/divfactor #start finding smaller a1 to look for more precise a1
            #maxstable = a1
            println("                 stab 1 and islandcheck")
            println("new a1 is $a1")
            println("before stable islandcheck master\n")
            #removed finding stability
            println("stable, islandcheck = true")
            println("precision is $precision")
            y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
            write(y,"Stable, islandcheck is true\n")
            close(y)
        
        ## Either stable but still checking island or unstable ##
        else
            #islandcheck = false
            #a1 = a1/divfactor
            #removed finding stability
            println("small count, before if/else")

            ## Stable, still checking islands ##
            if stability == 1 && counter <= 4 #stable system, but not enough consecutive steps yet
                println("test stable but small count")
                if counter == 0
                    maxstable = a1 #if this is the first stable a1 we've found, set maxstable
                end
                counter+=1 #increase counter for nbumber of checks
                a1 = a1/divfactor #decrease a1 to check next one is stable
                println("                    stab 1 counter <4")
                println("new a1 is $a1")
                y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
                write(y,"Stable, counter less than 4\n")
                close(y)

            ## Unstable ##    
            else #in this case, we have an unstable system
                println("test unstable")
                counter = 0 #reset counter in case we found an island
                minunstable = a1 #a1 will reset to minunstable after five checks, will this cause problem?
                if minunstable < maxstable #found an island of stability
                    maxstable = 0
                end
                a1 = a1/divfactor #decrease to look for next a1
                println("                       stab 0")
                println("new a1 is $a1")
                y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
                write(y,"Unstable\n")
                close(y)
            end
            println("precision is $precision")
        end

        y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
        write(y,"final\n")
        write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
        write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
        close(y)

        #=
        rm("Test_$([m[1],m[2],m[3],a1,a2]).txt") #deleting the input .txt file
        if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
            mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
        else
            rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
        end
        =#

        ### Calculate current precision ###
        diff = minunstable - maxstable
        println("minunstable = $minunstable, maxstable = $maxstable, diff is $diff")
        avg = (minunstable+maxstable)/2
        precision = diff/avg
        println("test end of while 1")
    end

    println("after while 1")
    ### Saving results in spreadsheet ###
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

function StabilityFinder2_noIslands(m, a2, theta, percent, t="1000P", hParam=0.01, fileSave="AutoSaveIn")
  #Note: Tester_Inner and Automatic tester have different meanings of percent and precision
  #here, m is an stringed array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
  #gotta get t here for the GR condition, even though NestedBinary will calculate it later
  #percent is such that 5 is 5%
  timelist = split(t, "P") #note: split removes the P
  if t[end] == 'P' #check  that this does not depend on theta
      t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
      notPeriods = timelist[1]*timelist[2]
  else
      t = parse(Float64,timelist[1])
      notPeriods = true
  end
  println("t is $t")
  println("G is $G")
  #gr = 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in Megayears
  a1 = .5*a2
  stability = 0 
  diff = a1 #the differences between consecutive guesses
  precision = 100 
  divfactor = 1.1 #for 10% decreases, making sizesteps logarithmic
  counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
  #i = 0 #counter to reset in while loops
  minunstable = a2
  maxstable = 0
  avg = 0
  islandcheck = true
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
      write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
      write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
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
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, just finished island check\n")
          close(y)
      elseif stability == 1 && islandcheck == true #we don't have to keep checking islands after the first time
          maxstable = a1
          divfactor = 1+(divfactor-1)/10
          #a1 = a1*divfactor
          a1 = minunstable
          a1 = a1/divfactor 
          #maxstable = a1
          println("                 stab 1 and islandcheck")
          println("new a1 is $a1")
          println("before stable islandcheck master\n")
          #removed finding stability
          println("stable, islandcheck = true")
          println("precision is $precision")
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, islandcheck is true\n")
          close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Stable, counter less than 4\n")
              close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Unstable\n")
              close(y)
          end
          println("precision is $precision")
      end
      y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
      write(y,"final\n")
      write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
      write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
      close(y)
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

function StabilityFinder2_noIslands_moveprecision(m, a2, theta, percent, t="1000P", hParam=0.01, fileSave="AutoSaveIn")
  #Note: Tester_Inner and Automatic tester have different meanings of percent and precision
  #here, m is an stringed array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
  #gotta get t here for the GR condition, even though NestedBinary will calculate it later
  #percent is such that 5 is 5%
  timelist = split(t, "P") #note: split removes the P
  if t[end] == 'P' #check  that this does not depend on theta
      t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
      notPeriods = timelist[1]*timelist[2]
  else
      t = parse(Float64,timelist[1])
      notPeriods = true
  end
  println("t is $t")
  println("G is $G")
  #gr = 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in Megayears
  a1 = .5*a2
  stability = 0 
  diff = a1 #the differences between consecutive guesses
  precision = 100 
  divfactor = 1.1 #for 10% decreases, making sizesteps logarithmic
  counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
  #i = 0 #counter to reset in while loops
  minunstable = a2
  maxstable = 0
  avg = 0
  islandcheck = true
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
      diff = minunstable - maxstable
      println("minunstable = $minunstable, maxstable = $maxstable, diff is $diff")
      avg = (minunstable+maxstable)/2
      precision = diff/avg
      y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
      write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
      write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
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
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, just finished island check\n")
          close(y)
      elseif stability == 1 && islandcheck == true #we don't have to keep checking islands after the first time
          maxstable = a1
          divfactor = 1+(divfactor-1)/10
          #a1 = a1*divfactor
          a1 = minunstable
          a1 = a1/divfactor 
          #maxstable = a1
          println("                 stab 1 and islandcheck")
          println("new a1 is $a1")
          println("before stable islandcheck master\n")
          #removed finding stability
          println("stable, islandcheck = true")
          println("precision is $precision")
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, islandcheck is true\n")
          close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Stable, counter less than 4\n")
              close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Unstable\n")
              close(y)
          end
          println("precision is $precision")
      end
      y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
      write(y,"final\n")
      write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
      write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
      close(y)
      #=
      rm("Test_$([m[1],m[2],m[3],a1,a2]).txt") #deleting the input .txt file
      if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
          mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
      else
          rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
      end
      =#
      println("test end of while 1")
  end
  println("after while 1")
  x = open("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
  write(x,"$m"[2:end-1],"\n")
  write(x,"$a1,0,$a2,0,$theta,0\n")
  write(x,"$t,$hParam, $percent")
  close(x)
  record, rowNumber, stability = Master("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta])")
  if stability == 0 #unstable
    a1 = maxstable
  end
  rm("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt")
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

function StabilityFinder2_stabfirst(m, a2, theta, percent, t="1000P", hParam=0.01, fileSave="AutoSaveIn")
  #Note: Tester_Inner and Automatic tester have different meanings of percent and precision
  #here, m is an stringed array of the masses [M1, M2, Mdonor] while a2 is a Float64 that is the outer separation
  #gotta get t here for the GR condition, even though NestedBinary will calculate it later
  #percent is such that 5 is 5%
  timelist = split(t, "P") #note: split removes the P
  if t[end] == 'P' #check  that this does not depend on theta
      t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
      notPeriods = timelist[1]*timelist[2]
  else
      t = parse(Float64,timelist[1])
      notPeriods = true
  end
  println("t is $t")
  println("G is $G")
  #gr = 1.135((t/(10^9))*m[1]*m[2]*(m[1]+m[2]))^(1/4) #GR CONDITION: assuming e = 0, t in Megayears
  a1 = .5*a2
  stability = 0 
  diff = a1 #the differences between consecutive guesses
  precision = 100 
  divfactor = 1.1 #for 10% decreases, making sizesteps logarithmic
  counter = 0 #this will help us test more cases, meaning that when a system is stable, we're more confident that we've stepped over that stability condition line
  #i = 0 #counter to reset in while loops
  minunstable = a2
  maxstable = 0
  avg = 0
  islandcheck = true
  #removed finding stability
  x = open("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
  write(x,"$m"[2:end-1],"\n")
  write(x,"$a1,0,$a2,0,$theta,0\n")
  write(x,"$t,$hParam, $percent")
  close(x)
  record, rowNumber, stability = Master("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta])") #run simulation, get stability. Note that we're saving stuff for the heck of it
  rm("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt") #deleting the input .txt file
  println("percent is $percent")
  while precision > percent || islandcheck == false || stability == 0 #islandcheck only evaluated if precision <= percent
  #while a1 > .5
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
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, just finished island check\n")
          close(y)
      elseif stability == 1 && islandcheck == true #we don't have to keep checking islands after the first time
          maxstable = a1
          divfactor = 1+(divfactor-1)/10
          #a1 = a1*divfactor
          a1 = minunstable
          a1 = a1/divfactor 
          #maxstable = a1
          println("                 stab 1 and islandcheck")
          println("new a1 is $a1")
          println("before stable islandcheck master\n")
          #removed finding stability
          println("stable, islandcheck = true")
          println("precision is $precision")
          y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
          write(y,"Stable, islandcheck is true\n")
          close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Stable, counter less than 4\n")
              close(y)
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
              y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
              write(y,"Unstable\n")
              close(y)
          end
          println("precision is $precision")
      end
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
      diff = minunstable - maxstable
      println("minunstable = $minunstable, maxstable = $maxstable, diff is $diff")
      avg = (minunstable+maxstable)/2
      precision = (diff/avg)*100
      println("precision is $precision")
      y = open("Total_$([m[1],m[2],m[3],a2,theta,percent]).txt","a") #creating file to record outputs, "a" means append
      write(y,"a1 is $a1, stability is $stability, counter is $counter, a2 is $a2\n")
      write(y,"maxstable is $maxstable, minunstable is $minunstable, divfactor is $divfactor, precision is $precision\n")
      close(y) #by opening and closing each time, should record even if pipe is broken
      #=
      rm("Test_$([m[1],m[2],m[3],a1,a2]).txt") #deleting the input .txt file
      if record #if the data was saved in the spreadsheet, we save the output .txt file with the corresponding row number. If it wasn't, we delete the output .txt file.
          mv("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt","h≈(r÷v) data files/$fileSave"*"_$rowNumber"*".txt", force=true)
      else
          rm("h≈(r÷v) data files/AutoSave_$([m[1],m[2],m[3],a1,a2]).txt")#deletes text file if data wasn't recorded in spreadsheet
      end
      =#
      println("test end of while 1")
  end
  println("after while 1")
  x = open("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt","w") #creating the input file
  write(x,"$m"[2:end-1],"\n")
  write(x,"$a1,0,$a2,0,$theta,0\n")
  write(x,"$t,$hParam, $percent")
  close(x)
  record, rowNumber, stability = Master("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt",true,"AutoSaveIn_$([m[1],m[2],m[3],a1,a2,theta])")
  rm("TestIn_$([m[1],m[2],m[3],a1,a2,theta]).txt")
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
