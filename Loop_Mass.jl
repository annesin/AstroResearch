include("Tester_Inner.jl")

function massloop(m, a2, theta=0, P=5, t="1000P", hParam=0.01, fileSave="AutoSave")
    #change to 10000P
    Minner = m[1]
    Md = m[2]
    if m[1] < 1.3
        error("Total inner mass must be higher than 1.3Msun to be greater than one neutron star")
    elseif m[2]>m[1]
        error("The donor star must be less massive than the inner system")
    end
    timelist = split(t, "P") #note: split removes the P
    if t[end] == 'P' #check  that this does not depend on theta
        t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
        notPeriods = timelist[1]*timelist[2]
    else
        t = parse(Float64,timelist[1])
        notPeriods = true
    end
    M1 = 1.3
    M2 = Minner - 1.3
    mstep = (.5*Minner-1.3)/5 #create a step to change the mass ratio by
    stability = 0
    stabmasslist = [] #will contain tuples of masses, a1 and stability
    while M1 <= M2
        M = [M1,M2,Md]
        gr = 1.135((t/(10^9))*M1*M2*(M1+M2))^(1/4) #GR CONDITION: assuming e = 0, t in gigayears?
        a1 = StabilityFinder2(M, a2, theta, P, t, hParam, fileSave)
        if a1 < gr
            stability = 0
        else
            stability = 1
        end
        stabletuple = (M1, M2, a1, stability)
        push!(stabmasslist, stabletuple)
        y = open("Loop_$([m[1],m[2],a2]).txt","a") #creating file to record outputs, "a" means append
        write(y,"a1 is $a1, stability is $stability\n")
        close(y)
        M1 += mstep
        M2 = Minner - M1
    end
    println(stabmasslist)
    #find the maximum a1
    maxa1 = 0
    maxindex = 0
    for i = 1:5
        if stabmasslist[i][3] == 1 #if it's stable
            if stabmasslist[i][4] > maxa1
                maxa1 = stabmasslist[i][4]
                maxindex = i
            end
        end
    end
    q = stabmasslist[maxindex][1]/stabmasslist[maxindex][2]
    return q
end
