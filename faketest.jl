function faketest(x)
    array1 = []
    for i = 1:x
        tup = (0,i)
        push!(array1, tup)
    end
    println("the thing is ", array1[1][2])
    return array1
end

function grradius(m, a2, t="1000P")
    G = 2945.49
    timelist = split(t, "P") #note: split removes the P
    if t[end] == 'P' #check  that this does not depend on theta
        t = sqrt((a2)^3*(4*pi^2)/(G*(m[1]+m[2]+m[3])))*parse(Float64, timelist[1])
        notPeriods = timelist[1]*timelist[2]
    else
        t = parse(Float64,timelist[1])
        notPeriods = true
    end
    gr = 1.135((t/(10^9))*(m[1]+m[2])*m[3]*(m[1]+m[2]+m[3]))^(1/4) #GR CONDITION: assuming e = 0, t in gigayears?
    return gr
end

function grtime(m, a1)
    #M1, M2 -> both inner
    G = 2945.49
    time = ((a1/1.135)^4)/(m[1]*m[2]*(m[1]+m[2])) #assume megayears for now
    return time #in gigayears or megayears?
end
    