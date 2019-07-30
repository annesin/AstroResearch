LlistCreator(List) = try
    parse.(Float64,split(List,","))
catch
    List = replace(List, "["=>"")
    firstArray = split(List,"],")
    firstArray[end]=firstArray[end][1:end-1]
    [parse.(Float64,split(x,",")) for x in firstArray]
end

file = "h≈(r÷v) data files/L.txt"
type = readlines(file)[end]
numBodies = parse.(Int64,readlines(file)[1])

datafile = file
Elist = []
L1list = []
L2list = []
L3list = []
E₁list = []
E₂list = []
L₁1list = []
L₁2list = []
L₁3list = []
L₂1list = []
L₂2list = []
L₂3list = []
lList = []
Tlist = []
X1 = []
X2 = []
X3 = []
Y1 = []
Y2 = []
Y3 = []
Z1 = []
Z2 = []
Z3 = []
if numBodies==4
    X4 = []
    Y4 = []
    Z4 = []
end
timesteps = parse(Int64,readlines(datafile)[end-2])
timeTaken = parse(Float64,readlines(datafile)[end-1])

page = readlines(datafile)
for i in 2:length(page)-5
    thisArray = parse.(Float64,split(page[i],","))
    push!(Elist, thisArray[1])
    push!(L1list, thisArray[2])
    push!(L2list, thisArray[3])
    push!(L3list, thisArray[4])
    push!(E₁list, thisArray[5])
    push!(E₂list, thisArray[6])
    push!(L₁1list, thisArray[7])
    push!(L₁2list, thisArray[8])
    push!(L₁3list, thisArray[9])
    push!(L₂1list, thisArray[10])
    push!(L₂2list, thisArray[11])
    push!(L₂3list, thisArray[12])
    push!(lList, thisArray[13])
    push!(Tlist, thisArray[14])
    push!(X1, thisArray[15])
    push!(X2, thisArray[16])
    push!(X3, thisArray[17])
    push!(Y1, thisArray[18])
    push!(Y2, thisArray[19])
    push!(Y3, thisArray[20])
    push!(Z1, thisArray[21])
    push!(Z2, thisArray[22])
    push!(Z3, thisArray[23])
    if numBodies==4
        push!(X4, thisArray[24])
        push!(Y4, thisArray[25])
        push!(Z4, thisArray[26])
    else
        Z4 = [0]
    end
end
Llist = []
L₁list = []
L₂list = []
for i in 1:length(L1list)
    push!(Llist,[L1list[i],L2list[i],L3list[i]])
    push!(L₁list,[L₁1list[i],L₁2list[i],L₁3list[i]])
    push!(L₂list,[L₂1list[i],L₂2list[i],L₂3list[i]])
end

file = "h≈(r÷v) data files/Master.txt"
type = readlines(file)[end]
numBodies = parse.(Int64,readlines(file)[1])
LlistM = LlistCreator(readlines(file)[2])
ElistM = parse.(Float64,split(readlines(file)[3],",")) 
TlistM = parse.(Float64,split(readlines(file)[4],",")) 
lListM = parse.(Float64,split(readlines(file)[5],",")) 
X1M = parse.(Float64,split(readlines(file)[6],",")) 
X2M = parse.(Float64,split(readlines(file)[7],",")) 
Y1M = parse.(Float64,split(readlines(file)[8],",")) 
Y2M = parse.(Float64,split(readlines(file)[9],",")) 
if numBodies > 2
    Z1M = parse.(Float64,split(readlines(file)[10],","))
    Z2M = parse.(Float64,split(readlines(file)[11],","))
    X3M = parse.(Float64,split(readlines(file)[12],","))
    Y3M = parse.(Float64,split(readlines(file)[13],","))
    Z3M = parse.(Float64,split(readlines(file)[14],","))
    if numBodies > 3
        X4 = parse.(Float64,split(readlines(file)[15],","))
        Y4 = parse.(Float64,split(readlines(file)[16],","))
        Z4 = parse.(Float64,split(readlines(file)[17],","))
    else  
        Z4 = [] #for deciding whether or not to plot 3D
    end
end
if type == "N"
    L₁listM = LlistCreator(readlines(file)[end-7])
    L₂listM = LlistCreator(readlines(file)[end-6])
    timeTakenM = parse.(Float64,split(readlines(file)[end-5],","))[1]
    E₁listM = parse.(Float64,split(readlines(file)[end-4],","))
    E₂listM = parse.(Float64,split(readlines(file)[end-3],","))
end