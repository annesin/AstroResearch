include("FileWriter.jl")

function runtimeTester(periods, maxPeriods=Inf, MemorySave="all")
    while periods<=maxPeriods
        fileWriter([8,8,1],[8,0,20,0,0,0],5,"$periods"*"P",0.01,"AutoSave",MemorySave)
        periods *=2
    end
end