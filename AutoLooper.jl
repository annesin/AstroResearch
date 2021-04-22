include("AutomaticTester.jl")

function Looper(m, a1, F, P, precision=2, t="1000P", MemorySave="all", hParam=0.01, fileSave="AutoSave", writeData=0)
    while a1 <= F
        StabilityFinder(m, a1, 0, P, precision, t, MemorySave, hParam, fileSave, writeData)
        a1 = 1.2*a1
    end
end