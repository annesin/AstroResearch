include("AutomaticTester.jl")

function Looper(m, a1, F, P, t="1000P", hParam=0.01, fileSave="AutoSave", writeData=0)
    while a1 <= F
        StabilityFinder(m, a1, P, t, hParam, fileSave, writeData)
        a1 = 1.2*a1
    end
end