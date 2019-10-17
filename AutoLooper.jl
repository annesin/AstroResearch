include("AutomaticTester.jl")

function Looper(m, a1, F, t="100P", hParam=0.01, fileSave="AutoSave")
    while a1 <= F
        StabilityFinder(m, a1, t, hParam, fileSave)
        a1 += 1
    end
end