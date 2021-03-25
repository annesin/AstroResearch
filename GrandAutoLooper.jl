include("AutoLooper.jl")

function grandLooper(a1, F, P, precision=2, t="1000P", hParam=0.01, fileSave="AutoSave")
    masses = [[1.5,1.5,1],[1.5,8,1],[8,8,1]]
    for i in masses
        Looper(i, a1, F, P, precision, t, hParam, fileSave)
    end
end