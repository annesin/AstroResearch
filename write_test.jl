function test1(x)
    i=0
    while i<10
        y = open("out_test_$x","a")
        write(y,"hello world $i")
        close(y)
        i+=1
    end
end
