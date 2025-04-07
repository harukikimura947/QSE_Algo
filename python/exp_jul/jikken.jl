module module_mys
    export MyStruct, sum_x

    struct MyStruct
        x
        y
    end

    function sum_x(x_1, x_2)
        
        y = x_1 + x_2
    
        return y
    
    end
    
end