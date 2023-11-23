function Prop_Dynamic(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
    for j in 2:101
        du[j] = (lambda + sigma - 1) * u[j] - (lambda + 3 * sigma) * (u[j] ^ 2 + 2 * sigma * u[j] ^ 3) 
            + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1] * (1- u[1])
    end
end

function Prop_order6(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
    #Orden 1
    for j in 2:1001
        du[j] = (lambda + sigma - 1) * u[j] - (lambda + 3 * sigma) * (u[j] ^ 2 + 2 * sigma * u[j] ^ 3) +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1] * (1- u[1])
    end
    ##Orden 2
    du[1002] = (lambda + sigma - 1) * u[1002] - (lambda + 3 * sigma) * (u[1002] ^ 2 + 2 * sigma * u[1002] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] * (1- u[3]) + sin(((3-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[1003] = (lambda + sigma - 1) * u[1003] - (lambda + 3 * sigma) * (u[1003] ^ 2 + 2 * sigma * u[1003] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[1004] = (lambda + sigma - 1) * u[1004] - (lambda + 3 * sigma) * (u[1004] ^ 2 + 2 * sigma * u[1004] ^ 3) + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[3] * (1- u[3]))
    for j in 2:4
        du[1003+j] = (lambda + sigma - 1) * u[1003+j] - (lambda + 3 * sigma) * (u[1003+j] ^ 2 + 2 * sigma * u[1003+j] ^ 3) + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j] * (1- u[j])
    end
    #Orden 3
    du[1008] = (lambda + sigma - 1) * u[1008] - (lambda + 3 * sigma) * (u[1008] ^ 2 + 2 * sigma * u[1008] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] * (1- u[1006]) + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002] * (1- u[1002]))
    du[1009] = (lambda + sigma - 1) * u[1009] - (lambda + 3 * sigma) * (u[1009] ^ 2 + 2 * sigma * u[1009] ^ 3) + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] * (1- u[1002]) + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005] * (1- u[1005]))
    du[1010] = (lambda + sigma - 1) * u[1010] - (lambda + 3 * sigma) * (u[1010] ^ 2 + 2 * sigma * u[1010] ^ 3) + sigma * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[1005] * (1- u[1005]))
    du[1011] = (lambda + sigma - 1) * u[1011] - (lambda + 3 * sigma) * (u[1011] ^ 2 + 2 * sigma * u[1011] ^ 3) + sigma * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[1006] * (1- u[1006]))
    #Orden 4
    du[1012] = (lambda + sigma - 1) * u[1012] - (lambda + 3 * sigma) * (u[1012] ^ 2 + 2 * sigma * u[1012] ^ 3) + sigma * sqrt(2/T) * (sqrt(4) * sin(((2-1) * pi * t)/ T) * u[1010] * (1- u[1010]))
    du[1013] = (lambda + sigma - 1) * u[1013] - (lambda + 3 * sigma) * (u[1013] ^ 2 + 2 * sigma * u[1013] ^ 3) + sigma * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[1011] * (1- u[1011]))
    #Orden 5
    du[1014] = (lambda + sigma - 1) * u[1014] - (lambda + 3 * sigma) * (u[1014] ^ 2 + 2 * sigma * u[1014] ^ 3) + sigma * sqrt(2/T) * (sqrt(5) * sin(((2-1) * pi * t)/ T) * u[1012] * (1- u[1012]))
    du[1015] = (lambda + sigma - 1) * u[1015] - (lambda + 3 * sigma) * (u[1015] ^ 2 + 2 * sigma * u[1015] ^ 3) + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013] * (1- u[1013]))
    #Orden 6
    du[1016] = (lambda + sigma - 1) * u[1016] - (lambda + 3 * sigma) * (u[1016] ^ 2 + 2 * sigma * u[1016] ^ 3) + sigma * sqrt(2/T) * (sqrt(6) * sin(((2-1) * pi * t)/ T) * u[1014] * (1- u[1014]))
    du[1017] = (lambda + sigma - 1) * u[1017] - (lambda + 3 * sigma) * (u[1017] ^ 2 + 2 * sigma * u[1017] ^ 3) + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015] * (1- u[1015]))
end
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
function Prop_order12_BM0(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
end

function Prop_order12_BM10(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
    #Orden 1
    for j in 2:11
        du[j] = (lambda + sigma - 1) * u[j] - (lambda + 3 * sigma) * (u[j] ^ 2 + 2 * sigma * u[j] ^ 3) + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1] * (1- u[1])
    end
    ##Orden 2 mixed
    du[12] = (lambda + sigma - 1) * u[12] - (lambda + 3 * sigma) * (u[12] ^ 2 + 2 * sigma * u[12] ^ 3) + sqrt(2/T) * sigma * (sin(((2-1) * pi * t)/ T) * u[3] * (1- u[3]) 
        + sin(((3-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[13] = (lambda + sigma - 1) * u[13] - (lambda + 3 * sigma) * (u[13] ^ 2 + 2 * sigma * u[13] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] * (1- u[4]) 
        + sin(((4-1) * pi * t)/ T) * u[2] * (1- u[2]) )
    du[14] = (lambda + sigma - 1) * u[14] - (lambda + 3 * sigma) * (u[14] ^ 2 + 2 * sigma * u[14] ^ 3) + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] * (1- u[4]) 
        + sin(((4-1) * pi * t)/ T) * u[3] * (1- u[3]) )
    ##Orden 2 diagonal
    for j in 2:11
        du[13+j] = (lambda + sigma - 1) * u[13+j] - (lambda + 3 * sigma) * (u[13+j] ^ 2 + 2 * sigma * u[13+j] ^ 3) 
            + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j] * (1- u[j])
    end
    #Orden 3 mixed
    du[25] = (lambda + sigma - 1) * u[25] - (lambda + 3 * sigma) * (u[25] ^ 2 + 2 * sigma * u[25] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[16] * (1- u[16]) 
        + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[12] * (1- u[12]))
    du[26] = (lambda + sigma - 1) * u[26] - (lambda + 3 * sigma) * (u[26] ^ 2 + 2 * sigma * u[26] ^ 3) + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[12] * (1- u[12]) 
        + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[15] * (1- u[15]))
    ##Orden 3 diagonal
    for j in 2:11
        du[25+j] = (lambda + sigma - 1) * u[25+j] - (lambda + 3 * sigma) * (u[25+j] ^ 2 + 2 * sigma * u[25+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[13+j] * (1- u[13+j]))
    end
    #Orden 4
    for j in 2:11
        du[35+j] = (lambda + sigma - 1) * u[35+j] - (lambda + 3 * sigma) * (u[35+j] ^ 2 + 2 * sigma * u[35+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[25+j] * (1- u[25+j]))
    end
    #Orden 5
    for j in 2:11
        du[45+j] = (lambda + sigma - 1) * u[45+j] - (lambda + 3 * sigma) * (u[45+j] ^ 2 + 2 * sigma * u[45+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[35+j] * (1- u[35+j]))
    end
    #Orden 6
    for j in 2:11
        du[55+j] = (lambda + sigma - 1) * u[55+j] - (lambda + 3 * sigma) * (u[55+j] ^ 2 + 2 * sigma * u[55+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[45+j] * (1- u[45+j]))
    end
    #Orden 7
    for j in 2:11
        du[65+j] = (lambda + sigma - 1) * u[65+j] - (lambda + 3 * sigma) * (u[65+j] ^ 2 + 2 * sigma * u[65+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[55+j] * (1- u[55+j]))
    end
    #Orden 8
    for j in 2:11
        du[75+j] = (lambda + sigma - 1) * u[75+j] - (lambda + 3 * sigma) * (u[75+j] ^ 2 + 2 * sigma * u[75+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[65+j] * (1- u[65+j]))
    end
    #Orden 9
    for j in 2:11
        du[85+j] = (lambda + sigma - 1) * u[85+j] - (lambda + 3 * sigma) * (u[85+j] ^ 2 + 2 * sigma * u[85+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[75+j] * (1- u[75+j]))
    end
    #Orden 10
    for j in 2:11
        du[95+j] = (lambda + sigma - 1) * u[95+j] - (lambda + 3 * sigma) * (u[95+j] ^ 2 + 2 * sigma * u[95+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[85+j] * (1- u[85+j]))
    end
    #Orden 11
    for j in 2:11
        du[105+j] = (lambda + sigma - 1) * u[105+j] - (lambda + 3 * sigma) * (u[105+j] ^ 2 + 2 * sigma * u[105+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[95+j] * (1- u[95+j]))
    end
    #Orden 12
    for j in 2:11
        du[115+j] = (lambda + sigma - 1) * u[115+j] - (lambda + 3 * sigma) * (u[115+j] ^ 2 + 2 * sigma * u[115+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[105+j] * (1- u[105+j]))
    end
end
#####
function Prop_order12_BM100(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
    #Orden 1
    for j in 2:101
        du[j] = (lambda + sigma - 1) * u[j] - (lambda + 3 * sigma) * (u[j] ^ 2 + 2 * sigma * u[j] ^ 3) + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1] * (1- u[1])
    end
    ##Orden 2 mixed
    du[102] = (lambda + sigma - 1) * u[102] - (lambda + 3 * sigma) * (u[102] ^ 2 + 2 * sigma * u[102] ^ 3) + sqrt(2/T) * sigma * (sin(((2-1) * pi * t)/ T) * u[3] * (1- u[3]) + sin(((3-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[103] = (lambda + sigma - 1) * u[103] - (lambda + 3 * sigma) * (u[103] ^ 2 + 2 * sigma * u[103] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[104] = (lambda + sigma - 1) * u[104] - (lambda + 3 * sigma) * (u[104] ^ 2 + 2 * sigma * u[104] ^ 3) + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[3] * (1- u[3]))
    ##Orden 2 diagonal
    for j in 2:101
        du[103+j] = (lambda + sigma - 1) * u[103+j] - (lambda + 3 * sigma) * (u[103+j] ^ 2 + 2 * sigma * u[103+j] ^ 3) + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j] * (1- u[j])
    end
    #Orden 3 mixed
    du[205] = (lambda + sigma - 1) * u[205] - (lambda + 3 * sigma) * (u[205] ^ 2 + 2 * sigma * u[205] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[106] * (1- u[106]) + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[102] * (1- u[102]))
    du[206] = (lambda + sigma - 1) * u[206] - (lambda + 3 * sigma) * (u[206] ^ 2 + 2 * sigma * u[206] ^ 3) + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[102] * (1- u[102]) + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[105] * (1- u[105]))
    ##Orden 3 diagonal
    for j in 2:101
        du[205+j] = (lambda + sigma - 1) * u[205+j] - (lambda + 3 * sigma) * (u[205+j] ^ 2 + 2 * sigma * u[205+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[103+j] * (1- u[103+j]))
    end
    #Orden 4
    for j in 2:101
        du[305+j] = (lambda + sigma - 1) * u[305+j] - (lambda + 3 * sigma) * (u[305+j] ^ 2 + 2 * sigma * u[305+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[205+j] * (1- u[205+j]))
    end
    #Orden 5
    for j in 2:101
        du[405+j] = (lambda + sigma - 1) * u[405+j] - (lambda + 3 * sigma) * (u[405+j] ^ 2 + 2 * sigma * u[405+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[305+j] * (1- u[305+j]))
    end
    #Orden 6
    for j in 2:101
        du[505+j] = (lambda + sigma - 1) * u[505+j] - (lambda + 3 * sigma) * (u[505+j] ^ 2 + 2 * sigma * u[505+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[405+j] * (1- u[405+j]))
    end
    #Orden 7
    for j in 2:101
        du[605+j] = (lambda + sigma - 1) * u[605+j] - (lambda + 3 * sigma) * (u[605+j] ^ 2 + 2 * sigma * u[605+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[505+j] * (1- u[505+j]))
    end
    #Orden 8
    for j in 2:101
        du[705+j] = (lambda + sigma - 1) * u[705+j] - (lambda + 3 * sigma) * (u[705+j] ^ 2 + 2 * sigma * u[705+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[605+j] * (1- u[605+j]))
    end
    #Orden 9
    for j in 2:101
        du[805+j] = (lambda + sigma - 1) * u[805+j] - (lambda + 3 * sigma) * (u[805+j] ^ 2 + 2 * sigma * u[805+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[705+j] * (1- u[705+j]))
    end
    #Orden 10
    for j in 2:101
        du[905+j] = (lambda + sigma - 1) * u[905+j] - (lambda + 3 * sigma) * (u[905+j] ^ 2 + 2 * sigma * u[905+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[805+j] * (1- u[805+j]))
    end
    #Orden 11
    for j in 2:101
        du[1005+j] = (lambda + sigma - 1) * u[1005+j] - (lambda + 3 * sigma) * (u[1005+j] ^ 2 + 2 * sigma * u[1005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[905+j] * (1- u[905+j]))
    end
    #Orden 12
    for j in 2:101
        du[1105+j] = (lambda + sigma - 1) * u[1105+j] - (lambda + 3 * sigma) * (u[1105+j] ^ 2 + 2 * sigma * u[1105+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[1005+j] * (1- u[1005+j]))
    end
end
##########
function Prop_order12_BM1000(du,u,p,t)
    du[1] = (lambda + sigma - 1) * u[1] - (lambda + 3 * sigma) * (u[1] ^ 2 + 2 * sigma * u[1] ^ 3)
    #Orden 1
    for j in 2:1001
        du[j] = (lambda + sigma - 1) * u[j] - (lambda + 3 * sigma) * (u[j] ^ 2 + 2 * sigma * u[j] ^ 3) + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1] * (1- u[1])
    end
    ##Orden 2 mixed
    du[1002] = (lambda + sigma - 1) * u[1002] - (lambda + 3 * sigma) * (u[1002] ^ 2 + 2 * sigma * u[1002] ^ 3) + sqrt(2/T) * sigma * (sin(((2-1) * pi * t)/ T) * u[3] * (1- u[3]) + sin(((3-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[1003] = (lambda + sigma - 1) * u[1003] - (lambda + 3 * sigma) * (u[1003] ^ 2 + 2 * sigma * u[1003] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[1004] = (lambda + sigma - 1) * u[1004] - (lambda + 3 * sigma) * (u[1004] ^ 2 + 2 * sigma * u[1004] ^ 3) + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[3] * (1- u[3]))
    ##Orden 2 diagonal
    for j in 2:1001
        du[1003+j] = (lambda + sigma - 1) * u[1003+j] - (lambda + 3 * sigma) * (u[1003+j] ^ 2 + 2 * sigma * u[1003+j] ^ 3) + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j] * (1- u[j])
    end
    #Orden 3 mixed
    du[2005] = (lambda + sigma - 1) * u[2005] - (lambda + 3 * sigma) * (u[2005] ^ 2 + 2 * sigma * u[2005] ^ 3) + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] * (1- u[1006]) + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002] * (1- u[1002]))
    du[2006] = (lambda + sigma - 1) * u[2006] - (lambda + 3 * sigma) * (u[2006] ^ 2 + 2 * sigma * u[2006] ^ 3) + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] * (1- u[1002]) + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005] * (1- u[1005]))
    ##Orden 3 diagonal
    for j in 2:1001
        du[2005+j] = (lambda + sigma - 1) * u[2005+j] - (lambda + 3 * sigma) * (u[2005+j] ^ 2 + 2 * sigma * u[2005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j] * (1- u[1003+j]))
    end
    #Orden 4
    for j in 2:1001
        du[3005+j] = (lambda + sigma - 1) * u[3005+j] - (lambda + 3 * sigma) * (u[3005+j] ^ 2 + 2 * sigma * u[3005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2005+j] * (1- u[2005+j]))
    end
    #Orden 5
    for j in 2:1001
        du[4005+j] = (lambda + sigma - 1) * u[4005+j] - (lambda + 3 * sigma) * (u[4005+j] ^ 2 + 2 * sigma * u[4005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j] * (1- u[3005+j]))
    end
    #Orden 6
    for j in 2:1001
        du[5005+j] = (lambda + sigma - 1) * u[5005+j] - (lambda + 3 * sigma) * (u[5005+j] ^ 2 + 2 * sigma * u[5005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4005+j] * (1- u[4005+j]))
    end
    #Orden 7
    for j in 2:1001
        du[6005+j] = (lambda + sigma - 1) * u[6005+j] - (lambda + 3 * sigma) * (u[6005+j] ^ 2 + 2 * sigma * u[6005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[5005+j] * (1- u[5005+j]))
    end
    #Orden 8
    for j in 2:1001
        du[7005+j] = (lambda + sigma - 1) * u[7005+j] - (lambda + 3 * sigma) * (u[7005+j] ^ 2 + 2 * sigma * u[7005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[6005+j] * (1- u[6005+j]))
    end
    #Orden 9
    for j in 2:1001
        du[8005+j] = (lambda + sigma - 1) * u[8005+j] - (lambda + 3 * sigma) * (u[8005+j] ^ 2 + 2 * sigma * u[8005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[7005+j] * (1- u[7005+j]))
    end
    #Orden 10
    for j in 2:1001
        du[9005+j] = (lambda + sigma - 1) * u[9005+j] - (lambda + 3 * sigma) * (u[9005+j] ^ 2 + 2 * sigma * u[9005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[8005+j] * (1- u[8005+j]))
    end
    #Orden 11
    for j in 2:1001
        du[10005+j] = (lambda + sigma - 1) * u[10005+j] - (lambda + 3 * sigma) * (u[10005+j] ^ 2 + 2 * sigma * u[10005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[9005+j] * (1- u[9005+j]))
    end
    #Orden 12
    for j in 2:1001
        du[11005+j] = (lambda + sigma - 1) * u[11005+j] - (lambda + 3 * sigma) * (u[11005+j] ^ 2 + 2 * sigma * u[11005+j] ^ 3) + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[10005+j] * (1- u[10005+j]))
    end
end
