function Prop_Dynamic(du,u,p,t)
    du[1] = alpha_0 * u[1]
    for j in 2:101
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
end

function Prop_order6(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    for j in 2:4
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[106] = alpha_0 * u[106] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[107] = alpha_0 * u[107] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3
    du[1008] = alpha_0 * u[1008] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[1009] = alpha_0 * u[1009] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    du[1010] = alpha_0 * u[1010] + sigma * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[1005])
    du[1011] = alpha_0 * u[1011] + sigma * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[1006])
    #Orden 4
    du[1012] = alpha_0 * u[1012] + sigma * sqrt(2/T) * (sqrt(4) * sin(((2-1) * pi * t)/ T) * u[1010])
    du[1013] = alpha_0 * u[1013] + sigma * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[1011])
    #Orden 5
    du[1014] = alpha_0 * u[1014] + sigma * sqrt(2/T) * (sqrt(5) * sin(((2-1) * pi * t)/ T) * u[1012])
    du[1015] = alpha_0 * u[1015] + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013])
    #Orden 6
    du[1016] = alpha_0 * u[1016] + sigma * sqrt(2/T) * (sqrt(6) * sin(((2-1) * pi * t)/ T) * u[1014])
    du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
end

function Prop_order6_1000(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:1001
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[106] = alpha_0 * u[106] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[107] = alpha_0 * u[107] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[2005] = alpha_0 * u[2005] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[2006] = alpha_0 * u[2006] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    ##Orden 3 diagonal
    for j in 2:1001
        du[2005+j] = alpha_0 * u[2005+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j])
        #du[1011] = alpha_0 * u[1011] + sigma * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[1006])
    end
    #Orden 4
    for j in 2:1001
        du[3005+j] = alpha_0 * u[3005+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2004+j])
        #du[1013] = alpha_0 * u[1013] + sigma * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[1011])
    end
    #Orden 5
    for j in 2:1001
     du[4005+j] = alpha_0 * u[4005+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j])
        #du[1015] = alpha_0 * u[1015] + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013])
    end
    #Orden 6
    for j in 2:1001
        du[5005+j] = alpha_0 * u[5005+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4006+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
end
#######################################################################################
########################## Propagador Orden 12 only 1000 BM in order 1
function Prop_order12(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    for j in 2:4
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[106] = alpha_0 * u[106] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    #du[107] = alpha_0 * u[107] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3
    du[1008] = alpha_0 * u[1008] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[1009] = alpha_0 * u[1009] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    du[1010] = alpha_0 * u[1010] + sigma * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[1005])
    du[1011] = alpha_0 * u[1011] + sigma * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[1006])
    #Orden 4
    du[1012] = alpha_0 * u[1012] + sigma * sqrt(2/T) * (sqrt(4) * sin(((2-1) * pi * t)/ T) * u[1010])
    du[1013] = alpha_0 * u[1013] + sigma * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[1011])
    #Orden 5
    du[1014] = alpha_0 * u[1014] + sigma * sqrt(2/T) * (sqrt(5) * sin(((2-1) * pi * t)/ T) * u[1012])
    du[1015] = alpha_0 * u[1015] + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013])
    #Orden 6
    du[1016] = alpha_0 * u[1016] + sigma * sqrt(2/T) * (sqrt(6) * sin(((2-1) * pi * t)/ T) * u[1014])
    du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    #Orden 7
    du[1018] = alpha_0 * u[1018] + sigma * sqrt(2/T) * (sqrt(7) * sin(((2-1) * pi * t)/ T) * u[1016])
    du[1019] = alpha_0 * u[1019] + sigma * sqrt(2/T) * (sqrt(7) * sin(((3-1) * pi * t)/ T) * u[1017])
    #Orden 8
    du[1020] = alpha_0 * u[1020] + sigma * sqrt(2/T) * (sqrt(8) * sin(((2-1) * pi * t)/ T) * u[1018])
    du[1021] = alpha_0 * u[1021] + sigma * sqrt(2/T) * (sqrt(8) * sin(((3-1) * pi * t)/ T) * u[1019])
    #Orden 9
    du[1022] = alpha_0 * u[1022] + sigma * sqrt(2/T) * (sqrt(9) * sin(((2-1) * pi * t)/ T) * u[1020])
    du[1023] = alpha_0 * u[1023] + sigma * sqrt(2/T) * (sqrt(9) * sin(((3-1) * pi * t)/ T) * u[1021])
    #Orden 10
    du[1024] = alpha_0 * u[1024] + sigma * sqrt(2/T) * (sqrt(10) * sin(((2-1) * pi * t)/ T) * u[1022])
    du[1025] = alpha_0 * u[1025] + sigma * sqrt(2/T) * (sqrt(10) * sin(((3-1) * pi * t)/ T) * u[1023])
    #Orden 11
    du[1026] = alpha_0 * u[1026] + sigma * sqrt(2/T) * (sqrt(11) * sin(((2-1) * pi * t)/ T) * u[1024])
    du[1027] = alpha_0 * u[1027] + sigma * sqrt(2/T) * (sqrt(11) * sin(((3-1) * pi * t)/ T) * u[1025])
    #Orden 12
    du[1028] = alpha_0 * u[1028] + sigma * sqrt(2/T) * (sqrt(12) * sin(((2-1) * pi * t)/ T) * u[1026])
    du[1029] = alpha_0 * u[1029] + sigma * sqrt(2/T) * (sqrt(12) * sin(((3-1) * pi * t)/ T) * u[1027])
end
##################################################################################################
#  First code order 12,  BM 1000 each order
##################################################################################################
# function Prop_order12_1000(du,u,p,t)
#     du[1] = alpha_0 * u[1]
#     #Orden 1
#     for j in 2:1001
#         du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
#     end
#     ##Orden 2 mixed
#     du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
#     + sin(((3-1) * pi * t)/ T) * u[2])
#     du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
#     + sin(((4-1) * pi * t)/ T) * u[2])
#     du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
#     + sin(((4-1) * pi * t)/ T) * u[3])
#     ##Orden 2 diagonal
#     for j in 2:1001
#         du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
#     #du[106] = alpha_0 * u[106] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
#     #du[107] = alpha_0 * u[107] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
#     end
#     #Orden 3 mixed
#     du[2005] = alpha_0 * u[2005] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
#     + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
#     du[2006] = alpha_0 * u[2006] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
#     + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
#     ##Orden 3 diagonal
#     for j in 2:1001
#         du[2005+j] = alpha_0 * u[2005+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j])
#         #du[1011] = alpha_0 * u[1011] + sigma * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[1006])
#     end
#     #Orden 4
#     for j in 2:1001
#         du[3005+j] = alpha_0 * u[3005+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2005+j])
#         #du[1013] = alpha_0 * u[1013] + sigma * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[1011])
#     end
#     #Orden 5
#     for j in 2:1001
#      du[4005+j] = alpha_0 * u[4005+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j])
#         #du[1015] = alpha_0 * u[1015] + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013])
#     end
#     #Orden 6
#     for j in 2:1001
#         du[5005+j] = alpha_0 * u[5005+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 7
#     for j in 2:1001
#         du[6005+j] = alpha_0 * u[6005+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[5005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 8
#     for j in 2:1001
#         du[7005+j] = alpha_0 * u[7005+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[6005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 9
#     for j in 2:1001
#         du[8005+j] = alpha_0 * u[8005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[7005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 10
#     for j in 2:1001
#         du[9005+j] = alpha_0 * u[9005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[8005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 11
#     for j in 2:1001
#         du[10005+j] = alpha_0 * u[10005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[9005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
#     #Orden 12
#     for j in 2:1001
#         du[11005+j] = alpha_0 * u[10005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[10005+j])
#         #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
#     end
# end
##########################################
function Prop_order12_BM0(du,u,p,t)
    du[1] = alpha_0 * u[1]
end

function Prop_order12_BM10(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:11
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[12] = alpha_0 * u[12] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[13] = alpha_0 * u[13] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[14] = alpha_0 * u[14] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:11
        du[13+j] = alpha_0 * u[13+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[25] = alpha_0 * u[25] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[16] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[12])
    du[26] = alpha_0 * u[26] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[12] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[15])
    ##Orden 3 diagonal
    for j in 2:11
        du[25+j] = alpha_0 * u[25+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[13+j])
    end
    #Orden 4
    for j in 2:11
        du[35+j] = alpha_0 * u[35+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[25+j])
    end
    #Orden 5
    for j in 2:11
     du[45+j] = alpha_0 * u[45+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[35+j])
    end
    #Orden 6
    for j in 2:11
        du[55+j] = alpha_0 * u[55+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[45+j])
    end
    #Orden 7
    for j in 2:11
        du[65+j] = alpha_0 * u[65+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[55+j])
    end
    #Orden 8
    for j in 2:11
        du[75+j] = alpha_0 * u[75+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[65+j])
    end
    #Orden 9
    for j in 2:11
        du[85+j] = alpha_0 * u[85+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[75+j])
    end
    #Orden 10
    for j in 2:11
        du[95+j] = alpha_0 * u[95+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[85+j])
    end
    #Orden 11
    for j in 2:11
        du[105+j] = alpha_0 * u[105+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[95+j])
    end
    #Orden 12
    for j in 2:11
        du[115+j] = alpha_0 * u[115+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[105+j])
    end
end

function Prop_order12_BM100(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:101
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[102] = alpha_0 * u[102] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[103] = alpha_0 * u[103] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[104] = alpha_0 * u[104] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:101
        du[103+j] = alpha_0 * u[103+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[205] = alpha_0 * u[205] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[106] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[102])
    du[206] = alpha_0 * u[206] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[102] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[105])
    ##Orden 3 diagonal
    for j in 2:101
        du[205+j] = alpha_0 * u[205+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[103+j])
    end
    #Orden 4
    for j in 2:101
        du[305+j] = alpha_0 * u[305+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[205+j])
    end
    #Orden 5
    for j in 2:101
     du[405+j] = alpha_0 * u[405+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[305+j])
        #du[1015] = alpha_0 * u[1015] + sigma * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[1013])
    end
    #Orden 6
    for j in 2:101
        du[505+j] = alpha_0 * u[505+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[405+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 7
    for j in 2:101
        du[605+j] = alpha_0 * u[605+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[505+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 8
    for j in 2:101
        du[705+j] = alpha_0 * u[705+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[605+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 9
    for j in 2:101
        du[805+j] = alpha_0 * u[805+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[705+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 10
    for j in 2:101
        du[905+j] = alpha_0 * u[905+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[805+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 11
    for j in 2:101
        du[1005+j] = alpha_0 * u[1005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[905+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
    #Orden 12
    for j in 2:101
        du[1105+j] = alpha_0 * u[1105+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[1005+j])
        #du[1017] = alpha_0 * u[1017] + sigma * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[1015])
    end
end

function Prop_order12_BM1000(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:1001
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[2005] = alpha_0 * u[2005] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[2006] = alpha_0 * u[2006] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    ##Orden 3 diagonal
    for j in 2:1001
        du[2005+j] = alpha_0 * u[2005+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j])
    end
    #Orden 4
    for j in 2:1001
        du[3005+j] = alpha_0 * u[3005+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2005+j])
    end
    #Orden 5
    for j in 2:1001
     du[4005+j] = alpha_0 * u[4005+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j])
    end
    #Orden 6
    for j in 2:1001
        du[5005+j] = alpha_0 * u[5005+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4005+j])
    end
    #Orden 7
    for j in 2:1001
        du[6005+j] = alpha_0 * u[6005+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[5005+j])
    end
    #Orden 8
    for j in 2:1001
        du[7005+j] = alpha_0 * u[7005+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[6005+j])
    end
    #Orden 9
    for j in 2:1001
        du[8005+j] = alpha_0 * u[8005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[7005+j])
    end
    #Orden 10
    for j in 2:1001
        du[9005+j] = alpha_0 * u[9005+j] + sigma * sqrt(2/T) * (sqrt(10) * sin(((j-1) * pi * t)/ T) * u[8005+j])
    end
    #Orden 11
    for j in 2:1001
        du[10005+j] = alpha_0 * u[10005+j] + sigma * sqrt(2/T) * (sqrt(11) * sin(((j-1) * pi * t)/ T) * u[9005+j])
    end
    #Orden 12
    for j in 2:1001
        du[11005+j] = alpha_0 * u[11005+j] + sigma * sqrt(2/T) * (sqrt(12) * sin(((j-1) * pi * t)/ T) * u[10005+j])
    end
end

function Prop_order15_BM1000(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:1001
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[2005] = alpha_0 * u[2005] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[2006] = alpha_0 * u[2006] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    ##Orden 3 diagonal
    for j in 2:1001
        du[2005+j] = alpha_0 * u[2005+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j])
    end
    #Orden 4
    for j in 2:1001
        du[3005+j] = alpha_0 * u[3005+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2005+j])
    end
    #Orden 5
    for j in 2:1001
     du[4005+j] = alpha_0 * u[4005+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j])
    end
    #Orden 6
    for j in 2:1001
        du[5005+j] = alpha_0 * u[5005+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4005+j])
    end
    #Orden 7
    for j in 2:1001
        du[6005+j] = alpha_0 * u[6005+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[5005+j])
    end
    #Orden 8
    for j in 2:1001
        du[7005+j] = alpha_0 * u[7005+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[6005+j])
    end
    #Orden 9
    for j in 2:1001
        du[8005+j] = alpha_0 * u[8005+j] + sigma * sqrt(2/T) * (sqrt(9) * sin(((j-1) * pi * t)/ T) * u[7005+j])
    end
    #Orden 10
    for j in 2:1001
        du[9005+j] = alpha_0 * u[9005+j] + sigma * sqrt(2/T) * (sqrt(10) * sin(((j-1) * pi * t)/ T) * u[8005+j])
    end
    #Orden 11
    for j in 2:1001
        du[10005+j] = alpha_0 * u[10005+j] + sigma * sqrt(2/T) * (sqrt(11) * sin(((j-1) * pi * t)/ T) * u[9005+j])
    end
    #Orden 12
    for j in 2:1001
        du[11005+j] = alpha_0 * u[11005+j] + sigma * sqrt(2/T) * (sqrt(12) * sin(((j-1) * pi * t)/ T) * u[10005+j])
    end
    #Orden 13
    for j in 2:1001
        du[12005+j] = alpha_0 * u[12005+j] + sigma * sqrt(2/T) * (sqrt(13) * sin(((j-1) * pi * t)/ T) * u[11005+j])
    end
    #Orden 14
    for j in 2:1001
        du[13005+j] = alpha_0 * u[13005+j] + sigma * sqrt(2/T) * (sqrt(14) * sin(((j-1) * pi * t)/ T) * u[12005+j])
    end
    #Orden 15
    for j in 2:1001
        du[14005+j] = alpha_0 * u[14005+j] + sigma * sqrt(2/T) * (sqrt(15) * sin(((j-1) * pi * t)/ T) * u[13005+j])
    end
end


function Prop_order8_BM1000(du,u,p,t)
    du[1] = alpha_0 * u[1]
    #Orden 1
    for j in 2:1001
        du[j] = alpha_0 * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * sigma * u[1]
    end
    ##Orden 2 mixed
    du[1002] = alpha_0 * u[1002] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[3] 
    + sin(((3-1) * pi * t)/ T) * u[2])
    du[1003] = alpha_0 * u[1003] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[2])
    du[1004] = alpha_0 * u[1004] + sigma * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] 
    + sin(((4-1) * pi * t)/ T) * u[3])
    ##Orden 2 diagonal
    for j in 2:1001
        du[1003+j] = alpha_0 * u[1003+j] + sigma * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[2005] = alpha_0 * u[2005] + sigma * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[1006] 
    + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[1002])
    du[2006] = alpha_0 * u[2006] + sigma * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[1002] 
    + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[1005])
    ##Orden 3 diagonal
    for j in 2:1001
        du[2005+j] = alpha_0 * u[2005+j] + sigma * sqrt(2/T) * (sqrt(3) * sin(((j-1) * pi * t)/ T) * u[1003+j])
    end
    #Orden 4
    for j in 2:1001
        du[3005+j] = alpha_0 * u[3005+j] + sigma * sqrt(2/T) * (sqrt(4) * sin(((j-1) * pi * t)/ T) * u[2005+j])
    end
    #Orden 5
    for j in 2:1001
     du[4005+j] = alpha_0 * u[4005+j] + sigma * sqrt(2/T) * (sqrt(5) * sin(((j-1) * pi * t)/ T) * u[3005+j])
    end
    #Orden 6
    for j in 2:1001
        du[5005+j] = alpha_0 * u[5005+j] + sigma * sqrt(2/T) * (sqrt(6) * sin(((j-1) * pi * t)/ T) * u[4005+j])
    end
    #Orden 7
    for j in 2:1001
        du[6005+j] = alpha_0 * u[6005+j] + sigma * sqrt(2/T) * (sqrt(7) * sin(((j-1) * pi * t)/ T) * u[5005+j])
    end
    #Orden 8
    for j in 2:1001
        du[7005+j] = alpha_0 * u[7005+j] + sigma * sqrt(2/T) * (sqrt(8) * sin(((j-1) * pi * t)/ T) * u[6005+j])
    end
end













