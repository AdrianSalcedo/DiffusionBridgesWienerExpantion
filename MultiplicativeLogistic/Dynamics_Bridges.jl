
F_m0_Drift(u,p,t) = alpha_0 * u

function F_m1_Drift(du,u,p,t)
 du[1] = alpha_0 * u[1] + ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[6]
 du[2] = alpha_0 * u[2] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[6]
 du[3] = alpha_0 * u[3] + ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[6]
 du[4] = alpha_0 * u[4] + ((2) ^ (1 / 2)) * sin(4 * pi * t) * sigma * u[6]
 du[5] = alpha_0 * u[5] + ((2) ^ (1 / 2)) * sin(5 * pi * t) * sigma * u[6]
 # u[6] is |m|=0, recalculate
 du[6] = alpha_0 * u[6]
end


function F_m2_Drift(du,u,p,t)
    # |m|=1
 du[1] = alpha_0 * u[1] + ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[12]
 du[2] = alpha_0 * u[2] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[12]
 du[3] = alpha_0 * u[3] + ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[12]
 du[4] = alpha_0 * u[4] + ((2) ^ (1 / 2)) * sin(4 * pi * t) * sigma * u[12]
 du[5] = alpha_0 * u[5] + ((2) ^ (1 / 2)) * sin(5 * pi * t) * sigma * u[12]
 # |m|=2
 du[6] = alpha_0 * u[6] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[2] +
        ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[1]
 du[7] = alpha_0 * u[7] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[1]
 du[8] = alpha_0 * u[8] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[2]
 du[9] = alpha_0 * u[9] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[1]
 du[10] = alpha_0 * u[10] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[2]
 du[11] = alpha_0 * u[11] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[3]
 # u[12] is |m|=0, recalculate
 du[12] = alpha_0 * u[12]
end

function F_m3_Drift(du,u,p,t)
    # |m|=1
 du[1] = alpha_0 * u[1] + ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[16]
 du[2] = alpha_0 * u[2] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[16]
 du[3] = alpha_0 * u[3] + ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[16]
 du[4] = alpha_0 * u[4] + ((2) ^ (1 / 2)) * sin(4 * pi * t) * sigma * u[16]
 du[5] = alpha_0 * u[5] + ((2) ^ (1 / 2)) * sin(5 * pi * t) * sigma * u[16]
 # |m|=2
 du[6] = alpha_0 * u[6] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[2] +
        ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[1]
 du[7] = alpha_0 * u[7] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[1]
 du[8] = alpha_0 * u[8] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[2]
 du[9] = alpha_0 * u[9] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[1]
 du[10] = alpha_0 * u[10] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[2]
 du[11] = alpha_0 * u[11] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[3]
 # |m| = 3
 du[12] = alpha_0 * u[12] +
    ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[10] +
        ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[6]
 du[13] = alpha_0 * u[13] +
    ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[6] +
        ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[9]
 du[14] = alpha_0 * u[14] +
    ((3) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[9]
 du[15] = alpha_0 * u[15] +
    ((3) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[10]
 # u[16] is |m|=0, recalculate
 du[16] = alpha_0 * u[16]
end

function F_m4_Drift(du,u,p,t)
    # |m|=1
 du[1] = alpha_0 * u[1] + ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[18]
 du[2] = alpha_0 * u[2] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[18]
 du[3] = alpha_0 * u[3] + ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[18]
 du[4] = alpha_0 * u[4] + ((2) ^ (1 / 2)) * sin(4 * pi * t) * sigma * u[18]
 du[5] = alpha_0 * u[5] + ((2) ^ (1 / 2)) * sin(5 * pi * t) * sigma * u[18]
 # |m|=2
 du[6] = alpha_0 * u[6] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[2] +
        ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[1]
 du[7] = alpha_0 * u[7] + ((2) ^ (1 / 2)) * sin( pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[1]
 du[8] = alpha_0 * u[8] + ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[3] +
        ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[2]
 du[9] = alpha_0 * u[9] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[1]
 du[10] = alpha_0 * u[10] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[2]
 du[11] = alpha_0 * u[11] + ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(3 * pi * t) * sigma * u[3]
 # |m| = 3
 du[12] = alpha_0 * u[12] +
    ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[10] +
        ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[6]
 du[13] = alpha_0 * u[13] +
    ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[6] +
        ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[9]
 du[14] = alpha_0 * u[14] +
    ((3) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[9]
 du[15] = alpha_0 * u[15] +
    ((3) ^ (1 / 2)) * ((2) ^ (1 / 2)) * ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[10]
 du[16] = alpha_0 * u[16] + 2 *  ((2) ^ (1 / 2)) * sin(pi * t) * sigma * u[14]
 du[17] = alpha_0 * u[17] + 2 *  ((2) ^ (1 / 2)) * sin(2 * pi * t) * sigma * u[15]

 # u[16] is |m|=0, recalculate
 du[18] = alpha_0 * u[18]
end
################################################################################
