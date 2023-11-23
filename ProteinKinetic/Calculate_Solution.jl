using Plots; plotly()
using CSV, DifferentialEquations
using IterableTables, DataTables
using DataFrames
using Statistics, Distributions, Random
path0 = "Dropbox/Artículos/Trabajos_Francisco/GeometricBrownianMotion/"
path1 = "/home/gabrielsalcedo/"
path2 = "C:/Users/adria/"
path = path1 * path0

Time = CSV.read(path * "Time_Brownian_motion.csv", DataFrame)
Propagator = CSV.read(path * "Propagator_Brownian_motion.csv", DataFrame)
Time = Matrix(Time)
Propagator = Matrix(Propagator)
m, n = size(Propagator)

function dW(number_max_of_steps)
    #input: step size, and number of iterations
    #output: Brownian motion
    #Random sample normal distribution
    Normal_sampler =  randn(number_max_of_steps - 1)
    W_t = zeros(number_max_of_steps)
    W_t[2:end] = cumsum(Normal_sampler)*sqrt(1/number_max_of_steps)
   return W_t
end

theta = 0.99
eta = 0.113
M = 1
Mx = zeros(m,M)
My = zeros(m,M)
Xz = zeros(m)
Vector_xi = zeros(m)
suma = Ym = Propagator - Propagator
X = 0
#w1 = w2 = w3 = w4 = w5 = 0
for k in 1:M
    Bt = dW(m)
    plot(Bt)
    sleep(1.0)
    Mx[1,k] = 0.1
    w1 = w2 = w3 = w4 = w5 = 0
    for t in 2:m
        Z = Bt[t]
        Zr = Bt[t-1]
        norm1 = 1#(1 - sin(2*pi)/(2*pi))
        norm2 = 1#(1 - sin(2*pi)/(2*pi))
        norm3 = 1#(1 - sin(2*pi)/(2*pi))
        norm4 = 1#(1 - sin(2*pi)/(2*pi))
        norm5 = 1#(1 - sin(2*pi)/(2*pi))
        w1 = w1 + sqrt(2)* sin(1*pi*Time[t-1])*(Z-Zr)
        w2 = w2 + sqrt(2)* sin(2*pi*Time[t-1])*(Z-Zr)
        w3 = w3 + sqrt(2)* sin(3*pi*Time[t-1])*(Z-Zr)
        w4 = w4 + sqrt(2)* sin(4*pi*Time[t-1])*(Z-Zr)
        w5 = w5 + sqrt(2)* sin(5*pi*Time[t-1])*(Z-Zr)
        W = [w1,w2,w3,w4,w5]

        xi1 = W[1]
        xi2 = W[2]
        xi3 = W[3]
        xi4 = W[4]
        xi5 = W[5]
        xi6 = (W[1]) * (W[2])
        xi7 = (W[1]) * (W[3])
        xi8 = (W[2]) * (W[3])
        xi9 = ((W[1]) ^ 2 - 1) / sqrt(2)
        xi10 = ((W[2]) ^ 2 - 1) / sqrt(2)
        xi11 = ((W[3]) ^ 2 - 1) / sqrt(2)
        xi12 = (W[1]) * ((W[2]) ^ 2 - 1) / sqrt(2)
        xi13 = (((W[1]) ^ 2 - 1)  * (W[2]) / sqrt(2))
        xi14 = (W[1]) * (((W[1]) ^ 2 - 3) / sqrt(6))
        xi15 = (W[2]) * (((W[2]) ^ 2 - 3) / sqrt(6))
        xi16 = ((W[1]) ^ 4 - 6 *(W[1]) ^ 2 + 3) / (2 * sqrt(6))
        xi17 = ((W[2]) ^ 4 - 6 *W[2] ^ 2 + 3) / (2 * sqrt(6))
        xi18  = 1

        Vector_xi = [xi1,xi2,xi3,xi4,xi5, xi6,xi7,xi8,xi9,xi10,xi11,xi12,
            xi13,xi14,xi15,xi16, xi17, xi18]

        X = sum(Propagator[t,:] .* Vector_xi)
        Mx[t,k] = X
        #Xz = mean(Mx,dims= 2)
        #plot(Time,Xz)
        suma = Ym = Propagator - Propagator

        for h in 1:17 #### for |m|>0
            suma[1,h] = Ym[1,h] = Propagator[1,h]
            for i in 2:m
                suma[i,h] = suma[i-1,h] +
                    (Propagator[i,h]-Propagator[i-1,h])/(1-Time[i-1])
                Ym[i,h] = (1-Time[i]) * suma[i,h]
            end
        end
        #### for |m|=0
        suma[1,18] = 0#Propagator[1,18]- Propagator[1,18]
        Ym[1,18]= eta
        for i in 2:m
            suma[i,18]= suma[i-1,18] +
                (Propagator[i,18] - Propagator[i-1,18]) / (1-Time[i-1])
            Ym[i,18] = eta + (theta-eta ) * (Time[i]/1) +
                (1-Time[i]) * suma[i-1,18]
        end
         My[t,k] = sum(Ym[t,:] .* Vector_xi)
    end # for t
    My[1,k] = sum(Ym[1,:] .* Vector_xi)
    println(k)

end # for k

Yz = mean(My,dims=2)
plot(Time, My)
My_datos = DataFrame(My,:auto)
CSV.write(path * "Solution_Means.csv", My_datos)
#CSV.write(path * "Solution_plot.csv", X_Solution)
