using DifferentialEquations
using Plots; plotly()
#using DifferentialEquations.EnsembleAnalysis
using CSV, IterableTables, DataTables
using DataFrames
using Random, Distributions
#D:\\DiffusionBridgesWienerExpantion\\DiffusionBridgesWienerExpantion\GeometricBrownian
#path0 = "Dropbox\\Artículos\\Trabajos_Francisco\\New_GeometricB\\"
#path2 = "C:\\Users\\Usuario1\\"
path0 = "DiffusionBridgesWienerExpantion\\GeometricBrownian\\"
path2 = "D:\\DiffusionBridgesWienerExpantion\\"
path = path2 * path0
pathd = "D:\\Propagador\\Geometric\\"
include(path * "Dynamic.jl")
sigma = 0.3
alpha_0 = 0.2
theta = 0.3

eta = 0.2
u_0 = eta
t0 = 0.0
T = 1.0
Time = (0.0,T)
dt = 0.001

u0 = zeros(8006)
u0[1]=u_0
step = dt
t_s = range(0.0,T, step=step)

################################################################################
problem = ODEProblem(Prop_order8_BM1000,u0,Time)
Solution = solve(problem,AutoTsit5(Rosenbrock23()),saveat = step)
df = DataFrame(Solution)
Data = Matrix(df)
CSV.write(path * "Propagator_order8_int0.2_alpha0.2_sigma0.3_MB1000.csv",df)

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
function BM(n,TiempoC)
      MBM = zeros(n,length(TiempoC))
      l=length(TiempoC)
      delta=TiempoC[2]-TiempoC[1]
      for i in 1:n
            for j in 2:l
                MBM[i,j] = MBM[i,j-1]+rand(Normal(0,1))#sqrt(delta)
            end
        end
      return(MBM)
end
#######################################################################
#######################################################################
M = 5
Ym = DataFrame()
M_Y = zeros(size(Data)[1],M)
w = zeros(1000)
vec = w
Xi = zeros(1017)


for i in 1:1000
    vec[i] = i
end

for k in 1:M # Loop for trajectories

    Bt = BM(1000,Data[:,1])
    w = zeros(1000)
    for i in 2:size(Data)[1]
        w = w + sqrt(2)* sin.(vec*pi*Data[i-1,1]) .* (Bt[:,i]- Bt[:,i-1])
        Xi[1] = 1
        Xi[2:1001] = w
        Xi[1002]= (w[1]) * (w[2])
        Xi[1003]= (w[1]) * (w[3])
        Xi[1004]= (w[2]) * (w[3])
        Xi[1005]= sqrt(2)*((w[1]) ^ 2 - 1) / sqrt(2)
        Xi[1006]= sqrt(2)*((w[2]) ^ 2 - 1) / sqrt(2)
        Xi[1007]= sqrt(2)*((w[3]) ^ 2 - 1) / sqrt(2)
        Xi[1008]= sqrt(2)*(w[1]) * ((w[2]) ^ 2 - 1) / sqrt(2)
        Xi[1009]= sqrt(2)*(((w[1]) ^ 2 - 1) * (w[2]) / sqrt(2))
        Xi[1010]= sqrt(6)*(w[1]) * (((w[1]) ^ 2 - 3) / sqrt(6))
        Xi[1011]= sqrt(6)*(w[2]) * (((w[2]) ^ 2 - 3) / sqrt(6))
        Xi[1012]= sqrt(24)* (w[1] ^ 4 - 6 * w[1] ^ 2 + 3) / (2 * sqrt(6))
        Xi[1013]= sqrt(24)* (w[2] ^ 4 - 6 * w[2] ^ 2 + 3) / (2 * sqrt(6))
        Xi[1014]= sqrt(120)*(w[1]*(15 - 10*w[1]^2 + w[1]^4))/(2*sqrt(30))
        Xi[1015]= sqrt(120)*(w[2]*(15 - 10*w[2]^2 + w[2]^4))/(2*sqrt(30))
        Xi[1016]= sqrt(720)*(-15 + 45*w[1]^2 - 15*w[1]^4 + w[1]^6)/(12*sqrt(5))
        Xi[1017]= sqrt(720)*(-15 + 45*w[2]^2 - 15*w[2]^4 + w[2]^6)/(12*sqrt(5))


        suma = Matrix(Data) - Matrix(Data)
        Ym = suma
        suma[1,2] = 0.0
        Ym[1,2] = eta
        for q in 2:size(Data)[1]  ########### Loop |m| =0
            suma[q,2] = suma[q-1,2] + (Data[q,2]-Data[q-1,2])/(T-Data[q-1,1])
            Ym[q,2] = eta + (theta-eta) * (Data[q,1]/T) + (T-Data[q,1])*suma[q-1,2]
        end
        for r in 3:118
            #Ym[1,r] = Matrix(Data[1,r])
            #suma[1,r] = Ym[1,r]
            for h in 2:size(Data)[1]
                suma[h,r] = suma[h-1,r] + (Data[h,r]-Data[h-1,r])/(T-Data[h-1,1])
                Ym[h,r] = (T-Data[h,1]) * suma[h,r]
            end
        end
        M_Y[i,k] = transpose(Ym[i,2:end])* Xi
    end
    M_Y[1,k] = transpose(Ym[1,2:end])*Xi
    println(k)
end
plot(Data[:,1],M_Y)
CSV.write(path * "Bridge_1000MB.csv",DataFrame(M_Y,:auto))