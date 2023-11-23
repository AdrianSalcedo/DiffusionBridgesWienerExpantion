using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV, IterableTables, DataTables
using DataFrames
using Random, Distributions

path0 = "Dropbox\\Artículos\\Trabajos_Francisco\\New_GeometricB\\"
path2 = "C:\\Users\\Usuario1\\"

path = path2 * path0
include(path * "Dynamic.jl")
alpha_0 = 0.2
sigma = 0.3
theta = 0.3
eta = 0.2
u_0 = eta
t0 = 0.0
T = 1.0
Time = (0.0,T)
dt = 0.001

u0 = [u_0,0.0]

################################################################################
################################################################################
################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################
j=1
Data = DataFrame()
step = 0.001
t_s = range(0.0,T, step=step)

for k in 1:100
    println(j)
    include(path * "Dynamic.jl")
    problem = ODEProblem(Prop_Dynamic,u0,Time)
    Solution = solve(problem,AutoTsit5(Rosenbrock23()),reltol=1e-8,saveat = step)
    df = DataFrame(Solution)  
    if (k==1)
        Data_aux =  DataFrame(t = df[:,1], x0 = df[:,2], x1 = df[:,3])
        Data = append!(Data, Data_aux)
    else 
        colname = "x$k"
        Data[!,colname] = df[:,3]
    end
    j+=1
end
CSV.write(path * "Propagator.csv",Data)

################################################################################
############################ PLot variables ####################################
D = Matrix(Data)
plot(D[:,1],D[:,2])
for i in 1:100
    Xi = W + sqrt(2/T) * sin((j * pi * t)/ T) * 
end



#######################################################################
#######################################################################
function BM(n,TiempoC)
      MBM = zeros(n,length(TiempoC))
      l=length(TiempoC)
      delta=TiempoC[2]-TiempoC[1]
      for i in 1:n
            for j in 2:l
                MBM[i,j] = MBM[i,j-1]+rand(Normal(0,sqrt(delta)))
            end
        end
      return(MBM)
end
#######################################################################
#######################################################################
Ym = DataFrame()
M_Y = DataFrame()
w = zeros(100)
vec = w
Xi = vec
for i in 1:100
    vec[i] = i
end

for k in 1:1 # Loop for trajectories
    Bt = BM(100,Data[:,1])
    for i in 2:size(Data)[1]
        w = w + sqrt(2)* sin.(vec*pi*Data[i-1,1]) * (Bt[:,i]- Bt[:,i-1])
        Xi = w # como son orden 1 no se ocupan más


        suma = Matrix(Data) - Matrix(Data)
        Ym = suma
        suma[1,2] = 0.0
        Ym[1,2] = eta
        for q in 2:size(Data)[1]  ########### Loop |m| =0
            suma[q,2] = suma[q-1,2] + (Data[q,2]-Data[q-1,2])/(T-Data[q-1,1])
            Ym[q,2] = eta + (theta-eta) * (Data[q,1]/T) + (T-Data[q,1])*suma[q-1,2]
        end
        for r in 3:102
            #Ym[1,r] = Matrix(Data[1,r])
            #suma[1,r] = Ym[1,r]
            for h in 2:size(Data)[1]
                suma[h,r] = suma[h-1,r] + (Data[h,r]-Data[h-1,r])/(T-Data[h-1,1])
                Ym[h,r] = (T-Data[h,1]) * suma[h,r]
            end
        end
        M_Y[t,k] = t(Ym[i,])*Xi
    end
    
end