
################################################################################
#|m|=0
Sol_Exc1 = eta * exp.(alpha_0*t_s)
#|m|=1
Sol_Exc2 = - (sqrt(2) * sigma * eta *exp.(alpha_0 * t_s) .* (cos.(pi * t_s).-1))/pi
Sol_Exc3 = - (eta * sigma * exp.(alpha_0 * t_s).*(cos.(2*pi*t_s).-1))./(sqrt(2)*pi)
#|m|=2 mixed
Sol_Exc1002 = -(exp.(alpha_0*t_s)*eta*sigma.*(.-3*sqrt(2)*pi.-8*sigma.+9*sigma*cos.(pi*t_s)+3*sqrt(2)*pi*cos.(2*pi*t_s)
.-sigma.*cos.(3*pi*t_s)))/(6*pi^2)
#|m|=2 diag
Sol_Exc1005 = -(exp.(alpha_0 * t_s)*eta*sigma^2 .*(.- 2 .+ cos.(pi*t_s).+2*cos.(2*pi*t_s).-cos.(3*pi*t_s)))/(sqrt(2)*pi^2)
#|m|=3 mixed
Sol_Exc2006 = -exp.(alpha_0*t_s)*eta*sigma^2 .*(.-160*sqrt(2)*pi .-187*sigma .+180*sqrt(2)*pi*cos.(pi*t_s) .+ 200*sigma*cos.(pi*t_s) 
.+20*sigma*cos.(2*pi*t_s) 
.-20*sqrt(2)*pi*cos.(3*pi*t_s) .-20*sigma*cos.(3*pi*t_s) .-25*sigma*cos.(4*pi*t_s) .+ 12*sigma*cos.(5*pi*t_s)
)/(120*pi^3)
#|m|=3 diag
Sol_Exc2007 = -exp.(alpha_0*t_s)*eta*sigma^3 .*(.-55 .+72*cos.(pi*t_s) .-12*cos.(2*pi*t_s) .-8*cos.(3*pi*t_s) .+3*cos.(4*pi*t_s))/(8*sqrt(3)*pi^3) 
#|m|=4 diag
Sol_Exc3007 = -exp.(alpha_0*t_s)*eta*sigma^4 .*(.-322 .+490*cos.(pi*t_s) .-200*cos.(2*pi*t_s) .+25*cos.(3*pi*t_s) .+10*cos.(4*pi*t_s) 
.-3*cos.(5*pi*t_s))/(20*sqrt(6)*pi^4)
#|m|=5 diag
Sol_Exc4007 = -exp.(alpha_0*t_s)*eta*sigma^5 .*(.-546 .+888*cos.(pi*t_s) .-465*cos.(2*pi*t_s) .+140*cos.(3*pi*t_s) .-14*cos.(4*pi*t_s) 
.-4*cos.(5*pi*t_s) .+cos.(6*pi*t_s))/(16*sqrt(15)*pi^5)
#|m|=6 diag
Sol_Exc5007 = -exp.(alpha_0*t_s)*eta*sigma^6 .*(.-7788 .+13167*cos.(pi*t_s) .-7854*cos.(2*pi*t_s) .+3157*cos.(3*pi*t_s) .-756*cos.(4*pi*t_s) 
.+63*cos.(5*pi*t_s) .+14*cos.(6*pi*t_s) .-3*cos.(7*pi*t_s))/(336*sqrt(5)*pi^6)
#|m|=7 diag
Sol_Exc6007 = -exp.(alpha_0*t_s)*eta*sigma^7 .*(.-35607 .+61776*cos.(pi*t_s) .-40040*cos.(2*pi*t_s) .+18928*cos.(3*pi*t_s) .-6188*cos.(4*pi*t_s) 
.+1232*cos.(5*pi*t_s) .-88*cos.(6*pi*t_s) .-16*cos.(7*pi*t_s).+3*cos.(8*pi*t_s))/(384*sqrt(70)*pi^7)
#|m|=8 diag
Sol_Exc7007 = -exp.(alpha_0*t_s)*eta*sigma^8 .*(.-52910 .+93522*cos.(pi*t_s) .-64272*cos.(2*pi*t_s) .+33852*cos.(3*pi*t_s) .-13272*cos.(4*pi*t_s) 
.+3660*cos.(5*pi*t_s) .-624*cos.(6*pi*t_s) .+39*cos.(7*pi*t_s).+6*cos.(8*pi*t_s) .-cos.(9*pi*t_s))/(576*sqrt(70)*pi^8)
#|m|=9 diag
Sol_Exc8007 = (1/(11520*sqrt(35)*pi^9)).*(exp.(alpha_0*t_s)*eta*sigma^9).*(
    695266 .-1246440*cos.(pi*t_s) .+895050*cos.(2*pi*t_s) .-510000*cos.(3*pi*t_s) .+226440*cos.(4*pi*t_s) 
    .-75888*cos.(5*pi*t_s) .+18105*cos.(6*pi*t_s) .-2700*cos.(7*pi*t_s) .+150*cos.(8*pi*t_s) .+20*cos.(9*pi*t_s) .-3*cos.(10*pi*t_s))
#|m|=10 diag imposible de escribir muy largo
#|m|=11 diag
Sol_Exc10007 = (1/(46080*sqrt(154)*pi^11)).*(exp.(alpha_0*t_s)*eta*sigma^11).*(
4291378 .-7850192*cos.(pi*t_s) .+5997464*cos.(2*pi*t_s) .-3808816*cos.(3*pi*t_s) .+1993233*cos.(4*pi*t_s) .-847704*cos.(5*pi*t_s) 
.+286748*cos.(6*pi*t_s) .-74536*cos.(7*pi*t_s) .+14014*cos.(8*pi*t_s) .-1672*cos.(9*pi*t_s) .+ 76*cos.(10*pi*t_s) .+ 8*cos.(11*pi*t_s) .-cos.(12*pi*t_s))
#|m|=12 diag
Sol_Exc11007 = .-(1/(599040*sqrt(231)*pi^12)).*(exp.(alpha_0*t_s)*eta*sigma^12).*(.-54707156 .+ 100826388*cos.(pi*t_s) .-78806832*cos.(2*pi*t_s)
.+52055003*cos.(3*pi*t_s).-28870842*cos.(4*pi*t_s).+13310583*cos.(5*pi*t_s).-5025592*cos.(6*pi*t_s).+1519518*cos.(7*pi*t_s).-355212*cos.(8*pi*t_s)
.+60398*cos.(9*pi*t_s) .-6552*cos.(10*pi*t_s).+273*cos.(11*pi*t_s).+26*cos.(12*pi*t_s)-3*cos.(13*pi*t_s))

plot(t_s, Sol_Exc1, c=:red)
plot!(t_s, Sol_Exc2, c=:red)
plot!(t_s, Sol_Exc3, c=:red)
plot!(t_s, Sol_Exc1002, c=:red)
plot!(t_s, Sol_Exc1005, c=:red)
plot!(t_s, Sol_Exc2006, c=:red)
plot!(t_s, Sol_Exc2007, c=:red)
plot!(t_s, Sol_Exc3007, c=:red)
plot!(t_s, Sol_Exc4007, c=:red)
plot!(t_s, Sol_Exc5007, c=:red)
plot!(t_s, Sol_Exc6007, c=:red)
plot!(t_s, Sol_Exc7007, c=:red)
plot!(t_s, Sol_Exc8007, c=:red)
plot!(t_s, Sol_Exc10007, c=:red)
plot!(t_s, Sol_Exc11007, c=:red)