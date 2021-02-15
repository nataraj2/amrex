from numpy import *
from matplotlib.pyplot import *

data_default = loadtxt('Tmax_vs_time_default.curve')
data_LS = loadtxt('Tmax_vs_time_LS.curve')

plot(data_default[:,0],data_default[:,1],'b',linewidth=2,label='Default')
plot(data_LS[:,0],data_LS[:,1],'k',linewidth=2,label='Least squares')

xlabel('Time (s)',fontsize=15)
ylabel('T (K)',fontsize=15)
legend()

title('T$_{max}$ vs time for plane',fontsize=20)

figure(2)

data_default = loadtxt('Tmax_vs_time_Cylinder_Default.curve')
data_LS = loadtxt('Tmax_vs_time_Cylinder_LS.curve')

plot(data_default[:,0],data_default[:,1],'b',linewidth=2,label='Default')
plot(data_LS[:,0],data_LS[:,1],'k',linewidth=2,label='Least squares')

legend()

title('T$_{max}$ vs time for cylinder',fontsize=20)
show()
