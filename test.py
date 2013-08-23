# Importing stuff
import os, sys, inspect
from pylab import *
import numpy as np
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"lib")))
if cmd_subfolder not in sys.path:
	sys.path.insert(0, cmd_subfolder)
import _PythonManager
#reload(_PythonManager)
from _PythonManager import PythonManager

matplotlib.pyplot.close('all')

# Load data
data = np.loadtxt('../LSETest2/data/measdata00588.dat')

# Initialize LSE
m = PythonManager("Parameters.xml")
m.resetEstimate(0)
m.clearMeas()
X = np.zeros((len(data),13))

r = range(len(data))

# Run LSE
for i in r:
	data_t = data[i,0]
	data_f = data[i,1:4]
	data_w = data[i,4:7]
	data_e = data[i,7:19]
	data_v = data[i,19:31]
	data_CF = data[i,31:35]
	data_r = data[i,35:38]
	data_q = data[i,38:42]
	m.addImuMeas(data_t,data_f,data_w)
	m.addEncMeas(data_t,data_e,data_v,data_CF)
	m.addPosMeas(data_t,data_r,data_q)
	m.update()
	m.getEst(X[i])

# Initialize LSE
m = PythonManager("ParametersOCEKF.xml")
m.resetEstimate(0)
m.clearMeas()
X_OCEKF = np.zeros((len(data),13))

# Run LSE
for i in r:
	data_t = data[i,0]
	data_f = data[i,1:4]
	data_w = data[i,4:7]
	data_e = data[i,7:19]
	data_v = data[i,19:31]
	data_CF = data[i,31:35]
	data_r = data[i,35:38]
	data_q = data[i,38:42]
	m.addImuMeas(data_t,data_f,data_w)
	m.addEncMeas(data_t,data_e,data_v,data_CF)
	m.addPosMeas(data_t,data_r,data_q)
	m.update()
	m.getEst(X_OCEKF[i])

# # Run delay estimator
# # m.delayIdentification(50,60)
# m.robotCalibration(35,30)
# print(m.getImuTD())
# print(m.getEncTD())
# print(m.getPosTD())
# N = m.getLengthOfBC()
# X = np.zeros((N,17))
# m.getBCData(X)
# 	
# # Plot Results
# figure(1)
# subplot(311)
# plot(X[:,0],X[:,1])
# subplot(312)
# plot(X[:,0],X[:,2])
# subplot(313)
# plot(X[:,0],X[:,3])
# 
# figure(2)
# subplot(311)
# plot(X[:,0],X[:,4])
# subplot(312)
# plot(X[:,0],X[:,5])
# subplot(313)
# plot(X[:,0],X[:,6])
# 
# figure(3)
# subplot(311)
# plot(X[:,0],X[:,7])
# subplot(312)
# plot(X[:,0],X[:,8])
# subplot(313)
# plot(X[:,0],X[:,9])
# 
# figure(4)
# subplot(411)
# plot(X[:,0],X[:,10])
# # plot(data[:,0],data[:,41])
# subplot(412)
# plot(X[:,0],X[:,11])
# # plot(data[:,0],data[:,38])
# subplot(413)
# plot(X[:,0],X[:,12])
# # plot(data[:,0],data[:,39])
# subplot(414)
# plot(X[:,0],X[:,13])
# # plot(data[:,0],data[:,40])
# 
# figure(5)
# subplot(311)
# plot(X[:,0],X[:,14])
# plot(data[:,0],data[:,4])
# subplot(312)
# plot(X[:,0],X[:,15])
# plot(data[:,0],data[:,5])
# subplot(313)
# plot(X[:,0],X[:,16])
# plot(data[:,0],data[:,6])

figure(1)
subplot(311)
plot(data[r,0],X[r,0]-X[0,0], 'b')
plot(data[r,0],X_OCEKF[r,0]-X_OCEKF[0,0], 'r')
plot(data[r,0],data[r,35]-data[0,35], 'k--')
subplot(312)
plot(data[r,0],X[r,1]-X[0,1], 'b')
plot(data[r,0],X_OCEKF[r,1]-X_OCEKF[0,1], 'r')
plot(data[r,0],data[r,36]-data[0,36], 'k--')
subplot(313)
plot(data[r,0],X[r,2]-X[0,2], 'b')
plot(data[r,0],X_OCEKF[r,2]-X_OCEKF[0,2], 'r')
plot(data[r,0],data[r,37]-data[0,37], 'k--')

figure(2)
subplot(311)
plot(data[r,0],X[r,3], 'b')
plot(data[r,0],X_OCEKF[r,3], 'r')
subplot(312)
plot(data[r,0],X[r,4], 'b')
plot(data[r,0],X_OCEKF[r,4], 'r')
subplot(313)
plot(data[r,0],X[r,5], 'b')
plot(data[r,0],X_OCEKF[r,5], 'r')
 
figure(3)
subplot(411)
plot(data[r,0],X[r,6], 'b')
plot(data[r,0],X_OCEKF[r,6], 'r')
plot(data[r,0],data[r,38], 'k--')
subplot(412)
plot(data[r,0],X[r,7], 'b')
plot(data[r,0],X_OCEKF[r,7], 'r')
plot(data[r,0],data[r,39], 'k--')
subplot(413)
plot(data[r,0],X[r,8], 'b')
plot(data[r,0],X_OCEKF[r,8], 'r')
plot(data[r,0],data[r,40], 'k--')
subplot(414)
plot(data[r,0],X[r,9], 'b')
plot(data[r,0],X_OCEKF[r,9], 'r')
plot(data[r,0],data[r,41], 'k--')
 
figure(4)
subplot(311)
plot(data[r,0],X[r,10], 'b')
subplot(312)
plot(data[r,0],X[r,11], 'b')
subplot(313)
plot(data[r,0],X[r,12], 'b')
 
show(block=False)

