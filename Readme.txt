The code includes fast algorithms for simulating Hodgkin-Huxley (HH) neuronal networks. It offers four methods: Runge-Kutta fourth-order scheme (RK4), adaptive scheme based on RK4, library-based scheme and exponential time differencing scheme with Runge-Kutta fourth-
order time stepping (ETD4RK). All the methods have an accuracy of fourth-order.

The RK4, adaptive and ETD4RK methods can obtain precise information of HH neurons like the precise spike shapes while the library method can only obtain statistical information like mean firing rate, fire patterns and chaotic dynamics. The ETD4RK and library methods are recommended to use because of their high efficiency. 

Before using this code, please add the path of the file to "addpath" in Matlab

Format of data
V:  [time; voltage of 1st neuron; voltage of 2nd neuron; ......]
Ras: [spike times; corresponding neuron index], like [2,9,10,15;0,0,1,0]
CM: Connection matirx (N*N), CM_ij = 1 means a connection from the i-th neuron to the j-th neuron
