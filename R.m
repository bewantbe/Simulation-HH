addpath('D:\code\HH\HH_Matlab_interface');    %path for this Run_HH_model.m file

pm=[];                % a new parameter set;
pm.method = 3;        % methods,0--RK4, 1--Adaptive,2--Library,3--ETD4RK
pm.t    = 2e3;        % Total run time
pm.dt   = 1.0/4;      % default: 1/4. Time step 
pm.NE = 100;          % default: 100. Number of Excitatory neurons.
pm.NI = 0;            % default: 0.   Number of Inhibitory neurons.
                      % Indexes are the later half
pm.SEE = 0.00;        % default: 0. Strength from Ex. to Ex.
pm.SIE = 0.00;        % default: 0. Strength from Ex. to In.
pm.SEI = 0.00;        % default: 0. Strength from In. to Ex.
pm.SII = 0.00;        % default: 0. Strength from In. to In.
pm.Pc = 0.1;          % default: 0.1. Connection probability  
pm.pS  = 0;           % distribution of Strength. default:0--fixed value
                      % 1--Uniform,2--Gaussion,3--Exponentional,4--Lognormal  
pm.I_const = 0;       % default: 0. Evolve model with constant input                         
pm.Nu = 0.1;          % default: 0.1. Poisson input rate(KHz)  
pm.f = 0.1;           % default: 0.1. Poisson input strength 
pm.pNu = 0;           % distribution of Poisson rate. default:0--fixed value
                      % 1--Uniform,2--Gaussion,3--Exponentional,4--Lognormal                        
pm.Lyapunov = 0;      % default: 0. Compute the largest lyapunov exponent
pm.record_Spike = 0;  % default: 0. Record the raster
pm.record_Voltage = 0;% default: 0. Record the voltages

tic;
[V, ras, CM] = gen_HH(pm);
toc
%V=[t;V1;V2;...], ras=[spike time;neuron index]

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 




















