



function [V, ras, CM] = gen_HH(pm)

pm.path = fileparts(mfilename('fullpath'));
pm.path = [pm.path,filesep];

if pm.method == 0
    pm.dt = 1.0/32;
end

if pm. Lyapunov
    pm.record_Spike = 0;
    pm.record_Voltage = 0;
end


parameter_filename = create_parameter_file(pm);
if isunix
    exe_name = 'HH';
else
    exe_name = 'HH.exe';
end
system([pm.path 'code' filesep exe_name ' ' parameter_filename]);
delete(parameter_filename);

[V,ras,CM] = get_info(pm);
end

function parameters_filename = create_parameter_file(pm)
aa = tempname;id=findstr(aa,filesep);aa=aa(id(end)+1:end);
parameters_filename=[pm.path,'data',filesep,aa,'_parameters.txt'];

fp=fopen(parameters_filename,'r');
if fp ~= -1
    fclose(fp);
    delete(parameters_filename);    
end

fp=fopen(parameters_filename,'a');
fprintf(fp,'NE:\t\t\t %d\r\n',pm.NE);
fprintf(fp,'NI:\t\t\t %d\r\n',pm.NI);
fprintf(fp,'T_Max(ms):\t\t %e\r\n',pm.t);
fprintf(fp,'T_step:\t\t\t %f\r\n',pm.dt);
fprintf(fp,'S:\t\t\t %f %f %f %f\r\n',pm.SEE,pm.SIE,pm.SEI,pm.SII);
fprintf(fp,'I_CONST:\t\t %d\r\n',pm.I_const);
fprintf(fp,'Nu(KHz):\t\t %f\r\n',pm.Nu);
fprintf(fp,'f:\t\t\t %f\r\n',pm.f);
fprintf(fp,'P_c:\t\t\t %f\r\n',pm.Pc);
fprintf(fp,'Random_S:\t\t %d   ',pm.pS);
fprintf(fp,'////(0-none 1-uniform 2-gauss 3-exponential 4-lognormal)\r\n');
fprintf(fp,'Random_Nu:\t\t %d   ',pm.pNu);
fprintf(fp,'////(0-none 1-uniform 2-gauss 3-exponential 4-lognormal)\r\n');
fprintf(fp,'method:\t\t\t %d   ',pm.method);
fprintf(fp,'////0-RK,1-Ad,2-Lib,3-ETD_RK \r\n');
fprintf(fp,'Lyapunov:\t\t %d\r\n',pm.Lyapunov);
fprintf(fp,'record_data(spike/V):\t %d %d\r\n',...
    pm.record_Spike,pm.record_Voltage);

str=[pm.path,'data',filesep];
fprintf(fp,'Record_path:\t\t %s\r\n',str);
str=[pm.path,'lib_data',filesep];
fprintf(fp,'Lib_path:\t\t %s\r\n',str);
fclose(fp);

end

function [V,ras,CM]=get_info(pm)

if pm.record_Spike == 0 && pm.record_Voltage == 0
    V=[];
    ras=[];
    CM=[];
    return;
end

path = [pm.path 'data' filesep];
N=pm.NE+pm.NI;
f = ['NE=',num2str(pm.NE),'NI=',num2str(pm.NI),'_'];

file_CM = [path f 'connect_matrix-p=',num2str(pm.Pc,'%0.3f'),'.dat'];
fp = fopen(file_CM,'rb');
CM = fread(fp,[N,N],'double')';
fclose(fp);

if pm.method == 0
   str=['RK4_t=' num2str(-log2(pm.dt),'%0.3f')];
elseif pm.method == 1
   str=['RK4_t_l=' num2str(-log2(pm.dt),'%0.3f'),'t_s=5.000'];     
elseif pm.method == 2
    str=['RK4_lib_t=' num2str(-log2(pm.dt),'%0.3f')];
else
    str=['ETD4RK_t=' num2str(-log2(pm.dt),'%0.3f')];
end

f1=['p=',num2str(pm.Pc,'%0.3f'),'SEE=',num2str(pm.SEE,'%0.3f'),...
    'SIE=',num2str(pm.SIE,'%0.3f'),'SEI=',num2str(pm.SEI,'%0.3f'),...
    'SII=',num2str(pm.SII,'%0.3f')];
f1=[f1 'f=',num2str(pm.f,'%0.3f'),'u=',num2str(pm.Nu,'%0.3f')];

if pm.record_Voltage
    fp = fopen([path f str f1 '_voltage.dat'],'rb');
    V = fread(fp,[N+1,inf],'double');
    fclose(fp);
else
    V=[];
end

if pm.record_Spike
    fp = fopen([path f str f1 '_spikes.dat'],'rb');
    ras = fread(fp,[2,inf],'double');
    fclose(fp);
else
    ras=[];
end

end
















