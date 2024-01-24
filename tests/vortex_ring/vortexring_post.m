clear; close all; clc;

addpath(fullfile('..','..','utils','parser_post'));
format long
time = 0:0.1:29.9; % seconds
re = 7500;
nu = 1.469387e-5;
gamma0 = re*nu;
R0 = 0.5;
delta0 = 0.2;
epsilon = delta0/R0;
C = -0.558; % gaussian
N_points = 100;
delta = sqrt(delta0^2 + 4*nu.*time);
a_e = 1.3607.*delta;
U_analytical = gamma0/(4*pi*R0)*(log(8*R0./a_e) + C); % analytical velocity
Z_analytical = U_analytical.*time;

%% matrix construction
for i = 1:numel(time)
    x = zeros(N_points,1);
    y = linspace(-(R0 + delta0)*1.33, (R0 + delta0)*1.33, N_points)';
    z = Z_analytical(i).*ones(N_points,1);
    fid = fopen(sprintf('points_%d.dat',i), 'wt');
    fprintf(fid, '%d\n', N_points);
    for j = 1:N_points
        fprintf(fid, '%16.16f %16.16f %16.16f\n', x(j), y(j), z(j));
    end
    fclose(fid);
end

%% write dust_post
fid = fopen('dust_post_probe.in', 'wt');

fprintf(fid, 'basename = Postpro/2DTest\n');
fprintf(fid, 'data_basename = Output/2DTest\n\n');

for i = 1:numel(time)
    fprintf(fid, 'analysis = {\n');
    fprintf(fid, '  type = probes\n');
    fprintf(fid, '  name = prb%d\n',i);
    fprintf(fid, '  start_res = %d\n',i);
    fprintf(fid, '  end_res = %d\n',i);
    fprintf(fid, '  step_res = 1\n');
    fprintf(fid, '  format = dat\n');
    fprintf(fid, '  variable = Velocity\n');
    fprintf(fid, '  input_type = from_file\n');
    fprintf(fid, '  file = points_%d.dat\n',i);
    fprintf(fid, '}\n\n');
end
fclose(fid);

% %% run dust_post assuming to be in the same folder
% if strcmp(computer,'GLNXA64') % linux
%     system('dust_post_test dust_post_probe.in');
% elseif strcmp(computer,'PCWIN64') % windows
%     system('bash -c dust_post_test dust_post_probe.in');
% end

%% read probes
meanvelz = zeros(numel(time),1);
for i = 1:numel(time)
    file = fullfile('Postpro',sprintf('2DTest_prb%d.dat',i));
    probe_data = probe_load(file);
    time(i) = probe_data.velocity.time;
    meanvelz(i) = mean(probe_data.velocity.value(3,:));
end


%% some comparisons
time_alvarez = [-0.0091	0.2252	0.4127	0.5533	0.7095	0.8970	1.0845	1.2407	1.3657	1.5532	1.7407	1.9750	2.1781	2.4281	2.5999	2.8655	3.0843	3.3655	3.5217	3.7248	3.9748	4.2248	4.5216	4.8497	5.1309	5.4902	5.8183	6.1308	6.5214	6.8182	7.1463	7.7243	8.1930	8.7555	9.3179	9.8335	10.1615	10.4428	10.8490	11.2239	11.6145	11.8957	12.1769	12.4738];
u_alvarez = [0.2265	0.2270	0.2281	0.2289	0.2298	0.2304	0.2304	0.2296	0.2289	0.2280	0.2270	0.2263	0.2262	0.2267	0.2272	0.2276	0.2277	0.2274	0.2271	0.2268	0.2268	0.2270	0.2271	0.2270	0.2267	0.2263	0.2263	0.2263	0.2263	0.2262	0.2259	0.2257	0.2255	0.2253	0.2252	0.2251	0.2251	0.2249	0.2248	0.2247	0.2247	0.2246	0.2244	0.2244];
time_dns = [0.8502	1.8969	2.7562	3.9279	5.0528	6.3495	7.6775	8.8804	10.0053	10.9583	12.0676	12.4425];
u_dns = [0.2296	0.2290	0.2286	0.2281	0.2274	0.2272	0.2263	0.2259	0.2253	0.2248	0.2244	0.2244];
time_dust  = time'*gamma0/R0^2;
u_dust = meanvelz*R0/gamma0;

figure();
set(gcf, 'Position', [0, 0, 700, 500])
plot(time_dust, u_dust, 'linewidth', 2)
hold on
plot(time_alvarez, u_alvarez,  'linewidth', 2)
plot(time_dns, u_dns,'--', 'linewidth', 2, 'color', 'k')
legend('DUST', 'Alvarez', 'DNS', 'interpreter', 'latex', 'fontsize', 16)
grid on
xlim([0 12.5])
ylim([0.21 0.24])
xlabel('$t \frac{\Gamma_0}{R_0^2}$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$U \frac{R_0}{\Gamma_0}$', 'interpreter', 'latex', 'fontsize', 16)
set(gca,'TickLabelInterpreter','latex', 'fontsize', 16)


%% export table

writematrix([time_dust, u_dust], 'table_xlx.dat', 'Delimiter', '\t')
