%./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\........
%.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////........
%..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\...........
%...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\..........
%....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\.........
%.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\........
%......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\.......
%.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
%........\///////////........\////////......\/////////..........\///......
%=========================================================================
%
% Copyright (C) 2018-2023 Politecnico di Milano,
%                           with support from A^3 from Airbus
%                    and  Davide   Montagnani,
%                         Matteo   Tugnoli,
%                         Federico Fonte
%
% This file is part of DUST, an aerodynamic solver for complex
% configurations.
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without
% restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following
% conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
%
% Author: Alessandro Cocco 
%
%=========================================================================

% This script serves to preview the airfoil mesh that will be created
% in dust_pre after having loaded a .dat profile in Selig format 

clear; close all; clc;

%% Input 

rr = readmatrix('naca4412.dat'); % file containing the profile coordinate 
n_elem = 30;                     % number of elements to discretize  

% discretization type (see DUST manual) 
%mesh_type = 'uniform';           
%mesh_type = 'cosineLE';
%mesh_type = 'cosineTE'; 
%mesh_type = 'cosine'; 
%-> geometric series discretization possibilities  
%mesh_type = 'geoseriesLE'; 
%mesh_type = 'geoseriesTE'; 
r = 1/5; 
%mesh_type = 'geoseries'; 
r_le = 1/5;
r_te = 1/10;
xh = 0.7;               % hinge chordwise position (adimensional)/position 
                        % from where the refinement start 
mesh_type = 'geoseriesHI';
theta = deg2rad(-30);   % hinge deflection (negative downwards) 
offset = 0.05;          % size of the blending region 
r_le_aft = 1/3;         % growth ratio at leading edge of the fixed part
r_te_aft = 1/3;         % growth ratio at trailing edge of the fixed part
r_le_fore = 1/3;        % growth ratio at leading edge of the hinged part
r_te_fore = 1/2;        % growth ratio at trailing edge of the hinged part

% analytic verication 
airfoil_str = 'NACA4412';
M = str2double(airfoil_str(5)  );
P = str2double(airfoil_str(6)  );
SS= str2double(airfoil_str(7:8));
[x_anal, y_anal] = naca4digit(M,P,SS,1,n_elem);

%% Interpolate dat file using a combination of Hermite Spline and a circle sector 

% build middle line 
id_up = floor(size(rr,1)/2);
rr_up = rr(1:id_up + 1,:); % assuming Seilig format 
rr_low = rr(id_up + 1:end,:);
rr_mean(:,1) = rr_low(:,1);
rr_mean(:,2) = (flipud(rr_up(:,2)) + rr_low(:,2))/2;

% get m and p
[m, idx_m] = max(rr_mean(:,2)); % y_mean line max thickness
p = rr_mean(idx_m,1);           % x_mean line max thickness
% The mean line is given by the expression: y_m = m/p^2(2px - x^2) 
% The expression is related NACA profile, but should be valid for all profile
% close to the leading edge. 

% Get analytical derivative of the mean line in x = 0  
m_line = 2*m/p; % need to get the center of the circle -> linearized approach

% calculate thickness in a simple manner 
t = max(rr_up(:,2)) - min(rr_low(:,2));  
%t = 0.12; 
% radius of the circle at the leading edge 
circ.radius = 1.10*t^2; 

% get center of the circle with the following conditions:
% 1. passing through [0 0]
% 2. center lying on the linearized middle line
% 3. radius given by the NACA formula
circ.x_c = sqrt(circ.radius^2/(m_line^2 + 1));
circ.y_c = sqrt(circ.radius^2 - circ.x_c^2);

% Define the region that is approximated by a circle arc  
% - get point at +/- 37.5 deg (angle that works the best) 
alpha = atan(circ.y_c/circ.x_c); % angle of m 
beta_up = (180-37.5)*pi/180;
beta_low = (180+37.5)*pi/180;
circ.x_plus = circ.radius*cos(beta_up + alpha) + circ.x_c;
circ.y_plus = circ.radius*sin(beta_up + alpha) + circ.y_c;
circ.x_minus = circ.radius*cos(beta_low + alpha) + circ.x_c;
circ.y_minus = circ.radius*sin(beta_low + alpha) + circ.y_c;

% get the derivative for tangent at leading edge (to use in the
% splining process 
circ.tan_plus = tan(atan(-cot(beta_up)) + alpha);
circ.tan_minus = tan(atan(-cot(beta_low)) + alpha);

% leading edge circle sector
theta_start = beta_up + alpha;
theta_end = beta_low + alpha;
circ.theta = linspace(theta_start, theta_end, 10);
circ.x = circ.radius.*cos(circ.theta) + circ.x_c;
circ.y = circ.radius.*sin(circ.theta) + circ.y_c;

% tangent at trailing edge: backward difference
tang_te_up = (rr_up(1,2) - rr_up(2, 2)) / ...
             (rr_up(1,1) - rr_up(2, 1));
tang_te_low = (rr_low(end-1, 2) - rr_low(end, 2)) / ...
              (rr_low(end-1, 1) - rr_low(end, 1));

% spline interpolation upper part
point_up = flipud(rr_up); % point from ,dat file 
% replace first point with circle point
point_up(1,:) = [circ.x_plus, circ.y_plus];
xq_up = linspace(circ.x_plus, 1, 10); % from circle_up to TE 
yq_up = hermite_spline(point_up(:,1), point_up(:,2), xq_up, ...
                       circ.tan_plus, tang_te_up);

% spline interpolation lower part
point_low = rr_low; % point from .dat file
point_low(1,:) = [circ.x_minus, circ.y_minus];
xq_low = linspace(circ.x_minus, 1, 10); 
yq_low = hermite_spline(point_low(:,1), point_low(:,2), xq_low, ...
                        circ.tan_minus, tang_te_low);

%% reorganize vectors to prepare mesh subdivision
[min_x, idx_circ_min] = min(circ.x); % get lowest x (in general is <= 0) 

circ.x_up = circ.x(1:idx_circ_min);
circ.y_up = circ.y(1:idx_circ_min);
circ.x_low = circ.x(idx_circ_min:end);
circ.y_low = circ.y(idx_circ_min:end);
circ.x_up = fliplr(circ.x_up);
circ.y_up = fliplr(circ.y_up);

% upper 
rr2mesh_up = [circ.x_up xq_up(2:end);
                circ.y_up yq_up(2:end)];
% lower part 
rr2mesh_low = [circ.x_low xq_low(2:end);
                circ.y_low yq_low(2:end)];

%% Mesh possibilites: 
% uniform 
% cosineLE
% cosineTE 
% cosine 
% geometric series 
%   - LE 
%   - TE 
%   - both 
%   - hinge 

switch mesh_type 
    case 'uniform'
        dcsi = linspace(0., 1., n_elem +1);
    case 'cosineLE'
        dcsi = linspace(0., 1., n_elem +1);
        dcsi = 1.0-cos(0.5*pi.*dcsi);   
    case 'cosineTE'
        dcsi = linspace(0., 1., n_elem +1);
        dcsi = sin(0.5*pi.*dcsi);       
    case 'cosine'
        dcsi = linspace(0., 1., n_elem +1);     
        dcsi = 1/2.0.*(1.0-cos(pi.*dcsi));   
    case 'geoseriesLE'
        [~, dcsi] = geoseries(0., 1., n_elem, r); 
    case 'geoseriesTE'
        [dcsi, ~] = geoseries(0., 1., n_elem, r); 
    case 'geoseries'
        dcsi = geoseries_both(0., 1., n_elem, xh, r_le, r_te);
    case 'geoseriesHI'
        dcsi = geoseries_hinge(0., 1., n_elem, xh, ...
                                r_le_aft, r_te_aft, r_le_fore, r_te_fore);
end 

% linear interpolation to get y coordinates
y_mesh_up = interp1(rr2mesh_up(1,:), rr2mesh_up(2,:), dcsi,'linear','extrap');
y_mesh_low = interp1(rr2mesh_low(1,:), rr2mesh_low(2,:), dcsi,'linear','extrap');

% mean line 
y_mesh_mean = (y_mesh_low + y_mesh_up)/2;
rr_mean = [dcsi; y_mesh_mean];

rr_final = [fliplr(dcsi) dcsi(2:end);
            fliplr(y_mesh_low) y_mesh_up(2:end)];

% hinge deflection
rr_deflected = hinge_deflection(rr_final, theta, xh, offset);

%% plotting results 
figure
plot(rr_final(1,:),rr_final(2,:),'r-+','LineWidth',2)
hold on 
plot(rr_mean(1,:),rr_mean(2,:))
plot(rr_deflected(1,:), rr_deflected(2,:), 'b-+', 'LineWidth',2) 
axis equal
grid on 

%% functions --------------------------------------------------------------
function yq = hermite_spline(x, y, xq, tang_start, tang_end)

    % Input data
    % x:  x database
    % y:  y database
    % xq: x query
    % tang_start: tangent coefficient at x(1)
    % tang_end:   tangent coefficient at x(end)
    %
    % Output data
    % yq: y interpolated
    
    % build tangent vector using centered differences
    m = zeros(1,numel(x));
    m(1) = tang_start;
    m(end) = tang_end;
    for i = 2:numel(x) - 1
        m(i) = 0.5*((y(i + 1) - y(i))/(x(i + 1) - x(i)) + ...
                    (y(i) - y(i - 1))/(x(i) - x(i - 1)));
    end
    
    yq = zeros(1,numel(xq));
    for i = 1:numel(xq)
        for j = 1:numel(x) - 1
            % two consecutive points x for which the point xq belongs
            if (xq(i) >= x(j)) && (xq(i) <= x(j + 1))
                t = (xq(i) - x(j))/(x(j + 1) - x(j));
    
                h00 = hermite_p1(t);
                h10 = hermite_p2(t);
                h01 = hermite_p3(t);
                h11 = hermite_p4(t);
    
                yq(i) = h00*y(j) + ...
                    h10*(x(j + 1) - x(j))*m(j) + ...
                    h01*y(j + 1) + ...
                    h11*(x(j + 1) - x(j))*m(j + 1);
            end
        end
    end
end

function h = hermite_p1(t)
    h = 2*t^3 - 3*t^2 + 1;
end

function h = hermite_p2(t)
    h = t^3 - 2*t^2 + t;
end

function h = hermite_p3(t)
    h = -2*t^3 + 3*t^2;
end

function h = hermite_p4(t)
    h = t^3 - t^2;
end

function [dcsi_te, dcsi_le] = geoseries(start_x, end_x, n_elem, r)
    
    dl = zeros(1,n_elem);
    dcsi_te = zeros(1,n_elem + 1);
    % check series convergence 
    if r >= 1 || r <= 0 
        error('insert a growth value between 0 and 1') 
    end

    r = 1-r; % growth ratio (intended as decreasing ratio) 
    % get size of the first element knowing the sum of the geometric series
    dl(1) = (1-r)/(1-r^(n_elem));
    for i = 2:n_elem
        dl(i) = r*dl(i-1); % geometric series-> get length of element  
    end
    for i = 1:n_elem
        dcsi_te(i + 1) = dcsi_te(i) + dl(i); %-> position  
    end

    dcsi_le = fliplr(1-dcsi_te); 
    % rescale in the interval
    dcsi_te = rescale_csi(dcsi_te, start_x, end_x);
    dcsi_le = rescale_csi(dcsi_le, start_x, end_x);
end


function dcsi = geoseries_both(start_x, end_x, n_elem, xh, r_le, r_te)

    n_elem_le = ceil(n_elem*xh);
    n_elem_te = n_elem - n_elem_le;
    mid_point = (start_x + end_x)*xh;
    start_x_le = start_x;
    end_x_le = mid_point;
    start_x_te = mid_point;
    end_x_te = end_x;
    [~, dcsi_le] = geoseries(start_x_le, end_x_le, n_elem_le, r_le);
    [dcsi_te, ~] = geoseries(start_x_te, end_x_te, n_elem_te, r_te);
    dcsi = [dcsi_le(1:end-1) dcsi_te];
end

function dcsi = geoseries_hinge(start_x, end_x, n_elem, xh, ...
                                   r_le_aft, r_te_aft, r_le_fore, r_te_fore)

    n_elem_aft = ceil(n_elem*xh);
    n_elem_fore = n_elem - n_elem_aft;
    start_x_aft = start_x;
    end_x_aft = xh;
    start_x_fore = xh;
    end_x_fore = end_x;
    dcsi_aft = geoseries_both(start_x_aft, end_x_aft, n_elem_aft, 0.5, r_le_aft, r_te_aft);
    dcsi_fore = geoseries_both(start_x_fore, end_x_fore, n_elem_fore, 0.5, r_le_fore, r_te_fore);
    dcsi = [dcsi_aft dcsi_fore];

end

function dcsi_rescaled = rescale_csi(dcsi, rmin, rmax)
    dcsi_rescaled = dcsi*(rmax - rmin) + rmin; 
end

% test hinge deflection 
function rr = hinge_deflection(rr, theta, xh, offset)
    position = [xh; 0];
    for ip = 1:size(rr,2)
        d = (rr(:,ip) - position)' * [1; 0] ;
        if d>= offset    
            Rot = [cos(theta), -sin(theta); ...
                   sin(theta),  cos(theta)];
            rr(:,ip) = position + Rot*(rr(:,ip) - position);
        elseif d >= -offset
            u  = offset;
            m = -cot(theta);
            yc = -m*u*(1 + cos(theta)) + u*sin(theta);
            thp = 0.5*(d+u)/u*theta;
            rr(:,ip) = [yc*sin(thp) - u + position(1) - rr(2,ip)*sin(thp); ...
                        yc*(1.-cos(thp))+ rr(2,ip)*cos(thp)             ];
        end
    end
end


%% %%-- just for verification --%%
function [x,y] = naca4digit(M,P,SS,c,n)

    m  = M  / 100;
    p  = P  / 10;
    t  = SS / 100;
    
    if ( m == 0 ) % p must be .ne. 0 
        p = 1 ;
    end
    
    %> === Chord discretization ===
    xv = linspace(0.0,c,n+1);            % uniform
    % xv = c/2.0 .*(1.0-cos(pi.*xv./c));   % cosine
    xv = c .*(1.0-cos(0.5*pi.*xv./c));   % half-cosine
    
    %> === Thickness ===
    ytfcn = @(x) 5.*t.*c.*(0.2969.*(x./c).^0.5 - 0.1260.*(x./c) ...
        - 0.3516.*(x./c).^2 + 0.2843.*(x./c).^3 - 0.1015.*(x./c).^4);
    % 1015: open TE
    % 1036: closed TE
    yt = ytfcn(xv);
    
    %> === Mean line ===
    yc = zeros(size(xv));

    for ii = 1 : n+1
        if xv(ii) <= p*c
            yc(ii) = c*(m/p^2 *(xv(ii)/c) * (2*p - (xv(ii)/c)));
        else
            yc(ii) = c*(m/(1-p)^2 * (1 + (2*p - (xv(ii)/c))*(xv(ii)/c) -2*p));
        end
    end
    
    %> === Mean line slope ===
    dyc = zeros(size(xv));
    
    for ii = 1 : n+1
        if xv(ii) <= p*c
            dyc(ii) = m/p^2 * 2*(p-xv(ii)/c);
        else
            dyc(ii) = m/(1-p)^2 * 2*(p-xv(ii)/c);
        end
    end
    
    %> === Upper and lower sides ===
    th = atan2(dyc,1);
    xU = xv - yt.*sin(th);
    yU = yc + yt.*cos(th);
    xL = xv + yt.*sin(th);
    yL = yc - yt.*cos(th);
    
    %> Sort points
    x = zeros(1,2*n+1);
    y = zeros(1,2*n+1);
    for ii = 1 : n
        x(ii) = xL(n+2-ii);
        y(ii) = yL(n+2-ii);
    end
    
    x(n+1:2*n+1) = xU;
    y(n+1:2*n+1) = yU;
end
