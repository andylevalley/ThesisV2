function [beta,beta_lower,beta_upper,output] = NMC_LightingConstraints_fun(problem)

%% Variable Extraction

UTC0 = problem.time.UTC0;
year = UTC0(1);
mo = UTC0(2);
d = UTC0(3);
h = UTC0(4);
mins = UTC0(5);
s = UTC0(6);

r_RSO_0 = problem.RSO.IC(1:3);
v_RSO_0 = problem.RSO.IC(4:6);
r_RSO_tf = problem.RSO.FC(1:3);
v_RSO_tf = problem.RSO.FC(4:6);
RSO_nu0 = problem.RSO.nu0;
sma = problem.RSO.sma;
w = problem.RSO.w;
period_s = problem.RSO.period_s;
p = problem.RSO.p;
ecc = problem.RSO.ecc;
incl = problem.RSO.incl;
omega = problem.RSO.RAAN;
argp = problem.RSO.argp;
arglat = problem.RSO.arglat;
truelon = problem.RSO.truelon;
lonper = problem.RSO.lonper;

traj.margin = problem.traj.margin;
beta_margin = problem.traj.beta_margin;
ae = problem.traj.ae;
x_d = problem.traj.x_d;
y_d0 = problem.traj.y_d0;
z_max = problem.traj.z_max;
gamma = problem.traj.gamma;

tf = problem.time.tf;
JD_0 = problem.time.JD_0;
UTC_tf = problem.time.UTC_tf;

moonangle = problem.FOV.Moon.angle;

%% Given tf, Find Beta and Applicable Margins

if problem.suncon.status == 1
    
    r_sun2RSO_unit_tf = sun2RSO(UTC_tf(1),UTC_tf(2),UTC_tf(3),UTC_tf(4),UTC_tf(5),UTC_tf(6),r_RSO_tf,v_RSO_tf);
    r_sun2RSO_unit_plane_tf = r_sun2RSO_unit_tf(1:2)/norm(r_sun2RSO_unit_tf(1:2));
    alpha2entry = atan2(r_sun2RSO_unit_plane_tf(2),r_sun2RSO_unit_plane_tf(1));
    % xt_sun_1 = (-2*tan(alpha2entry)*y_d0+sqrt(4*tan(alpha2entry)^2*y_d0^2-4*(tan(alpha2entry)^2+1)*(y_d0^2-ae^2)))/(2*(tan(alpha2entry)^2+1));
    xt_sun_1 = (-tan(alpha2entry)*y_d0+sqrt(tan(alpha2entry)^2*ae^2-4*y_d0^2+4*ae^2))/(tan(alpha2entry)^2+4);
    % xt_sun_2 = (-2*tan(alpha2entry)*y_d0-sqrt(4*tan(alpha2entry)^2*y_d0^2-4*(tan(alpha2entry)^2+1)*(y_d0^2-ae^2)))/(2*(tan(alpha2entry)^2+1));
    xt_sun_2 = (-tan(alpha2entry)*y_d0-sqrt(tan(alpha2entry)^2*ae^2-4*y_d0^2+4*ae^2))/(tan(alpha2entry)^2+4);
    % Choose correct sign
    if r_sun2RSO_unit_plane_tf(1) > 0
        xt_sun_choose_pos = 0;
    else
        xt_sun_choose_pos = 1;
    end
    if xt_sun_1 > 0 && xt_sun_choose_pos == 1
        xt_sun = xt_sun_1;
    else
        xt_sun = xt_sun_2;
    end
    yt_sun_1 = y_d0+sqrt(ae^2-4*xt_sun^2);
    yt_sun_2 = y_d0-sqrt(ae^2-4*xt_sun^2);
    if r_sun2RSO_unit_plane_tf(2) > 0 % THIS LOGIC ONLY WORKS FOR y_d0=0
        yt_sun_choose_pos = 0;
    else
        yt_sun_choose_pos = 1;
    end
    if yt_sun_1 > 0 && yt_sun_choose_pos == 1
        yt_sun = yt_sun_1;
    else
        yt_sun = yt_sun_2;
    end
    xdt_sun = (yt_sun-y_d0)*w/2;
    ydt_sun = -2*xt_sun*w;
    
    beta = atan2(xdt_sun,3*w*xt_sun+2*ydt_sun); % SAME order compared to the way Lovell writes it
    beta_lower = beta-beta_margin;
    beta_upper = beta+beta_margin;
    
else
    
    beta = problem.traj.beta;
    beta_lower = beta-problem.traj.beta_margin;
    beta_upper = beta+problem.traj.beta_margin;
    
end

%% Check Feasibility and Modify if Needed

% Calculate the minimum time solution and ensure that the given tf is
% greater than the minimum time solution to inject into the NMC

%% Analyze Moon Conflicts

% Find vector from deputy to chief during region of interest
t_nomoon = period_s;
JD_nomoon_start = JD_0+tf/3600/24;
JD_nomoon_vec = linspace(JD_nomoon_start,JD_nomoon_start+t_nomoon/3600/24);
for i = 1:length(JD_nomoon_vec)
    RSO_nu = RSO_nu0 + w*tf + w*(JD_nomoon_vec(i)-JD_nomoon_start)*3600*24; % deg % THIS ASSUMES CIRCULAR ORBIT
    [r_RSO,v_RSO] = coe2rvh(p,ecc,incl,omega,argp,RSO_nu,RSO_nu,RSO_nu,lonper,problem.mu);
    r_RSO2moon_unit_nomoon = RSO2moon(JD_nomoon_vec(i),r_RSO,v_RSO);
    
    beta_nomoon = beta + w*(JD_nomoon_vec(i)-JD_nomoon_start)*3600*24;
    X_dep = LROE2X(ae,x_d,y_d0,z_max,gamma,beta_nomoon,w);
    X_dep_pos_unit_NEG = -X_dep(1:3)/norm(X_dep(1:3));
    
    alpha_M(i) = acos(dot(r_RSO2moon_unit_nomoon,X_dep_pos_unit_NEG));
end

output.Moon.p1x = (JD_nomoon_vec-JD_nomoon_start)*24;
output.Moon.p1y = rad2deg(alpha_M);
output.Moon.t_nomoon = t_nomoon;

%% Modify Outputs if Moon FOV Constraint is Enabled

if problem.traj.margin == 1
    
    beta_range = linspace(beta-beta_margin,beta+beta_margin);
    clear alpha_M
    for i = 1:length(beta_range)
        
        for j = 1:length(JD_nomoon_vec)
            RSO_nu = RSO_nu0 + w*tf + w*(JD_nomoon_vec(j)-JD_nomoon_start)*3600*24; % deg % THIS ASSUMES CIRCULAR ORBIT
            [r_RSO,v_RSO] = coe2rvh(p,ecc,incl,omega,argp,RSO_nu,RSO_nu,RSO_nu,lonper,problem.mu);
            r_RSO2moon_unit_nomoon = RSO2moon(JD_nomoon_vec(j),r_RSO,v_RSO);
            
            beta_nomoon = beta_range(i) + w*(JD_nomoon_vec(j)-JD_nomoon_start)*3600*24;
            X_dep = LROE2X(ae,x_d,y_d0,z_max,gamma,beta_nomoon,w);
            X_dep_pos_unit_NEG = -X_dep(1:3)/norm(X_dep(1:3));
            
            alpha_M(i,j) = acos(dot(r_RSO2moon_unit_nomoon,X_dep_pos_unit_NEG));
            
        end
        
        alpha_M_min(i) = min(alpha_M(i,:));
        
    end
    
    % Collect output
    output.Moon.p2x = rad2deg(beta_range);
    output.Moon.p2y = rad2deg(alpha_M_min);
    
    if problem.FOV.Moon.status == 1
    
    % Find largest range of beta where there would be no conflict
    goodranges = find(alpha_M_min >= moonangle);
    if isempty(goodranges) == 1
        disp('**Not possible to satisfy Moon FOV constraint. Continuing with violation**');
    else
        marker = 1;
        counter = 1;
        for i = 2:length(goodranges)
            if goodranges(i)-goodranges(i-1) > 1
                goodranges_cell{counter} = goodranges(marker:i-1);
                marker = i;
                counter = counter + 1;
            elseif goodranges(i)-goodranges(i-1) == 1 && i == length(goodranges)
                goodranges_cell{counter} = goodranges(marker:i);
            end
            
        end
        
        for i = 1:size(goodranges_cell,2)
            if i == 1
                largestrange = length(goodranges_cell{i});
                largestrangecol = i;
            elseif length(goodranges_cell{i}) > largestrange
                largestrangecol = i;
            end
        end
        
        largestrangeindices = goodranges_cell{largestrangecol};
        beta_range_largest = beta_range(largestrangeindices);
        beta_lower = beta_range_largest(1);
        beta_upper = beta_range_largest(end);
        if beta < beta_lower
            beta = beta_lower;
        elseif beta > beta_upper
            beta = beta_upper;
        end
        
    end
    
    end
    
end

