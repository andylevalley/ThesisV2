function [beta_new,beta_lower,beta_upper,output] = Teardrop_LightingConstraints_fun(problem)

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

beta_cutoff = problem.traj.pogo.beta_cutoff;
T_p = problem.traj.pogo.T_p;


%% Calculate Time Left Until Desired Sun Lighting

if problem.suncon.status == 1

    r_sun2RSO_unit_0 = sun2RSO(year,mo,d,h,mins,s,r_RSO_0,v_RSO_0);

    % Find the in-plane angle
    alpha_s_proj = atan2(r_sun2RSO_unit_0(2),r_sun2RSO_unit_0(1));

    % alpha_s_proj is how far the angle is from pointing along the radial axis,
    % which is where we want the Sun when at the top of the (lower) teardrop

    % atan2 returns a value between -pi and pi.  The sun vector can be
    % considered inertial if we are dealing with just several days

    % Change the angle to be between 0 and 2*pi
    if alpha_s_proj < 0
        alpha_s_proj = alpha_s_proj + 2*pi;
    end

    % Change the angle to a time, assuming that it is mostly inertial and will
    % rotate in the orbital plane of the RSO at a rate equal to the mean motion
    % of the RSO
    tau = alpha_s_proj/w;
    if traj.margin == 1
        tau_max = (alpha_s_proj+beta_margin)/w;
        tau_min = (alpha_s_proj-beta_margin)/w;
    end

    %% Check Feasibility and Modify if Needed

%     % See if there's even time to do half of the teardrop first
%     if traj.margin == 0 && tau < 1/2*T_p
%         proceedflag = 0;
%         tflowflag = -1;
%         disp('1/2*T_p greater than time left until desired Sun vector.  Either reduce T_p or wait until next Sun vector pass.  Calculating solution for next Sun vector pass');
%         tau = tau+2*pi/w; % wait another period, assuming that the Sun moves 360 deg in that time (which is an approximation)
%     elseif traj.margin == 1 && tau_max < 1/2*T_p
%         proceedflag = 0;
%         tflowflag = -1;
%         disp('1/2*T_p greater than time left until desired Sun vector.  Either reduce T_p or wait until next Sun vector pass.  Calculating solution for next Sun vector pass');
%         tau = tau+2*pi/w; % wait another period, assuming that the Sun moves 360 deg in that time (which is an approximation)
%     end

%% MIN TIME STUFF

%     % Find the minimum time solution to inject into the first half of the teardrop
%     problem.time.tau = tau;
%     if traj.margin == 1
%         problem.time.tau_max = tau_max;
%         problem.time.tau_min = tau_min;
%     end
%     minTime_solution = RPO_minTime_fun(problem); % the bounds on beta should allow for the second sun pass if needed - just subtract 2*pi from original lower bound, right?, yes, but leaving it for now...
%     if traj.margin == 0
%         problem.Xguess = minTime_solution(2:end-1); % skip tf and don't include the beta from min time solution
%     elseif traj.margin == 1
%         problem.Xguess = minTime_solution(2:end); % keep beta for sunconsoft initial guess
%         %     problem.bounds.lb = [problem.bounds.lb,];
%         %     problem.bounds.ub = [problem.bounds.ub,];
%     end
% 
%     tf_min = minTime_solution(1)
%     problem.time.tf_min = tf_min;
% 
%     if traj.margin == 0 && (tau < 1/2*T_p + tf_min)
%         disp('Not enough time to inject (even with minimum time solution).  This should have been corrected with solution from above.  If not, need to recalculate minimum time solution based on new timing requirement for tau = tau + 2*pi/w.');
%         tau = tau+2*pi/w;
%     elseif traj.margin == 1 && (tau_max < 1/2*T_p + tf_min)
%         disp('Not enough time to inject (even with minimum time solution).  This should have been corrected with solution from above.  If not, need to recalculate minimum time solution based on new timing requirement for tau = tau + 2*pi/w.');
%         tau = tau+2*pi/w;
%     end
% 
%     % For an individual case
%     if problem.time.tf < tf_min
%         tflowflag = 1;
%         disp('Chosen tf is too low.  It is lower than the minimum time solution.  Increasing tf to minimum time solution...');
%         tf = tf_min;
%         problem.time.tf = tf;
%     end
% 
%     %     if tau < 1/2*T_p + tf_min
%     %         proceedflag = 0;
%     %         disp('Not enough time (with minimum time solution) to obtain desired Sun vector.  Either reduce T_p or wait until next Sun vector pass.');
%     %     else
%     %         proceedflag = 1;
%     %     end

%%
nextpass = 0;
% THIS IS FOR THE INDIVIDUAL CASE
if traj.margin == 0 && (tau < 1/2*T_p + problem.time.tf) % assuming it checks the first argument first...
    proceedflag = 0;
    tflowflag = 0; % meaning it is too high
    disp('Does not arrive in time to obtain desired lighting.  Reduce t_f, or reduce T_p, or wait until next Sun vector pass.  Calculating solution for next Sun vector pass');
%     pause;
%     if problem.options.varytf_pogo == 0
        tau = tau+2*pi/w;
        nextpass = 1;
%     end
elseif traj.margin == 1 && (tau_max < 1/2*T_p + problem.time.tf)
    proceedflag = 0;
    tflowflag = 0;
    disp('Does not arrive in time to obtain desired lighting.  Reduce t_f, or reduce T_p, or wait until next Sun vector pass.  Calculating solution for next Sun vector pass');
%     pause;
%     if problem.options.varytf_pogo == 0
        tau = tau+2*pi/w;
        nextpass = 1;
%     end
end
output.nextpass = nextpass;

% Calculate the amount of time needed to coast up to the teardrop
% if problem.options.varytf_pogo == 1
%     coasttime = tau-1/2*T_p-tf_min;
% elseif problem.options.varytf_pogo == 0
    coasttime = tau-1/2*T_p-problem.time.tf;
% end
% Cacluate the amount of beta to shift the entry point
betashift = w*coasttime;
proceedflag = 1;
tflowflag = -1; % meaning it is meaningless


beta_new = beta_cutoff-betashift;

% sunconsoft = problem.options.sunconsoft;
% if sunconsoft == 1
beta_upper = min(beta_cutoff,beta_new + beta_margin);
beta_lower = beta_new - beta_margin;

else
    coasttime = (beta_cutoff-problem.traj.beta)/w;
    beta_new = problem.traj.beta;
    beta_lower = beta_new-problem.traj.beta_margin;
    beta_upper = beta_new+problem.traj.beta_margin;

end
% else
%     beta_upper = beta_new;
%     beta_lower = beta_new;
% end

% % Check bounds here - but don't check beta with sunconsoft...
% if traj.margin == 1
%     problem.bounds.lb = [problem.bounds.lb,beta_lower];
%     problem.bounds.ub = [problem.bounds.ub,beta_upper];
% end
% problem.Xguess
% for i = 1:length(problem.Xguess)
%     if problem.Xguess(i) == problem.bounds.lb(i)
%         problem.Xguess(i) = problem.Xguess(i) + eps;
%     elseif problem.Xguess(i) == problem.bounds.ub(i)
%         problem.Xguess(i) = problem.Xguess(i) - eps;
%     end
% end

%% Varytf stuff

% varytf_pogo = problem.options.varytf_pogo;
% varytf_max = tf_min + 2*pi/w;
% 
% a_e = problem.LROEs.ae;
% x_d = problem.LROEs.x_d;
% y_d0 = problem.LROEs.y_d0;
% z_max = problem.LROEs.z_max;
% gamma = problem.LROEs.gamma;
% 
% problem.varytf.newIGflag = 1;
% IGchangeflag = 0;
% if varytf_pogo == 1
%     tf_range = linspace(tf_min,varytf_max,5);
%     varytftic = tic;
%     hwaitbar = waitbar(0);
%     for i = 1:length(tf_range)
%         
%         if traj.margin == 0
%             if tau-1/2*T_p-tf_range(i) >= 0
%                 tauvec(i) = tau;
%             elseif tau-1/2*T_p-tf_range(i) < 0
%                 tauvec(i) = tau+2*pi/w;
%                 if IGchangeflag == 0
%                     %                 problem.Xguess = Xguess_minTime;
%                     problem.varytf.newIGflag = 1;
%                     IGchangeflag = 1;
%                     nextpassindex = i;
%                     
%                     
%                 end
%             end
%         elseif traj.margin == 1
%             if tau_max-1/2*T_p-tf_range(i) >= 0
%                 tauvec(i) = tau;
%             elseif tau_max-1/2*T_p-tf_range(i) < 0
%                 tauvec(i) = tau+2*pi/w;
%                 if IGchangeflag == 0
%                     %                 problem.Xguess = Xguess_minTime;
%                     problem.varytf.newIGflag = 1;
%                     IGchangeflag = 1;
%                     nextpassindex = i;
%                     
%                     
%                 end
%             end
%         end
%         
%         
%         
%         problem.time.tf = tf_range(i);
%         problem.time.tau = tauvec(i);
%         coasttime = tauvec(i)-1/2*T_p-tf_range(i);
%         betashift = w*coasttime;
%         beta_new = beta_cutoff-betashift;
%         if traj.margin == 1
%             problem.beta_upper = min(beta_cutoff,beta_new + beta_margin);
%             problem.beta_lower = beta_new - beta_margin;
%             problem.bounds.lb(8) = problem.beta_lower;
%             problem.bounds.ub(8) = problem.beta_upper;
%         else
%             problem.beta_upper = beta_new;
%             problem.beta_lower = beta_new;
%             % and bounds remain as they were
%         end
%         problem.Xt = LROE2X(a_e,x_d,y_d0,z_max,gamma,beta_new,w);
%         varytf_RPO_minFuel_solution(i,:) = RPO_minFuel_fun(problem);
%         %         if newIGflag == 0
%         problem.Xguess = varytf_RPO_minFuel_solution(i,4:end);
%         problem.varytf.newIGflag = 1;
%         %         end
%         waitbar(i/length(tf_range));
%     end
%     varytftime = toc(varytftic)
%     close(hwaitbar);
%     
%     t1f = varytf_RPO_minFuel_solution(:,4);
%     alpha1 = varytf_RPO_minFuel_solution(:,5);
%     phi1 = varytf_RPO_minFuel_solution(:,6);
%     tcf = varytf_RPO_minFuel_solution(:,7);
%     alpha2 = varytf_RPO_minFuel_solution(:,8);
%     phi2 = varytf_RPO_minFuel_solution(:,9);
%     t2f = varytf_RPO_minFuel_solution(:,10);
%     
%     exitflags = varytf_RPO_minFuel_solution(:,1);
%     
%     t1dur = tcf.*t2f.*tf_range'-t1f.*tcf.*t2f.*tf_range';
%     t2dur = tf_range'-t2f.*tf_range';
%     engine_on = t1dur+t2dur;
%     
%     nextpassline = linspace(max(engine_on),min(engine_on));
%     
%     figure;
%     sp1 = subplot(2,1,1);
%     sp11 = plot(tf_range/3600,engine_on/60);
%     set(sp11,'LineWidth',1.25);
%     xl = xlabel('$t_f$, Allowed Maneuver Time (hr)','Interpreter','LaTeX');
%     yl = ylabel('Engine-on Time, $t_{on}$ (min)','Interpreter','LaTeX');
%     set(xl,'FontSize',18);
%     set(yl,'FontSize',18);
%     hold on;
%     grid on;
%     sp12 = plot(linspace(tf_range(nextpassindex)/3600,tf_range(nextpassindex)/3600),nextpassline/60,'--k');
%     set(sp12,'LineWidth',1.25);
%     set(gca,'FontSize',16,'FontName','Times');
%     ll = legend([sp12],{'Next Pass'});
%     set(ll,'FontSize',18,'Interpreter','LaTeX');
%     sp2 = subplot(2,1,2);
%     sp21 = plot(tf_range/3600,exitflags,'-o');
%     set(sp21,'LineWidth',1.25);
%     hold on;
%     grid on;
%     sp22 = plot(linspace(tf_range(nextpassindex)/3600,tf_range(nextpassindex)/3600),linspace(min(exitflags),max(exitflags)),'--k');
%     set(sp22,'LineWidth',1.25);
%     ll = legend([sp22],{'Next Pass'});
%     set(ll,'FontSize',18,'Interpreter','LaTeX');
%     set(gca,'FontSize',16,'FontName','Times');
%     linkaxes([sp1,sp2],'x');
%     xl = xlabel('$t_f$, Allowed Maneuver Time (hr)','Interpreter','LaTeX');
%     yl = ylabel('\emph{fmincon} Exit Flag','Interpreter','LaTeX');
%     set(xl,'FontSize',18);
%     set(yl,'FontSize',18);
%     set(gca,'ytick',min(exitflags):max(exitflags));
%     
%     plotinputs.makenewfig = 1;
%     plotinputs.IC = problem.dep.IC;
%     plotinputs.coast = problem.options.coast;
%     plotinputs.rpo_type = problem.options.rpo_type;
%     plotinputs.w = problem.RSO.w;
%     plotinputs.dep.engine.a0 = problem.dep.engine.a0;
%     plotinputs.dep.engine.c = problem.dep.engine.c;
%     for i = 1:length(tf_range)
%         if i == 2
%             plotinputs.makenewfig = 0;
%             plotinputs.fig2D = plotreturn.fig2D;
%             plotinputs.fig3D = plotreturn.fig3D;
%         end
%         plotinputs.x = varytf_RPO_minFuel_solution(i,4:end);
%         plotinputs.time.tf = tf_range(i);
%         plotreturn = cbcb_plot(plotinputs);
%     end
%     figure(plotreturn.fig2D);
%     xl = xlabel('$y$ (km)','Interpreter','LaTeX');
%     yl = ylabel('$x$ (km)','Interpreter','LaTeX');
%     set(xl,'FontSize',18);
%     set(yl,'FontSize',18);
% %     ll = legend('Coast 1','','Burn 1','Coast 2','','Burn 2','Resulting Path');
%     ll = legend([plotreturn.D2.c,plotreturn.D2.b1,plotreturn.D2.b2,plotreturn.D2.nm],{'Coast','Burn 1','Burn 2','Result'});
%     set(ll,'FontSize',18,'Interpreter','LaTeX');
%     set(gca,'FontSize',16,'FontName','Times');
%     figure(plotreturn.fig3D);
%     xl = xlabel('$y$ (km)','Interpreter','LaTeX');
%     yl = ylabel('$x$ (km)','Interpreter','LaTeX');
%     zl = zlabel('$z$ (km)','Interpreter','LaTeX');
%     set(xl,'FontSize',18);
%     set(yl,'FontSize',18);
%     set(zl,'FontSize',18);
% %     ll = legend('Coast 1','','Burn 1','Coast 2','','Burn 2','Resulting Path');
%     ll = legend([plotreturn.D3.c,plotreturn.D3.b1,plotreturn.D3.b2,plotreturn.D3.nm],{'Coast','Burn 1','Burn 2','Result'});
%     set(ll,'FontSize',18,'Interpreter','LaTeX');
%     set(gca,'FontSize',16,'FontName','Times');
%     if problem.PCA.switch == 1
%         [xell, yell, zell] = ellipsoid(0,0,0,problem.PCA.ellipse(1)/1000,problem.PCA.ellipse(2)/1000,problem.PCA.ellipse(3)/1000,20);
%         ell = surf(yell, xell, zell);
%         ell.FaceAlpha = .5;
%     end
%     pause;
% end

%% Analyze Moon Conflicts

% Find vector from deputy to chief during region of interest
deltabeta = 2*(pi-acos(3*x_d/(2*ae)));
t_nomoon = deltabeta/w;
coast2Wtime = coasttime+1/2*T_p-1/2*t_nomoon;
JD_nomoon_start = JD_0+(problem.time.tf+coast2Wtime)/3600/24;
JD_nomoon_vec = linspace(JD_nomoon_start,JD_nomoon_start+t_nomoon/3600/24);
for i = 1:length(JD_nomoon_vec)
    RSO_nu = RSO_nu0 + w*(problem.time.tf+coast2Wtime) + w*(JD_nomoon_vec(i)-JD_nomoon_start)*3600*24; % deg % THIS ASSUMES CIRCULAR ORBIT
    [r_RSO,v_RSO] = coe2rvh(p,ecc,incl,omega,argp,RSO_nu,RSO_nu,RSO_nu,lonper,problem.mu);
%     r_RSO = [sma*cos(RSO_nu);sma*sin(RSO_nu);0];
%     v_RSO = [-sma*w*sin(RSO_nu);sma*w*cos(RSO_nu);0];
    r_RSO2moon_unit_nomoon = RSO2moon(JD_nomoon_vec(i),r_RSO,v_RSO);
    
    beta_nomoon = beta_new+(coast2Wtime)*w + w*(JD_nomoon_vec(i)-JD_nomoon_start)*3600*24;
    X_dep = LROE2X(ae,x_d,y_d0,z_max,gamma,beta_nomoon,w);
    X_dep_pos_unit_NEG = -X_dep(1:3)/norm(X_dep(1:3));
    
    alpha_M(i) = acos(dot(r_RSO2moon_unit_nomoon,X_dep_pos_unit_NEG));
end

output.Moon.p1x = (JD_nomoon_vec-JD_nomoon_start)*24;
output.Moon.p1y = rad2deg(alpha_M);
output.Moon.t_nomoon = t_nomoon;

%% Modify Outputs if Moon FOV Constraint is Enabled

if problem.traj.margin == 1
    
    beta_range = linspace(beta_lower,beta_upper);
    clear alpha_M
    for i = 1:length(beta_range)
        
        coast2Wtime = (pi-t_nomoon/2*w-beta_range(i))/w;
        for j = 1:length(JD_nomoon_vec)
            RSO_nu = RSO_nu0 + w*(problem.time.tf+coast2Wtime) + w*(JD_nomoon_vec(j)-JD_nomoon_start)*3600*24; % deg % THIS ASSUMES CIRCULAR ORBIT
            [r_RSO,v_RSO] = coe2rvh(p,ecc,incl,omega,argp,RSO_nu,RSO_nu,RSO_nu,lonper,problem.mu);
%             r_RSO = [sma*cos(RSO_nu);sma*sin(RSO_nu);0]; % Assumes equatorial
%             v_RSO = [-sma*w*sin(RSO_nu);sma*w*cos(RSO_nu);0];
            r_RSO2moon_unit_nomoon = RSO2moon(JD_nomoon_vec(j),r_RSO,v_RSO);
            
            beta_nomoon = beta_range(i) + (coast2Wtime)*w+w*(JD_nomoon_vec(j)-JD_nomoon_start)*3600*24;
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
            if beta_new < beta_lower
                beta_new = beta_lower;
            elseif beta_new > beta_upper
                beta_new = beta_upper;
            end

        end
        
    end
    
end
