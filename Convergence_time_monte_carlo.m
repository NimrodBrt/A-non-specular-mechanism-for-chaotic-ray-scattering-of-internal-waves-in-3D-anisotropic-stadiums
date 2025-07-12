c = fix(clock);
c(4:end)

tYpe = 'stadium';

tau = 1;
tic
if strcmp(tYpe,'cone') == 1
    cone_specs
elseif strcmp(tYpe,'channel') == 1
    channel_specs
elseif strcmp(tYpe,'par_channel') == 1
    par_channel_specs
elseif strcmp(tYpe,'paraboloid') == 1
    paraboloid_specs
elseif strcmp(tYpe,'par_stadium') == 1
    [H,l,a] = par_stadium_specs(tau);
elseif strcmp(tYpe,'stadium') == 1
    [s,H,fl,l] = stadium_specs(tau);
else
    disp('Non valid tYpe')
end


dt = 0.1;
Njumps = 8000;
Nrays = 100;

Ntau = 100;
TAU = linspace(pi/4 + 0.001,pi/2 - 0.001,Ntau);
mean_conv_time = zeros(1,Ntau);
% delta = 1e-6;


gamma = 1;

for t = 1:Ntau
    R = cell(Nrays,1);
    LR = cell(Nrays,1);


    tau = TAU(t)

    
    % tophits = cell(Nrays,1);
    % bothits = cell(Nrays,1);
    % conehits = cell(Nrays,1);
    
    attractor_length = 20;
        attractorR = cell(Nrays,1);
        attractorLR = cell(Nrays,1);
    %     attractorTH = cell(Nrays,1);
    %     attractorBH = cell(Nrays,1);
    %     attractorCH = cell(Nrays,1);
    % end
    Y = zeros(1,Nrays);
    parfor n = 1:Nrays
        
        r0 = [100,100,100];
        
        phi0 = PHI0(n);


            while InBound(tYpe,r0,tau) == 0
                r0 = [2*(2*rand - 1), 2*(2*rand - 1),0.999*H];
                phi0 = 2*pi*rand;

            end

        % deltaVec = delta*[sin(phi0),-cos(phi0),0];
        % Lr0 = r0 + deltaVec;
            % r0 = [0.3,0,1];
            % phi0 = pi/3;
        % end
        % end
    
    
    
        v0 = [cos(phi0)/sqrt(1 + gamma^2), sin(phi0)/sqrt(1 + gamma^2), -gamma/sqrt(1 + gamma^2)];
        [R{n},] = singleRay(tYpe,r0,v0,dt,Njumps,gamma,tau);
        Y(n) = R{n}(end,2);
        % [LR{n},] = singleRay(tYpe,Lr0,v0,dt,Njumps,gamma,tau);

        % tophits{n} = R{n}(R{n}(:,3) > H - 0.0001,:);
        % bothits{n} = R{n}(R{n}(:,3) < 0.0001,:);
        % conehits{n} = R{n}(R{n}(:,3) > 0.0001,:);
        % conehits{n} = conehits{n}(conehits{n}(:,3) < H - 0.0001,:);
        if attractor_length < size(R{n},1)
            attractorR{n} = R{n}(end-attractor_length:end,:);
            % attractorLR{n} = LR{n}(end-attractor_length:end,:);

            % attractorTH{n} = attractorR{n}(attractorR{n}(:,3) > H - 0.0001,:);
            % attractorBH{n} = attractorR{n}(attractorR{n}(:,3) < 0.0001,:);
            % attractorCH{n} = attractorR{n}(attractorR{n}(:,3) > 0.0001,:);
            % attractorCH{n} = attractorCH{n}(attractorCH{n}(:,3) < H - 0.0001,:);
        end
    
    end
    conv_times = cell2mat(cellfun(@size,R,'UniformOutput',false));
    CONV_TIMES = conv_times(:,1);

    ang_ind = find(CONV_TIMES < 7900);

    CONV_TIMES = CONV_TIMES(ang_ind);
    ANG = PHI0(ang_ind);
    mean_conv_time(t) = mean(conv_times(:,1));
end
toc
Y = Y(ang_ind);


% % figure(3)  % y_infty plot
% fntSZ = 12;
% plot(ANG,Y,'.','MarkerSize',3.5)
% hold on
% plot(-[0.5,0.5,0.383588,0.383588,0.5],-0.3*[1,0.5,0.5,1,1],'Color','r','LineWidth',0.8)
% % xlabel('launching angle $\phi_0$',Interpreter='latex',FontSize= fntSZ-1)
% ylabel('$y_\infty$',Interpreter='latex',FontSize= fntSZ,Rotation=0)
% hold off
% xlim([PHI0(1),PHI0(end)])
% % title(strcat('$\mu = ', num2str(tau),'$'),Interpreter='latex',FontSize= fntSZ-2)
% grid on
% ylim(0.5*l*[-1,1])
% hold off

figure(4) %3D convergence time plot
plot(TAU,mean_conv_time)  
grid on
fntSZ = 12;
xlabel('$\mu$',Interpreter='latex',FontSize= fntSZ)
ylabel('mean 3D convergence time',Interpreter='latex',FontSize= fntSZ)
title('3D convergence plot $\qquad W = 2 \sqrt{2} \qquad \tau = 1 \qquad l=0.6$',Interpreter='latex',FontSize= fntSZ)
xlim([TAU(1),TAU(end)])
ylim([0,Njumps])

% fig2 = figure(2);  %full paths
% % nexttile
% plot_bound(tYpe,gamma,tau)
% grid on
% hold on
% title('full path')
% for n = 1:Nrays
%     plot3(R{n}(:,1), R{n}(:,2), R{n}(:,3),'LineWidth',1.2)
%     plot3(LR{n}(:,1), LR{n}(:,2), LR{n}(:,3),'LineWidth',1.2)
%     scatter3(R{n}(1,1), R{n}(1,2), R{n}(1,3),'red','filled')
% end
% axis equal
% zlim([0,H])
% hold off
% xlim([-1.5,2.5])
% ylim([-3,3])
% % 
% fig1 = figure(1); %full paths from above
% plot_bound_2D(tYpe,gamma,tau)
% for n = 1:Nrays
%     plot(R{n}(:,1), R{n}(:,2),'LineWidth',1)
%     scatter(tophits{n}(:,1),tophits{n}(:,2),'o','filled','cyan')
%     scatter(bothits{n}(:,1),bothits{n}(:,2),'o','filled','magenta')
%     scatter(conehits{n}(:,1),conehits{n}(:,2),'o','filled','green')
% end
% scatter(0,0,'o','filled','k')
% grid on
% axis equal
% xlabel('x')
% ylabel('y')
% title('full paths from above')
% hold off
% % 
% 
% fig6 = figure(6); %single ray convergence plot
% % nexttile
% max_conv_time = max(cell2mat(cellfun(@size,R,'UniformOutput',false)),[],"all");
% 
% plot(nan,nan)
% plot([0,max_conv_time],0.5*[l,l],'LineStyle','--','Color','k')
% hold on
% plot([0,max_conv_time],-0.5*[l,l],'LineStyle','--','Color','k')
% for n = 1:Nrays
%     plot(R{n}(:,2))
%     % plot(LR{n}(:,2))
% end
% grid on
% xlabel('step #')
% ylabel('y')
% title('convergence plot')
% hold off
% % if attractor_length < min(conv_times)
% % 
    % fig3 = figure(3); %attractors from above
    % plot_bound_2D(tYpe,gamma,tau)
    % for n = 1:Nrays
    %     plot(attractorR{n}(:,1),attractorR{n}(:,2),'LineWidth',1.2)
    %     % plot(attractorLR{n}(:,1),attractorLR{n}(:,2),'LineWidth',1.2)
    %     % % scatter(attractorTH{n}(:,1),attractorTH{n}(:,2),'o','filled','cyan')
    %     % % scatter(attractorBH{n}(:,1),attractorBH{n}(:,2),'o','filled','magenta')
    %     % % scatter(attractorCH{n}(:,1),attractorCH{n}(:,2),'o','filled','green')
    % end
    % grid on
    % axis equal
    % xlabel('x')
    % ylabel('y')
    % title('attractors from above')
    % hold off
% % 
% %     % 
% %     % fig5 = figure(5);  %attractors only
% %     % plot_bound(tYpe,gamma,tau)
% %     % grid on
% %     % hold on
% %     % title('attractors only')
% %     % for n = 1:Nrays
% %     %     plot3(attractorR{n}(:,1),attractorR{n}(:,2),attractorR{n}(:,3),'LineWidth',1.2)
% %     % end
% %     % 
% %     % hold off
% %     % axis equal
% %     % zlim([0,H])
% % end