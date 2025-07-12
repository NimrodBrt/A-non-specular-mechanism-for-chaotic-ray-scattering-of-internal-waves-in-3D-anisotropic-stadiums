% c = fix(clock);
% c(4:end)

tYpe = 'stadium';

tau = 0.9 + 4*0.025;
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
Njumps = 10000;
Nrays = 10000;

Nzoom = 5;
height = zeros(Nzoom,1);
left = zeros(Nzoom,1);
right = zeros(Nzoom,1);
width = zeros(Nzoom,1);

mean_conv_time = zeros(1,Nzoom);
delta = 1e-6;


% PHI0 = linspace(-pi,pi,Nrays);
% PHI0 = linspace(-0.42,-0.261369,Nrays); %for mu = 1.1
PHI0 = linspace(-0.5,-0.38359175,Nrays); %for mu = 1
% PHI0 = linspace(-0.545,-0.4234618,Nrays); %for mu = 0.975
% PHI0 = linspace(-0.565,-0.468563,Nrays); %for mu = 0.95
% PHI0 = linspace(-0.6,-0.52001,Nrays); %for mu = 0.925
% PHI0 = linspace(-0.65,-0.579238,Nrays); %for mu = 0.9

gamma = 1;
    % tophits = cell(Nrays,1);
    % bothits = cell(Nrays,1);
    % conehits = cell(Nrays,1);

% figure(1)
% TL = tiledlayout('vertical');
fntSZ = 12;

for t = 1:Nzoom
    t
    R = cell(Nrays,1);
    LR = cell(Nrays,1);
    LyapunovEXP = cell(Nrays,1);


    % tau = TAU(t)

    
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
        % n
        
        r0 = [0.7,0.1,H];
        phi0 = PHI0(n);

    
    
    
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
    Y = Y(ang_ind);

    [pk, lc] = findpeaks(Y,ANG,'MinPeakWidth',abs(PHI0(1) - PHI0(end))/17);
    [edges, edgsLC] = findpeaks(-Y,ANG,'MinPeakWidth',abs(PHI0(1) - PHI0(end))/100);
    
    height(t,1) = pk(1);
    left(t,1) = edgsLC(1);
    right(t,1) = edgsLC(2);
    width(t,1) = abs(edgsLC(1) - edgsLC(2));

    nexttile
    plot(ANG,Y,'.','MarkerSize',3.5)
    hold on
    % plot(PHI0,height(t,1)*ones(size(PHI0)))
    % plot([left(t),left(t)],[-0.5,0])
    % plot([right(t),right(t)],[-0.5,0])
    RedBoxBoundY = (-0.3 + 0.5*(0.3 + height(t,1)));
    RedBoxBoundX = right(t,1) - width(t,1)/4;
    if t == Nzoom
        xlabel('launching angle $\phi_0$',Interpreter='latex',FontSize= fntSZ)
    else
        plot(RedBoxBoundX*[1,1],[-0.5*l, RedBoxBoundY],'Color','r','LineWidth',0.8)
        plot([RedBoxBoundX,PHI0(end)],RedBoxBoundY*[1,1],'Color','r','LineWidth',0.8)
    end
    grid on
    xlim([PHI0(1),PHI0(end)])
    yupperlim = -0.3 + 2*(0.3 + height(t,1));
    ylim([-0.3, yupperlim])
    ylabel('$y_\infty$',Interpreter='latex',FontSize= fntSZ,Rotation=0)
    hold off

    % phi0L = (right(t,1) + lc(1))/2;
    phi0L = RedBoxBoundX;
    % phi0R = right(t,1) + width(t,1)/2;
    phi0R = PHI0(end);
    PHI0 = linspace(phi0L,phi0R,Nrays);
end
toc
title(TL,'$W = 2\sqrt{2} \qquad \tau = 1 \qquad l = 0.6 \qquad \mu = 1$',Interpreter='latex',FontSize=fntSZ + 2,FontWeight='bold')

% logWidth = log10(width);
% logHeight = log10(height + 0.3);
% 
% regX = [ones(Nzoom,1),logWidth];
% regMatrix = regX\logHeight;
% logHeightCalc = regX*regMatrix;
% Rsq = 1 - sum((logHeight - logHeightCalc).^2)/sum((logHeight - mean(logHeight)).^2);
% RsqSTR = num2str(Rsq);
% 
% % regX = [ones(Nzoom,1),width];
% % regMatrix = regX\(height + 0.3);
% % HeightCalc = regX*regMatrix;
% % Rsq = 1 - sum((logHeight - HeightCalc).^2)/sum((logHeight - mean(logHeight)).^2);
% % RsqSTR = num2str(Rsq);
% 
% 
% 
% figure(2)
% TL2 = tiledlayout('horizontal');
% 
% nexttile([1,3])
% plot(logWidth,logHeight,'o')
% hold on
% plot(logWidth,logHeightCalc,'--')
% hold off
% grid on 
% xlabel('$\log_{10} (peak \, width)$',Interpreter='latex',FontSize=fntSZ)
% ylabel('$\log_{10} (peak \, height)$',Interpreter='latex',FontSize=fntSZ)
% 
% title('linear regression',Interpreter='latex',FontSize=fntSZ)
% lgd = legend('simulated data',strcat('linear regeression $y = ', num2str(regMatrix(2)), '\cdot x + ', num2str(regMatrix(1)), '$'));
% lgd.Location = "northwest";
% lgd.Interpreter = "latex";
% dim = [.2 .5 .3 .3];
% dim = [.2 .5 .25 .25];
% str = strcat('R^2 = ',' ', RsqSTR);
% annotation('textbox',dim,'String',str,'FitBoxToText','on',FontSize=7.5);

% title('$-0.385<\phi_0<-0.383$',Interpreter='latex',FontSize= fntSZ-2)
% grid on
% ylim(0.5*l*[-1,1])

% feigen = width(1:end-1)./width(2:end);
% % figure(3)
% nexttile([1,2])
% plot(feigen,'-*','MarkerSize',8,'LineWidth',1)
% xlabel('zoom in number',Interpreter='latex',FontSize=fntSZ)
% ylabel('$\frac{width_i}{width_{i+1}}$',Interpreter='latex',FontSize=fntSZ+4)
% title('Feigenbaum convergence',Interpreter='latex',FontSize=fntSZ)
% grid on
% hold on
% 
% 
% title(TL2,'$W = 2\sqrt{2} \qquad \tau = 1 \qquad l = 0.6 \qquad \mu = 1$',Interpreter='latex',FontSize=fntSZ + 2,FontWeight='bold')
% plot([1,Nzoom-1],4.669201609102*[1,1],'--','Color','r','LineWidth',1.2)
% hold off

% figure(4)%3D convergence plot
% plot(TAU,mean_conv_time)  
% grid on
% fntSZ = 12;
% xlabel('$\mu$',Interpreter='latex',FontSize= fntSZ)
% ylabel('mean 3D convergence time',Interpreter='latex',FontSize= fntSZ)
% title('3D convergence plot $\qquad W = 2.5 \qquad \tau = 1 \qquad l=0.6$',Interpreter='latex',FontSize= fntSZ)
% xlim([TAU(1),TAU(end)])
% ylim([0,Njumps])

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