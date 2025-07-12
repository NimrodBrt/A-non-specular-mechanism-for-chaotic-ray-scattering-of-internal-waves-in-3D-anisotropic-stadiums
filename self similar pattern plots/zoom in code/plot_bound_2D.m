function [] = plot_bound_2D(tYpe,gamma,tau)
    if strcmp(tYpe,'cone') == 1
        cone_specs
        A = (H + fl)/tanAlpha;
        X = linspace(-A,A,10000);
        
        a = fl/tanAlpha;
        x = linspace(-a,a,10000);
        
        
        plot(X,b*sqrt(A^2 - X.^2),'Color','k',LineWidth=1)
        hold on
        plot(X,-b*sqrt(A^2 - X.^2),'Color','k',LineWidth=1)
        plot(x,b*sqrt(a^2 - x.^2),'Color','k',LineWidth=1)
        plot(x,-b*sqrt(a^2 - x.^2),'Color','k',LineWidth=1)
        % plot(X,b*sqrt((tanAlpha^2 - gamma^2)/(gamma^2 - tanAlpha^2/b^2))*X,'Color','m',LineWidth=1)
        % plot(X,-b*sqrt((tanAlpha^2 - gamma^2)/(gamma^2 - tanAlpha^2/b^2))*X,'Color','m',LineWidth=1)
        grid on
        
        % title('ful')
        % legend('upper perimeter','','bottom perimeter','','critical line','','Location','north')
        
        % ylim([-2,2.5
    elseif strcmp(tYpe,'channel')
        channel_specs
        b = H/s;
        six_points = [0, 0; 0, L; a, L; a+b, L; a+b, 0; a, 0];
        sixX = six_points(:,1);
        sixY = six_points(:,2);
        plot(sixX(1:2), sixY(1:2), 'Color', 'k', LineWidth=1)
        hold on
        plot(sixX(2:4), sixY(2:4), 'Color', 'k', LineWidth=1)
        plot(sixX(4:5), sixY(4:5), 'Color', 'k', LineWidth=1)
        plot([sixX(1), sixX(5)], [sixY(1), sixY(5)], 'Color', 'k', LineWidth=1)
        plot([sixX(3), sixX(6)], [sixY(3), sixY(6)], '--', 'Color', 'k', LineWidth=1)
        axis equal
    elseif strcmp(tYpe,'par_channel')
        par_channel_specs
        xcrit = 0.5*gamma/a;
        plot([d,d],[0,L], 'Color', 'k', LineWidth=1)
        hold on
        plot(-[d,d],[0,L], 'Color', 'k', LineWidth=1)
        plot([0,0],[0,L], '-.', 'Color', 'k', LineWidth=1)
        plot([xcrit,xcrit],[0,L], '-.', 'Color', 'm', LineWidth=1)
        plot(-[xcrit,xcrit],[0,L], '-.', 'Color', 'm', LineWidth=1)
    elseif strcmp(tYpe,'paraboloid')
        paraboloid_specs
        dx = sqrt(H/a);
        MXC = 0.5*gamma/a;
        
        X = linspace(-dx,dx,1000);
        Y = sqrt((H/a - X.^2)/b);
        critx = linspace(-MXC,MXC,1000);
        crity = sqrt(0.25*gamma^2/a^2 - critx.^2)/b;

        plot(X,Y,'Color','k',LineWidth=1)
        hold on
        plot(X,-Y,'Color','k',LineWidth=1)
        plot(critx,crity,'--','Color','m',LineWidth=1)
        plot(critx,-crity,'--','Color','m',LineWidth=1)
    elseif strcmp(tYpe,'par_stadium') == 1
        par_stadium_specs
        d = sqrt(H/a);
        circN = 150;
        theta = linspace(0,pi,circN);
        xcrit = 0.5*gamma/a;

        plot([d,d],[-l/2,l/2], 'Color', 'k', LineWidth=1)
        hold on
        plot(-[d,d],[-l/2,l/2], 'Color', 'k', LineWidth=1)
        plot(d*cos(theta),d*sin(theta) + l/2, 'Color', 'k', LineWidth=1)
        plot(d*cos(theta),-d*sin(theta) - l/2, 'Color', 'k', LineWidth=1)
        plot([xcrit,xcrit],[-l/2,l/2],'--', 'Color', 'm', LineWidth=1)
        plot(-[xcrit,xcrit],[-l/2,l/2],'--', 'Color', 'm', LineWidth=1)
        plot(xcrit*cos(theta),xcrit*sin(theta) + l/2,'--', 'Color', 'm', LineWidth=1)
        plot(xcrit*cos(theta),-xcrit*sin(theta) - l/2,'--', 'Color', 'm', LineWidth=1)
    elseif strcmp(tYpe,'stadium') == 1
        [s,H,fl,l] = stadium_specs(tau);
        circN = 150;
        theta = linspace(0,pi,circN);
        d = (H+fl)/s;
        xcrit = fl/s;

        plot([d,d],[-l/2,l/2], 'Color', 'k', LineWidth=1)
        hold on
        plot(-[d,d],[-l/2,l/2], 'Color', 'k', LineWidth=1)
        plot(d*cos(theta),d*sin(theta) + l/2, 'Color', 'k', LineWidth=1)
        plot(d*cos(theta),-d*sin(theta) - l/2, 'Color', 'k', LineWidth=1)
        plot([xcrit,xcrit],[-l/2,l/2],'--', 'Color', 'k', LineWidth=1)
        plot(-[xcrit,xcrit],[-l/2,l/2],'--', 'Color', 'k', LineWidth=1)
        plot(xcrit*cos(theta),xcrit*sin(theta) + l/2,'--', 'Color', 'k', LineWidth=1)
        plot(xcrit*cos(theta),-xcrit*sin(theta) - l/2,'--', 'Color', 'k', LineWidth=1)
    end

end