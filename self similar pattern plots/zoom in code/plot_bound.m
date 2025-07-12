function plot_bound(tYpe,gamma,tau)
    if strcmp(tYpe,'cone') == 1
        cone_specs
        xmax = (H+fl)/tanAlpha;
        ymax = b*xmax;
        [X,Y] = meshgrid(-xmax:0.005:xmax,-ymax:0.005:ymax);
        Z = tanAlpha*(X.^2 + Y.^2/b^2).^0.5 - fl;
%         axis equal
        zlim([0,H])
        surf(X, Y, Z,EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
    elseif strcmp(tYpe,'channel') == 1
        channel_specs
        b = H/s;
        six_points = [0, 0; 0, L; a, L; a+b, L; a+b, 0; a, 0];
        sixX = six_points(:,1);
        sixY = six_points(:,2);
        six_points_bot = [0, 0, 0; 0, L, 0; a, L, 0; a+b, L, 0; a+b, 0, 0; a, 0, 0];
        six_points_top = [0, 0, H; 0, L, H; a, L, H; a+b, L, H; a+b, 0, H; a, 0, H];
        [Xflat,Yflat] = meshgrid([sixX(1),sixX(6)],sixY(1:2));
        [Xtilt,Ytilt] = meshgrid([sixX(6),sixX(5)],[sixY(6),sixY(3)]);
        [Ywall,Zwall] = meshgrid([0,L],[0,H]);
        surf(Xflat, Yflat, 0*Xflat,EdgeColor="none",FaceAlpha=0.15,FaceColor='k')
        hold on
        surf(Xtilt,Ytilt,s*(Xtilt - a),EdgeColor="none",FaceAlpha=0.15,FaceColor='k')
        surf(0*Ywall, Ywall, Zwall,EdgeColor="none",FaceAlpha=0.15,FaceColor='k')
    elseif strcmp(tYpe,'par_channel') == 1
        par_channel_specs
        xcrit = 0.5*gamma/a;
        zcrit = a*xcrit^2;
        [X,Y] = meshgrid(-d:0.005:d,0:0.005:L);
        Z = a*X.^2;
        surf(X, Y, Z,EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        hold on
        plot3([xcrit,xcrit],[0,L],[zcrit,zcrit], '--', 'Color', 'm', LineWidth=1.5)
        plot3(-[xcrit,xcrit],[0,L],[zcrit,zcrit], '--', 'Color', 'm', LineWidth=1.5)
        plot3([0,0],[0,L],[0,0], '--', 'Color', 'k', LineWidth=1.5)
    elseif strcmp(tYpe,'paraboloid') == 1
        paraboloid_specs
        dx = sqrt(H/a);
        dy = sqrt(H/(a*b));
        [X,Y] = meshgrid(-dx:0.005:dx,-dy:0.005:dy);
        Z = a*(X.^2 + b*Y.^2);
        MXC = 0.5*gamma/a;
        critx = linspace(-MXC,MXC,1000);
        crity = sqrt(0.25*gamma^2/a^2 - critx.^2)/b;
        critz = a*(critx.^2 + b*crity.^2);
        surf(X, Y, Z,EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        hold on
        plot3(critx,crity,critz,'--','Color','m', LineWidth=1.5)
        plot3(critx,-crity,critz,'--','Color','m', LineWidth=1.5)
    elseif strcmp(tYpe,'par_stadium') == 1
        [H,l,a] = par_stadium_specs(tau);
        d = sqrt(H/a);
        [Xchan, Ychan] = meshgrid(-d:0.005:d,-l/2:0.005:l/2);
        Zchan = a*Xchan.^2;
        [XcircR,YcircR] = meshgrid(-d:0.005:d,l/2:0.005:(l/2+d));
        ZcircR = a*(XcircR.^2 + (YcircR - l/2).^2);
        [XcircL,YcircL] = meshgrid(-d:0.005:d,-(l/2+d):0.005:-l/2);
        ZcircL = a*(XcircL.^2 + (YcircL + l/2).^2);

        xcrit = 0.5*gamma/a;
        zcrit = a*xcrit^2;
        critcircN = 150;
        theta = linspace(0,pi,critcircN);

        surf(Xchan, Ychan, Zchan, EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        hold on
        surf(XcircR, YcircR, ZcircR, EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        surf(XcircL, YcircL, ZcircL ,EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        plot3([xcrit,xcrit],[-l/2,l/2],[zcrit,zcrit], '--', 'Color', 'm', LineWidth=1.5)
        plot3(-[xcrit,xcrit],[-l/2,l/2],[zcrit,zcrit], '--', 'Color', 'm', LineWidth=1.5)
        plot3(xcrit*cos(theta),xcrit*sin(theta) + l/2,zcrit*ones(1,critcircN), '--', 'Color', 'm', LineWidth=1.5)
        plot3(xcrit*cos(theta),-xcrit*sin(theta) - l/2,zcrit*ones(1,critcircN), '--', 'Color', 'm', LineWidth=1.5)
    elseif strcmp(tYpe,'stadium') == 1
        [s,H,fl,l] = stadium_specs(tau);
        d = s*(H+fl);
        [Xchan, Ychan] = meshgrid(-d:0.005:d,-l/2:0.005:l/2);
        Zchan = s*abs(Xchan) - fl;
        [XcircR,YcircR] = meshgrid(-d:0.005:d,l/2:0.005:(l/2+d));
        ZcircR = s*(XcircR.^2 + (YcircR - l/2).^2).^0.5 - fl;
        [XcircL,YcircL] = meshgrid(-d:0.005:d,-(l/2+d):0.005:-l/2);
        ZcircL = s*(XcircL.^2 + (YcircL + l/2).^2).^0.5 - fl;

        surf(Xchan, Ychan, Zchan, EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        hold on
        surf(XcircR, YcircR, ZcircR, EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
        surf(XcircL, YcircL, ZcircL ,EdgeColor = 'none',FaceAlpha=0.15,FaceColor='k')
    end
        
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

