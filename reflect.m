function [rBound,vRef] = reflect(tYpe,rIn,rOut,v,gamma,tau)

    if strcmp(tYpe,'cone') == 1
        cone_specs
    elseif strcmp(tYpe,'channel') == 1
        channel_specs
    elseif strcmp(tYpe,'paraboloid') == 1
        paraboloid_specs
    elseif strcmp(tYpe,'par_channel')
        par_channel_specs
    elseif strcmp(tYpe,'par_stadium') == 1
    [H,l,a] = par_stadium_specs(tau);
    elseif strcmp(tYpe,'stadium')
        [s,H,fl,l] = stadium_specs(tau);
    end

    rMid = (rIn + rOut)/2;
        while norm(rOut - rIn) > 1e-14
            % log10(norm(rOut - rIn))
            % rMid
            % InBound(tYpe,rIn)
            if InBound(tYpe,rMid,tau) == 1
                rIn = rMid;
                rMid = (rIn + rOut)/2;
            elseif InBound(tYpe,rMid,tau) == 0
                rOut = rMid;
                rMid = (rIn + rOut)/2;
            end
            if strcmp(tYpe,'channel')
                if (rMid(3) == 0)|| (rMid(3) == H) || (rMid(2) == 0) || (rMid(2) == L) || (rMid(1) == 0) || (rMid(1) == (rMid(3) + a*s)/s)
                    break
                end
            end
        end
    rBound = rMid; 

    x = rBound(1);
    y = rBound(2);
    z = rBound(3);
    wall = 0;
    if (rOut(3) >= H) || (rOut(3) <= 0) %floor or ceiling
        vRef = [v(1),v(2),-v(3)];
    else
        if strcmp(tYpe,'cone')
            gradZ = (tanAlpha^2/(z + fl))*[x,y/b^2];
            s = norm(gradZ);
        elseif strcmp(tYpe,'channel')
            channel_specs
            if (rOut(1) <= 0) || (rOut(2) <= 0) || (rOut(2) >= L)
                wall = 1;
            end
            if (rOut(2) <= 0) || (rOut(2) >= L)
                gradZ = sign(y - L/2)*[0,1];
            else
                gradZ = sign(x - a/2)*[1,0];
            end
        elseif strcmp(tYpe,'par_channel')
            par_channel_specs
            if (rOut(2) <= 0) || (rOut(2) >= L)
                wall = 1;
                gradZ = sign(y - L/2)*[0,1];
            else
                gradZ = 2*a*x*[1,0];
                s = abs(2*a*x);
            end
        elseif strcmp(tYpe,'paraboloid')
            paraboloid_specs
            gradZ = 2*a*[x,b*y];
            s = norm(gradZ);
        elseif strcmp(tYpe,'par_stadium') == 1
            [H,l,a] = par_stadium_specs(tau);
            if abs(y) < l/2
                gradZ = 2*a*x*[1,0];
                s = abs(2*a*x);
            else
                gradZ = 2*a*[x,y - 0.5*l*sign(y)];
                s = norm(gradZ);
            end
        elseif strcmp(tYpe,'stadium')
            [s,H,fl,l] = stadium_specs(tau);
            if abs(y) < l/2
                gradZ = s*sign(x)*[1,0];
            else
                gradZ = s*[x,y - 0.5*l*sign(y)];
            end
        end

        
        if wall == 1
            h = 0; % (s -> inf) ==> (h = 1/(gamma*s) -> 0)
            s = gamma + 1; % this way we make sure sign(gamma - s) = -1, which we want for a supercritical reflection
        else
            h = gamma/s;
        end

        gradZ = real(gradZ);
        phiGZ = atan2(gradZ(2),gradZ(1));
        v = real(v);
        phiIn = atan2(v(2),v(1)) - phiGZ; %relative to gradH
        if phiIn > pi
            phiIn = phiIn - 2*pi;
        end
        phiR = sign(phiIn)*acos(-sign(1-h)*((1+h^2)*cos(phiIn) - sign(v(3))*2*h)/((1+h^2) - sign(v(3))*2*h*cos(phiIn)));
        % phiR = acos(-sign(1-h)*((1+h^2)*cos(phiIn) - sign(v(3))*2*h)/((1+h^2) - sign(v(3))*2*h*cos(phiIn)));
        % if abs(phiIn) > acos(-2*h/(1+h^2))
        %     phiR = -phiR;
        % end



        % phiR  = pi - asin(sign(1-h)*(1-h^2)*sin(phiIn)/((1+h^2) - sign(v(3))*2*h*cos(phiIn))); %relative to gradH
        % if (h < 1) %can probably be shortenned by a few lines
        %     abPhiI = abs(phiIn);
        %     if abPhiI > acos(sign(v(3))*2*h/(1+h^2))
        %         phiR  = asin(sign(1-h)*(1-h^2)*sin(phiIn)/((1+h^2) - sign(v(3))*2*h*cos(phiIn))); %relative to gradH
        %     end
        % elseif (h > 1)
        %     abPhiI = abs(phiIn);
        %     % acos(-2*h/(1+h^2));
        %     if abPhiI < acos(sign(v(3))*2*h/(1+h^2))
        %         phiR  = asin(sign(1-h)*(1-h^2)*sin(phiIn)/((1+h^2) - sign(v(3))*2*h*cos(phiIn))); %relative to gradH
        %     end
        % end

        if wall == 1
            phiR  = pi - asin(sign(1-h)*(1-h^2)*sin(phiIn)/((1+h^2) - sign(v(3))*2*h*cos(phiIn))); %relative to gradH
        end
        
        if abs(h - 1) < 0.0001
            phiR = pi - phiIn;
        end
        phiR = phiR + phiGZ; %relative to x axis

        vRef = [cos(phiR)/sqrt(1+gamma^2), sin(phiR)/sqrt(1+gamma^2), -sign(gamma - s)*v(3)];
        % if (vRef(3) == -v(3)) && (vRef(1)*v(1) + vRef(2)*v(2) < 0)
        %     vRef(1) = -vRef(1);
        %     vRef(2) = -vRef(2);
        % end
    end
end