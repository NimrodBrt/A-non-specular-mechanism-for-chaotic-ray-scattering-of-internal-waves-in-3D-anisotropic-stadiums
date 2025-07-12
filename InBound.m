function InBound = InBound(tYpe,r,tau)
    if strcmp(tYpe,'cone') == 1
        cone_specs
    elseif strcmp(tYpe,'channel') == 1
        channel_specs
    elseif strcmp(tYpe,'paraboloid') == 1
        paraboloid_specs
    elseif strcmp(tYpe,'par_channel') == 1
        par_channel_specs
    elseif strcmp(tYpe,'par_stadium') == 1
        [H,l,a] = par_stadium_specs(tau);
    elseif strcmp(tYpe,'stadium') == 1
        [s,H,fl,l] = stadium_specs(tau);
    end
%     specs(tYpe)
    x = r(1);
    y = r(2);
    z = r(3);
    InBound = true;
    if strcmp(tYpe,'cone') == 1
        if (x^2 + (y/b)^2 >= ((z + fl)/tanAlpha)^2) || (z<=0) || (z>=H)
            InBound = false;
        end
    elseif strcmp(tYpe,'channel') == 1
        if (z <= 0)|| (z >= H) || (y <= 0) || (y >= L) || (x <= 0) || (x >= (z + a*s)/s)
            InBound = false;
        end
    elseif strcmp(tYpe,'paraboloid') == 1
        if (a*(x^2 + b*y^2) >= z) || (z >= H)
            InBound = false;
        end
    elseif strcmp(tYpe, 'par_channel') == 1
        if (z <= a*x^2) || (z >= H) || (y <= 0) || (y >= L)
            InBound = false;
        end
    elseif strcmp(tYpe,'par_stadium') == 1
        if (z >= H) || ((abs(y) >= l/2) && (z <= a*(x^2 + (y - 0.5*l*sign(y))^2))) || (((abs(y) <= l/2)) && (z <= a*x^2))
            InBound = false;
        end
    elseif strcmp(tYpe,'stadium') == 1
        if (z >= H) || (z <= 0) || ((abs(y) >= l/2) && (z + fl <= s*(x^2 + (y - 0.5*l*sign(y))^2)^0.5)) || (((abs(y) <= l/2)) && (z + fl <= s*abs(x)))
            InBound = false;
        end
    end