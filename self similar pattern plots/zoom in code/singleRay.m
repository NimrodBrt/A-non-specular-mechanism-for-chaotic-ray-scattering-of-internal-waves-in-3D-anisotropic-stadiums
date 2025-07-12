function [singleRay,vend] = singleRay(tYpe,r0,v0,dt,Njumps,gamma,tau)

    conv_window = 10;

    r = zeros(Njumps,3);
    R = r0;
    r(1,:) = R;
    v = v0;
    R = R + v0*dt;
    for k = 2:Njumps
        while InBound(tYpe,R,tau) == 1
            R = R + v*dt;
        end
        [R,v] = reflect(tYpe,R - v*dt,R,v,gamma,tau);
        r(k,:) = R;
        R = R + v*dt;
        % SRtest = InBound(tYpe,R)
        if (k > conv_window + 1) && (abs(min(r(k-conv_window:k,2) - max(r(k-conv_window:k,2)))) < 10^-8)
            break
        end
    end
    singleRay = r(1:k,:);
    vend = v;
end