function [alpha, delta, lon, lat] = groundTrack(Y, lambda0, Tperiod, omega_E, N) %Y = (r,v)

L = length(Y(:,1));

% Creating lists to avoid changing the size in the for loops
r = zeros(L,1);
delta = zeros(L*floor(N),1);
alpha = zeros(L*floor(N),1);
lon = zeros(L*floor(N),1);
Lon = zeros(L*floor(N),1);
lat = zeros(L*floor(N),1);

deltaLambda = Tperiod * omega_E; 

for k = 0:floor(N)-1
    for i = 1:L
      r(i) = sqrt(Y(i,1)^2 + Y(i,2)^2 + Y(i,3)^2);
      n = Y(i,3)/r(i); %rad
      delta(i + k*L) = asin(n); %start again shifting the groundtrack
      m = Y(i,2)/r(i);
      l = Y(i,1)/r(i);
      if m>0
         alpha(i + k*L) = acos(l/cos(delta(i)));
      else
         alpha(i + k*L) = 2*pi() - acos(l/cos(delta(i)));
      end
      Lon(i + k*L) = alpha(i)/pi()*180 + lambda0 - k*deltaLambda ; %- k*2*pi() %longitude = lambda   %lambda0 = lon of Greenwhich meridian
    end
    [~,imini] = min(Lon(1+k*L : L+k*L));
    lon_k = [Lon(imini +k*L : L +k*L) ; Lon(1 +k*L : imini-1 +k*L)];
    Q = floor(abs(Lon(1 +k*L)/360));
    lon(1 +k*L : L +k*L) = lon_k - sign(Lon(1 +k*L))*Q*360;
    lat_k = [delta(imini +k*L : L +k*L)*180/pi() ; delta(1 +k*L : imini-1 +k*L)*180/pi()];
    lat(1 +k*L : L +k*L) = lat_k;
    
    
    % Shifting the groundtracks to get the longitude between -180 and 180
    if lon(1 +k*L)<-180
        j=1;
        while lon(j +k*L)<-180
            j = j+1;
            if j==L
                break
            end
        end
        lon(1 +k*L : L +k*L) = [lon(j +k*L : L+k*L) ; lon(1 +k*L : j-1 +k*L)+360];
        lat(1 +k*L : L +k*L) = [lat(j +k*L : L+k*L) ; lat(1 +k*L : j-1 +k*L)];

    elseif lon(L +k*L)<-180
        j=1;
        while lon(j +k*L)>-180
            j = j+1;
            if j==length(lon(1:L +k*L))
                break
            end
        end
        lon(1 +k*L : L +k*L) = [lon(j +k*L : L+k*L)+360 ; lon(1 +k*L : j-1 +k*L)];
        lat(1 +k*L : L +k*L) = [lat(j +k*L : L+k*L) ; lat(1 +k*L : j-1 +k*L)];
    
    elseif lon(1 +k*L)>180
        j=1;
        while lon(j +k*L)>180
            j = j+1;
            if j==length(lon(1:L +k*L))
                break
            end
        end
        lon(1 +k*L : L +k*L) = [lon(j +k*L : L+k*L) ; lon(1 +k*L : j-1 +k*L)-360];
        lat(1 +k*L : L +k*L) = [lat(j +k*L : L+k*L) ; lat(1 +k*L : j-1 +k*L)];
        
    elseif lon(L +k*L)>180 
        j = 1;
        while lon(j +k*L)<180
            j = j+1;
            if j==length(lon(1:L +k*L))
                break
            end
        end
        lon(1 +k*L : L +k*L) = [lon(j +k*L : L+k*L)-360 ; lon(1 +k*L : j-1 +k*L)];
        lat(1 +k*L : L +k*L) = [lat(j +k*L : L+k*L) ; lat(1 +k*L : j-1 +k*L)];
    end
end

