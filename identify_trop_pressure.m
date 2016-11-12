function trp = identify_trop_pressure(t, p)

% from top to bottom
% http://www.inscc.utah.edu/~reichler/research/projects/TROPO/code.txt

% implicit none
% integer,intent(in)                  :: level
% real,intent(in),dimension(level)    :: t, p
% real,intent(in)                     :: plimu, pliml, gamma
% real,intent(out)                    :: trp
%
% real,parameter                      :: kap=0.286
% real,parameter                      :: faktor = -9.81/287.0
% real,parameter                      :: deltaz = 2000.0
% real,parameter                      :: ka1=kap-1.
%
% real                                :: pmk, pm, a, b, tm, dtdp, dtdz
% real                                :: ag, bg, ptph
% real                                :: pm0, pmk0, dtdz0
% real                                :: p2km, asum, aquer
% real                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
% integer                             :: icount, jj
% integer                             :: j

kap    = 0.286; 
faktor = -9.81/287.0;
deltaz = 2000.0; 
ka1    = kap-1.;

gamma = -0.002;  % K/m
plimu = 45000;   % upper limit for tropopause pressure in Pa
pliml = 7500;    % lower limit for tropopause pressure in Pa
  
level = length(p);

dtdz0=0;
trp=-99.0;                           % negative means not valid
ptph = -99;

next88=0;
next_l=0;

for j=level:-1:2
    
    % dt/dz
    pmk= .5 * (p(j-1).^kap+p(j).^kap);
    pm = pmk.^(1/kap);
    a = (t(j-1)-t(j))/(p(j-1).^kap-p(j).^kap);
    b = t(j)-(a*p(j).^kap);
    %tm = a * pmk + b   ;
    %tm - (t(j)+t(j-1))/2
    tm = (t(j)+t(j-1))/2 ;
    dtdp = a * kap * (pm.^ka1);
    dtdz = faktor*dtdp*pm/tm;
    
    % dt/dz valid?
    if (j==level);    next_l=1;  end % no, start level, initialize first
    if (dtdz<gamma);  next_l=1;  end %go to 999     % no, dt/dz < -2 K/km
    if (pm>plimu);    next_l=1;  end %   go to 999     % no, too low
    
    % dtdz is valid, calculate tropopause pressure
    if(next_l ~=1)
        if (dtdz0 < gamma)
            ag = (dtdz-dtdz0) / (pmk-pmk0)   ;
            bg = dtdz0 - (ag * pmk0)        ;
            %ptph1 = exp(log((gamma-bg)/ag)/kap);
            pi_tmp = (gamma-dtdz0)/ag+pmk0;
            ptph = exp(log(pi_tmp)/kap);
            %[ptph1-ptph]/ptph            
        else
            ptph = pm;
        end
    end
    
    %if (ptph < pliml); next_l=1;  end % go to 999
    %if (ptph > plimu); next_l=1;  end % go to 999
     
    %[j p(j)/1e2 dtdz*1e3 t(j)];
    
    if (ptph < pliml || ptph > plimu)
        
        pm0 = pm;
        pmk0 = pmk;
        dtdz0  = dtdz;
        next_l=0;
        continue;
    else
        
        
        
        % 2nd test: dtdz above 2 km must not exceed gamma
        p2km = ptph + deltaz*(pm/tm)*faktor ;         % p at ptph + 2km
        asum = 0.0    ;                               % dtdz above
        icount = 0    ;                               % number of levels above
        
        % test until apm < p2km
        for jj=j:-1:2
            
            pmk2 = .5 * (p(jj-1).^kap+p(jj).^kap);    % p mean ^kappa
            pm2 = pmk2.^(1/kap)      ;                % p mean
            if(pm2 > ptph); continue; end;            % go to 110                % doesn't happen
            if(pm2 < p2km); next88 = 1; break; end;   % go to 888                % ptropo is valid
            
            a2 = (t(jj-1)-t(jj))  ;                   % a
            a2 = a2/(p(jj-1).^kap-p(jj).^kap);
            b2 = t(jj)-(a2*p(jj).^kap) ;              % b
            tm2 = a2 * pmk2 + b2 ;                    % T mean
            dtdp2 = a2 * kap * (pm2.^(kap-1)) ;       % dt/dp
            dtdz2 = faktor*dtdp2*pm2/tm2;
            asum = asum+dtdz2;
            icount = icount+1;
            aquer = asum/(icount)  ;                  % dt/dz mean
            
            % discard ptropo ?
            if (aquer <= gamma );
                next_l=1;  break; 
            end           % go to 999; dt/dz above < gamma
            
            % 110 continue
        end                           % test next level
        
        
        if(next88==1) %888 continue                        % ptph is valid
            trp = ptph;            
            return
        else
            next88=0;
            % continue
        end
        
    end
    
    % 999 continue                        % continue search at next higher level
    if( next_l ==1)
        pm0 = pm;
        pmk0 = pmk;
        dtdz0  = dtdz;
        next_l=0;
    end
    
end

% no tropopouse found
return

