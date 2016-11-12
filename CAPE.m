function [cape, cin, plcl, tlcl] =  cape_bryan(  p_in , t_in , td_in, adiabat)
% [cape, cin, lcl, nbl] =  getcape(  p_in , t_in , td_in)
%
%  Input:     nk - number of levels in the sounding (integer)
%
%           p_in - one-dimensional array of pressure (mb) (real)
%
%           t_in - one-dimensional array of temperature (C) (real)
%
%          td_in - one-dimensional array of dewpoint temperature (C) (real)
%
%  Output:  cape - Convective Available Potential Energy (J/kg) (real)
%
%            cin - Convective Inhibition (J/kg) (real)
%
%-----------------------------------------------------------------------
%
%  getcape - a fortran90 subroutine to calculate Convective Available
%            Potential Energy (CAPE) from a sounding.
%
%  Version 1.02                           Last modified:  10 October 2008
%
%  Author:  George H. Bryan
%           Mesoscale and Microscale Meteorology Division
%           National Center for Atmospheric Research
%           Boulder, Colorado, USA
%           gbryan@ucar.edu
%
%  Disclaimer:  This code is made available WITHOUT WARRANTY.
%
%  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
%               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
%               Equation (38)
%
%-----------------------------------------------------------------------
%  User options:
nk = numel(p_in);

pinc = 100.0;   % Pressure increment (Pa)
pinc = 100.0;   % Pressure increment (Pa)

% (smaller number yields more accurate
%  results,larger number makes code
%  go faster)

source = 1;    % Source parcel:
% 1 = surface
% 2 = most unstable (max theta-e)
% 3 = mixed-layer (specify ml_depth)

ml_depth =  200.0 ; % depth (m) of mixed layer
% for source=3

if(nargin<4)
    adiabat = 3;   % Formulation of moist adiabat:
    % 1 = pseudoadiabatic, liquid only
    % 2 = reversible, liquid only
    % 3 = pseudoadiabatic, with ice
    % 4 = reversible, with ice
end

%-----------------------------------------------------------------------
%            No need to modify anything below here:
%-----------------------------------------------------------------------
%
% logical :: doit,ice,cloud,not_converged
% integer :: k,kmax,n,nloop,i,orec
% real, dimension(nk) :: p,t,td,pi,q,th,thv,z,pt,pb,pc,pn,ptv
%
% real :: the,maxthe,parea,narea,lfc
% real :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
% real :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
% real :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
% real*8 :: avgth,avgqv
% real :: getqvs,getqvi,getthe
p=zeros(size(p_in)); t=p; td=p; pi=p; q=p; th=p; thv=p; z=p; pt=p; pb=p; pc=p; pn=p; ptv=p;
lcl=0; nbl=0;
%-----------------------------------------------------------------------

g     = 9.81 ;
p00   = 100000.0;
cp    = 1005.7;
rd    = 287.04;
rv    = 461.5;
xlv   = 2501000.0;
xls   = 2836017.0;
t0    = 273.15;
cpv   = 1875.0;
cpl   = 4190.0;
cpi   = 2118.636;
lv1   = xlv+(cpl-cpv)*t0;
lv2   = cpl-cpv;
ls1   = xls+(cpi-cpv)*t0;
ls2   = cpi-cpv;

rp00  = 1.0/p00;
eps   = rd/rv;
reps  = rv/rd;
rddcp = rd/cp;
cpdrd = cp/rd;
cpdg  = cp/g;

converge = 0.0002;

debug_level =   400;
debug_level =   0;

%-----------------------------------------------------------------------

%---- convert p,t,td to mks units; get pi,q,th,thv ----%

for k=1:nk
    p(k) = 100.0*p_in(k);
    t(k) = 273.15+t_in(k);
    td(k) = 273.15+td_in(k);
    pi(k) = (p(k)*rp00)^rddcp;
    q(k) = getqvs(p(k),td(k));
    th(k) = t(k)/pi(k);
    thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k));
end

%---- get height using the hydrostatic equation ----%

z(1) = 0.0;
for k=2:nk
    dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1));
    z(k) = z(k-1) + dz;
end

%---- find source parcel ----%

if(source==1)
    % use surface parcel
    kmax = 1;
    
elseif(source==2)
    % use most unstable parcel (max theta-e)
    
    if(p(1)<50000.0)
        % first report is above 500 mb ... just use the first level reported
        kmax = 1;
        maxthe = getthe(p(1),t(1),td(1),q(1));
    else
        % find max thetae below 500 mb
        maxthe = 0.0;
        for k=1:nk
            if(p(k)>=50000.0)
                the = getthe(p(k),t(k),td(k),q(k));
                if( the>maxthe )
                    maxthe = the;
                    kmax = k;
                end
            end
        end
    end
    if(~exist('kmax')); kmax=1; end
    %if(debug_level>=100) %print *,'  kmax,maxthe = ',kmax,maxthe
    
elseif(source==3)
    % use mixed layer
    
    if( (z(2)-z(1))> ml_depth )
        % the second level is above the mixed-layer depth:  just use the
        % lowest level
        
        avgth = th(1);
        avgqv = q(1);
        kmax = 1;
        
    elseif( z(nk)< ml_depth )
        % the top-most level is within the mixed layer:  just use the
        % upper-most level
        
        avgth = th(nk);
        avgqv = q(nk);
        kmax = nk;
        
    else
        % calculate the mixed-layer properties:
        
        avgth = 0.0;
        avgqv = 0.0;
        k = 2;
        %if(debug_level>=100) %print *,'  ml_depth = ',ml_depth
        %if(debug_level>=100) %print *,'  k,z,th,q:'
        % if(debug_level>=100) %print *,1,z(1),th(1),q(1)
        
        while( (z(k)<=ml_depth) && (k<=nk) )
            
            %if(debug_level>=100) %print *,k,z(k),th(k),q(k)
            
            avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1));
            avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1));
            
            k = k + 1;
            
        end
        
        th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1));
        qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1));
        
        %if(debug_level>=100) %print *,999,ml_depth,th2,qv2
        
        avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1));
        avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1));
        
        %if(debug_level>=100) %print *,k,z(k),th(k),q(k);
        
        avgth = avgth/ml_depth;
        avgqv = avgqv/ml_depth;
        
        kmax = 1;
    end
    
    %if(debug_level>=100) %print *,avgth,avgqv
    
end

%---- define parcel properties at initial location ----%
narea = 0.0;

if( (source==1)||(source==2) )
    k    = kmax;
    th2  = th(kmax);
    pi2  = pi(kmax);
    p2   = p(kmax);
    t2   = t(kmax);
    thv2 = thv(kmax);
    qv2  = q(kmax);    
    b2   = 0.0;
elseif( source==3 )
    k    = kmax;
    th2  = avgth;
    qv2  = avgqv;
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2);
    pi2  = pi(kmax);
    p2   = p(kmax);
    t2   = th2*pi2;
    b2   = g*( thv2-thv(kmax) )/thv(kmax);
end

ql2 = 0.0;
qi2 = 0.0;
qt  = qv2;

cape = 0.0;
cin  = 0.0;
lfc  = 0.0;

doit = 1;
cloud = 0;
if(adiabat==1||adiabat==2)
    ice = 0;
else
    ice = 1;
end

the = getthe(p2,t2,t2,qv2);
if(debug_level>=100); disp(['  the = ',num2str(the)]); end

%---- begin ascent of parcel ----%

if(debug_level>=100)
    disp('  Start loop:');
    disp(['  p2,th2,qv2 = ',num2str([p2,th2,qv2])]);
end

lcl_tmp = [];
while( doit && (k<nk) )
    
    k = k+1;
    b1 =  b2;
    
    dp = p(k-1)-p(k);
    
    if( dp<pinc )
        nloop = 1;
    else
        nloop = 1 +  ( dp/pinc );
        dp = dp/ (nloop);
    end
    
    for n=1:nloop
        
        p1 =  p2;
        t1 =  t2;
        pi1 = pi2;
        th1 = th2;
        qv1 = qv2;
        ql1 = ql2;
        qi1 = qi2;
        thv1 = thv2;
        
        p2 = p2 - dp;
        pi2 = (p2*rp00)^rddcp;
        
        thlast = th1;
        i = 0;
        not_converged = 1;
        
        while( not_converged )
            i = i + 1;
            t2 = thlast*pi2;
            if(ice)
                fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0);
                fice = 1.0-fliq;
            else
                fliq = 1.0;
                fice = 0.0;
            end
            qv2 = min( qt , fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) );
            qi2 = max( fice*(qt-qv2) , 0.0 );
            ql2 = max( qt-qv2-qi2 , 0.0 );
            
            tbar  = 0.5*(t1+t2);
            qvbar = 0.5*(qv1+qv2);
            qlbar = 0.5*(ql1+ql2);
            qibar = 0.5*(qi1+qi2);
            
            lhv = lv1-lv2*tbar;
            lhs = ls1-ls2*tbar;
            lhf = lhs-lhv;
            
            rm=rd+rv*qvbar;
            cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar;
            th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar) +lhs*(qi2-qi1)/(cpm*tbar)     ...
                +(rm/cpm-rd/cp)*log(p2/p1) );
            % Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
            % Equation (38)
            
            if(i>90); disp([i,th2,thlast,th2-thlast]); end
            if(i>100)
                 error('Lack of convergence')
            end
            if( abs(th2-thlast)>converge )
                thlast=thlast+0.3*(th2-thlast);
            else
                not_converged = 0;
            end
        end
        
        lcl_tmp = [lcl_tmp; [p2 t2 getqvs(p2, t2) qt]] ;
         
        % Latest pressure increment is complete.  Calculate some important stuff:
        
        if( ql2>=1.0e-10 ); cloud = 1; end
        
        if(adiabat==1||adiabat==3)
            % pseudoadiabat
            qt  = qv2;
            ql2 = 0.0;
            qi2 = 0.0;
        elseif(adiabat<=0||adiabat>=5)            
            error('  Undefined adiabat');            
            return
        end
        
    end
    
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2);
    b2 = g*( thv2-thv(k) )/thv(k);
    dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1));
    
    the = getthe(p2,t2,t2,qv2);
    
    % Get contributions to CAPE and CIN:
    
    if( (b2>=0.0) && (b1<0.0) )
        % first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1);
        frac = b2/(b2-b1);
        parea =  0.5*b2*dz*frac;
        narea = narea-0.5*b1*dz*(1.0-frac);
        if(debug_level>=200)
            %print *,'      b1,b2 = ',b1,b2
            %print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
            %print *,'      frac = ',frac
            %print *,'      parea = ',parea
            %print *,'      narea = ',narea
        end
        cin  = cin  + narea;
        narea = 0.0;
        lfc=ps; % level of free convection
        
    elseif( (b2<0.0) && (b1>0.0) )
        % first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1);
        nbl =ps;
        frac = b1/(b1-b2);
        parea =  0.5*b1*dz*frac;
        narea = -0.5*b2*dz*(1.0-frac);
        if(debug_level>=200);
            %print *,'      b1,b2 = ',b1,b2
            %print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
            %print *,'      frac = ',frac
            %print *,'      parea = ',parea
            %print *,'      narea = ',narea
        end
    elseif( b2<0.0 )
        % still collecting negative buoyancy
        parea =  0.0;
        narea = narea-0.5*dz*(b1+b2);
    else
        % still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2);
        narea =  0.0;
    end
    
    cape = cape + max(0.0,parea);
    
    % if(debug_level>=200)
    %   write(6,102) p2,b1,b2,cape,cin,cloud;
    %    102     format(5(f13.4),2x,l1)
    %end
    
    if( (p(k)<=10000.0)&&(b2<0.0) )
        % stop if b < 0 and p < 100 mb
        doit = 0;
    end
    
end


[tmp,mid] = min(abs(lcl_tmp(:,4)-q(1)));
plcl = lcl_tmp(mid,1)/1e2;
tlcl = lcl_tmp(mid,2)-273.15;


end % end if function
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

function [out] = getqvs(p,t)

eps = 287.04/461.5;

es = 611.2*exp(17.67*(t-273.15)/(t-29.65));
out = eps*es/(p-es);

end

%-----------------------------------------------------------------------

function [out] = getqvi(p,t)


eps = 287.04/461.5;

es1 = 611.2*exp(21.8745584*(t-273.15)/(t-7.66));
es2 = 611.2*exp(17.67*(t-273.15)/(243.5+(t-273.15))); % used by emanuel
es = es2;
out = eps*es/(p-es);

end


%-----------------------------------------------------------------------

function [out] = getthe(p, t, td, q)


if( (td-t)>=-0.1 )
    tlcl = t;
else
    tlcl = 56.0 + ( (td-56.0)^(-1) + 0.00125*log(t/td) )^(-1);
end

out=t*( (100000.0/p)^(0.2854*(1.0-0.28*q)) )*exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) );

end
%-----------------------------------------------------------------------
