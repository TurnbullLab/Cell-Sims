%ProbMatrix - computed as per Charlie's notes 1/23/2018
clear;              % clear parameters from other programs

% Parameter structure to pass to functions...
par.Td      = 20;   % gcp doubling time (in h)
par.tau     = 8*24; % clone differentiation time (in h)
par.lambda  = 32;    % clone "tuning" constant (dimensionless)

% Time is specified in hours, from t0 to t1 with time step, dt
t0 = 2*24;          % start time(0=E15=P-4; 2 days = E17) 
t1 = 20*24;         % stop time (E15 + 20 days = P16)
dt = 0.1;             % time step (in h)

ts=t0:dt:t1;        % vector to hold times
tMax=length(ts);    % number of time steps

kmax=2000;          % k=n+m; kmax = max(k)

%Previous terminology in Charlie's notes...
%tmax=40
%clockmax=2000;
%dt=tmax/clockmax

P=zeros(1+kmax,1+kmax);
P(1+1,1+0)=1;
%for clock=1:clockmax
%    t=clock*dt;
for l=1:tMax
    t=ts(l);                % time is now indexed by l
    for k=1:kmax
        for m=0:k
            n=k-m;
            if(n>0)
                P(1+n,1+m)=P(1+n,1+m)+dt*(n-1)*alpha(t,par)*P(1+n-1,1+m);
            end
            if(m>0)
                P(1+n,1+m)=P(1+n,1+m)+dt*(n+1)*beta(t,par)*P(1+n+1,1+m-1);
            end
            P(1+n,1+m)=P(1+n,1+m)/(1+dt*n*(alpha(t,par)+beta(t,par)));
        end
    end
end
%Mean (clone size = 2*M) and variance for long times, when N=0.
%For the general case, we need to count both proliferating and
%differentiated cells...
mean=0;
var=0;
for m=1:kmax
    mean=mean+2*m*P(1,m+1);
    var=var+((2*m)^2)*P(1,m+1);
end
mean
var=var-(mean)^2
sd=sqrt(var)
