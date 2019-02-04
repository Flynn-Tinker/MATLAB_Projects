function [U10,cd,ust,cd10,tau] = mf_dragNC35(Zvec,U,Ts)
% FUNCTION MF_DRAGNC35
% Compute the neutral transfer coefficients using the TC 3.5 bulk algorithm.
% [U10,cd,ust,cd10,tau] = mf_dragNC35(Zvec,U,Ts); where
% Zvec is the height of the wind measurement (m) (may be a single number or
%   a function of time)
% U is the wind measuremnt (m/s)
% Ts is the sea temperature to compute the kinematic viscosity, visa (oC)
%
% from Jim Edson 10/2012 as dragNC35.m
% modified by Melanie Fewings 12/10/12 to allow Z and Ts to vary with time, 
%   and calculate wind stress, 
% 4/17/13 and remove factor of 1000 that Edson had multiplying Cd


k=0;
ug=0.2;  %Neutral gustiness value
von=0.4; %von Karman constant
         
umax=19;   %Charnock is fixed after this value
a1=0.0017; %Coefficients used in parameterization of Charnock variable
a2=-0.005;

% if Zvec is a scalar, make it a vector (constant in time)
if sum(size(Zvec))==2
    Zvec = repmat(Zvec,size(U));
end
% if Ts is a scalar, make it a vector (constant in time)
if sum(size(Ts))==2
    Ts = repmat(Ts,size(U));
end

for k=1:length(U)
    uZ=U(k);               %wind speed
    t=Ts(k);
    Z = Zvec(k);
    visa = 1.326e-5*(1+6.542e-3.*t+8.301e-6.*t.^2-4.84e-9.*t.^3);
    ut=sqrt(uZ*uZ+ug*ug);  %wind speed plus gustiness
    us=0.035*ut;           %initial guess of u*
    charn=a1*ut+a2;        %wind speed dependent Charnock variable
    if ut>umax
        charn=a1*umax+a2;  %maximum value
    end
    for i=1:10
        zo=visa/us/9+charn*us*us/9.8;  %roughness length
        us=ut*von/log(Z/zo);           %friction velocity
        u10=us/von*log(10/zo);         %U at 10-m
        charn=a1*u10+a2;               %parameterization requires u10
        if u10>umax
            charn=a1*umax+a2;
        end
    end
    cd(k)=von*von/log(Z/zo)/log(Z/zo);  %Drag coefficient
    ust(k)=us;                               %Friction velocity 
    U10(k)=u10;                              %U at 10-m
    cd10(k)=von*von/log(10/zo)/log(10/zo);  %Drag coefficient at 10-m
end
rhoa = 1.2200; % from stresstc
tau = rhoa*(cd10.*U10.^2);
