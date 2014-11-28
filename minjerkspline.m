%%
% Copyright (C) 2013 Gennaro Raiola, ENSTA-ParisTech
%
% minjerkspline is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% minjerkspline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with minjerkspline. If not, see <http://www.gnu.org/licenses/>.
%%
function [y yd ydd yddd] = minjerkspline(dt,x_time,x,velocity_conditions,acceleration_condtions,figure_handle)
% This algorithm generates a spline using minjerk equations.
%%
% Input:
%  - dt: sample time
%  - x_time: time instants
%  - x: points to interpolate (At least 3 points)
%  - velocity_conditions: bound conditions in velocity
%  - acceleration_conditions: bound conditions in acceleration
%  - figure_handle: do you want to plot something?
%%
% Ouput:
%  - y: output points
%  - yd: ouput velocities
%  - ydd: output accelerations
%  - yddd: output jerk (third derivative)
%
% Author: Gennaro Raiola
%%

if (nargin < 6), figure_handle = 0; end

n = length(x); % Number of points 
t = x_time;
T = diff(t);

% Initial and final conditions
v1 = velocity_conditions(1);
vn = velocity_conditions(end);
a1 = acceleration_condtions(1);
an = acceleration_condtions(end);

% Useful initiliazions for the loops
c_v1_eq4 = zeros(1,n-2);
c_a1_eq4 = zeros(1,n-2);
c_v1_eq5 = zeros(1,n-2);
c_a1_eq5 = zeros(1,n-2);
c_v2_eq4 = zeros(1,n-2);
c_a2_eq4 = zeros(1,n-2);
c_v2_eq5 = zeros(1,n-2);
c_a2_eq5 = zeros(1,n-2);
c_v3_eq4 = zeros(1,n-2);
c_a3_eq4 = zeros(1,n-2);
c_v3_eq5 = zeros(1,n-2);
c_a3_eq5 = zeros(1,n-2);
ck_eq4 = zeros(1,n-2);
ck_eq5 = zeros(1,n-2);
alpha0 = zeros(1,n-1);
alpha1 = zeros(1,n-1);
alpha2 = zeros(1,n-1);
alpha3 = zeros(1,n-1);
alpha4 = zeros(1,n-1);
alpha5 = zeros(1,n-1);

for k = 1:n-2
  c_v1_eq4(k) =  -8/T(k)^2;
  c_a1_eq4(k) = -1/T(k);
  c_v1_eq5(k) = -42/T(k)^3;
  c_a1_eq5(k) = -6/T(k)^2; 
  c_v2_eq4(k) = (-12/T(k)^2 + 12/T(k+1)^2);
  c_a2_eq4(k) = 3/T(k+1) + 3/T(k);
  c_v2_eq5(k) = (-48/T(k)^3 - 48/T(k+1)^3);
  c_a2_eq5(k) = (-9/T(k+1)^2 + 9/T(k)^2);
  c_v3_eq4(k) = 8/T(k+1)^2;
  c_a3_eq4(k) = -1/T(k+1);
  c_v3_eq5(k) = -42/T(k+1)^3;
  c_a3_eq5(k) = 6/T(k+1)^2;
  ck_eq4(k) = 20/T(k)^3 * x(k) + -(20/T(k+1)^3 + 20/T(k)^3) * x(k+1) + 20/T(k+1)^3 * x(k+2);
  ck_eq5(k) = 90/T(k)^4 * x(k) - (-90/T(k+1)^4 + 90/T(k)^4) * x(k+1) - 90/T(k+1)^4 * x(k+2);
end

c_v1_eq4 = [c_v1_eq4 0];
c_a1_eq4 = [c_a1_eq4 0];
c_v1_eq5 = [c_v1_eq5 0];
c_a1_eq5 = [c_a1_eq5 0];

c_v2_eq4 = [0 c_v2_eq4 0];
c_a2_eq4 = [0 c_a2_eq4 0];
c_v2_eq5 = [0 c_v2_eq5 0];
c_a2_eq5 = [0 c_a2_eq5 0];

c_v3_eq4 = [0 c_v3_eq4];
c_a3_eq4 = [0 c_a3_eq4];
c_v3_eq5 = [0 c_v3_eq5];
c_a3_eq5 = [0 c_a3_eq5];

ck_eq4(1) = ck_eq4(1) - c_v1_eq4(1) * v1 - c_a1_eq4(1) * a1;
ck_eq5(1) = ck_eq5(1) - c_v1_eq5(1) * v1 - c_a1_eq5(1) * a1;
ck_eq4(end) = ck_eq4(end) - c_v3_eq4(end) * vn - c_a3_eq4(end) * an;
ck_eq5(end) = ck_eq5(end) - c_v3_eq5(end) * vn - c_a3_eq5(end) * an;

ck_eq4 = ck_eq4(:);
ck_eq5 = ck_eq5(:);

c = [ck_eq4;ck_eq5];

% Construction of the linear system
Av_eq4 = diag(c_v1_eq4,-1) + diag(c_v2_eq4) + diag(c_v3_eq4,1);
Av_eq4 = Av_eq4(2:end-1,:);
Aa_eq4 = diag(c_a1_eq4,-1) + diag(c_a2_eq4) + diag(c_a3_eq4,1);
Aa_eq4 = Aa_eq4(2:end-1,:);

Av_eq5 = diag(c_v1_eq5,-1) + diag(c_v2_eq5) + diag(c_v3_eq5,1);
Av_eq5 = Av_eq5(2:end-1,:);
Aa_eq5 = diag(c_a1_eq5,-1) + diag(c_a2_eq5) + diag(c_a3_eq5,1);
Aa_eq5 = Aa_eq5(2:end-1,:);

A = [Av_eq4 Aa_eq4 ; Av_eq5 Aa_eq5];

indices = [1 n n+1 2*n];

A(:,indices) = [];

va = A\c;
 
v = va(1:n-2);
a = va(n-1:end);

v = [v1 ; v ; vn];
a = [a1 ; a ; an];
 
for k = 1:n-1
  alpha0(k) = x(k);
  alpha1(k) = v(k);
  alpha2(k) = a(k)/2;
  alpha3(k) = 1/2*(-20*x(k)-12*v(k)*T(k)-3*a(k)*T(k)^2+20*x(k+1)-8*v(k+1)*T(k)+a(k+1)*T(k)^2)/T(k)^3;
  alpha4(k) = -1/2*(-16*v(k)*T(k)-3*a(k)*T(k)^2-30*x(k)+30*x(k+1)-14*v(k+1)*T(k)+2*a(k+1)*T(k)^2)/T(k)^4;
  alpha5(k) = 1/2*(-a(k)*T(k)^2-12*x(k)-6*v(k)*T(k)+12*x(k+1)-6*v(k+1)*T(k)+a(k+1)*T(k)^2)/T(k)^5;
end

y = [];
yd = [];
ydd = [];
yddd = [];

for k = 1:n-1
  tau = 0:dt:T(k);
  N(k) = length(tau)-1; %#ok<AGROW>
  for i = 1:N(k)
    cur_y(i) = alpha0(k) + alpha1(k) * tau(i) + alpha2(k) * tau(i)^2 + alpha3(k) * tau(i)^3 + ...
      alpha4(k) * tau(i)^4 + alpha5(k) * tau(i)^5; %#ok<AGROW>
    cur_yd(i) = alpha1(k) + 2 * alpha2(k) * tau(i) + 3 * alpha3(k) * tau(i)^2 + ...
      4 * alpha4(k) * tau(i)^3 + 5 * alpha5(k) * tau(i)^4; %#ok<AGROW>
    cur_ydd(i) = 2 * alpha2(k) + 6 * alpha3(k) * tau(i) + ...
      12 * alpha4(k) * tau(i)^2 + 20 * alpha5(k) * tau(i)^3; %#ok<AGROW>
    cur_yddd(i) = 6 * alpha3(k) + ...
      24 * alpha4(k) * tau(i) + 60 * alpha5(k) * tau(i)^2; %#ok<AGROW>
  end
  y = [y cur_y]; %#ok<AGROW>
  yd = [yd cur_yd]; %#ok<AGROW>
  ydd = [ydd cur_ydd]; %#ok<AGROW>
  yddd = [yddd cur_yddd]; %#ok<AGROW>
  cur_y = [];
  cur_yd = [];
  cur_ydd = [];
  cur_yddd = []; 
end
y = [y x(end)];
yd = [yd v(end)];
ydd = [ydd a(end)];

if (figure_handle)
  tt = t(1):dt:t(end);
 
  figure(figure_handle)
  subplot(2,2,1)
  plot(tt,y,'b')
  title('position')
  hold on
  plot(t,x,'ro')

  subplot(2,2,2)
  plot(tt,yd,'b')
  title('velocity')

  subplot(2,2,3)
  plot(tt,ydd,'b')
  title('acceleration')

  subplot(2,2,4)
  plot(tt(1:end-1),yddd,'b')
  title('jerk')
end

end
