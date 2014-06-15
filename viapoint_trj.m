%%
% Copyright (C) 2013 Gennaro Raiola, ENSTA-ParisTech
%
% viapoint_trj is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% viapoint_trj is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with viapoint_trj. If not, see <http://www.gnu.org/licenses/>.
%%
function [ pos vel acc N_viapoint ] = viapoint_trj(dt,time,x0,xd0,xdd0,xvia,xg,viapoint_time_ratio,figure_handle)
% This function generates a trajectory passing for a defined viapoint. The trajectory is generated using a spline interpolating minjerk equations.
% The script provides a test function; run the scipt without args to call
% the test.
%%
% Input:
%  - dt: sample time
%  - time: trajectory length in secs
%  - x0: initial condition
%  - xd0: first derivate initial condition
%  - xdd0: second derivate initial condition
%  - xvia: viapoint
%  - viapoint_time_ratio: define the time of the viapoint
%  - figure_handle: do you want to plot something?
%%
% Ouput:
%  - pos: output points
%  - vel: ouput velocities
%  - acc: output accelerations
%  - N_viapoint: vector index for the viapoint
%
% Author: Gennaro Raiola
%%
if (nargin==0),test_viapoint_trj; return; end
if (nargin<8), viapoint_time_ratio = 0.5; end
if (nargin<9), figure_handle = 0; end  

n_dim = length(x0);
N = ceil(time/dt)+1;

time_viapoint = viapoint_time_ratio*time;
N_viapoint = ceil(time_viapoint/dt)+1;
x_time = [0 time_viapoint time];

for i_dim = 1:n_dim
  x_values = [x0(i_dim) xvia(i_dim) xg(i_dim)];
  velocity_conditions = [xd0(i_dim) 0]; % 0 = final value
  acceleration_condtions = [xdd0(i_dim) 0];
  [ pos(:,i_dim) vel(:,i_dim) acc(:,i_dim) ] = minjerkspline(dt,x_time,x_values,velocity_conditions,acceleration_condtions);
end

if (figure_handle)
  figure(figure_handle)
  subplot(3,2,1)
  plot(pos,'-k');
  hold on;
  plot(1,x0,'or')
  plot(N_viapoint,xvia,'ob')
  plot(N,xg,'og')
  hold off
  xlabel('time')
  ylabel('x')

  subplot(3,2,3)
  plot(vel,'-k');
  xlabel('time')
  ylabel('v')

  subplot(3,2,5)
  plot(acc,'-k');
  xlabel('time')
  ylabel('a')

  subplot(3,2,2:2:6)
  plot(pos(:,1),pos(:,2),'-k');
  hold on;
  plot(x0(1),x0(2),'or')
  plot(xvia(1),xvia(2),'ob')
  plot(xg(1),xg(2),'og')
  hold off
  xlabel('x_1')
  ylabel('x_2')
  axis equal
end

  function test_viapoint_trj
    dt = 0.01;
    time = 6;
    x0 = [1 4];
    xd0 = [0 0];
    xdd0 = [0 0];
    xvia = [4 1];
    xg = [5 3];
    figure_handle = 1;
    viapoint_time_ratio = 0.5;
    [pos vel acc N_viapoint] = viapoint_trj(dt,time,x0,xd0,xdd0,xvia,xg,viapoint_time_ratio,figure_handle);
  end

end