function g = nearFieldChan_withAMplitude(d_t,azimuth,elevation,U,lambda)

  

%%
% The nearFieldChan function computes the exact near-field phase shifts for each RIS
% element without using far-field approximations.
% It calculates how far each RIS antenna element is from the user.
% It then computes the phase shift required for each RIS element to correctly focus the
% signal towards the user.
% The function accounts for near-field effects, meaning the distance to each RIS element is different 
% (unlike in far-field models where all elements are assumed equidistan).
%%
% To generate the real near field phase array without any approximation
% Input:
%   - d_t: the distance with respect to the reference element(center)
%   - U: is the relative location of each element from center of the array
% Output: 
%   - phase: is the phase vector when the full expression of near field
%            applied
U = -U;
%   If we don't flip the coordinate system (U = -U), we are measuring distances as if
%  the RIS is the one observing the user instead of the user observing the RIS.
M = size(U,2);
T = length(d_t);
d_m = zeros(M,T);
% d_m is an M × T matrix, which stores the distance between each RIS element 
% and the user.

parfor t = 1:T
    d_m(:,t) = d_t(t)*sqrt(1+2*sin(azimuth(t))*cos(elevation(t))*U(2,:)/d_t(t)...
        + 2*sin(elevation(t))*U(3,:)/d_t(t) ...
        + (U(2,:).^2+U(3,:).^2) /d_t(t)^2 ); 

    % dm_t = dt *√(1+2sin(φ)cos(θ)(Uy/dt)+2sin(θ)(Uz/dt)+ (Uy^2+Uz^2)/dt^2

%%
%d_t: Distance of the user from the center of the RIS.
% U(2,:) = U_y: Y-coordinates of the RIS elements.
% U(3,:) = U_z: Z-coordinates of the RIS elements.
% azimuth = ϕ: User's azimuth angle.
% elevation = θ: User's elevation angle
%%
end


beta_n = (lambda ./ (4*pi .* d_m)).^2;



 
phase = exp(-1i*2*pi*d_m/lambda);
% Forward-propagating wave (default convention).
g = sqrt(beta_n) .* phase;  % [M × T]
end

%%
% this function calculates the distance of all the antenna elements based
% upon the users location and then apply a phase shift to each of them/
% So that each RIS element's reflection is aligned to maximize power at the receiver.
