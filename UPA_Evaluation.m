function [h, U] = UPA_Evaluation(lambda,M_V,M_H,Azimuth,Elevation,vd,hd,Rate)

%   Inputs:
%       lambda: Wavelength
%       M_V:    Number of vertical antenna elements
%       M_H:    Number of horizental antenna elements
%       Azimuth
%       Elevation
%
%   Outputs:
%       UPA_response: Antenna response with respect to Azimuth and
%       Elevation WHICH IS THE ARRAY RESPONSE
%       U:  Antenna Elements positions
%           Arrays are on Y-Z Axis. The index of the array is as follows:
%           6 7 8 
%           3 4 5
%           0 1 2
%

%% Parameters Initializations

d_H = hd*lambda; % Horizontal antenna spacing   

d_V = vd*lambda; % Vertical antenna spacing

M = M_H*M_V; % Total number of antennas
U = zeros(3,M); % Matrix containing the position of the antennas.   
%  3 for 3d position of RIS x,y,z

i = @(m) mod(m-1,M_H);% Horizontal index for RIS index

j = @(m) floor((m-1)/M_H);% Vertical index for RIS index



%% Element position evaluation

% This one is the normal convention to evaluate position of the elements
for m = 1:M
    
        U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the m-th element
%U is a 3×M matrix where each column represents the 3D coordinates of an antenna element.

end


% The evaluation is based on the following book:
% Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency
% by Emil Björnson, Jakob Hoydis, Luca Sanguinetti


N = numel(Azimuth); 
        if isscalar(Rate), Rate = repmat(Rate,1,N); end
h = zeros(M,N); a = zeros(M,N);

for n = 1:N
  az = Azimuth(n); el = Elevation(n);
  k  = (2*pi/lambda)*[cos(az)*cos(el); sin(az)*cos(el); sin(el)];
  a_n =  transpose(exp(1i* k'*U)); 
  beta_n = (lambda/(4*pi*Rate(n)))^2;      % Friis with Gtx=Grx=1
  h(:,n) = sqrt(beta_n) * a_n;
  a(:,n) = a_n;
end
end