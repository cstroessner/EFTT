function [y] = borehole(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOREHOLE FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = water flow rate
% xx = [rw, r, Tu, Hu, Tl, Hl, L, Kw]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rw ∈ [0.05, 0.15]	radius of borehole (m)
% r ∈ [100, 50000]	radius of influence (m)
% Tu ∈ [63070, 115600]   	transmissivity of upper aquifer (m2/yr)
% Hu ∈ [990, 1110]	potentiometric head of upper aquifer (m)
% Tl ∈ [63.1, 116]	transmissivity of lower aquifer (m2/yr)
% Hl ∈ [700, 820]	potentiometric head of lower aquifer (m)
% L ∈ [1120, 1680]	length of borehole (m)
% Kw ∈ [9855, 12045]	hydraulic conductivity of borehole (m/yr)
left = [0.05;100;63070;990;63.1;700;1120;9855]';
right = [0.15;50000;115600;1110;116;820;1680;12045]';

xx = xx.*(right -left)/2 + (right +left)/2;

rw = xx(1);
r  = xx(2);
Tu = xx(3);
Hu = xx(4);
Tl = xx(5);
Hl = xx(6);
L  = xx(7);
Kw = xx(8);

frac1 = 2 * pi * Tu * (Hu-Hl);

frac2a = 2*L*Tu / (log(r/rw)*rw^2*Kw);
frac2b = Tu / Tl;
frac2 = log(r/rw) * (1+frac2a+frac2b);

y = frac1 / frac2;

end
