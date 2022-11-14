function [y] = wingweight(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WING WEIGHT FUNCTION
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
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = wing weight
% xx = [Sw, Wfw, A, LamCaps, q, lam, tc, Nz, Wdg, Wp]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sw ∈ [150, 200]	wing area (ft2)
% Wfw ∈ [220, 300]	weight of fuel in the wing (lb)
% A ∈ [6, 10]	aspect ratio
% Λ ∈ [-10, 10]	quarter-chord sweep (degrees)
% q ∈ [16, 45]	dynamic pressure at cruise (lb/ft2)
% λ ∈ [0.5, 1]	taper ratio
% tc ∈ [0.08, 0.18]	aerofoil thickness to chord ratio
% Nz ∈ [2.5, 6]	ultimate load factor
% Wdg ∈ [1700, 2500]   	flight design gross weight (lb)
% Wp ∈ [0.025, 0.08]	paint weight (lb/ft2)

left = [150;220;6;-10;16;0.5;0.08;2.5;1700;0.025]';
right = [200;300;10;10;45;1;0.18;6;2500;0.08]';

xx = xx.*(right -left)/2 + (right +left)/2;


Sw      = xx(1);
Wfw     = xx(2);
A       = xx(3);
LamCaps = xx(4) * (pi/180);
q       = xx(5);
lam     = xx(6);
tc      = xx(7);
Nz      = xx(8);
Wdg     = xx(9);
Wp      = xx(10);

fact1 = 0.036 * Sw^0.758 * Wfw^0.0035;
fact2 = (A / ((cos(LamCaps))^2))^0.6;
fact3 = q^0.006 * lam^0.04;
fact4 = (100*tc / cos(LamCaps))^(-0.3);
fact5 = (Nz*Wdg)^0.49;

term1 = Sw * Wp;

y = fact1*fact2*fact3*fact4*fact5 + term1;

end
