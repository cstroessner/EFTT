function [Vm] = otlcircuit(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OTL CIRCUIT FUNCTION
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
% Vm = midpoint voltage
% xx = [Rb1, Rb2, Rf, Rc1, Rc2, beta]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rb1 ∈ [50, 150]	resistance b1 (K-Ohms)
% Rb2 ∈ [25, 70]	resistance b2 (K-Ohms)
% Rf ∈ [0.5, 3]	resistance f (K-Ohms)
% Rc1 ∈ [1.2, 2.5]	resistance c1 (K-Ohms)
% Rc2 ∈ [0.25, 1.2]   	resistance c2 (K-Ohms)
% β ∈ [50, 300]	current gain (Amperes)
left = [50;25;0.5;1.2;0.25;50]';
right = [150;70;3;2.5;1.2;300]';
xx = xx.*(right -left)/2 + (right +left)/2;

Rb1  = xx(1);
Rb2  = xx(2);
Rf   = xx(3);
Rc1  = xx(4);
Rc2  = xx(5);
beta = xx(6);

Vb1 = 12*Rb2 / (Rb1+Rb2);
term1a = (Vb1+0.74) * beta * (Rc2+9);
term1b = beta*(Rc2+9) + Rf;
term1 = term1a / term1b;

term2a = 11.35 * Rf;
term2b = beta*(Rc2+9) + Rf;
term2 = term2a / term2b;

term3a = 0.74 * Rf * beta * (Rc2+9);
term3b = (beta*(Rc2+9)+Rf) * Rc1;
term3 = term3a / term3b;

Vm = term1 + term2 + term3;

end
