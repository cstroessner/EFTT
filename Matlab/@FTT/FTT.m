classdef FTT

properties
    Ahat
    
    basis = "chebyshev"
    numAllEvals
    numTTswp
end

methods
    function ftt = FTT(varargin)
        ftt = constructor(ftt, varargin{:});
    end
end

methods (Access = public)
    varargout = degree(ftt);
    varargout = ttrank(ftt);
    varargout = dof(ftt);
    
    varargout = fevaluate(ftt,x);
    varargout = fintegrate(ftt);
end
end

