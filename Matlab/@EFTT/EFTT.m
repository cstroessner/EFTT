classdef EFTT

properties
    chat
    TF
    
    basis
    numAllEvals
    numTTEvals
    numTTswp
end

methods
    function eftt = EFTT(varargin)
        eftt = constructor(eftt, varargin{:});
    end
end

methods (Access = public)
    varargout = degree(eftt);
    varargout = ttrank(eftt);
    varargout = tuckerrank(eftt);
    varargout = dof(eftt);
    
    varargout = fevaluate(eftt,x);
    varargout = fintegrate(eftt);
end
end

