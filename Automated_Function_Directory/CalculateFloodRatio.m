function [varargout]=CalculateFloodRatio(varargin)
% include_Optimisation
global CalculateFloodRatio_Handle
try
nOut=nargout(CalculateFloodRatio_Handle);
catch
include_Optimisation
nOut=nargout(CalculateFloodRatio_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateFloodRatio_Handle(varargin{:});
end
