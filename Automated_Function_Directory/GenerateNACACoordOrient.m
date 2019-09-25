function [varargout]=GenerateNACACoordOrient(varargin)
% include_Optimisation
global GenerateNACACoordOrient_Handle
try
nOut=nargout(GenerateNACACoordOrient_Handle);
catch
include_Optimisation
nOut=nargout(GenerateNACACoordOrient_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNACACoordOrient_Handle(varargin{:});
end
