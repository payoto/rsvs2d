function [varargout]=IdentifyCorners(varargin)
% include_Optimisation
global IdentifyCorners_Handle
try
nOut=nargout(IdentifyCorners_Handle);
catch
include_Optimisation
nOut=nargout(IdentifyCorners_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifyCorners_Handle(varargin{:});
end
