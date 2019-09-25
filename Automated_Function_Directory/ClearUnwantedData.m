function [varargout]=ClearUnwantedData(varargin)
% include_Optimisation
global ClearUnwantedData_Handle
try
nOut=nargout(ClearUnwantedData_Handle);
catch
include_Optimisation
nOut=nargout(ClearUnwantedData_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ClearUnwantedData_Handle(varargin{:});
end
