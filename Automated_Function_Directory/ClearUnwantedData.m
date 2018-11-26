function [varargout]=ClearUnwantedData(varargin)
% include_Optimisation
global ClearUnwantedData_Handle
nOut=nargout(ClearUnwantedData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ClearUnwantedData_Handle(varargin{:});
end
