function [varargout]=ReadAirfoilData(varargin)
% include_Optimisation
global ReadAirfoilData_Handle
nOut=nargout(ReadAirfoilData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadAirfoilData_Handle(varargin{:});
end
