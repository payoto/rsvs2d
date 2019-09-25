function [varargout]=ReadAirfoilData(varargin)
% include_Optimisation
global ReadAirfoilData_Handle
try
nOut=nargout(ReadAirfoilData_Handle);
catch
include_Optimisation
nOut=nargout(ReadAirfoilData_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadAirfoilData_Handle(varargin{:});
end
