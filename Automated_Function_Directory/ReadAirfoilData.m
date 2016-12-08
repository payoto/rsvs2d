function [varargout]=ReadAirfoilData(varargin)
global ReadAirfoilData_Handle
nOut=nargout(ReadAirfoilData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReadAirfoilData_Handle(varargin{:});
end
