function [varargout]=GetStructureData(varargin)
% include_Utilities
global GetStructureData_Handle
try
nOut=nargout(GetStructureData_Handle);
catch
include_Utilities
nOut=nargout(GetStructureData_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GetStructureData_Handle(varargin{:});
end
