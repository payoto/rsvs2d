function [varargout]=GetStructureData(varargin)
% include_Utilities
global GetStructureData_Handle
nOut=nargout(GetStructureData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GetStructureData_Handle(varargin{:});
end
