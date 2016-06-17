function [varargout]=DataArrayToStruct(varargin)
global DataArrayToStruct_Handle
nOut=nargout(DataArrayToStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DataArrayToStruct_Handle(varargin{:});
end
