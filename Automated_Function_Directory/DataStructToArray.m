function [varargout]=DataStructToArray(varargin)
global DataStructToArray_Handle
nOut=nargout(DataStructToArray_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DataStructToArray_Handle(varargin{:});
end
