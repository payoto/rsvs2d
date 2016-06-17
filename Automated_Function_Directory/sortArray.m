function [varargout]=sortArray(varargin)
global sortArray_Handle
nOut=nargout(sortArray_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=sortArray_Handle(varargin{:});
end
