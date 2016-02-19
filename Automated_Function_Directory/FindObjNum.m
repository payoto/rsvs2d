function [varargout]=FindObjNum(varargin)
global FindObjNum_Handle
nOut=nargout(FindObjNum_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindObjNum_Handle(varargin{:});
end
