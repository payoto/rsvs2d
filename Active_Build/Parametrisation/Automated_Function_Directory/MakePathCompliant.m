function [varargout]=MakePathCompliant(varargin)
global MakePathCompliant_Handle
nOut=nargout(MakePathCompliant_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakePathCompliant_Handle(varargin{:});
end
