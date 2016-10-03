function [varargout]=ReduceCornerFrac(varargin)
global ReduceCornerFrac_Handle
nOut=nargout(ReduceCornerFrac_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReduceCornerFrac_Handle(varargin{:});
end
