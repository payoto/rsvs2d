function [varargout]=OutputFreefempp(varargin)
global OutputFreefempp_Handle
nOut=nargout(OutputFreefempp_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OutputFreefempp_Handle(varargin{:});
end
