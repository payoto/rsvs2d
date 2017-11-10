function [varargout]=StartFromFill(varargin)
global StartFromFill_Handle
nOut=nargout(StartFromFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=StartFromFill_Handle(varargin{:});
end
