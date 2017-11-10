function [varargout]=REFINE_desvargradadvanced(varargin)
global REFINE_desvargradadvanced_Handle
nOut=nargout(REFINE_desvargradadvanced_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_desvargradadvanced_Handle(varargin{:});
end
