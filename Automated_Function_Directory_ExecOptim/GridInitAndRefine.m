function [varargout]=GridInitAndRefine(varargin)
global GridInitAndRefine_Handle
nOut=nargout(GridInitAndRefine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridInitAndRefine_Handle(varargin{:});
end
