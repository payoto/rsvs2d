function [varargout]=OptimisationOutput_Init(varargin)
% OptimisationOutput
global OptimisationOutput_Init_Handle
nOut=nargout(OptimisationOutput_Init_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_Init_Handle(varargin{:});
end
