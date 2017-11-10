function [varargout]=OptimisationOutput_Final(varargin)
% OptimisationOutput
global OptimisationOutput_Final_Handle
nOut=nargout(OptimisationOutput_Final_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_Final_Handle(varargin{:});
end
