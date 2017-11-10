function [varargout]=OptimisationOutput_Final_Post(varargin)
% OptimisationOutput
global OptimisationOutput_Final_Post_Handle
nOut=nargout(OptimisationOutput_Final_Post_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_Final_Post_Handle(varargin{:});
end
