function [varargout]=PrepareCFDPostProcessing(varargin)
% OptimisationOutput
global PrepareCFDPostProcessing_Handle
nOut=nargout(PrepareCFDPostProcessing_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PrepareCFDPostProcessing_Handle(varargin{:});
end
