function [varargout]=RebuildFromIter(varargin)
% OptimisationOutput
global RebuildFromIter_Handle
nOut=nargout(RebuildFromIter_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RebuildFromIter_Handle(varargin{:});
end
