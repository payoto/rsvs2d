function [varargout]=ExtractOptimalSnake(varargin)
% OptimisationOutput
global ExtractOptimalSnake_Handle
nOut=nargout(ExtractOptimalSnake_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractOptimalSnake_Handle(varargin{:});
end
