function [varargout]=CalculateFloodRatio(varargin)
global CalculateFloodRatio_Handle
nOut=nargout(CalculateFloodRatio_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateFloodRatio_Handle(varargin{:});
end
