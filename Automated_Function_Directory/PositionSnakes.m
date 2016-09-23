function [varargout]=PositionSnakes(varargin)
global PositionSnakes_Handle
nOut=nargout(PositionSnakes_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PositionSnakes_Handle(varargin{:});
end
