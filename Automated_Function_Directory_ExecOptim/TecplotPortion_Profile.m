function [varargout]=TecplotPortion_Profile(varargin)
% OptimisationOutput
global TecplotPortion_Profile_Handle
nOut=nargout(TecplotPortion_Profile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TecplotPortion_Profile_Handle(varargin{:});
end
