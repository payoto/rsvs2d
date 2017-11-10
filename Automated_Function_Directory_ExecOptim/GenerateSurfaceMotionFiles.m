function [varargout]=GenerateSurfaceMotionFiles(varargin)
global GenerateSurfaceMotionFiles_Handle
nOut=nargout(GenerateSurfaceMotionFiles_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSurfaceMotionFiles_Handle(varargin{:});
end
