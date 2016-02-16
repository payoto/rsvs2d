function [varargout]=PointGeneration(varargin)
global PointGeneration_Handle
nOut=nargout(PointGeneration_Handle);
[varargout{1:nOut}]=PointGeneration_Handle(varargin{:});
end
