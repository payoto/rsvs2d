function [varargout]=FindIdenticalVector(varargin)
global FindIdenticalVector_Handle
nOut=nargout(FindIdenticalVector_Handle);
[varargout{1:nOut}]=FindIdenticalVector_Handle(varargin{:});
end
