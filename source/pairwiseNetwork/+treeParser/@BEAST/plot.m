function [ht, phyData]=plot(treeParserObj,varargin)
% pass through to parent class plotter

inputs=varargin;
[ht, phyData]=treeParser.plot(treeParserObj,inputs{:});

end