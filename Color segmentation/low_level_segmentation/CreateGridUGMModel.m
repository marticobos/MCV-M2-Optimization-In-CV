function [edgePot,edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K, lambda, data_term)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor

tic
nRows = NumFils;
nCols = NumCols;
nStates = K;

nNodes = nRows*nCols;

 
adj = sparse(nNodes,nNodes);
 
% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;
 
% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;
 
% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
   n1 = edgeStruct.edgeEnds(e,1);
   n2 = edgeStruct.edgeEnds(e,2);

   %pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
   pot_same = exp(1.8 + .3*1/(1+abs(data_term(n1)-data_term(n2))));
   pot = ones(K,K);
   for i = 1:K
    pot(i,i) = pot_same;
   end
   edgePot(:,:,e) = pot;
   %edgePot(:,:,e) = [pot_same 1;1 pot_same];
end


toc;