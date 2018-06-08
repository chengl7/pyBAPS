function Z = myLinakge(X)
% design an n^2log(n) linkage algorithm that is space effient, so not many
% elements need to be moved
% Lu Cheng
% 1.6.2018

global HED_VAL
global END_VAL
global DEL_VAL

HED_VAL=-1;
END_VAL=-2;
DEL_VAL=NaN;

n = size(X,1);

% calculate distance matrix, only upper triangle is used
dist= zeros(n,n)+DEL_VAL;
for i=1:n
    for j=i+1:n
        dist(i,j) = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
    end
end

nodeFlag = true(1,n);

% sort each row, and generate the double direction link list
hedInd = zeros(n-1,1);
hedVal = zeros(n-1,1);
prev = zeros(n,n)+DEL_VAL;
next = zeros(n,n)+DEL_VAL;
for i=1:n-1
    tmpinds = i+find(nodeFlag(i+1:end));
    [prev(i,tmpinds), next(i,tmpinds), hedInd(i), hedVal(i)] = genPointers(dist(i,i+1:end), nodeFlag(i+1:end), i);
end

% printHedInfo(hedInd,prev,next)

% complete linkage core algorithm
treeNodeArr=1:n;
Z = zeros(n-1,3);
for i=1:n-1
    
    % find the pair to be merged
    [minval, minind] = min(hedVal);
    ii = minind;
    jj = hedInd(ii);
    
%     fprintf('%dth step, merger index:node %d:%d and %d:%d.\n',i,ii,treeNodeArr(ii),jj,treeNodeArr(jj));
    
    Z(i,:) = [sort([treeNodeArr(ii) treeNodeArr(jj)]) minval];
    
    nodeFlag([ii jj])=false;
    newdist = calPairDist(dist,nodeFlag,ii,jj);
    
    % update distance matrix
    % update ii row
    tmpInds = find(nodeFlag(1:ii-1));
    dist(tmpInds,ii) = newdist(tmpInds)';
    tmpInds = ii+find(nodeFlag(ii+1:end));
    dist(ii,tmpInds) = newdist(tmpInds);
    
%     % delete jj row
%     dist(jj,jj+1:end) = DEL_VAL;
%     dist(1:jj-1,jj) = DEL_VAL;
    
    % update tree node array
    treeNodeArr(ii) = i+n;
    treeNodeArr(jj) = 0;
    
%     printHedInfo(hedInd,prev,next)
    
    % update the sorted priority arrays
    [prev, next, hedInd, hedVal] = delPointers(dist, prev, next, hedInd, hedVal, nodeFlag, ii);
    [prev, next, hedInd, hedVal] = delPointers(dist, prev, next, hedInd, hedVal, nodeFlag, jj);
    
    if jj<n
        hedInd(jj)=NaN;
        hedVal(jj)=NaN;
    end
%     printHedInfo(hedInd,prev,next)

%     [hedInd hedVal]'
    
    
    nodeFlag(ii) = true;
%     for t=1:ii
%         if ~nodeFlag(t)
%             continue;
%         end
%         [prev, next, hedInd, hedVal] = insertPointers1(dist,prev,next,hedInd,hedVal, nodeFlag,t);
%     end

    [prev, next, hedInd, hedVal] = inserPointers(dist, prev, next, hedInd, hedVal, nodeFlag, ii);
    
%     plotMatrix(dist,prev,next);
    
%     printHedInfo(hedInd,prev,next)
        
end

function printHedInfo(hedInd,prev,next)
n=size(prev,1);
for t=1:n-1
    fprintf('%d %d %d\n',hedInd(t),prev(t,hedInd(t)),next(t,hedInd(t)));
end
fprintf('\n\n');

function y = calPairDist(distMat,nodeFlag,i,j)
global DEL_VAL
veci = extractRow(distMat,i);
vecj = extractRow(distMat,j);
y = zeros(size(veci));
y(nodeFlag) = max(veci(nodeFlag),vecj(nodeFlag));
y(~nodeFlag) = DEL_VAL;


function y = extractRow(distMat, i)
y = distMat(i,:);
y(1:i-1) = distMat(1:i-1,i)';


function [prev, next, hedInd, hedVal] = delPointers(dist, prev, next, hedInd, hedVal, nodeFlag, i)
% update the previous and next node of for all elements in row i

global HED_VAL
global END_VAL
global DEL_VAL

for ii=1:i-1
    
%     if i==4 && ii==3
%         keyboard
%     end
    
    if ~nodeFlag(ii)
        continue;
    end
    
    prevind = prev(ii,i);
    nextind = next(ii,i);
    
%     % elements that has been deleted    
%     if isnan(prevind) && isnan(nextind)
%         continue;
%     end
    
    % elements is the single left element
    if prevind==HED_VAL && nextind==END_VAL
        hedInd(ii)=DEL_VAL;
        hedVal(ii)=DEL_VAL;
        prev(ii,i)=DEL_VAL;
        next(ii,i)=DEL_VAL;
        continue;        
    end
    
    % remove the element from the list
    if prevind==HED_VAL
        hedInd(ii) = nextind;
        hedVal(ii) = dist(ii,nextind);
        prev(ii, nextind) = prevind;
    elseif nextind==END_VAL
        next(ii, prevind) = nextind;
    else
        prev(ii, nextind) = prevind;
        next(ii, prevind) = nextind;
    end
    
%     prev(ii,i)=DEL_VAL;
%     next(ii,i)=DEL_VAL;

end

hedInd(i)=DEL_VAL;
hedVal(i)=DEL_VAL;

% % not necessary
% prev(i,i+1:end) = DEL_VAL;
% next(i,i+1:end) = DEL_VAL;


function [prev, next, hedInd, hedVal] = inserPointers(dist, prev, next, hedInd, hedVal, nodeFlag, i)
% update the previous and next node of for all elements in row i

global HED_VAL
global END_VAL

for ii=1:i-1
    
    % when a node is deleted
    if ~nodeFlag(ii)
        continue;
    end
    
    targetVal = dist(ii,i);
    curNodeInd = hedInd(ii);
    
    % in case of all elements deleted, insert one new 
    if isnan(curNodeInd) 
        hedInd(ii)=i;
        hedVal(ii)=targetVal;
        prev(ii,i) = HED_VAL;
        next(ii,i) = END_VAL;
        continue;
    end
    
    % insert in the head
    if dist(ii, curNodeInd)>=targetVal
        hedInd(ii)=i;
        hedVal(ii)=targetVal;
        prev(ii,i) = HED_VAL;
        next(ii,i) = curNodeInd;     
        prev(ii,curNodeInd) = i;
        continue;
    end
    
    prevNodeInd=prev(ii,curNodeInd);
    assert(prevNodeInd==HED_VAL);
    while curNodeInd~=END_VAL && dist(ii, curNodeInd)<targetVal
        prevNodeInd = curNodeInd;
        curNodeInd = next(ii,curNodeInd);
    end

    if curNodeInd==END_VAL
        % add to the end
        next(ii,prevNodeInd) = i;
        next(ii,i) = END_VAL; %-1;
        prev(ii,i) = prevNodeInd;
    else
        
        % add to the middle
        next(ii,prevNodeInd) = i;
        prev(ii,curNodeInd) = i;
        
        
        next(ii,i) = curNodeInd;
        prev(ii,i) = prevNodeInd;
    end
    
%     plotMatrix(dist, prev, next)
    
end

if ~any(nodeFlag(i+1:end))
    return
end

[prevVec, nextVec, hedInd(i), hedVal(i)] = genPointers(dist(i,i+1:end),nodeFlag(i+1:end),i);
tmpinds = i+find(nodeFlag(i+1:end));
prev(i,tmpinds) = prevVec;
next(i,tmpinds) = nextVec;


function [prevVec, nextVec, hedInd, hedVal] = genPointers(arr, flagArr, offset)

global HED_VAL
global END_VAL

flagInds = find(flagArr);

[~, idx] = sort(arr(flagInds));
nidx = [-offset-1 flagInds(idx) -offset-1];

prevVec = zeros(size(idx));
prevVec(idx) = offset + nidx(1:end-2);
prevVec(idx(1)) = HED_VAL; %-1;

nextVec=zeros(size(idx));
nextVec(idx) = offset + nidx(3:end);
nextVec(idx(end)) = END_VAL; %-1;

hedInd = flagInds(idx(1)) + offset;
hedVal = arr(flagInds(idx(1)));

function plotMatrix(dist, prev, next)
clf
subplot(2,2,1)
myplot(prev)
title('prev')
subplot(2,2,3)
myplot(next)
title('next')
subplot(2,2,2)
myplot(dist)
title('dist')

function myplot(mat)
mat(isnan(mat))=-5;
imagesc(mat);

function [prev, next, hedInd, hedVal] = insertPointers1(dist,prev,next,hedInd,hedVal, nodeFlag,i)
arr = dist(i,i+1:end);
flagArr = nodeFlag(i+1:end);

if ~any(flagArr)
    hedInd(i)=NaN;
    hedVal(i)=NaN;
    return
end

prevVec = zeros(size(arr))-1;
nextVec = zeros(size(arr))-1;
[prevVec(flagArr), nextVec(flagArr), hedInd(i), hedVal(i)] = genPointers(arr, flagArr, i);
prev(i,i+1:end) = prevVec;
next(i,i+1:end) = nextVec;



