function Z = myLinakge1(X)
% design an n^2log(n) linkage algorithm that is space effient, so not many
% elements need to be moved

% when modeling priority queues
% merge delete and insertion opertations line by line

% Lu Cheng
% 6.6.2018

global HED_VAL
global END_VAL
global DEL_VAL

HED_VAL=-1;
END_VAL=-2;
DEL_VAL=-3;

n = size(X,1);

% calculate distance matrix, only upper triangle is used
% dist= zeros(n,n)+DEL_VAL;
dist= zeros(n,n);
for i=1:n
    for j=i+1:n
        dist(i,j) = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
    end
end

nodeFlag = true(1,n);

% sort each row, and generate the double direction link list
hedInd = zeros(n-1,1);
hedVal = zeros(n-1,1);
% prev = zeros(n,n)+DEL_VAL;
% next = zeros(n,n)+DEL_VAL;
prev = zeros(n,n);
next = zeros(n,n);
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
    [minval, minind] = mymin(hedVal);
    ii = minind;
    jj = hedInd(ii);
    
%     [hedInd hedVal]
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
    
    % update tree node array
    treeNodeArr(ii) = i+n;
    treeNodeArr(jj) = 0;
    
    nodeFlag([ii jj])=true;
    
    for k=1:n-1
        
        if ~nodeFlag(k)
            continue;
        end
        
        if k>=1 && k<ii
            [prev(k,:), next(k,:), tmpHedInd] = delPointers1(prev(k,:), next(k,:), hedInd(k), ii);
            [prev(k,:), next(k,:), tmpHedInd] = delPointers1(prev(k,:), next(k,:), tmpHedInd, jj);
            
            [prev(k,:), next(k,:), tmpHedInd, tmpHedVal] = insertPointers1(dist(k,:), prev(k,:), next(k,:), tmpHedInd, ii);
            hedInd(k) = tmpHedInd;
            hedVal(k) = tmpHedVal;
            
        elseif k==ii
            tmpNodeFlag = nodeFlag;
            tmpNodeFlag(jj)=false;
            
            [prevVec, nextVec, hedInd(ii), hedVal(ii)] = genPointers(dist(ii,ii+1:end),tmpNodeFlag(ii+1:end),ii);
            tmpinds = ii+find(tmpNodeFlag(ii+1:end));
            prev(ii,tmpinds) = prevVec;
            next(ii,tmpinds) = nextVec;
            
        elseif k>ii && k<jj
            [prev(k,:), next(k,:), tmpHedInd] = delPointers1(prev(k,:), next(k,:), hedInd(k), jj);
            if tmpHedInd~=hedInd(k)
                if tmpHedInd==DEL_VAL
                    hedInd(k) = DEL_VAL;
                    hedVal(k) = DEL_VAL;
                else
                    hedInd(k) = tmpHedInd;
                    hedVal(k) = dist(k,tmpHedInd);
                end
            end 
        else
            break;
        end
    end
    
    nodeFlag(ii)=true;
    nodeFlag(jj)=false;
    
    if jj<n
        hedInd(jj)=DEL_VAL;
        hedVal(jj)=DEL_VAL;
    end
    
end


function [minVal, minInd] = mymin(vec)
global DEL_VAL
inds = find(vec~=DEL_VAL);
[tmpVal, tmpInd]=min(vec(inds));
minVal = tmpVal;
minInd = inds(tmpInd);


function y = calPairDist(distMat,nodeFlag,i,j)
% global DEL_VAL
veci = extractRow(distMat,i);
vecj = extractRow(distMat,j);
y = zeros(size(veci));
y(nodeFlag) = max(veci(nodeFlag),vecj(nodeFlag));
% y(~nodeFlag) = DEL_VAL;


function y = extractRow(distMat, i)
y = distMat(i,:);
y(1:i-1) = distMat(1:i-1,i)';


function [prevVec, nextVec, hedInd] = delPointers1(prevVec, nextVec, hedInd, i)
% delete the ith element from the double linked array
global HED_VAL
global END_VAL
global DEL_VAL

prevind = prevVec(i);
nextind = nextVec(i);

% not necessary, since deleted values only occur on rows with nodeFlag(k)==false
% elements that has been deleted    
% if prevind==DEL_VAL && nextind==DEL_VAL
%     return;
% end

% elements is the single left element
if prevind==HED_VAL && nextind==END_VAL
    hedInd=DEL_VAL;
    prevVec(i)=DEL_VAL;
    nextVec(i)=DEL_VAL;
    return;        
end

% remove the element from the list
if prevind==HED_VAL
    hedInd = nextind;
    prevVec(nextind) = prevind;
elseif nextind==END_VAL
    nextVec(prevind) = nextind;
else
    prevVec(nextind) = prevind;
    nextVec(prevind) = nextind;
end


function [prevVec, nextVec, hedInd, hedVal] = insertPointers1(distVec, prevVec, nextVec, hedInd, i)
% insert the updated value of distVec(i) into the double linked list
global HED_VAL
global END_VAL
global DEL_VAL

targetVal = distVec(i);
curNodeInd = hedInd;

% in case of all elements deleted, insert one new 
if curNodeInd==DEL_VAL
    hedInd=i;
    hedVal=targetVal;
    prevVec(i) = HED_VAL;
    nextVec(i) = END_VAL;
    return;
end

 % insert in the head
if distVec(curNodeInd)>=targetVal
    hedInd=i;
    hedVal=targetVal;
    prevVec(i) = HED_VAL;
    nextVec(i) = curNodeInd;     
    prevVec(curNodeInd) = i;
    return;
end

% insert in the middle or end
prevNodeInd=prevVec(curNodeInd);
assert(prevNodeInd==HED_VAL);
while curNodeInd~=END_VAL && distVec(curNodeInd)<targetVal
    prevNodeInd = curNodeInd;
    curNodeInd = nextVec(curNodeInd);
end

if curNodeInd==END_VAL
    % add to the end
    nextVec(prevNodeInd) = i;
    nextVec(i) = END_VAL;
    prevVec(i) = prevNodeInd;
else
    % add to the middle
    nextVec(prevNodeInd) = i;
    prevVec(curNodeInd) = i;

    nextVec(i) = curNodeInd;
    prevVec(i) = prevNodeInd;
end

% retrieve the distance of the hed element
hedVal = distVec(hedInd);


function [prevVec, nextVec, hedInd, hedVal] = genPointers(arr, flagArr, offset)

global HED_VAL
global END_VAL
global DEL_VAL

% devoted to handle the following case
% e.g. in total 8 nodes, 4 merge with 8 and 5,6,7 are deleted
% 8 will be delete so the input flagArr is all false
if ~any(flagArr)
    hedInd=DEL_VAL;
    hedVal=DEL_VAL;
    prevVec=[];
    nextVec=[];
    return;
end

flagInds = find(flagArr);

[~, idx] = sort(arr(flagInds));
nidx = [-offset-1 flagInds(idx) -offset-1];

prevVec = zeros(size(idx));
prevVec(idx) = offset + nidx(1:end-2);
try
    prevVec(idx(1)) = HED_VAL;
catch
    keyboard
end

nextVec=zeros(size(idx));
nextVec(idx) = offset + nidx(3:end);
nextVec(idx(end)) = END_VAL;

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



