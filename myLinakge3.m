function Z = myLinakge3(X)
% design an n^2log(n) linkage algorithm that is space effient, so not many
% elements need to be moved

% when modeling priority queues
% merge delete and insertion opertations line by line

% now update by chunk and record all modifications

% Lu Cheng
% 8.6.2018

global HED_VAL
global END_VAL
global DEL_VAL

global BLOCK_SIZE
global N_BLOCK
global N_NODE

HED_VAL=-1;
END_VAL=-2;
DEL_VAL=-3;

n = size(X,1);
N_NODE = n;

N_BLOCK = ceil(sqrt(n));
BLOCK_SIZE = ceil(n/N_BLOCK);
% nn = N_BLOCK*BLOCK_SIZE;

nodeFlag = true(1,n);

% sort each row, and generate the double direction link list
hedInd = zeros(n-1,1);
hedVal = zeros(n-1,1);

bprev = cell(N_BLOCK,N_BLOCK);
bnext = cell(N_BLOCK,N_BLOCK);
bdist = cell(N_BLOCK,N_BLOCK);
% bedit = cell(N_BLOCK,N_BLOCK);
for bi=1:N_BLOCK
    for bj=bi:N_BLOCK
        bdist{bi,bj} = calDist(X, bi, bj);
    end
end

for bi=1:N_BLOCK
    [resMat, offset] = getMatBlock(bdist,bi);
    
    tmpPrev = zeros(size(resMat));
    tmpNext = zeros(size(resMat));
    for j=1:BLOCK_SIZE
        tmpinds = (j+1):(n-offset);
        tmprowind = BLOCK_SIZE*(bi-1)+j;
        
        if tmprowind>n-1
            continue
        end
        
        tmprowflag = nodeFlag(tmprowind+1:n);
        assert(length(tmpinds)==length(tmprowflag));
        [tmpPrev(j,tmpinds(tmprowflag)), tmpNext(j,tmpinds(tmprowflag)), hedInd(tmprowind), hedVal(tmprowind)] = genPointers(resMat(j,tmpinds), tmprowflag, offset+j);
    end
    
    for j=bi:N_BLOCK
        bprev(bi,bi:end) = distMatBlock(tmpPrev);
        bnext(bi,bi:end) = distMatBlock(tmpNext);
    end
    
end

% complete linkage core algorithm
treeNodeArr = 1:n;
Z = zeros(n-1,3);
for i=1:n-1
    
    % find the pair to be merged
    [minval, minind] = mymin(hedVal);
    ii = minind;
    jj = hedInd(ii);
    
%     [hedInd hedVal]
    fprintf('%dth step, merger index:node %d:%d and %d:%d.\n',i,ii,treeNodeArr(ii),jj,treeNodeArr(jj));
    
    Z(i,:) = [sort([treeNodeArr(ii) treeNodeArr(jj)]) minval];
    
    nodeFlag([ii jj])=false;
    newdist = calPairDist1(bdist,nodeFlag,ii,jj);

    % update distance matrix
    % update ii row
    tmpInds = find(nodeFlag(1:ii-1));
    bdist = updateMatCol(bdist, tmpInds, ii, newdist(tmpInds));
    
    tmpInds = ii+find(nodeFlag(ii+1:end));
    bdist = updateMatRow(bdist, ii, tmpInds, newdist(tmpInds));
        
    % update tree node array
    treeNodeArr(ii) = i+n;
    treeNodeArr(jj) = 0;
    
    nodeFlag([ii jj])=true;
    
    [bii, iii] = getbi(ii);
    [bjj, jjj] = getbi(jj);
    
    for bk = 1:bii-1
        % retrieve data blocks
        [resDist, resPrev, resNext, offset] = retrieveBlocks(bdist, bprev, bnext, bk);
        
        for kk=1:BLOCK_SIZE
            mk = getmi(bk,kk);
            if ~nodeFlag(mk)
                continue;
            else
                [resPrev(kk,:), resNext(kk,:), hedInd(mk), hedVal(mk)] = del2ins1(resDist(kk,:), resPrev(kk,:), resNext(kk,:), offset, hedInd(mk), ii, jj);
            end
        end
        
        bprev(bk,bk:end) = distMatBlock(resPrev);
        bnext(bk,bk:end) = distMatBlock(resNext);
    end
    
    for bk = bii
        [resDist, resPrev, resNext, offset] = retrieveBlocks(bdist, bprev, bnext, bk);
        for kk=1:iii-1
            mk = getmi(bk,kk);
            if ~nodeFlag(mk)
                continue;
            else
                [resPrev(kk,:), resNext(kk,:), hedInd(mk), hedVal(mk)] = del2ins1(resDist(kk,:), resPrev(kk,:), resNext(kk,:), offset, hedInd(mk), ii, jj);
            end
        end
        
        tmpNodeFlag = nodeFlag;
        tmpNodeFlag(jj)=false;
        tmpInds = iii+find(tmpNodeFlag(ii+1:N_NODE));
        % note that distvec, prevvec, nextvec are 1*(n-ii) vector
        [resPrev(iii,tmpInds), resNext(iii,tmpInds), hedInd(ii), hedVal(ii)] = genPointers(resDist(iii,iii+1:end),tmpNodeFlag(ii+1:N_NODE),ii);
        
        if bii==bjj
            endRow=jjj-1;
        else
            endRow=BLOCK_SIZE;
        end
        
        for kk=iii+1:endRow
            mk = getmi(bk,kk);
            if ~nodeFlag(mk)
                continue;
            else
                [resPrev(kk,:), resNext(kk,:), hedInd(mk), hedVal(mk)] = delPointers(resDist(kk,:), resPrev(kk,:), resNext(kk,:), offset, hedInd(mk), jj);
            end
        end
        
        bprev(bk,bk:end) = distMatBlock(resPrev);
        bnext(bk,bk:end) = distMatBlock(resNext);
    end
    
    for bk=bii+1:bjj
        [resDist, resPrev, resNext, offset] = retrieveBlocks(bdist, bprev, bnext, bk);
        
        if bk==bjj
            endRow=jjj-1;
        else
            endRow=BLOCK_SIZE;
        end
        
        for kk=1:endRow
            mk = getmi(bk,kk);
            if ~nodeFlag(mk)
                continue;
            else
                [resPrev(kk,:), resNext(kk,:), hedInd(mk), hedVal(mk)] = delPointers(resDist(kk,:), resPrev(kk,:), resNext(kk,:), offset, hedInd(mk), jj);
            end
        end
        
        bprev(bk,bk:end) = distMatBlock(resPrev);
        bnext(bk,bk:end) = distMatBlock(resNext);
        
    end
    
    nodeFlag(ii)=true;
    nodeFlag(jj)=false;
    
    if jj<n
        hedInd(jj)=DEL_VAL;
        hedVal(jj)=DEL_VAL;
    end
    
%     [hedInd hedVal]
end

function [prevvec, nextvec, tmpHedInd, tmpHedVal] = del2ins1(distvec, prevvec, nextvec, offset, hedInd, ii, jj)
[prevvec, nextvec, tmpHedInd] = delPointers(distvec, prevvec, nextvec, offset, hedInd, ii);
[prevvec, nextvec, tmpHedInd] = delPointers(distvec, prevvec, nextvec, offset, tmpHedInd, jj);
[prevvec, nextvec, tmpHedInd, tmpHedVal] = insertPointers(distvec, prevvec, nextvec, offset, tmpHedInd, ii);

function y = calDist(X, bi, bj)
global BLOCK_SIZE
global N_NODE
assert(bi<=bj);
indsi = getmi(bi,1):getmi(bi,BLOCK_SIZE);
indsj = getmi(bj,1):getmi(bj,BLOCK_SIZE);
y = zeros(BLOCK_SIZE,BLOCK_SIZE);
for i=1:BLOCK_SIZE
    for j=1:BLOCK_SIZE
        ii = indsi(i);
        jj = indsj(j);
        if ii>=jj || ii>N_NODE || jj>N_NODE
            continue
        else
            y(i,j) = sqrt((X(ii,:)-X(jj,:))*(X(ii,:)-X(jj,:))');
        end
    end
end

function [minVal, minInd] = mymin(vec)
global DEL_VAL
inds = find(vec~=DEL_VAL);
[tmpVal, tmpInd]=min(vec(inds));
minVal = tmpVal;
minInd = inds(tmpInd);

function y = calPairDist1(bdist,nodeFlag,i,j)
% global DEL_VAL
veci = extractRow1(bdist,i);
vecj = extractRow1(bdist,j);
y = zeros(size(veci));
y(nodeFlag) = max(veci(nodeFlag),vecj(nodeFlag));

function y = extractRow1(bdist,i)
% global DEL_VAL

global BLOCK_SIZE
global N_BLOCK
global N_NODE

n=N_NODE;
[bi, ii] = getbi(i);
y = zeros(1,n);
% y = zeros(1,n)+DEL_VAL;

k=0;
for bj=1:bi-1
    y(k+1:k+BLOCK_SIZE)=bdist{bj,bi}(:,ii);
    k=k+BLOCK_SIZE;
end

y(k+1:k+ii-1) = bdist{bi,bi}(1:ii-1,ii);
y(k+ii+1:k+BLOCK_SIZE) = bdist{bi,bi}(ii,ii+1:end);
k=k+BLOCK_SIZE;

for bj=bi+1:N_BLOCK
    y(k+1:k+BLOCK_SIZE)=bdist{bi,bj}(ii,:);
    k=k+BLOCK_SIZE;
end

y = y(1:n);

function bmat = updateMatRow(bmat, ri, colInds, vals)
% updates row ri of the distance matrix
assert(length(colInds)==length(vals))
for i=1:length(colInds)
    [bi, bj, ii, jj] = mi2bi(ri,colInds(i));
    bmat{bi,bj}(ii,jj) = vals(i);
end

function bmat = updateMatCol(bmat, rowInds, ci, vals)
% update column ci of the distance matrix
assert(length(rowInds)==length(vals))
for i=1:length(rowInds)
    [bi, bj, ii, jj] = mi2bi(rowInds(i),ci);
    bmat{bi,bj}(ii,jj) = vals(i);
end

function [prevVec, nextVec, hedInd, hedVal] = delPointers(distVec, prevVec, nextVec, offset, hedInd, i)
% delete the ith element from the double linked array
global HED_VAL
global END_VAL
global DEL_VAL

prevind = prevVec(i-offset);
nextind = nextVec(i-offset);

hedVal = distVec(hedInd-offset);

% not necessary, since deleted values only occur on rows with nodeFlag(k)==false
% elements that has been deleted    
% if prevind==DEL_VAL && nextind==DEL_VAL
%     return;
% end

% elements is the single left element
if prevind==HED_VAL && nextind==END_VAL
    hedInd=DEL_VAL;
    hedVal=DEL_VAL;
    prevVec(i-offset)=DEL_VAL;
    nextVec(i-offset)=DEL_VAL;
    return;        
end

% remove the element from the list
if prevind==HED_VAL
    hedInd = nextind;
    hedVal = distVec(hedInd-offset);
    prevVec(nextind-offset) = prevind;
elseif nextind==END_VAL
    nextVec(prevind-offset) = nextind;
else
    prevVec(nextind-offset) = prevind;
    nextVec(prevind-offset) = nextind;
end

function [prevVec, nextVec, hedInd, hedVal] = insertPointers(distVec, prevVec, nextVec, offset, hedInd, i)
% insert the updated value of distVec(i) into the double linked list
global HED_VAL
global END_VAL
global DEL_VAL

targetVal = distVec(i-offset);
curNodeInd = hedInd;

% in case of all elements deleted, insert one new 
if curNodeInd==DEL_VAL
    hedInd=i;
    hedVal=targetVal;
    prevVec(i-offset) = HED_VAL;
    nextVec(i-offset) = END_VAL;
    return;
end

 % insert in the head
if distVec(curNodeInd-offset)>=targetVal
    hedInd=i;
    hedVal=targetVal;
    prevVec(i-offset) = HED_VAL;
    nextVec(i-offset) = curNodeInd;     
    prevVec(curNodeInd-offset) = i;
    return;
end

% insert in the middle or end
prevNodeInd=prevVec(curNodeInd-offset);
assert(prevNodeInd==HED_VAL);
while curNodeInd~=END_VAL && distVec(curNodeInd-offset)<targetVal
    prevNodeInd = curNodeInd;
    curNodeInd = nextVec(curNodeInd-offset);
end

if curNodeInd==END_VAL
    % add to the end
    nextVec(prevNodeInd-offset) = i;
    nextVec(i-offset) = END_VAL;
    prevVec(i-offset) = prevNodeInd;
else
    % add to the middle
    nextVec(prevNodeInd-offset) = i;
    prevVec(curNodeInd-offset) = i;

    nextVec(i-offset) = curNodeInd;
    prevVec(i-offset) = prevNodeInd;
end

% retrieve the distance of the hed element
hedVal = distVec(hedInd-offset);

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

prevVec = zeros(size(idx))+DEL_VAL;
prevVec(idx) = offset + nidx(1:end-2);
try
    prevVec(idx(1)) = HED_VAL;
catch
    keyboard
end

nextVec=zeros(size(idx))+DEL_VAL;
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

function [bi, ii] = getbi(i)
% get block index of original index i
global BLOCK_SIZE
bi = floor(i/BLOCK_SIZE)+1;
ii = mod(i,BLOCK_SIZE);
if ii==0
    ii=BLOCK_SIZE;
    bi = bi -1;
end

function i = getmi(bi, ii)
global BLOCK_SIZE
i = (bi-1)*BLOCK_SIZE + ii;

function [bi, bj, ii, jj] = mi2bi(i,j)
% maxtrix index to block index
% transform the index of the whole matrix into two level index
% block index and index within block

[bi, ii] = getbi(i);
[bj, jj] = getbi(j);

function [i,j] = bi2mi(bi, bj, ii, jj)
% block index to matrix index
% bi, bj block index
% ii, jj, index within block
i = getmi(bi,ii);
j = getmi(bj,jj);

function [resMat, offset] = getMatBlock(bmat,bi)
% concatenate the blocks of bi row
global DEL_VAL
global BLOCK_SIZE
global N_BLOCK

offset = (bi-1)*BLOCK_SIZE;
resMat = zeros(BLOCK_SIZE,BLOCK_SIZE*(N_BLOCK-bi+1))+DEL_VAL;

for i=bi:N_BLOCK
    ii = i-bi+1;
    tmpinds = ((ii-1)*BLOCK_SIZE+1):(ii*BLOCK_SIZE);
    resMat(:,tmpinds) = bmat{bi,i};
end

function resArr = distMatBlock(resMat)
% distribute the combined matrix of bi row into blocks 
global BLOCK_SIZE
% global DEL_VAL

len = size(resMat,2);

nb = len/BLOCK_SIZE;

% if mod(len,BLOCK_SIZE)==0
%     nb = len/BLOCK_SIZE;
% else
%     nb = ceil(len/BLOCK_SIZE);
%     newMat = zeros(BLOCK_SIZE,nb*BLOCK_SIZE)+DEL_VAL;
%     newMat(:,1:len) = resMat;
%     resMat = newMat;
% end

resArr = cell(1,nb);

for i=1:nb
    tmpinds = (i-1)*BLOCK_SIZE+(1:BLOCK_SIZE);
    resArr{i} = resMat(:,tmpinds);
end

function [resDist, resPrev, resNext, offsetD] = retrieveBlocks(bdist, bprev, bnext, bk)
% global N_NODE
% retrieve row bk for the block matrix

[resDist, offsetD] = getMatBlock(bdist,bk);
[resPrev, offsetP] = getMatBlock(bprev,bk);
[resNext, offsetN] = getMatBlock(bnext,bk);
% assert(offsetP==offsetN);
% assert(offsetP==offsetD);

% resDist = resDist(:,1:N_NODE-offsetD);
% resPrev = resPrev(:,1:N_NODE-offsetD);
% resNext = resNext(:,1:N_NODE-offsetD);

