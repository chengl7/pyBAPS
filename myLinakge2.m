function Z = myLinakge2(X)
% design an n^2log(n) linkage algorithm that is space effient, so not many
% elements need to be moved

% when modeling priority queues
% merge delete and insertion opertations line by line

% Lu Cheng
% 6.6.2018

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
nn = N_BLOCK*BLOCK_SIZE;

% calculate distance matrix, only upper triangle is used
% dist= zeros(n,n)+DEL_VAL;
% dist= zeros(nn,nn)+DEL_VAL;
dist = zeros(nn,nn);
for i=1:n
    for j=i+1:n
        dist(i,j) = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
    end
end
nodeFlag = true(1,n);

% sort each row, and generate the double direction link list
hedInd = zeros(n-1,1);
hedVal = zeros(n-1,1);

bprev = cell(N_BLOCK,N_BLOCK);
bnext = cell(N_BLOCK,N_BLOCK);
bdist = cell(N_BLOCK,N_BLOCK);
for bi=1:N_BLOCK
    for bj=bi:N_BLOCK
        [tmpul1, tmpul2]=bi2mi(bi,bj,1,1);
        [tmplr1, tmplr2]=bi2mi(bi,bj,BLOCK_SIZE,BLOCK_SIZE);
        
        bdist{bi,bj} = dist(tmpul1:tmplr1,tmpul2:tmplr2);
    end
end

nodeBlockIndex = zeros(1,n);
for i=1:n
    nodeBlockIndex(i) = getbi(i);
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

% % check if the block style data agree with original method
% for bi=1:N_BLOCK
%     for bj=bi:N_BLOCK
%         [tmpul1, tmpul2]=bi2mi(bi,bj,1,1);
%         [tmplr1, tmplr2]=bi2mi(bi,bj,BLOCK_SIZE,BLOCK_SIZE);
%         try
%         assert(all(all(bprev{bi,bj} == prev(tmpul1:tmplr1,tmpul2:tmplr2))));
%         assert(all(all(bnext{bi,bj} == next(tmpul1:tmplr1,tmpul2:tmplr2))));
%         catch
%             [bi bj]
%             keyboard
%         end
%     end
% end

% printHedInfo(hedInd,prev,next)

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
    
    for k=1:n-1
        
        if ~nodeFlag(k)
            continue;
        end
        
        if k>=1 && k<ii
%             tmpInds = k+find(nodeFlag(k+1:n));
            tmpInds = k+1:n;
            tmpvec = exMatRow(bprev, k, tmpInds);
            prevvec = zeros(1,n);
            prevvec(tmpInds) = tmpvec;
            tmpvec = exMatRow(bnext, k, tmpInds);
            nextvec = zeros(1,n);
            nextvec(tmpInds) = tmpvec;
            tmpvec = exMatRow(bdist, k, tmpInds);
            distvec = zeros(1,n);
            distvec(tmpInds) = tmpvec;
            assert(all(dist(k,tmpInds)-distvec(tmpInds)<1e-4))


            [prevvec, nextvec, tmpHedInd] = delPointers(prevvec, nextvec, hedInd(k), ii);
            [prevvec, nextvec, tmpHedInd] = delPointers(prevvec, nextvec, tmpHedInd, jj);
            [prevvec, nextvec, tmpHedInd, tmpHedVal] = insertPointers(distvec, prevvec, nextvec, tmpHedInd, ii);
            
            bprev = updateMatRow(bprev, k, tmpInds, prevvec(tmpInds));
            bnext = updateMatRow(bnext, k, tmpInds, nextvec(tmpInds));
                        
            hedInd(k) = tmpHedInd;
            hedVal(k) = tmpHedVal;
            
        elseif k==ii
            tmpNodeFlag = nodeFlag;
            tmpNodeFlag(jj)=false;

            distvec = exMatRow(bdist, k, ii+1:n);
            % note that distvec, prevvec, nextvec are 1*(n-ii) vector
            [prevvec, nextvec, hedInd(ii), hedVal(ii)] = genPointers(distvec,tmpNodeFlag(ii+1:n),ii);
            tmpInds = ii+find(tmpNodeFlag(ii+1:end));
            
            bprev = updateMatRow(bprev, k, tmpInds, prevvec);
            bnext = updateMatRow(bnext, k, tmpInds, nextvec);
            
        elseif k>ii && k<jj
%             tmpInds = k+find(nodeFlag(k+1:n));
            tmpInds = k+1:n;
            tmpvec = exMatRow(bprev, k, tmpInds);
            prevvec = zeros(1,n);
            prevvec(tmpInds) = tmpvec;
            tmpvec = exMatRow(bnext, k, tmpInds);
            nextvec = zeros(1,n);
            nextvec(tmpInds) = tmpvec;
            
            [prevvec, nextvec, tmpHedInd] = delPointers(prevvec, nextvec, hedInd(k), jj);
            bprev = updateMatRow(bprev, k, tmpInds, prevvec(tmpInds));
            bnext = updateMatRow(bnext, k, tmpInds, nextvec(tmpInds));
            
            if tmpHedInd~=hedInd(k)
                if tmpHedInd==DEL_VAL
                    hedInd(k) = DEL_VAL;
                    hedVal(k) = DEL_VAL;
                else
                    hedInd(k) = tmpHedInd;
                    [tbi, tbj, tii, tjj] = mi2bi(k,hedInd(k));
                    hedVal(k) = bdist{tbi,tbj}(tii,tjj);
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
    
%     [hedInd hedVal]
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

function vals = exMatRow(bmat, ri, colInds)
% extract row ri of the matrix
vals = zeros(1,length(colInds));
for i=1:length(colInds)
    [bi, bj, ii, jj] = mi2bi(ri,colInds(i));
    vals(i) = bmat{bi,bj}(ii,jj);
end

function bmat = updateMatRow(bmat, ri, colInds, vals)
% updates row ri of the distance matrix
assert(length(colInds)==length(vals))
for i=1:length(colInds)
    [bi, bj, ii, jj] = mi2bi(ri,colInds(i));
    bmat{bi,bj}(ii,jj) = vals(i);
end

function vals = exMatCol(bmat, rowInds, ci)
% extract column ci of the distance matrix
vals = zeros(1,length(rowInds));
for i=1:length(rowInds)
    [bi, bj, ii, jj] = mi2bi(rowInds(i),ci);
    vals(i) = bmat{bi,bj}(ii,jj);
end

function bmat = updateMatCol(bmat, rowInds, ci, vals)
% update column ci of the distance matrix
assert(length(rowInds)==length(vals))
for i=1:length(rowInds)
    [bi, bj, ii, jj] = mi2bi(rowInds(i),ci);
    bmat{bi,bj}(ii,jj) = vals(i);
end


function y = calPairDist(distMat,nodeFlag,i,j)
% global DEL_VAL
veci = extractRow(distMat,i);
vecj = extractRow(distMat,j);
% y = zeros(size(veci))+DEL_VAL;
y = zeros(size(veci));
y(nodeFlag) = max(veci(nodeFlag),vecj(nodeFlag));
% y(~nodeFlag) = DEL_VAL;


function y = extractRow(distMat, i)
global N_NODE
y = distMat(i,1:N_NODE);
% y = distMat(i,:);
y(1:i-1) = distMat(1:i-1,i)';

function [prevVec, nextVec, hedInd] = delPointers(prevVec, nextVec, hedInd, i)
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


function [prevVec, nextVec, hedInd, hedVal] = insertPointers(distVec, prevVec, nextVec, hedInd, i)
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

len = size(resMat,2);
nb = len/BLOCK_SIZE;

resArr = cell(1,nb);

for i=1:nb
    tmpinds = (i-1)*BLOCK_SIZE+(1:BLOCK_SIZE);
    resArr{i} = resMat(:,tmpinds);
end
