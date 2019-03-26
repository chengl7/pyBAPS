# -*- coding: utf-8 -*-
"""
Dendrogram operations
"""
import numpy as np
class TreeNode:    
    def __init__(self, id=-1):
        self.leftChild = None
        self.rightChild = None
        self.id = id
        self.height = 0
        self.nLeafNode = 1
        self.bleft = 0
        self.bmid = 0
        self.bright = 0
        self.depth = 0
    
        
def build_tree(Z):
    n = Z.shape[0]+1
    tmplist = []
    for i in range(n-1):
        currNode = TreeNode()
        lnode,rnode,dist,_ = Z[i,:]
        lnode = int(lnode)
        rnode = int(rnode)
        currNode.leftChild = get_node(lnode, n, tmplist)
        currNode.rightChild = get_node(rnode, n, tmplist)
        currNode.id = i+n
        currNode.height = dist
        currNode.nLeafNode = currNode.leftChild.nLeafNode + currNode.rightChild.nLeafNode
        tmplist.append(currNode)
#    global arr, idxarr
    arr = np.zeros(n,dtype='i')
    idxarr = np.zeros(n,dtype='i')
    set_boarder(currNode, arr, idxarr, 0, n)
    set_depth(currNode,0)
    return currNode,arr,idxarr

def get_node(inode, n, innerList):
    if inode<n:
        return TreeNode(inode)
    first = 0
    last = len(innerList)
    while first<last:
        mid = (first+last) // 2
        item = innerList[mid].id
        if item == inode:
            retNode = innerList[mid]
            del innerList[mid]
            return retNode
        elif inode < item:
            last = mid
        else:
            first = mid

def bin_search(mylist,inode):
    first = 0
    last = len(mylist)
    while first<last:
        mid = (first+last) // 2
        item = mylist[mid]
        if item == inode:
            print(f'{inode} is {mylist[mid]}')
            return
        elif inode < item:
            last = mid
        else:
            first = mid  
            
def set_depth(root,depth):
    if root:
        root.depth = depth
        set_depth(root.leftChild, depth+1)
        set_depth(root.rightChild, depth+1)
        return root

def leaf_generator(root):
    if root.leftChild is None or root.rightChild is None:
        yield root
    if root.leftChild:
        yield from leaf_generator(root.leftChild)
    if root.rightChild:
        yield from leaf_generator(root.rightChild)
    
        
def inner_generator(root):
    if root.leftChild is None and root.rightChild is None:
        return root
    else:
        yield root
    if root.leftChild:
        yield from inner_generator(root.leftChild)
    if root.rightChild:
        yield from inner_generator(root.rightChild)   

# generate the leaf nodes of all inner nodes of height lower than cutoff        
def inner_generator_leaf(root,cutoff):
    if root.height<=cutoff:
        res = np.array([n.id for n in leaf_generator(root)])
        yield res
        return
    if root.leftChild:
        yield from inner_generator_leaf(root.leftChild,cutoff)
    if root.rightChild:
        yield from inner_generator_leaf(root.rightChild,cutoff)

        
def disp_tree(root):
    print('  '*root.depth+'|-',end='')
    print(root.id)
    if root.leftChild:
        disp_tree(root.leftChild)
    if root.rightChild:
        disp_tree(root.rightChild)     
            
def set_boarder(root, arr, idxarr, left, right):
    if root.leftChild is None or root.rightChild is None:
        root.bleft = left
        root.bright = right
        root.bmid = right
        
        
        arr[left] = root.id
        idxarr[root.id] = left
#        print('setting arr[%d]=%d => %d '%(left, arr[left], root.id))
#        print('setting idxarr[%d]=%d => %d '%(root.id, idxarr[left], left))
#        print(arr)
#        print(idxarr)
#        print()
        return
    else:
        root.bleft = left
        root.bright = right
        root.bmid = left + root.leftChild.nLeafNode
        
    set_boarder(root.leftChild, arr, idxarr, left, root.bmid)
    set_boarder(root.rightChild, arr, idxarr, root.bmid, right)            
            
# e.g. st = subtree(r,idxarr[subarr],arr)
def subtree(root, subidx, arr):
    # subidx = idxarr[subarr]
    assert len(np.unique(subidx)) == len(subidx), f'input subarr not unique {arr[subidx]}'
    if len(subidx)==0:
        return None
    if len(subidx)==1:
        return TreeNode(int(arr[subidx]))
        
    flag = subidx<root.bmid
    leftIdx = subidx[flag]
    rightIdx = subidx[~flag]    
    
    if len(leftIdx)==0:
        return subtree(root.rightChild, rightIdx, arr)
    elif len(rightIdx)==0:
        return subtree(root.leftChild, leftIdx, arr)
    else:
        retNode = TreeNode()
        retNode.height = root.height
        retNode.id = root.id
        retNode.nLeafNode = len(subidx)
        
        retNode.leftChild = subtree(root.leftChild, leftIdx, arr)
        retNode.rightChild = subtree(root.rightChild, rightIdx, arr)
        
#        if len(leftIdx)<=len(rightIdx):
#            retNode.leftChild = subtree(root.leftChild, leftIdx, arr)
#            retNode.rightChild = subtree(root.rightChild, rightIdx, arr)
#        else:
#            retNode.rightChild = subtree(root.leftChild, leftIdx, arr)
#            retNode.leftChild = subtree(root.rightChild, rightIdx, arr)
        return retNode
    
# return the heights of all inner nodes in ascending order
def get_inner_heights(root):
    reslist = []
    for inode in inner_generator(root):
        reslist.append(inode.height)
    res = np.sort(reslist)
    return res

def tree_split_cutoff(root, cutoffPerctile=None):
    if cutoffPerctile:
        cutoff = cutoffPerctile
    else:
        cutoff = np.random.rand(1)*30+20
    heightArr = get_inner_heights(root)
    heightCutoff = np.percentile(heightArr,cutoff)
    return [i for i in inner_generator_leaf(root,heightCutoff)]
    
def tree_split_two(root):
    lres = np.array([n.id for n in leaf_generator(root.leftChild)])
    rres = np.array([n.id for n in leaf_generator(root.rightChild)])
    return (lres,rres)    

# count the number of outliers given the cutoff in cutoffArr
# pass cutoffArr, cntArr by reference, flagArr by copy
def count_outlier(root, cutoffArr, cntArr, flagArr):
    if not any(flagArr):
        return
    if root.leftChild is None and root.rightChild is None:
        cntArr[flagArr] += 1
#        print('rootid=',root.id,'cntArr=',cntArr)
        return
    currFlagArr = root.height>=cutoffArr
    if any(currFlagArr):
        newFlagArr = np.logical_and(currFlagArr,flagArr)
        count_outlier(root.leftChild, cutoffArr, cntArr, newFlagArr)
        count_outlier(root.rightChild, cutoffArr, cntArr, newFlagArr)
        
def outlier_gen(root, cutoff):
    if root.leftChild is None and root.rightChild is None:
        yield root.id
        return
    if root.height>cutoff:
        yield from outlier_gen(root.leftChild, cutoff)
        yield from outlier_gen(root.rightChild, cutoff)
        
def get_tree_outliers(root):
    heightArr = get_inner_heights(root)
    # 0.5(1-x)n < i < (1-x)n, where x is the expected outlier percentage
    # set the range as (0.5*(1-x)*100, (1-x)*100+1)
    heightCutoffArr = np.percentile(heightArr,np.arange(45,91,5))  # get around 10% outliers
    cntArr = np.zeros(len(heightCutoffArr),dtype='i')
    flagArr = np.ones(len(heightCutoffArr),dtype=bool)
    count_outlier(root, heightCutoffArr, cntArr, flagArr)
    mind = np.argmin(np.abs(cntArr - 0.1*root.nLeafNode))  # checked that subtree root nLeafNode is correct
    return np.sort([i for i in outlier_gen(root, heightCutoffArr[mind])])

def tree_split_k(root, k):
    assert not root is None
    n = root.nLeafNode
    assert k<n
    heightArr = get_inner_heights(root)
    cutoff = heightArr[n-k-1]
    return [i for i in inner_generator_leaf(root,cutoff)]

# think about trees that are not monotonic, e.g. https://uk.mathworks.com/help/stats/linkage.html#Tips
# does not seem to cause any code to break down, 25.03.2019
# tested many trees constructed using centroid, ward, median and did not see any brokendown case

if __name__=="__main__":                
    
    from scipy.cluster.hierarchy import linkage
#    np.random.seed(0)
    n=10
    #n = 100
    #arr = list()
    #idxarr = list()
    #arr = np.zeros(n,dtype='i')
    #idxarr = np.zeros(n,dtype='i')
#    np.random.seed(0)
    x = np.random.randn(n,4)
#    Z = linkage(x,'complete')
    Z = linkage(x,method='median')
    r,arr,idxarr = build_tree(Z)
    
    print('Constructed tree matrix')
    print(Z)
    print()
    
    print('display constucted tree')
    disp_tree(r)
    print()
    
#    subarr=[1,2,3,4]
    k = round(n*0.3)
    subarr=  np.sort(np.random.choice(n,k,replace=False))
    assert len(np.unique(subarr)==len(subarr))
    assert len(subarr)<=n
    st = subtree(r,idxarr[subarr],arr)
    set_depth(st,0)
    print(f'display subtree for subarr={subarr}')
    disp_tree(st)
    print()
    
    print('display leaf nodes')
    for ln in leaf_generator(r):
        print(ln.id, end=' ')
    print()
    
    print('display inner nodes')
    for ln in inner_generator(r):
        print(ln.id, end=' ')
    print()
    
    print('display inner nodes (subtree)')
    for ln in inner_generator(st):
        print(ln.id,'-',ln.height)
    print()
    
    print('cut the tree at 30% percentile')    
    print(tree_split_cutoff(r,30))
    print()
    
    print('Split the tree into two')
    print(tree_split_two(r))
    print()
    
    k = 9
    print(f'Split the tree into k={k} clusters')
    print(tree_split_k(r, k))
    print()
    
    print('testing cutting the tree at a given set of cutoffs')
    cutoffArr = np.arange(0.9,5.5,0.5)
    cntArr = np.zeros(len(cutoffArr))
    flagArr = np.ones(len(cutoffArr),dtype=bool)
    count_outlier(r, cutoffArr, cntArr, flagArr)
    print(cntArr)
    print([i for i in outlier_gen(r, 1)])
    print([i for i in outlier_gen(r, 1.4)])
    print()
    
    print('getting tree outliers')
    print(get_tree_outliers(r))
    
    print('getting tree outliers (substree)')
    print(get_tree_outliers(st))