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
        
#    def __iter__(self):
#        #print(f'id={self.id}')
#        if self.leftChild:
#            print('-',end='')
#            yield from self.leftChild
#        yield self
#        if self.rightChild:
#            print('-',end='')
#            yield from self.rightChild
#     
#    def __next__(self):
#        if self is None:
#            raise StopIteration
#        else:
#            print(f'id={self.id}')
             
    
        
def build_tree(Z, n):
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
    global arr, idxarr
    arr = np.arange(n)
    idxarr = np.arange(n)
    set_boarder(currNode, arr, idxarr, 0, n)
    set_depth(currNode,0)
    return currNode

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
        
def disp_tree(root):
    print('  '*root.depth+'-',end='')
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
        print('setting arr[%d] %d => %d '%(left, arr[left], root.id))
        print('setting idxarr[%d] %d => %d '%(root.id, idxarr[left], left))
        print(arr)
        print(idxarr)
        print()
        return
    else:
        root.bleft = left
        root.bright = right
        root.bmid = left + root.leftChild.nLeafNode
        
    set_boarder(root.leftChild, arr, idxarr, left, root.bmid)
    set_boarder(root.rightChild, arr, idxarr, root.bmid, right)            
            
def subtree(root, subidx, arr):
    # subidx = idxarr[subarr]
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
        
        if len(leftIdx)<=len(rightIdx):
            retNode.leftChild = subtree(root.leftChild, leftIdx, arr)
            retNode.rightChild = subtree(root.rightChild, rightIdx, arr)
        else:
            retNode.rightChild = subtree(root.leftChild, leftIdx, arr)
            retNode.leftChild = subtree(root.rightChild, rightIdx, arr)
        return retNode

                
    
from scipy.cluster.hierarchy import linkage
np.random.seed(0)
n = 6
arr = list()
idxarr = list()
x = np.random.randn(n,4)
Z = linkage(x,'complete')
r = build_tree(Z,n)

disp_tree(r)

subarr=[1,2,3,4]
st = subtree(r,idxarr[subarr],arr)
set_depth(st,0)
disp_tree(st)
    