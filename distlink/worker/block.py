import numpy as np
from distlink.common.constants import Constants
from distlink.common.misc import MinTurple
from distlink.common.server import Server

class Block:
    #! JS: should the block functions be taking global coordinates and translating them in every function?
    """Represents a block matrix and associated variables.

    Attributes:
        bmat (np.array): bs x bs, block matrix, bs stands for block size
        browflag (np.array): block row flag, indicate if an element deleted or not
        bcolflag (np.array): block column flag, indicate if an element deleted or not
        hedInd (np.array): bs x 1, index of min value in each row
        hedVal (np.array): bs x 1, minimum value of each row
        bi (int): block row index
        bj (int): block column index 
        minInd (tuple): (ri, ci) row and col index of the minimum value in bmat
        minHedVal (constants.DATA_TYPE): minimum value in bmat 
        count (int): number of elements in this block
    """
    
    def __init__(self, bi, bj):
        """Initialize Block with two indices, loading block data from file.

        Args:
            bi (int): first index
            bj (int): second index
        """
        assert(bi<=bj)
        self.bmat = np.load(Constants.get_dist_block_file(bi,bj))
        print(self.bmat)
        
        self.browflag = np.ones(Constants.BLOCK_SIZE, dtype=bool)
        self.bcolflag = np.ones(Constants.BLOCK_SIZE, dtype=bool)
        if bi==Constants.N_BLOCK-1 and (Constants.N_NODE % Constants.BLOCK_SIZE)!=0:
            ii = Constants.N_NODE % Constants.BLOCK_SIZE
            self.browflag[ii:]=False
        if bj==Constants.N_BLOCK-1 and (Constants.N_NODE % Constants.BLOCK_SIZE)!=0:
            jj = Constants.N_NODE % Constants.BLOCK_SIZE
            self.bcolflag[jj:]=False
        if bi==bj:
            self.browflag[-1]=False
            self.bcolflag[0]=False
        self.count = int(np.sum(self.browflag)) 
        
        self.minHedVal = Constants.DEL_VAL_DIST
        self.minInd = (Constants.DEL_VAL_IND,Constants.DEL_VAL_IND)
        
        self.bi = bi
        self.bj = bj
        
        self.hedInd = np.zeros(Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE)
        self.hedVal = np.zeros(Constants.BLOCK_SIZE, dtype=Constants.DATA_TYPE_DIST)
        self.update_batch_head(range(Constants.BLOCK_SIZE))   # minInd, minHedVal updated
    
    # return a tuple that is the minimum, should be call when minHedVal is not DEL_VAL
    def get_min_tuple(self):
        """Calculate and return the tuple with minimum value."""
        if self.minHedVal==Constants.DEL_VAL_DIST:
#            delval = Constants.DEL_VAL
            return MinTurple(Constants.DEL_VAL_IND,Constants.DEL_VAL_IND,Constants.DEL_VAL_DIST)
        ii,jj=self.minInd
        mi = Constants.getmi(self.bi, ii)
        mj = Constants.getmi(self.bj, jj)
        return MinTurple(mi,mj,self.minHedVal)
        
    def __le__(self,obj):
        if self.minHedVal==Constants.DEL_VAL_DIST:
            return False
        elif obj.minHedVal==Constants.DEL_VAL_DIST:
            return True
        else:
            return self.minHedVal <= obj.minHedVal
             
    # extract the block portion of global row xi, delete row xi at the same time
    def extract_row(self,xi):
        """Extract a segment of global row xi from this block."""
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            return self.bmat[:,xii]
        elif self.bi==xbi and self.bj==xbi:
            tmparr = np.copy(self.bmat[:,xii])
            tmparr[xii+1:] = self.bmat[xii,xii+1:]
            return tmparr
        elif self.bi==xbi and self.bj>xbi:
            return self.bmat[xii,:]
        else:
            return None
        
    # assign the block portion of global row xi
    def assign_row(self,xi,arr):
        """Assign a row of this block corresponding to global row xi."""
        assert len(arr)==Constants.BLOCK_SIZE
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.bmat[:,xii]=arr
        elif self.bi==xbi and self.bj==xbi:
            self.bmat[:xii,xii] = arr[:xii]
            self.bmat[xii,xii+1:] = arr[xii+1:]
        elif self.bi==xbi and self.bj>xbi:
            self.bmat[xii,:] = arr
        else:
            return None
        
    # update hedInd and hedVal for block row i    
    def update_head(self,i):
        """Find the minimum j and value for row (point) i.

        Args:
            i: local index specifying data element.
        """
        if not self.browflag[i]:
            self.hedInd[i] = Constants.DEL_VAL_IND
            self.hedVal[i] = Constants.DEL_VAL_DIST
            return
        if self.bi!=self.bj:
            valinds = np.where(self.bcolflag)[0]
        elif self.bi==self.bj:
            valinds = i+1+np.where(self.bcolflag[i+1:])[0]
        
        if valinds.size==0:
            self.hedInd[i] = Constants.DEL_VAL_IND
            self.hedVal[i] = Constants.DEL_VAL_DIST
        else:
            tmpminind = np.argmin(self.bmat[i,valinds])
            self.hedInd[i] = valinds[tmpminind]
            self.hedVal[i] = self.bmat[i,self.hedInd[i]]
    
    # update hedInd and hedVal for index (block) in batch, batch is iterable
    def update_batch_head(self,batch,minFlag=True):
        """Call update_head for a list of indices.
        
        Args:
            batch: list of indices.
            minFlag (bool): whether or not to update_min_head.    
        """
        for i in batch:
            self.update_head(i)
        if minFlag:
            self.update_min_head()    
                
    # update minInd, minHedVal
    def update_min_head(self):
        """Find the row i with the closest neighbor j."""
        if self.count<=0:
            assert not any(self.browflag)
            self.minHedVal = Constants.DEL_VAL_DIST
            self.minInd = (Constants.DEL_VAL_IND,Constants.DEL_VAL_IND)
        else:
            valinds = np.where(self.browflag)[0]
            rowind = np.argmin(self.hedVal[valinds])
            rowind = valinds[rowind]
            rowind = getattr(np,Constants.DATA_TYPE)(rowind)  # convert int to numpy int
            self.minHedVal = self.hedVal[rowind]
            self.minInd = (rowind,self.hedInd[rowind])
                    
    # delete row xi in the global matrix, update browflag, bcolflag, 
    def delete_row(self,xi):
        """Delete global row xi from the block and update_batch_head."""
        #! JS: both deletes and updates. Is this too much responsibility?
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.bcolflag[xii]=False
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            if self.browflag[xii]:
                self.browflag[xii]=False
                self.count -= 1
            self.bcolflag[xii]=False
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.browflag[xii]=False
            self.count -= 1
            self.update_batch_head(range(xii,xii+1))
        else:
            pass
            
    # insert row xi (global matrix), browflag: from global flag matrix
    def insert_row(self,xi):
        #! JS: does this insert anything? Seems to update not insert.
        #! JS: perhaps different language would improve clarity.
        """Call update_batch_head for row xi."""
        xbi,xii = Constants.getbi(xi)
        if self.bi<xbi and self.bj==xbi:
            self.update_batch_head(range(Constants.BLOCK_SIZE))
        elif self.bi==xbi and self.bj==xbi:
            self.update_batch_head(range(xii+1))
        elif self.bi==xbi and self.bj>xbi:
            self.update_batch_head(range(xii,xii+1))
        else:
            pass

class BlockProcess(Server):
    """BlockProcess server class for operating directly on blocks.

    Attributes:
        block: Block object
    """
    def __init__(self, parentConn=None, serverName=None, blockIndex=None, logFile=None):
        """Initialize with optional server connection to parent.

        Keyword args:
            parentConn: connection to parent.
            serverName: name of this server.
            blockIndex: index of contained block.
            logFile: file to which logs are written.
        """
        super().__init__(parentConn=parentConn, serverName=serverName, logFile=logFile)
        
        ######### special for block process ###########
        self.block=Block(blockIndex[0],blockIndex[1])
        
    def update_child_block_list(self):
        """Fetches a list of one element, the (bi,bj) tuple."""
        # JS: does not update
        return [(self.block.bi,self.block.bj)]
    
    def del_blocks(self,bi):
        """Log that block bi has been deleted, report None to parent."""
        # JS: see above comments.
        self.log('debug','block (%d,%d) deleted.' % (self.block.bi,self.block.bj))
        self.parentConn.send(None)
        self.log('debug','Block process (%d,%d) closed.' % (self.block.bi,self.block.bj))
        self.close()
    
    def get_min(self):
        """Gets the stored block's min tuple (mi,mj,d)."""
        # JS: see above comments.
        self.log('debug','minimal value in block (%d,%d) is %s.' % (self.block.bi, self.block.bj, str(self.block.get_min_tuple()) ))
        return self.block.get_min_tuple()
    
    def ex_row(self, xi, delFlag=False):
        """Extracts row xi from stored block.

        Args:
            xi (int): global row index xi.
            delFlag (bool): whether to also delete the row xi !!!
        """
        arr = self.block.extract_row(xi)
        assert arr is not None
        if delFlag:
            self.block.delete_row(xi)
        xbi,xii = Constants.getbi(xi)
        k = self.block.bi if self.block.bi!=xbi else self.block.bj
        self.log('debug','block (%d,%d) for xi=%d (bi=%d) k=%d extracted.' % (self.block.bi,self.block.bj, xi, xbi, k))
        return [(k,arr)]
    
    def ins_row(self, xi, segList):
        """Inserts row xi with data contained in segList.
        
        Args:
            xi (int): global row index xi.
            segList: list containing a single element (index, array)!!!
        """
        assert(len(segList)==1)
        ind,arr = segList[0]
        assert(self.block.bi==ind[0] and self.block.bj==ind[1])
        
        self.block.assign_row(xi,arr)
        self.block.insert_row(xi)
        return None


