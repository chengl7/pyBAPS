from numpy import array,zeros,sort,arange,argsort,append,unique
from numpy.random import randint

class HierPartitionNode:
    """
    This class implement a node class of hierarchical clustering class
    """
    type='int32'
    
    def __init__(self, sorted_index, group_boarder_index=None, previous_node=None):
        self.sorted_index = array(sorted_index,dtype=self.type)
        if group_boarder_index is None:
            self.group_boarder_index = append(sorted_index,len(self.sorted_index))
        else:
            self.group_boarder_index = array(group_boarder_index,dtype=self.type)
        self.previous_node = previous_node
        self.n_group = len(self.group_boarder_index)-1
        
    def is_base_partition(self):
        return self.previous_node is None
    
    def rand_sample(self,group_arr):
        # sample a random index from each group given in group_arr
        ngroup = len(group_arr)
        rand_ind_arr = zeros(ngroup,dtype=self.type)
        for i in range(0,ngroup):
            tmps = self.group_boarder_index[group_arr[i]]
            tmpe = self.group_boarder_index[group_arr[i]+1]
            tmpi = randint(tmps,tmpe);
            rand_ind_arr[i]=self.sorted_index[tmpi]
        return sort(rand_ind_arr,kind='mergesort')    
    
    def gen_partition(self,hash_key_arr):
        # generate a partition of the current groups
        # hash_key_arr is the whole hash key array for all input sequences, numpy array
        
        # extract a random key for each group
        curr_part = self
        tmpinds = arange(0,self.n_group)
        while not curr_part.is_base_partition():
            tmpinds = curr_part.rand_sample(tmpinds)
            curr_part = curr_part.previous_node
        
        # this operations is not needed since the base partition is one to one map
        # tmpinds = curr_part.rand_sample(tmpinds)
        
        key_arr = hash_key_arr[tmpinds]
        sorted_index = argsort(key_arr,kind='mergesort')
        unival,uniind=unique(key_arr[sorted_index],return_index=True)
        group_boarder_index = append(uniind,len(sorted_index))
        #print(key_arr[sorted_index],unival)
        
        return HierPartitionNode(sorted_index, group_boarder_index, self)
        
        
if __name__ == "__main__":
   nseq=10000
   index = arange(0,nseq)

   hashkey = randint(1,20,size=nseq)
   hashkey1 = randint(1,10,size=nseq)
   hashkey2 = randint(1,10,size=nseq)
   
   base_part=HierPartitionNode(index)
   a1 = base_part.gen_partition(hash_key_arr=hashkey)
   a2 = a1.gen_partition(hash_key_arr=hashkey1)
   a3 = a2.gen_partition(hash_key_arr=hashkey1)




        
            
        
        
        
        
                
            