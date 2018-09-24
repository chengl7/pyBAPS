classdef editPool
   properties
      N_BLOCK
      BLOCK_SIZE
      DEL_VAL
      
      normEdit % cell array, each element is edit queue for del, insertion 
      rowEdit % cell array, each element is edit queue for whole row insertion
      
      editFlag % whether a block is editted
      insertRow % index of the insertion row (within block)
   end
   
   methods
       function obj = editPool(nblock, blocksize, delval)
           obj.N_BLOCK = nblock;
           obj.BLOCK_SIZE = blocksize;
           obj.DEL_VAL = delval;
           
           obj.normEdit = cell(1,nblock);
           obj.editFlag = false(1,nblock); % 0-no update for this block, 1-update for this block
           
           obj.insertRow = delval;
           obj.rowEdit = cell(1,nblock);
           
           for i=1:nblock
               obj.normEdit{i} = editQueue(blocksize,4);
               obj.rowEdit{i} = editQueue(1,blocksize);
           end
           
       end
       
       function obj = clear(obj,bi)
           obj.insertRow = obj.DEL_VAL;
           for i=bi:obj.N_BLOCK
               obj.normEdit{i} = clear(obj.normEdit{i});
               obj.editFlag(i) = false;
               
               obj.rowEdit{i} = clear(obj.rowEdit{i});
           end
       end
       
       function obj = insertEditRep(obj, rowInd, ind, val)
           [bi, ii] = getbi(obj,ind);
           obj.normEdit{bi} = insertRep(obj.normEdit{bi}, rowInd, ii, val);
           obj.editFlag(bi) = true;
       end
       
       function obj = insertRowEdit(obj, rowInd, ind, val)
           [bi, ii] = getbi(obj,ind);
           if obj.insertRow==obj.DEL_VAL
               obj.insertRow = rowInd;
           end
           if ~obj.editFlag(bi)
               obj.editFlag(bi)=true;
           end
           obj.rowEdit{bi} = insert(obj.rowEdit{bi}, 1, ii, val);
       end
       
       function obj = sortEdit(obj,bi)
           for bk=bi:obj.N_BLOCK
              if obj.editFlag(bk)
%                   continue
%               else
                  obj.normEdit{bk} = sort(obj.normEdit{bk}); 
              end
           end
       end
       
       function [bi, ii] = getbi(obj,i)
           % get block index of original index i
            n = obj.BLOCK_SIZE;
            bi = floor(i/n)+1;
            ii = mod(i,n);
            if ii==0
                ii=n;
                bi = bi -1;
            end
       end
   end
end