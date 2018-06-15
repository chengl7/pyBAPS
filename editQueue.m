classdef editQueue
   properties
      % NUM * DIM matrix, each row is a queue
      NUM
      DIM
      
      index
      value
      pointer
   end
   methods
       function obj = editQueue(num, dim)
           obj.NUM = num;
           obj.DIM = dim;
           
           obj.index = zeros(num,dim);
           obj.value = zeros(num,dim);
           obj.pointer = zeros(num,1);
       end
       
       function obj = clear(obj)
           obj.pointer(:)=0;
       end
       
       function obj = insert(obj, rowind, ind, val)
           tmppointer = obj.pointer(rowind);
           tmppointer = tmppointer+1;
           
           obj.pointer(rowind) = tmppointer;
           obj.index(rowind,tmppointer) = ind;
           obj.value(rowind,tmppointer) = val;
       end
       
       function obj = insertRep(obj, rowind, ind, val)
           tmppointer = obj.pointer(rowind);
           for i=1:tmppointer
               if obj.index(rowind,i)==ind
                   obj.value(rowind,i) = val;
                   return;
               end
           end
           obj = insert(obj, rowind, ind, val);
       end
       
       function obj = sort(obj)
           for ri=1:obj.NUM
               tmppointer = obj.pointer(ri);
               if tmppointer==0
                   continue
               else
                   [tmparr, tmpidx] = sort(obj.index(ri,1:tmppointer));
                   obj.index(ri,1:tmppointer) = tmparr;
                   obj.value(ri,1:tmppointer) = obj.value(ri,tmpidx);
               end
           end
       end

   end
end