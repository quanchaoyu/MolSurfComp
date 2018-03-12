classdef treenode %%treenode denote the node of the hierarchical tree
    properties
        activenode;
        set;        
        activeelement;
        n1;
        n2;
        leftnode;
        rightnode;
    end
    
     methods
       function obj=treenode(activenode,set,activeelement,n1,n2)
           obj.activenode=activenode;
           obj.set=set;
           obj.activeelement=activeelement;
           obj.n1=n1;
           obj.n2=n2;
           obj.leftnode=0;
           obj.rightnode=0;
       end
   end
    
end