classdef activefront

   properties
        meshed;
        Ae;
        Nae;
       
   end 
   
   methods
       function obj=activefront(meshed,Ae,Nae)
           obj.meshed=meshed;
           obj.Ae=Ae;
           obj.Nae=Nae;
       end
   end
  
end