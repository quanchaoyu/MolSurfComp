function call_pdb2pqr(filename,option)
%
%   Call the pdb2pqr package to obtain position, charge, radius information
%
if (nargin == 1)
    option = 'peoepb';
end

command = ['python Package_PDB2PQR/pdb2pqr/pdb2pqr.py --ff=',option,...
    ' PDB_Files/',filename,'.pdb PQR_Files/',filename,'.pqr'];
[status,cmdout] = system(command);
cmdout

end