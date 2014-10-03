There are multiple directories to process the 3d time migrated data from 
Parihaka New Zealand.  There are multiple volumes on the distributed by Dennis 
Cooke in the Summer of 2014.  Some of these data are: 
PSTM_full_angle.sgy               a 5.1 gbyte subset selected by Dennis Cooke
parihaka\ inline\ cross\ line.png an inline/xline image of the subset     

There are 47-49 gybte data volumes.  These volumes are not on a server for you
to download.  I am trying to extract a subset to make available on a public
server. 
C40502_Parihaka_M3D_Final_PSTM_Far_Stack_IL825-2815.SEGY 
C40502_Parihaka_M3D_Final_PSTM_Mid_Stack_IL824-2813.SEGY 
C40502_Parihaka_M3D_Final_PSTM_Near_Stack_IL824-2811.SEGY 

The subdirectories to load the data in rsf format and apply simple processing:
fetch- fetch from amazon server, convert to rsf
displayin-put fetch file in 3d volume, display an il/xl, display as 3d cube
loadavo - convert farstack to rsf format, extract subset, display il/xl/3dcube
  
