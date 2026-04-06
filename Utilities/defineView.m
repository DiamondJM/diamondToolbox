function viewStruct = defineView

viewStruct = struct; 

[az, el] = view; 
viewStruct.az = az; viewStruct.el = el; 
tempXlim = xlim; tempYlim = ylim; tempZlim = zlim;


viewStruct.xlim = tempXlim; viewStruct.ylim = tempYlim; viewStruct.zlim = tempZlim; 

end