function assignView(viewStruct)

az = viewStruct.az; el = viewStruct.el; 
view(az,el); 

tempXlim = viewStruct.xlim; tempYlim = viewStruct.ylim; tempZlim = viewStruct.zlim;

xlim(tempXlim); ylim(tempYlim); zlim(tempZlim);

end

