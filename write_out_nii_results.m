function write_out_nii_results(nii_hdr,Data_array,filename,processed_dir)


fname=[processed_dir filename];

Data_array(isnan(Data_array))=0;
Data_array(le(Data_array,0))=0;

nii.img=Data_array;
nii.fileprefix=fname;

nii_hdr.dime.datatype=64;
nii_hdr.dime.bitpix=64;
nii_hdr.dime.dim(1)=3;
nii_hdr.dime.dim(5)=1;

nii.hdr=nii_hdr;

save_nii(nii, fname);


end