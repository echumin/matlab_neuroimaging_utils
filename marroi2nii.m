%function marroi2nii(Mroi,Nout,Rspace)
%example: marroi2nii('_-10_36_20_roi.mat','fROI_TEST.nii','../All_fROIs.nii')
%{  
Jenya 11.03.2016
Quick and dirty.
This will take an input marsbar _roi.mat (Mroi) file and convert it into
a nifti volume in the space of the reference image (Nspace). The reference
image needs to be a T1 for which the ROI were made or at least in the space
as the _roi.mat. The output nifti ROI is given the name Nout.
---------------------------------------------------------------------------
This needs marsbar and nifti toolbox in path in order to run.
%}  

addpath('/usr/local/spm12/toolbox/marsbar-0.44')
addpath('/usr/local/nifti_toolbox')

%These are the _roi.mat input files
Mroi = {
'_-10_36_20_roi.mat';
'_10_40_28_roi.mat';
'_-14_-10_46_roi.mat'
}

%These are the names of the output files, in register w/ Mroi
Nout = {
'fROI1.nii';
'fROI2.nii';
'fROI3.nii'
}

%This is the .nii that contains the image space for the ROIs to be saved
%out to (note: must be .nii, not .hdr or .img)
Rspace = ('../All_fROIs.nii')
for x = 1:size(Mroi)

% read in the marsbar roi file as object
load(Mroi{x});

% read in the image in the same space as (reference image) _roi.mat
% this is generally the T1 for which roi's were made
Nspace=MRIread(Rspace);

% pull out the image data from the structure.
%---------------------------- This is new to convert point_place roi
roi_space=mars_space(Rspace);
mo=maroi_matrix(roi,roi_space);
%-------------------------------
dat=matrixdata(mo);

% change the coordinate system from marsbar to nifti format.
%transpose the x and y dimensions of the matrix.
niidat=zeros(size(Nspace.vol));
for i=1:size(dat,3);
    niidat(:,:,i)=dat(:,:,i)';
end

% replace the roi image data into the reference image.
Nspace.vol=niidat;

%write out the new roi image.
MRIwrite(Nspace,Nout{x});
end