function [CRSPD,Brelabeled] = compare_modularity_partitions(A,B,nii,flag,outnii)
%                  COMPARE_MODULARITY_PARTITIONS   
%    Compare and relabel a modularity partition B, with respect
%                   to a second partition A
%
%   [CRSPD,Brelabeled] = compare_modularity_partitions(A,B,nii,flag,outnii)
%
%INPUTS:
%
%A:     Reference partition (REF) deviations from which are identified.
%B:     Partition of interest (POI) that is compared to reference.
%nii:   Path to nifti parcellation volume of nodes. (REQUIRES NIfTI TOOLBOX)
%flag:  1 or 0; write out a NIfTI volume of REF partition. (REQUIRES NIfTI TOOLBOX)
%outnii:name of the relabeled outout partition nifti
%
%OUTPUTS:
%
%CRSPD:     Is a 2x number of communities in POI matrix, which matches the
%           POI communities (column2) to REF communities (column1). 
%Brelabeld: A vector with POI communities relabeld to match REF. If
%           multiple POI communities matched to a single REF community,
%           they were relabeld in decimal (e.g if 2,3 in POI belonged to 2
%           in REF, they would be relabeld as 2 and 2.5).
%
%Example Usage:
%[CRSPD,Brelabeled]=compare_modularity_partitions(REF,POI,'/usr/path/parcellation.nii.gz',1);
%
%   Evgeny Chumin, Indiana University School of Medicine, 2018

if nargin<2 % check that two partitions are input
    error('Missing second partition input!')
end
if nargin<3 % check that string path to parcellation nii exists
    warning('No path to nifti parcellation! No images will be written.')
end
if length(A)~= length(B) % length A and B must be equal
    error('Reference and interest partitions must be the same length!')
end
if ~exist('flag','var')
    flag=0;
end 

%% Initialization
REFcomm = cell(max(A),1);
POIcomm = cell(max(B),1);
simi = nan(max(A),max(B));
CRSPD = nan(max(B),2);
surp = cell(max(B),1);
defi = cell(max(B),1);

%% Load Parcellation
if nargin>=3
Parc=MRIread(nii);
end

%% Match Partitions
% Parse REF into a cell by community 
for i=1:max(A)
    [REFcomm{i,1},~]=find(A==i);
end
% Parse POI into a cell by community
for i=1:max(B)
    [POIcomm{i,1},~]=find(B==i);
end

for i=1:max(B)
    for j=1:max(A)
        % proportion POI comm regions found in REF comms
        simi(j,i)=(sum(ismember(POIcomm{i,1},REFcomm{j,1})))/max(size(POIcomm{i,1}));
    end
        % find highest proportion partition
        [~,match]=max(simi(:,i));
        CRSPD(i,2)= i; %POI community number
        CRSPD(i,1)= match; % matched REF community number
end

%% Identify regions of difference
for i=1:size(CRSPD,1)
    % Found in POI but not REF
    surp{i,1}=setdiff(POIcomm{CRSPD(i,2),1},REFcomm{CRSPD(i,1),1});
    % Absent in POI but not in REF
    defi{i,1}=setdiff(REFcomm{CRSPD(i,1),1},POIcomm{CRSPD(i,2),1});
end

%% Relabel POI partition to match REF
%initialize
Brelabeled=B;
isrelabeld=[];
k=0; % initiate counter
% for every community in POI
for i=1:max(Brelabeled)
    % if it has not been relabeld
    if isempty(find(isrelabeld==i)) %#ok<*EFIND>
        % find location of REF community corresponding to POI community i
        newlabelpos=find(CRSPD(:,1)==CRSPD(find(CRSPD(:,2)==i),1)); %#ok<*FNDSB>
        % if several POI communities match to one REF community
        if length(newlabelpos)>1
            % for every matched POI community j 
            for j=1:length(newlabelpos)
                % find the number of regions in each community that matched
                numRegInRef(j)=simi(CRSPD(i,1),newlabelpos(j))*max(size(POIcomm{newlabelpos(j),1})); %#ok<*AGROW>
            end
            % community with the most matched regions gets the label
            [~,I]=max(numRegInRef);
            % relabel POI community with matched REF decreasing label by a
            % magnitude of 100 (avoids overlap after relabeling)
            Brelabeled(Brelabeled==CRSPD(newlabelpos(I),2))=(CRSPD(newlabelpos(I),1)/100);
            % add relabeld POI community to completed list
            isrelabeld(end+1,1)=CRSPD(newlabelpos(I),2);
            % for every matched community
            for j=1:length(numRegInRef) 
                if j~=I % excluding the one that won the label above
                    k=k+1;  % increase counter
                    % set new label to be 1+ maximum label in REF
                    Brelabeled(Brelabeled==CRSPD(newlabelpos(j),2))=(max(CRSPD(:,1))+k)/100;
                    % add relabeld POI community to completed list
                    isrelabeld(end+1,1)=CRSPD(newlabelpos(j),2);
                end
            end
                    clear numRegInRef % clear region counter
        % if only one POI community matched REF
        else
            % relabel POI community with matched REF decreasing label by a
            % magnitude of 100 (avoids overlap after relabeling)
            Brelabeled(Brelabeled==CRSPD(newlabelpos,2))=CRSPD(newlabelpos,1)/100;
            % add relabeld POI community to completed list
            isrelabeld(end+1,1)=CRSPD(newlabelpos,2);
        end
    end
end
% increase by magnitude of 100 to return to original scale
Brelabeled=Brelabeled*100; 
Brelabeled=round(Brelabeled); % make sure labels are all integers

%% Write NIfTI volume of relabeled POI
if nargin>=3
BrelabeledNii=Parc;
% For every region in volume, relabel according to partition
for i=1:max(max(max(BrelabeledNii.vol)))
    BrelabeledNii.vol(BrelabeledNii.vol==i)=Brelabeled(i,1);
end
% Write out the .nii.gz partition image
MRIwrite(BrelabeledNii,outnii)
end

%% Write NIfTI volume of REF
if flag==1
    REFnii=Parc;
    % For every region in volume, relabel according to partition
    for i=1:max(max(max(REFnii.vol)))
        REFnii.vol(REFnii.vol==i)=A(i,1);
    end
    % Write out the .nii.gz partition image
    MRIwrite(REFnii,'REF_partition.nii.gz')
end

end