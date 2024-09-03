% This code takes the input of the RD and RS files from one treatment plan
% for a patient and applys the dose-distance metric on the outermost voxels
% of the body contour. 

% Load in the structure set
RSpath = %your path to RS file here
[structureSet] = dicominfo(RSpath);
% Extract contours from structure set
[Structures,structureNames] = ExtractContours(structureSet);
% Extract contour points using index of the 'Body' Structure 
bodyIndex = %index of Body Structure 
contourPts = Structures(bodyIndex).contour;
% Process your DICOM RD file
RDpath = %your path to RD file here
%Get dicom info on the treatment plan
PlanInfo = dicominfo(RDpath);
% Get position for center of each dose voxel
planPtsZ = PlanInfo.ImagePositionPatient(3) + PlanInfo.GridFrameOffsetVector(:);
planPtsX = (1:1:double(PlanInfo.Width))*PlanInfo.PixelSpacing(1)+PlanInfo.ImagePositionPatient(1)-1.5;
planPtsY = (1:1:double(PlanInfo.Height))*PlanInfo.PixelSpacing(2)+PlanInfo.ImagePositionPatient(2)-1.5;
% Get voxel data for RD file
RawPlan = dicomread(RDpath);
DoseScaling = PlanInfo.DoseGridScaling;
plan = double(squeeze(RawPlan)*DoseScaling);    % Plan array accounts for scaling
% Get dose to the body points
bodyDoses = interp3(planPtsX,planPtsY,planPtsZ,plan,contourPts(1,:),contourPts(2,:),contourPts(3,:));
max(bodyDoses,[],'omitnan')
% DoseMatrix contains the coodinates x,y,z, the corresponding dose, and the
% EQD2 dose
dose_matrix(:,1) = contourPts(1,:);
dose_matrix(:,2) =contourPts(2,:);
dose_matrix(:,3) =contourPts(3,:);
dose_matrix(:,4) = bodyDoses;
f= %number of fractions for treatment here
ab=11; %alpha beta ratio for moist desquamation
dose_matrix(:,5) = dose_matrix(:,4).*(dose_matrix(:,4)/f+ab)/(2+ab);
% Find Euclidean Distance between coordinates
dist = squareform(pdist(dose_matrix(:,1:3)));

% Define Dose-Distance Parameters:
a = % a parameter
d0 = %d0 parameter
ar = % surface area of voxel

% Here is the Dose-Distance Equation.
dist_thresh = % your distance threshold here
for i = 1:size(dose_matrix(:,1),1)
    for j = 1:size(dose_matrix(:,1),1)
        if dist_mat(i,j) <= dist_thresh
                    temp(j) = dose_matrix(i,5)*dose_matrix(j,5)*ar/(exp(a*(dist_mat(i,j)/d0 -1))+1);
        else
            temp(j) = NaN;
        end
    end
    dose_matrix(i,6) = mean(mean(temp(:),'omitnan'),'omitnan');
    temp = [];
end

% Dose_matrix(:,6) now contains the dose-distance value corresponding to
% each coordinate.

