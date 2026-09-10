function acpcXform = ssf_rsfMRIAlignToT1(acpcT1, T1preproc, t1MaskFile, useStdXformFlag, figNum, unwarpDti, sepParam)
% ssf_rsfMRIAlignToT1
%
% Align the T1-weighted anatomical image used for electrode localization
% to the preprocessed T1-weighted anatomical image associated with
% resting-state fMRI.
%
% INPUTS
%   acpcT1          - T1 image used for electrode localization
%   T1preproc       - Preprocessed T1 image associated with rs-fMRI
%   t1MaskFile      - Optional target T1 brain mask
%   useStdXformFlag - Use canonical transform for initial alignment
%   figNum          - Figure number; >0 enables SPM coregistration display
%   unwarpDti       - Legacy option controlling 6- vs 12-parameter fit
%   sepParam        - SPM coregistration sampling separation
%
% OUTPUT
%   acpcXform       - Transformation from electrode-localization T1 space
%                     to rs-fMRI T1 space
%
% Adapted from the Vistasoft function dtiRawAlignToT1.
% Original implementation by Robert F. Dougherty (RFD), Vista Lab.
%
% Adaptation for resting-state fMRI:
% Maria Guadalupe Yanez Ramos
% Developed with scientific and technical guidance from Dora Hermes
% and the Multimodal Neuroimaging Lab (MNL) team.
% 2025-2026


%% Set defaults

if(~exist('t1MaskFile','var')), t1MaskFile = []; end
if(~exist('unwarpDti','var') || isempty(unwarpDti))
  unwarpDti = false;
end
if(~exist('sepParam','var') || isempty(sepParam))
  sepParam = [16 8 4 2];
end


% Initialize SPM default params for the coregistration.
estParams        = spm_get_defaults('coreg.estimate');
estParams.params = [0 0 0 0 0 0];% Rigid-body (6-params)
estParams.sep    = sepParam; 

if (unwarpDti)
  estParams.params = [0 0 0 0 0 0 1 1 1 0 0 0]; % 12-param affine
else
  estParams.params = [0 0 0 0 0 0]; % 6-param Rigid body
end


%% Align electrode-localization T1 to rs-fMRI T1

fprintf('[%s] Aligning electrode-localization T1 to rs-fMRI T1...\n', ...
    mfilename);

source.uint8 = uint8(round(mrAnatHistogramClip(double(acpcT1.data),0.4,0.99)*255));
% % Blur images given the specified highest-resolution sampling density
fwhm = sqrt(max([1 1 1]*estParams.sep(end)^2 - acpcT1.pixdim.^2, [0 0 0]))./acpcT1.pixdim;
source.uint8 = mrAnatSmoothUint8(source.uint8,fwhm);
% We also need a reasonable starting guess at the mnB0 ac-pc xform.
if (acpcT1.qform_code>0)
  source.mat = acpcT1.qto_xyz;
elseif (acpcT1.sform_code>0)
  source.mat = acpcT1.sto_xyz;
else
  error('Source T1 requires qform_code>0 or sform_code>0.');
end

% if(ischar(T1preproc))
%     T1preproc = niftiRead(T1preproc);
% end

target.uint8 = uint8(round(mrAnatHistogramClip(double(T1preproc.data),0.50,0.995)*255));
fwhm = sqrt(max([1 1 1]*estParams.sep(end)^2 - T1preproc.pixdim(1:3).^2, [0 0 0]))./T1preproc.pixdim(1:3);
target.uint8 = mrAnatSmoothUint8(target.uint8,fwhm);
if(T1preproc.qform_code>0)
  target.mat = T1preproc.qto_xyz;
elseif(T1preproc.sform_code>0)
  target.mat = T1preproc.sto_xyz;
else
  error('Requires that t1 qform_code>0 OR sform_code>0.');
end
if(~isempty(t1MaskFile))
  brainMask = niftiRead(t1MaskFile);
  target.uint8(~brainMask.data) = 0;
  clear brainMask;
end
% Sanity check for rough alignment xforms from image header. If the center
% coordinate falls outside the image volume, then we replace the
% translations with a default to the image center. This WILL NOT fix bad
% rotations or bad scales.
centerCoord = inv(target.mat)*[0 0 0 1]'; centerCoord = centerCoord(1:3)';
if (any(centerCoord<1) || any(centerCoord>size(target.uint8)))
    [t,r,s] = affineDecompose(target.mat);
    t = size(target.uint8)./2.*s;
    im2std = mrAnatComputeCannonicalXformFromDicomXform(target.mat,size(target.uint8));
    im2std(:,4) = [0 0 0 1]';
    target.mat = im2std*affineBuild(-t,[0 0 0],s);
    %t = t.*-sign(target.mat(1:3,1:3)*[size(target.uint8)./2]')';
    %target.mat = affineBuild(t,[0 0 0],s);
end
centerCoord = inv(source.mat)*[0 0 0 1]'; centerCoord = centerCoord(1:3)';
if (useStdXformFlag || any(centerCoord<1) || any(centerCoord > size(source.uint8)))
    [t,r,s] = affineDecompose(source.mat);
    t = size(source.uint8)./2.*s;
    im2std = mrAnatComputeCannonicalXformFromDicomXform(source.mat,size(source.uint8));
    im2std(:,4) = [0 0 0 1]';
    source.mat = im2std*affineBuild(-t,[0 0 0],s);
    %t = t.*-sign(source.mat(1:3,1:3)*[size(source.uint8)./2]')';
    %source.mat = affineBuild(t,[0 0 0],s);
end
%headerMI = mrAnatComputeMutualInfo([0 0 10 0 0 0],source,target,[1 1 1],estParams.cost_fun,estParams.fwhm);
% sc = fileparts(acpcT1.fname); if(length(sc)>25) sc = ['...' sc(end-18:end)]; end
% if (figNum>0)
%   figure(figNum);
%   dtiShowAlignFigure(figNum, target, source, [], [], ['initial align ( ' sc ' )']);
% end


if(figNum>0)
    transRot = spm_coreg(source,target,estParams);
else
    % suppress the verbose optimizer output.
    msg = evalc('transRot = spm_coreg(source,target,estParams);');
end
acpcXform = spm_matrix(transRot(end,:))*source.mat;
source.mat = acpcXform;
% if(figNum>0)
%   figure(figNum);
%   %dtiShowAlignFigure(figNum, target, source, [], [], ['final align ( ' sc ' )']);
% end
%disp(['Saving b0 acpcXform to ' outAcpcXform '...']);
%save(outAcpcXform, 'acpcXform');
%dlmwrite([outAcpcXform '.txt'],acpcXform,'delimiter',' ','precision',6);

return;
