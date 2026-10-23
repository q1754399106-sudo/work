function varargout=p1_synthesize(action,c,baseModels,suite,seed,scenario)
% SYNTHESIS ONLY. This file may read Q/fpeak/amplitudes and clean references.
% Its outputs handed to inversion are D, offsets, and uq_public_config ONLY.
% CoherentNoiseAVO is a scene/geometry name; coherent contamination is enabled
% exclusively by StressCoherent, keeping Main/LowerSNR white-noise suites.
public=uq_public_config(c,c.nominal_frequency);
if strcmp(suite.perturbation,'velocity')
 public.vp=public.vp.*(1+suite.magnitude); % deliberately misspecified kinematics
end
trueModels=baseModels;
if strcmp(suite.perturbation,'q_error')
 key=sprintf('%s_%g',c.model_name,suite.magnitude);
 trueModels=q_stress_models(key,c,suite.magnitude);
end
if strcmp(action,'calibration')
 scene=build_scene('SmoothAVO',c,suite);
 [D,sigma]=observation(scene,c,trueModels,suite,seed,true);
 varargout={public,D,c.offset,struct('noise_sigma',sigma,'synthesis_only',true)};
elseif strcmp(action,'observation')
 scene=build_scene(scenario,c,suite);
 [D,sigma]=observation(scene,c,trueModels,suite,seed,false);
 varargout={D,scene,sigma};
else, error('Unknown synthesis action'); end
end

function models=q_stress_models(key,c,magnitude)
% Reused only within this synthesis file; never accessible to an inverse.
persistent keys values
if isempty(keys), keys={}; values={}; end
ix=find(strcmp(keys,key),1);
if isempty(ix)
 truthConfig=c; truthConfig.Q=c.Q.*[1-magnitude,1+magnitude,1-magnitude,1+magnitude];
 models=p1_truth_models(truthConfig);
 % Keep only one model in this function to bound memory across geology.
 keys={key}; values={models};
else, models=values{ix}; end
end

function scene=build_scene(name,c,suite)
scene=make_stability_scene(name,c);
if strcmp(suite.perturbation,'source_phase') || strcmp(suite.perturbation,'source_shape') || strcmp(suite.perturbation,'residual_nmo')
 for j=1:size(scene.truth,2)
  times=scene.events{j};
  if strcmp(suite.perturbation,'residual_nmo')
   % A deterministic residual quadratic moveout, with no stretch. Its sign
   % alternates by event and reaches magnitude seconds at the far offset.
   moveout=suite.magnitude*(scene.offsets(j)/max(c.offset))^2;
   times=times+moveout*[1;-1;1]; scene.events{j}=times;
  end
  W=zeros(c.nt,3);
  for k=1:3
   x=c.t-times(k); z=pi*c.fpeak*x; w=(1-2*z.^2).*exp(-z.^2);
   if strcmp(suite.perturbation,'source_phase')
    w=cosd(suite.magnitude)*w+sind(suite.magnitude)*quadrature(w,c.nfft);
   elseif strcmp(suite.perturbation,'source_shape')
    % Fixed doublet, outside the single-frequency phase-rotated Ricker family.
    delay=suite.magnitude; z2=pi*(.68*c.fpeak)*(x-delay);
    w=w+.32*(1-2*z2.^2).*exp(-z2.^2);
   end
   W(:,k)=w;
  end
  scene.templates{j}=W; scene.truth(:,j)=W*scene.coeff{j};
 end
end
end

function [D,sigma]=observation(scene,c,models,suite,seed,isCalibration)
n=size(scene.truth,1); nx=size(scene.truth,2); clean=zeros(n,nx);
for j=1:nx
 ix=find(c.offset==scene.offsets(j),1); assert(~isempty(ix));
 clean(:,j)=models{ix}.G*scene.truth(:,j);
end
rng(seed,'twister'); noise=randn(n,nx);
if strcmp(suite.perturbation,'colored_noise')
 % AR(1) temporal correlation; normalized per trace so nominal RMS remains
 % comparable. Filter realization/burn-in uses only this independent seed.
 draw=randn(n+256,nx); draw=filter(1,[1 -suite.magnitude],draw,[],1);
 noise=draw(257:end,:);
end
noise=noise./max(sqrt(mean(noise.^2,1)),realmin);
sigma=sqrt(mean(clean.^2,1))*10^(-suite.snr_db/20);
D=clean+noise.*sigma;
if strcmp(suite.perturbation,'coherent_noise') && ~isCalibration
 coherent=zeros(n,nx);
 for j=1:nx
  z=pi*35*(c.t-(mean(c.event_t(1:2))+.004*(j-1)));
  coherent(:,j)=(1-2*z.^2).*exp(-z.^2);
 end
 D=D+suite.magnitude*norm(clean,'fro')*coherent/max(norm(coherent,'fro'),realmin);
elseif strcmp(suite.perturbation,'impulse_noise') && ~isCalibration
 % Contamination is additional to nominal white SNR and reported separately.
 count=max(1,round(suite.magnitude*numel(D))); indices=randperm(numel(D),count); indices=indices(:);
 scales=repmat(sigma,n,1); signs=2*(rand(count,1)>.5)-1;
 D(indices)=D(indices)+10*scales(indices).*signs;
end
end

function q=quadrature(w,N)
phase=[0;-1i*ones(N/2-1,1);0;1i*ones(N/2-1,1)];
u=real(ifft(fft(w,N).*phase)); q=u(1:numel(w));
end
