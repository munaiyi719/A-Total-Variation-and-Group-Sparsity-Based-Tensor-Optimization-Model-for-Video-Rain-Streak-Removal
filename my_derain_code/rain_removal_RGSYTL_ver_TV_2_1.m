function [B,iter]=rain_removal_RGSYTL_ver_TV_2_1(O,opts)
%%%%log det




%%%--- Preparation ---%%%

Dx = def3Dx;	DxT = def3DxT;    % up-down direction
Dy = def3Dy;	DyT = def3DyT;
Dz = def3Dz;	DzT = def3DzT;



% filter.x(1,:,:) = gpuArray(1);      filter.x(2,:,:) = gpuArray(-1);
% filter.y(:,1,:) = gpuArray(1);      filter.y(:,2,:) = gpuArray(-1);
% filter.z(:,:,1) = gpuArray(1);      filter.z(:,:,2) = gpuArray(-1);
filter.x(1,:,:) = 1;      filter.x(2,:,:) = -1;
filter.y(:,1,:) = 1;      filter.y(:,2,:) = -1;
filter.z(:,:,1) = 1;      filter.z(:,:,2) = -1;
Size = size(O);
Dim = length(Size);
eigsDxTDx = abs(psf2otf(filter.x,Size)).^2;
eigsDyTDy = abs(psf2otf(filter.y,Size)).^2;
eigsDzTDz = abs(psf2otf(filter.z,Size)).^2;

%%%--- Parameters ---%%%
maxit=opts.maxit;
tol=opts.tol;
alpha1=opts.alpha1;     beta1=opts.beta;
if alpha1==0
    beta1=0;
end
alpha2=opts.alpha2;     beta2=opts.beta;
if alpha2==0
    beta2=0;
end
alpha3=opts.alpha3;     beta3=opts.beta;
if alpha3==0
    beta3=0;
end
alpha4=opts.alpha4;     beta4=opts.beta;
if alpha4==0
    beta4=0;
end





%weight=opts.weight;
%%%--- Initialize ---%%%
R=zeros(Size);
Lambda1 = (zeros(Size));      Lambda2 = (zeros(Size));
Lambda3 = (zeros(Size));      Lambda4 = (zeros(Size));


%%%--- Main loop---%%%
Denom=beta1*eigsDxTDx+beta2*eigsDyTDy+beta3*eigsDzTDz+beta4*ones(Size);

iter=0;relcha=1;relcha_all=zeros(maxit);
PSNR_all=[];
B=O;
%fprintf('\nIteration:     ');
s=1;
epsilon=1e-4;
dontknow=0;
maxit=maxit+dontknow;
G=zeros(Size);
Y=zeros(Size);
T=zeros(Size);
S=zeros(Size);
while relcha>tol &&  iter<maxit
 
    %--- G-subproblem---%
    if beta1~=0
        G = wthresh(Dx(R)+Lambda1/beta1,'s',alpha1/beta1);
    end
    
    
    %--- Y-subproblem---%
    if beta2~=0
        Y = wthresh(Dy(O-R)+Lambda2/beta2,'s',alpha2/beta2);
    end
    
    
    %--- T-subproblem---%
    if beta3~=0
        T = wthresh(Dz(O-R)+Lambda3/beta3,'s',alpha3/beta3);
    end
    
    
    %--- S-subproblem---%
    if beta4~=0
        tS=Unfold(R,Size,1);
        tLam=Unfold(Lambda4,Size,1);
        r=tS+tLam/beta4;
        tnorm=sqrt(sum(r.^2));
        tS=r.*repmat(max( tnorm -alpha4/beta4,0)./(tnorm+eps),Size(1),1);
        S=Fold(tS,Size,1);
    end
    
    
    %--- R-subproblem---%
    R_k=R;
    B_k=B;
    R=real(ifftn(fftn(DxT(beta1*G-Lambda1)+DyT(beta2*Dy(O)-beta2*Y+Lambda2)+DzT(beta3*Dz(O)-beta3*T+Lambda3)+beta4*S-Lambda4)./Denom));

    R(R<0) = 0;
    R(R>1) = 1;
    %B=Fold(b,Size,1);
    B=O-R;
    B(B<0) = 0;
    B(B>1) = 1;
    R=O-B;
     relcha=norm(R_k(:)-R(:),'fro')/norm(R_k(:),'fro');
    
    %--- Multipliers updating---%
    Lambda1 = Lambda1 + 1.618*beta1*(Dx(R)-G);
    Lambda2 = Lambda2 + 1.618*beta2*(Dy(O-R)-Y);
    Lambda3 = Lambda3 + 1.618*beta3*(Dz(O-R)-T);
    Lambda4 = Lambda4 + 1.618*beta4*(R-S);
    
    
    
    iter=iter+1;

    %--- Relative change reporting---%
    fprintf('Iteration: %3i   Relative chang: %5.5f\n',iter,relcha);

end
