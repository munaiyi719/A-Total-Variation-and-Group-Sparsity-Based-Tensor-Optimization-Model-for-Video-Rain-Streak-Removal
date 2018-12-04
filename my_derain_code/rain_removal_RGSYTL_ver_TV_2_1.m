function [B,iter]=rain_removal_RGSYTL_ver_TV_2_1(O,opts)
%%%%log det


%%% date: 31/08/2016
%%% This demo is for the rain steak model, i.e., the $\mathcal{R}$ model (3 Directions, with Low-Rank)
%%%
%%% $\min\limits_{\mathcal{G,Y,T,S,L}}\quad\alpha1||\mathcal{G}||_1+\alpha2||\mathcal{Y}||_1+\alpha3||\mathcal{T}||_1+\alpha4||\mathcal{S}||_1+\alpha5||\mathcal{L}||_*$
%%%         s.t.
%%%              $\mathcal{G}=D_x(\mathcal{R})$
%%%              $\mathcal{Y}=D_y(\mathcal{O}-\mathcal{R})$
%%%              $\mathcal{T}=D_t(\mathcal{O}-\mathcal{R})$
%%%              $\mathcal{S}=\mathcal{R}$
%%%              $\mathcal{L}=\mathcal{O}-\mathcal{R}$
%%%   input:      the original rainy video $\mathcal{O}$
%%%   output:   the rainy steak $\mathcal{R}$
%%%                  the rain-free video $\mathcal{B}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---date: 2016/08/26

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
% alpha5_1=opts.alpha5;     beta5_1=opts.beta;
% alpha5_2=opts.alpha5;     beta5_2=opts.beta;
% alpha5_3=opts.alpha5;     beta5_3=opts.beta;
%%%%%%  new
%alpha6=opts.alpha6;     beta6=opts.beta;
%alpha7=opts.alpha7;     beta7=opts.beta;





%weight=opts.weight;
%%%--- Initialize ---%%%
R=zeros(Size);
Lambda1 = (zeros(Size));      Lambda2 = (zeros(Size));
Lambda3 = (zeros(Size));      Lambda4 = (zeros(Size));
% Lambda5_1 = (zeros(Size));
% Lambda5_2 = (zeros(Size));
% Lambda5_3 = (zeros(Size));
%Lambda6 = (zeros(Size));
%Lambda7 = (zeros(Size));


%%%--- Main loop---%%%
Denom=beta1*eigsDxTDx+beta2*eigsDyTDy+beta3*eigsDzTDz+beta4*ones(Size);
%     +beta5_2*ones(Size)+beta5_3*ones(Size);

%Denom=Denom+beta6*ones(Size)+beta7*ones(Size);
iter=0;relcha=1;relcha_all=zeros(maxit);
PSNR_all=[];
B=O;
%fprintf('\nIteration:     ');
s=1;
epsilon=1e-4;
% lastV_1=ones(min(size(Unfold(O,Size,1))),1);
% lastV_2=ones(min(size(Unfold(O,Size,2))),1);
% lastV_3=ones(min(size(Unfold(O,Size,3))),1);
dontknow=0;
maxit=maxit+dontknow;
G=zeros(Size);
Y=zeros(Size);
T=zeros(Size);
S=zeros(Size);
while relcha>tol &&  iter<maxit
    %fprintf('\b\b\b\b\b%5i',iter)
    %--- L-subproblem---%
%     L=zeros(Size);
%      [tL_1,lastV_1]=Pro2TraceNorm_logdet(Unfold(O-R+Lambda5_1/beta5_1,Size,1),alpha5_1/beta5_1,lastV_1+epsilon);
%      L_1=Fold(tL_1,Size,1);
%      [tL_2,lastV_2]=Pro2TraceNorm_logdet(Unfold(O-R+Lambda5_2/beta5_2,Size,2),alpha5_2/beta5_2,lastV_2+epsilon);
%      L_2=Fold(tL_2,Size,2);
%      [tL_3,lastV_3]=Pro2TraceNorm_logdet(Unfold(O-R+Lambda5_3/beta5_3,Size,3),alpha5_3/beta5_3,lastV_3+epsilon);
%      L_3=Fold(tL_3,Size,3);
    
%    a = Unfold(O-R+Lambda5_1/beta5_1,Size,1);
%     [a,~] = Pro2TraceNorm(a,alpha5_1/beta5_1);
%     L_1=Fold(a,Size,1);
     
     
%     a = Unfold(O-R+Lambda5_2/beta5_2,Size,2);
%     [a,~] = Pro2TraceNorm(a,alpha5_2/beta5_2);
%     L_2=Fold(a,Size,2);
     
     
%     a = Unfold(O-R+Lambda5_3/beta5_3,Size,3);
%     [a,~] = Pro2TraceNorm(a,alpha5_3/beta5_3);
%     L_3=Fold(a,Size,3);
    
    
%     for i = 1:Dim
%         a = Unfold(O-R+Lambda5/beta5,Size,i);
%         if s==1
%             [a,~] = Pro2TraceNorm(a,alpha5/beta5);
%         else
%             a = Normker(a,alpha5/beta5,delta);
%         end
%         L = L+(1/Dim)*weight(i)*Fold(a,Size,i);
%     end
    
    
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
%     if s==1
%     S = wthresh(R +Lambda4/beta4,'s',alpha4/beta4);
%     else
%     S = Norm0(R +Lambda4/beta4,alpha4/beta4,delta);
%     end
%       2_1
    if beta4~=0
        tS=Unfold(R,Size,1);
        tLam=Unfold(Lambda4,Size,1);
        r=tS+tLam/beta4;
        tnorm=sqrt(sum(r.^2));
        tS=r.*repmat(max( tnorm -alpha4/beta4,0)./(tnorm+eps),Size(1),1);

    %     for i=1:Size(2)*Size(3)
    %          r=tG(:,i)+tLam(:,i)/beta1;
    %         tG(:,i)=r.*max(norm(r)-alpha1/beta1,0)/(norm(r)+eps);
    %     end
        S=Fold(tS,Size,1);
    end
%    new item    
%    Z=Normk(Unfold(R+Lambda6/beta6,size,2),alpha6/beta6,delta);
%    K=Normk(Unfold(R+Lambda7/beta7,size,3),alpha7/beta7,delta);
    
    
    %--- R-subproblem---%
    R_k=R;
    B_k=B;
%     R=real(ifftn(fftn(DxT(beta1*G-Lambda1)+DyT(beta2*Dy(O)-beta2*Y+Lambda2)+DzT(beta3*Dz(O)-beta3*T+Lambda3)+beta4*S-Lambda4+beta5_1*(O-L_1)+Lambda5_1...
%     +beta5_2*(O-L_2)+Lambda5_2+beta5_3*(O-L_3)+Lambda5_3)./Denom));
    R=real(ifftn(fftn(DxT(beta1*G-Lambda1)+DyT(beta2*Dy(O)-beta2*Y+Lambda2)+DzT(beta3*Dz(O)-beta3*T+Lambda3)+beta4*S-Lambda4)./Denom));

    %R=real(ifftn(fftn(DxT(beta1*G-Lambda1)+DyT(beta2*Dy(O)-beta2*Y+Lambda2)+DzT(beta3*Dz(O)-beta3*T+Lambda3)+beta4*S-Lambda4+beta5*(O-L)+Lambda5...
    %+beta6*Z-Lambda6+beta7*K-Lambda7)./Denom));
    R(R<0) = 0;
    R(R>1) = 1;
    % b(b<0) = 0;
    % b(b>1) = 1;
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
%     Lambda5_1 = Lambda5_1 + 1.618*beta5_1*(O-R-L_1);
%     Lambda5_2 = Lambda5_2 + 1.618*beta5_2*(O-R-L_2);
%     Lambda5_3 = Lambda5_3 + 1.618*beta5_3*(O-R-L_3);
    %Lambda6 = Lambda6 + 1.618*beta4*(R-Z);
    %Lambda7 = Lambda7 + 1.618*beta4*(R-K);
    
    
    
    iter=iter+1;
%     if iter==125+dontknow
%         B125=B;
%     end
%     if iter==50+dontknow
%         B50=B;
%     end
%     if iter==75+dontknow
%         B75=B;
%     end
%     if iter==100+dontknow
%         B100=B;
%     end
    %--- Relative change reporting---%
    fprintf('Iteration: %3i   Relative chang: %5.5f\n',iter,relcha);
    %relcha_all(iter)=relcha;
%     [psnr ~]=getdata(B,padsize,O_Rainy,O_clean,Rain);
%     [relcha iter psnr]
    %PSNR_all(iter)=PSNR3D(gather(B),gather(O_clean));
end