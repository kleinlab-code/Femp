clear, format compact,close all
misspec=1;
extNoise=0; %  =1?
for N=5:100 % number of independent electrodes and time points
    F=randn(3,N); % the true forward model for two sources
    T=randn(3,N); % the two true temporal response
    V=F'*T + extNoise*randn(N,N); % electrode voltages plus noise
    Ffmri=F+misspec*randn(3,N); % 3-shell estimate of F (note the noisy misspecification)
    T1=Ffmri'\V; % temporal response from algorithm of Section 2
    F1=V/T1; % improved estimate of forward model
    T2=F1\V; % improved estimate of time function
    F2=V/T2; % Do another iteration of forward model
    Cor=corrcoef([T(1,:)' T1(1,:)' T2(1,:)']); % correlations of time functions
    CorallT(N,:)=[Cor(1,2) Cor(1,3)]; % correls of true T to two estimates
    FF=[F(1,:)' Ffmri(1,:)' F1(:,1) F2(:,1) ];
    Cor=corrcoef(FF); % correlation of true F to the three estimates
    Corall(N,:)=[Cor(1,2) Cor(1,3) Cor(1,4)];   
    
end

N=1:100;
subplot(1,2,1); plot(N,CorallT(N,1),'b.',N,CorallT(N,2),'r','linewidth',1.5)
axis([1 100 0 1]); hold on
xlabel('number of electrode/patches and time points')
ylabel('correlation')
title('Correlations of estimated time function to true')
subplot(1,2,2); plot(N,Corall(N,:));hold on
xlabel('number of electrode/patches and time points')
ylabel('correlation');axis([1 100 0 1]);
title('Correlations of estimated forward model to true')
%% Use the empirical head model for a single patch (smallish N)
for N=1:30,   %NOTE THAT THIS IS BASED ON THE SAME NOISE!
    T1Patch=F2(1:N,:)\V(1:N,:);  %this is for a small subset of patch/elec
    Cor=corrcoef([T(1,:)' T1Patch(1,:)']);
    CorPatchT(N)=Cor(1,2);
end
N=1:30;
subplot(1,2,1);plot(N,CorPatchT(N),'k','linewidth',3);

%% Use the empirical head model for a single patch (smallish N)
Nelec = 1; % Let's just say there are 25 electrodes
Npatch = 4; % And there are 4 patches we think have same time function

% Get Tcorr for empirical head model
for iPatch = 1:Npatch
    ep = (1:Nelec)+(Nelec*(iPatch-1));
    T1Patch=F2(ep,:)\V(ep,:);  % Goes through all the patches
    Cor=corrcoef([T(1,:)' T1Patch(1,:)']);
    CorPatchTemp(iPatch,:)=Cor(1,2);
end

% Get Tcorr for BEM head model
for iPatch = 1:Npatch
    ep = (1:Nelec)+(Nelec*(iPatch-1));
    T1Patch=Ffmri(:,ep)'\V(ep,:);
    Cor=corrcoef([T(1,:)' T1Patch(1,:)']);
    CorPatchTbem(iPatch,:)=Cor(1,2);
end
figure(2)
plot(1:Npatch,CorPatchTemp,'k*-','linewidth',3);hold on;
plot(1:Npatch,CorPatchTbem,'r*-','linewidth',3);
xlabel('Patch #'); ylabel('Corr');
legend('T_{corr} - Empirical', 'T_{corr} - BEM', 'Location', 'SouthEast')
axis([1 Npatch .6 1])

