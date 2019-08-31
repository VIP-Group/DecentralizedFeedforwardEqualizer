% ----------------------------------------------------------
% -- Simple Massive MU-MIMO Simulator with Decentralized
% -- Feedforward Architectures (v0.1)
% -- 2019 (c) studer@cornell.edu
% ----------------------------------------------------------
% -- Code to reproduce parts of Fig. 3 in the paper :
% -- C. Jeon, K. Li, J. Cavallaro, C. Studer, "Decentralized
% -- Equalization with Feedforward Architectures for Massive
% -- MU-MIMO. IEEE Transactions on Signal Processing,
% -- pp. 4418 -4432, Vol. 67, No. 17, July 2019
% ----------------------------------------------------------

function MIMOsim_DBP(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.simName = 'DBP_128x8_16QAM'; % simulation name (used for saving results)
    par.runId = 0; % simulation ID (used to reproduce results)
    par.save = true; % save results?
    par.B = 256; % receive antennas
    par.U = 16; % transmit antennas (set not larger than MR!)
    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 10000; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = 0:2:10; % list of SNR [dB] values to be simulated
    par.detector = {'MRC','ZF','LMMSE','PD_LMMSE','FD_LMMSE'}; % define detector(s) to be simulated
    par.C = 8; % number of antenna clusters
    
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);

% check whether antenna cluster size is valid
if round(par.B/par.C)~=par.B/par.C
    error('Invalid cluster size par.C')
end

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK',
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% initialize result arrays (detector x SNR)
res.VER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.Q,par.trials);

%initialize parameters for ADMM

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.B,1)+1i*randn(par.B,1));
    H = sqrt(0.5)*(randn(par.B,par.U)+1i*randn(par.B,par.U));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
        N0 = par.U*par.Es*10^(-par.SNRdB_list(k)/10);
        
        % transmit data over noisy channel
        y = x+sqrt(N0)*n;
        
        % algorithm loop
        for d=1:length(par.detector)
            
            switch (par.detector{d}) % select algorithms
                case 'MRC', % maximum ratio combining detection
                    [idxhat,bithat] = MRC(par,H,y);
                case 'ZF', % zero-forcing detection
                    [idxhat,bithat] = ZF(par,H,y);
                case 'LMMSE', % linear MMSE detector
                    [idxhat,bithat] = LMMSE(par,H,y,N0);
                case 'PD_LMMSE'
                    [idxhat,bithat] = PD_LMMSE(par,H,y,N0);
                case 'FD_LMMSE'
                    [idxhat,bithat] = FD_LMMSE(par,H,y,N0);
                otherwise,
                    error('par.detector type not defined.')
            end
            
            % -- compute error metrics
            err = (idx~=idxhat);
            res.VER(d,k) = res.VER(d,k) + any(err);
            res.SER(d,k) = res.SER(d,k) + sum(err)/par.U;
            res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.U*par.Q);
            
        end % algorithm loop
        
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structure)
if par.save
    save([ par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plot)

marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
figure(1)
for d=1:length(par.detector)
    if d==1
        semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
        hold on
    else
        semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
    end
end
hold off
grid on
xlabel('average SNR per receive antenna [dB]','FontSize',12)
ylabel('symbol error rate (SER)','FontSize',12)
axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)

end

% -- detectord and feedforward architecture functions

%% maximum ratio combining (MRC) detector (unbiased)
function [idxhat,bithat] = MRC(par,H,y)
shat = (1./diag(H'*H)).*(H'*y);
[~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
shat = H\y;
[~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% unbiased linear MMSE detector (LMMSE)
function [idxhat,bithat] = LMMSE(par,H,y,N0)
% preprocessing
Ainv = inv(H'*H+(N0/par.Es)*eye(par.U));
WH = Ainv*H';
% compensate for bias
% magic trick is to use matrix inversion lemma to compute bias
% diag(WH_c*H_c') = 1 - N0/par.Es*diag(Ainv)
debias_c = 1./(1-N0/par.Es*diag(Ainv));
shat = debias_c.*(WH*y);
% detection
[~,idxhat] = min(abs(shat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% Partially-decentralized LMMSE (PD-LMMSE) detector (feedforward)
function [idxhat,bithat] = PD_LMMSE(par,H,y,N0)

% -- decentralized preprocessing
C_l = par.B/par.C;
G_c = zeros(par.U,par.U,par.C);
y_MRC_c = zeros(par.U,par.U,par.C);
for cc=1:par.C
    H_c = H(C_l*(cc-1)+1:C_l*cc,:);
    G_c(:,:,cc) = H_c'*H_c;
    y_c = y(C_l*(cc-1)+1:C_l*cc);
    y_MRC_c(:,cc) = H_c'*y_c;
end

% -- fuse Gram matrix and matched filter
G = zeros(par.U,par.U);
y_MRC = zeros(par.U,1);
for cc=1:par.C
    G = G + G_c(:,:,cc);
    y_MRC = y_MRC + y_MRC_c(:,cc);
end

% -- equalization
s_hat = (G + eye(par.U)*N0/par.Es)\y_MRC;

% -- perform detection
[~,idxhat] = min(abs(s_hat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end


%% Fully-decentralized LMMSE (FD-LMMSE) detector (feedforward)
function [idxhat,bithat] = FD_LMMSE(par,H,y,N0)

% decentralized preprocessing and equalization
C_l = par.B/par.C;
shat_c = zeros(par.U,par.C);
sigma2_c = zeros(par.U,par.C);
for cc=1:par.C
    % extract local data
    H_c = H(C_l*(cc-1)+1:C_l*cc,:);
    y_c = y(C_l*(cc-1)+1:C_l*cc);
    % local preprocessing
    G_c = H_c'*H_c;
    Ainv = inv(G_c + eye(par.U)*N0/par.Es);
    WH_c = Ainv*H_c';
    % compensate for bias
    % magic trick is to use matrix inversion lemma to compute bias
    % diag(WH_c*H_c') = 1 - N0/par.Es*diag(Ainv)
    debias_c = 1./(1-N0/par.Es*diag(Ainv));
    % local detection
    shat_c(:,cc) = debias_c.*(WH_c*y_c);
    sigma2_c(:,cc) = par.Es*real(diag((WH_c*H_c-eye(par.U))*(WH_c*H_c-eye(par.U))'))+N0*real(diag(WH_c*WH_c'));
end

% fuse equalized signals by weighting
nu_c = diag(1./sum(1./sigma2_c,2))*(1./sigma2_c);
s_hat = sum(nu_c.*shat_c,2);

% centralized equalization and detection
[~,idxhat] = min(abs(s_hat*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end
