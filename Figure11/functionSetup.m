function [gainOverNoisedB_new,pilotIndex] = functionSetup(K,L,tau_p)

%Size of the coverage area (as a square with wrap-around)
squareLength = 500; %meter

%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 9;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Parameters for the shadow fading from [15]
sigma_sf = 8;
delta = 0.5;
decorr = 100;

%Prepare to save results
distancesUE = zeros(L,K);

gainOverNoisedB = zeros(L,K);
gainOverNoisedB_new = zeros(K,L);

%% Go through all setups

    
    %Random AP locations with uniform distribution
    APpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Random UE locations with uniform distribution
    UEpositions = (rand(K,1) + 1i*rand(K,1)) * squareLength;
       
    %Compute alternative AP locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    UEpositionsWrapped = repmat(UEpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[K 1]);
        
    %Compute the correlation matrices for the shadow fading
    shadowCorrMatrix_APs = zeros(L,L);
    shadowCorrMatrix_UEs = zeros(K,K);
    
    for l = 1:L
        distancetoAP = min(abs(APpositionsWrapped - repmat(APpositions(l),size(APpositionsWrapped))),[],2);
        shadowCorrMatrix_APs(:,l) = 2.^(-distancetoAP/decorr);
    end
    
    for k = 1:K
        distancetoUE = min(abs(UEpositionsWrapped - repmat(UEpositions(k),size(UEpositionsWrapped))),[],2);
        shadowCorrMatrix_UEs(:,k) = 2.^(-distancetoUE/decorr);
    end
      
    %Generate shadow fading realizations
    a = sigma_sf*sqrtm(shadowCorrMatrix_APs)*randn(L,1);
    b = sigma_sf*sqrtm(shadowCorrMatrix_UEs)*randn(K,1);
    
    for k = 1:K
        
        %Compute distances between each of the APs to UE k
        [distancetoUE,~] = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
        distancesUE(:,k) = distancetoUE;
        
        %Compute the channel gain divided by the noise power
        gainOverNoisedB(:,k) = pathloss_threeslope(distancesUE(:,k)) - noiseVariancedBm;
        
        %Add shadow fading to all channels from APs to UE k that have a
        %distance that is larger than 50 meters
        gainOverNoisedB(distancetoUE>50,k) = gainOverNoisedB(distancetoUE>50,k) + sqrt(delta)*a(distancetoUE>50) + sqrt(1-delta)*b(k);
        
    end

for k=1:K
   for l=1:L
    gainOverNoisedB_new(k,l)=gainOverNoisedB(l,k);
   end
end

    %Assign random pilots while guaranteeing that each pilot is used
    %equally many times, as an initiation to the greedy algorithm from [15]
    pilotIndex = mod(randperm(K),tau_p)+1;






