% 5. Variance decomposition
clear
clc

for Model = 2:6
switch(Model)
    case 1
        disp('Benchmark')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Benchmark_2019.mat';
    case 2
        disp('Primiceri')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Primiceri_2019.mat';
    case 3
        disp('Benchmark covariance constant')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Benchmark_covariance_constant_2019.mat';
    case 4
        disp('Benchmark variance covariance constant')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Benchmark_variance_covariance_constant_2019.mat';
    case 5
        disp('Benchmark coefficient constant')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Benchmark_coefficients_constant_2019.mat';
    case 6
        disp('Standard VAR')
        load 'C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\Standard_Invariant_VAR_2019.mat';
end 
% Shock
shocks=1:M;     
% Responses
responses=1:M;
        
Ndraws=10000; 
        
% Select horizon
nhor=20;    
    
% Select final date
% -----------------
yy1=2019;   %year
qq1=3;      %quarter     
date1=yy1+round2(qq1/4,0.0001); %Con 2017 + q=3 tienes al final 2017q4
    
impresp = zeros(Ndraws,M,nhor*M,t);
vd = zeros (Ndraws,M,nhor*M,t);
    
%Matrices previas
bigj = zeros(M,M*p);
bigj(1:M,1:M) = eye(M);

it_print = 100;  % Print in the screen every "it_print"-th iteration
tic;
disp('Number of iterations')

for k=1:Ndraws
    
    % Print iterations
    if mod(k,it_print) == 0
        disp(k);toc;
    end
    
    R=round(unifrnd(1,nrep/nthin));
    
    biga = zeros(M*p,M*p);
    
    for j = 1:p-1
        biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
    end
    
    At_posttemp=At_post{R};
    Bt_posttemp=Bt_post{R};
    Sigt_temp=Sigt_post{R};    
    
    impresp_t = zeros(M,M*nhor,t);
    vd_t= zeros(M,nhor*M,t);
for i = 1:t %Get impulses recurssively for each time period
    
    capatemp=reshape(S_A*At_posttemp(:,i)+s_A,M,M);          
    stem = diag(Sigt_temp(:,i));   
    Hsd = capatemp\stem;
    Ht = Hsd*Hsd'; 
                  
    bbtemp = Bt_posttemp(M+1:K,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
   % ------------Identification code:                
   % St dev matrix for structural VAR
    diagonal = diag(diag(Hsd));
    Hsd = inv(diagonal)*Hsd;    % Unit initial shock
    impresp_temp = zeros(M,M*nhor); % First shock is the Cholesky of the VAR covariance
    impresp_temp(1:M,1:M) = Hsd;
    bigai = biga;
    for j = 1:nhor-1
        impresp_temp(:,j*M+1:(j+1)*M) = bigj*bigai*bigj'*Hsd;    
        bigai = bigai*biga;
    end
    impresp_t(:,:,i) = impresp_temp;
    irf_t = impresp_temp;
    
    %FEVD
    mse=zeros(M,nhor);
    contr=zeros(M,M,nhor);
    for j=1:nhor
        temp2=eye(M);
        for f=1:M
            if j==1
                mse(f,j)=temp2(:,f)'*irf_t(:,((j-1)*M+1):(j*M))*Ht*irf_t(:,((j-1)*M+1):(j*M))'*temp2(:,f);
                for l=1:M
                    contr(f,l,j)=(temp2(:,f)'*irf_t(:,((j-1)*M+1):(j*M))*chol(Ht)'*temp2(:,l))^2;
                end
            else
                mse(f,j)=mse(f,j-1)+temp2(:,f)'*irf_t(:,((j-1)*M+1):(j*M))*Ht*irf_t(:,((j-1)*M+1):(j*M))'*temp2(:,f);
                for l=1:M
                    contr(f,l,j)=contr(f,l,j-1)+(temp2(:,f)'*irf_t(:,((j-1)*M+1):(j*M))*chol(Ht)'*temp2(:,l))^2;
                end
            end
        end
    end
    vd_temp=zeros(M,nhor*M);
    for ll=1:M
        vd_temp(:,(ll-1)*nhor+1:(ll*nhor))=squeeze(contr(:,ll,:))./mse;
    end
    vd_t(:,:,i)=vd_temp;           
end %END geting impulses for each time period
    impresp(k,:,:,:)=impresp_t;
    vd(k,:,:,:)=vd_t;
end

toc; % Stop timer and print total time

vd_mean=squeeze(mean(vd));

switch(Model)
    case 1
        disp('FEVD Benchmark')
        vd_mean_benchmark = vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_benchmark','vd_mean_benchmark')
    case 2
        disp('FEVD Primiceri')
        vd_mean_primiceri = vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_primiceri','vd_mean_primiceri')
    case 3
        disp('FEVD Benchmark covariance constant')
        vd_mean_covariance = vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_benchmark_covariance_constant','vd_mean_covariance')
    case 4
        disp('FEVD Benchmark variance covariance constant')
        vd_mean_variance_covariance = vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_benchmark_variance_covariance_constant','vd_mean_variance_covariance')
    case 5
        disp('FEVD Benchmark coefficient constant')
        vd_mean_coefficient = vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_benchmark_coefficient_constant','vd_mean_coefficient')
    case 6
        disp('FEVD Standard VAR')
        vd_mean_standard_var =vd_mean;
        save('C:\Users\DELL\Desktop\Tesis-Paper\Resultados_2021\FEVD_standard_var','vd_mean_standard_var')
end
clear
close
clc
warning off all
end