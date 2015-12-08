

%% this is to simulate power as a function of number of phenotypes combined.
printfig=true;
%% specify parameters %%
% n: sample size
% maf: testing snp minor allele freqency
% pct_explained: total variances explained by a major snp.
% m: # of phenotypes
% r_g: genotypic correlation
% r_p: phenotypic correlation
% h2: heritability of the phenotypes
% fix h2=0.4 and r_g=0.9
h2_1=0.4;
h2_2=0;

n=2996;
sigma_x=1;
b1 = -0.15:0.01:0.15;
b2 = -0.15:0.01:0.15;
m=1:5;
ColOrd = jet(6);

r_g=0;
r_p=0.8;
g=r_g*sqrt(h2_1*h2_1);

p_max = zeros(length(b1),length(b2));
p_pca = zeros(length(b1),length(b2));
p_single = zeros(length(b1),length(b2));

m=2;

%% Noncentral parameters %%%%

for i=1:length(b1)
    for j=1:length(b2)
        %ncp and power
        a2_1 = b1(i)^2;
        a2_2 = b2(j)^2;
        %delta = max(n*a2_1*sigma_x/(1-a2_1*sigma_x),n*a2_2*sigma_x/(1-a2_2*sigma_x));
        delta = n*a2_1/(1-a2_1);
 
        %first is the single phenotype
        p_single(i,j)=1-ncx2cdf(17.764,1,delta);

        V_g=zeros(m,m);
        V_p=eye(m,m);
        V_g(logical(eye(size(V_g))))=[h2_1 h2_2];
        V_p(~logical(eye(size(V_p))))=r_p;
        V_g(~logical(eye(size(V_g))))=g;
        A=V_p\V_g;
        [L,H] = eig(A);
        [Lpca,Hpca]=eig(V_p);
        [g2_max,I]=max(diag(H));
        [g2_pca,Ipca]=max(diag(Hpca));
        
        a2_1 = b1(i);
        a2_2 = b2(j);
        %if(g2_max>1)
        %    p_max(i,j) = 0;
        %    continue
        %else
            v=L(:,I);
            a_max = (v'*[a2_1;a2_2])^2;
            delta = n*a_max/(v'*V_p*v-a_max);
            %delta = n*a_max*sigma_x/(1-a_max*sigma_x);
            p_max(i,j) = 1-ncx2cdf(16.448,1,delta);
        %end
        vpca=Lpca(:,Ipca);
        a_max = (vpca'*[a2_1;a2_2])^2;
        %delta = n*a_max*sigma_x/(1-a_max*sigma_x);
        delta = n*a_max/(vpca'*V_p*vpca-a_max);
        p_pca(i,j) = 1-ncx2cdf(16.448,1,delta);
    end
end

figure;
fs = 10;
set(gca,'FontSize',fs);
pcolor(b1,b2,p_max)
colorbar
xlabel('b_1');
ylabel('b_2');
if printfig
    filename = 'power_b1_vs_b2_maxh';
    print('-depsc2',filename);
end

figure;
fs = 10;
set(gca,'FontSize',fs);
pcolor(b1,b2,p_pca)
colorbar
xlabel('b_1');
ylabel('b_2');
if printfig
    filename = 'power_b1_vs_b2_pca';
    print('-depsc2',filename);
end

figure;
fs = 10;
set(gca,'FontSize',fs);
pcolor(b1,b2,p_single)
colorbar
xlabel('b_1');
ylabel('b_2');
if printfig
    filename = 'power_b1_vs_b2_single';
    print('-depsc2',filename);
end


%% Plot %%
figure;
fs = 10;
set(gca,'FontSize',fs);
plot(m,p_max,':.','Color',ColOrd(1,:),'MarkerSize',30);

%% this is to simulate power as a function of number of phenotypes combined.
printfig=true;
%% specify parameters %%
% n: sample size
% maf: testing snp minor allele freqency
% pct_explained: total variances explained by a major snp.
% m: # of phenotypes
% r_g: genotypic correlation
% r_p: phenotypic correlation
% h2: heritability of the phenotypes
% fix h2=0.4 and r_g=0.9
g_cor=[0.9 0.8 0.7];
p_cor=0.4;
h2=0.4;

n=2996;
sigma_x=1;
pct_explained = 0.03;
m=1:5;
ColOrd = jet(6);

r_g=g_cor(1);
r_p=p_cor;
g=r_g*sqrt(h2*h2);


%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);

%ncp and power
p_max = zeros(1,length(m));
%first is the single phenotype
p_max(1)=1-ncx2cdf(17.764,1,b);

%ncp and power
p_max_pca = zeros(1,length(m));
%first is the single phenotype
p_max_pca(1)=1-ncx2cdf(17.764,1,b);


for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);
    [Lpca,Hpca]=eig(V_p);
    [g2_max,I]=max(diag(H));
    [g2_pca,Ipca]=max(diag(Hpca));
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(16.448,1,s_max);
    end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(16.448,1,s_max);
end


%% Plot %%
figure;
fs = 10;
set(gca,'FontSize',fs);
plot(m,p_max,':.','Color',ColOrd(1,:),'MarkerSize',30);
hold on

%%
r_g=g_cor(2);
g=r_g*sqrt(h2*h2);

%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
end


%% Plot %%
hold on 
plot(m,p_max,':.','Color',ColOrd(3,:),'MarkerSize',30);

%%
r_g=g_cor(3);
g=r_g*sqrt(h2*h2);

%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    %if(g2_max>1)
    %    p_max(i) = 0;
    %else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    %end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained;
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
end


%% Plot %%
hold on 
plot(m,p_max,':.','Color',ColOrd(5,:),'MarkerSize',30);

xlabel('# of phenotypes');
ylabel('Power');
xlim([1 max(m)]);
ylim([0 1]);
set(gca,'Xtick',m);
name={};
for k = 1:3
   name{k} = ['r_g= ' num2str(g_cor(int8(k/2)))];
end
legend(gca, name,'Location','NorthWest');

grid on;
figure(gcf);

if printfig
    filename = 'power_positive_corr_gbp_g_vs_m';
    print('-depsc2',filename);
end

%% specify parameters %%
% n: sample size
% maf: testing snp minor allele freqency
% pct_explained: total variances explained by a major snp.
% m: # of phenotypes
% r_g: genotypic correlation
% r_p: phenotypic correlation
% h2: heritability of the phenotypes
% fix h2=0.4 and r_g=0.9
g_cor=[0.85 0.8 0.75];
p_cor=0.9;
h2=0.4;

n=3000;
maf=0.2;
%sigma_x = 2*maf*(1-maf);
sigma_x=1;
pct_explained = 0.01;
m=1:5;
ColOrd = jet(6);

r_g=g_cor(1);
r_p=p_cor;
g=r_g*sqrt(h2*h2);


%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);

%ncp and power
p_max = zeros(1,length(m));
%first is the single phenotype
p_max(1)=1-ncx2cdf(16.448,1,b);

%ncp and power
p_max_pca = zeros(1,length(m));
%first is the single phenotype
p_max_pca(1)=1-ncx2cdf(16.448,1,b);


for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);
    [Lpca,Hpca]=eig(V_p);
    [g2_max,I]=max(diag(H));
    [g2_pca,Ipca]=max(diag(Hpca));
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
end


%% Plot %%
figure;
fs = 10;
set(gca,'FontSize',fs);
plot(m,p_max,':.','Color',ColOrd(1,:),'MarkerSize',30);
hold on
plot(m,p_max_pca,':.','Color',ColOrd(2,:),'MarkerSize',30);
%sim_power=[1 0.2442;2 0.3432;3 0.3655;4 0.3656];
%plot(sim_power(:,1),sim_power(:,2),'bo');

%%
r_g=g_cor(2);
g=r_g*sqrt(h2*h2);

%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
end


%% Plot %%
hold on 
plot(m,p_max,':.','Color',ColOrd(3,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(4,:),'MarkerSize',30);
%sim_power=[1 0.2478;2 0.3152;3 0.3331;4 0.3559];
%plot(sim_power(:,1),sim_power(:,2),'ro');

%%
r_g=g_cor(3);
g=r_g*sqrt(h2*h2);

%% Noncentral parameters %%%%

a2 = h2*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_g=zeros(m(i),m(i));
    V_p=eye(m(i),m(i));
    V_g(logical(eye(size(V_g))))=h2;
    V_p(~logical(eye(size(V_p))))=r_p;
    V_g(~logical(eye(size(V_g))))=g;
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    %if(g2_max>1)
    %    p_max(i) = 0;
    %else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    %end
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained;
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
end


%% Plot %%
hold on 
plot(m,p_max,':.','Color',ColOrd(5,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(6,:),'MarkerSize',30);
%sim_power=[1 0.2432;2 0.3012;3 0.3047;4 0.3164];
%plot(sim_power(:,1),sim_power(:,2),'ko');

xlabel('# of phenotypes');
ylabel('Power');
xlim([1 max(m)]);
ylim([0 1]);
set(gca,'Xtick',m);
name={};
for k = 1:2:6
   name{k} = ['r_g= ' num2str(g_cor(int8(k/2)))];
   name{k+1} = ['PCA r_g= ' num2str(g_cor(int8(k/2)))];
end
legend(gca, name,'Location','NorthWest');
grid on;
figure(gcf);

if printfig
    filename = 'power_positive_corr_pbg_g_vs_m';
    print('-depsc2',filename);
end

%% specify parameters %%
%% this is to simulate the cases where k^2 is dropping.
% k=0.98, k=0.95, k=0.9
% m: # of phenotypes
% r_g: genotypic correlation
% r_p: phenotypic correlation
% fix r_p=0.9 and r_g=0.8
% fix highest het h2_1=0.4
% h2=k^i*h2_1
prop=[0.95 0.9 0.8];

m=1:5;
r_g=0.9;
r_p=0.4;
k=prop(1);
h2_1=0.4;

ColOrd = jet(6);

%% Noncentral parameters %%%%

a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
 
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);

p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%
figure;
fs = 10;
set(gca,'FontSize',fs);
plot(m,p_max,':.','Color',ColOrd(1,:),'MarkerSize',30);
hold on;
plot(m,p_max_pca,':.','Color',ColOrd(2,:),'MarkerSize',30);

%% specify parameters %%
% m: # of phenotypes
% g: genotypic correlation
% r: phenotypic correlation

k=prop(2);
ColOrd = jet(6);

%% Noncentral parameters %%%%

a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);

p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%

hold on
plot(m,p_max,':.','Color',ColOrd(3,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(4,:),'MarkerSize',30);

%% specify parameters %%
% m: # of phenotypes
% g: genotypic correlation
% r: phenotypic correlation

k=prop(3);

%% Noncentral parameters %%%%


a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%

hold on
plot(m,p_max,':.','Color',ColOrd(5,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(6,:),'MarkerSize',30);

figure(gcf);
xlabel('# of phenotypes');
ylabel('Power');
ylim([0 1]);
xlim([1 max(m)]);
set(gca,'Xtick',m);
name={};
for k = 1:2:6
   name{k} = ['k= ' num2str(prop(int8(k/2)))];
   name{k+1} = ['PCA k= ' num2str(prop(int8(k/2)))];
end
legend(gca, name,'Location','NorthWest');
grid on;

if printfig
    filename = 'power_positive_corr_gbp_unequal_g_vs_m';
    print('-depsc2',filename);
end

%% specify parameters %%
%% this is to simulate the cases whhere k^2 is dropping.
% k=0.98, k=0.95, k=0.9
% m: # of phenotypes
% r_g: genotypic correlation
% r_p: phenotypic correlation
% fix r_p=0.9 and r_g=0.8
% fix highest het h2_1=0.4
% h2=k^i*h2_1
prop=[0.95 0.9 0.8];

m=1:5;
r_g=0.8;
r_p=0.9;
k=prop(1);
h2_1=0.4;

ColOrd = jet(6);

%% Noncentral parameters %%%%

a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
 
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);

p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%
figure;
fs = 10;
set(gca,'FontSize',fs);
plot(m,p_max,':.','Color',ColOrd(1,:),'MarkerSize',30);
hold on;
plot(m,p_max_pca,':.','Color',ColOrd(2,:),'MarkerSize',30);

%% specify parameters %%
% m: # of phenotypes
% g: genotypic correlation
% r: phenotypic correlation

k=prop(2);
h2=zeros(1,length(m));
ColOrd = jet(6);

%% Noncentral parameters %%%%

a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);

p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%

hold on
plot(m,p_max,':.','Color',ColOrd(3,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(4,:),'MarkerSize',30);

%% specify parameters %%
% m: # of phenotypes
% g: genotypic correlation
% r: phenotypic correlation

k=prop(3);

%% Noncentral parameters %%%%


a2 = h2_1*pct_explained;
b = n*a2*sigma_x/(1-a2*sigma_x);
p_max = zeros(1,length(m));
p_max(1)=1-ncx2cdf(16.448,1,b);
p_max_pca = zeros(1,length(m));
p_max_pca(1)=1-ncx2cdf(16.448,1,b);

for i=2:length(m)
    V_p=eye(m(i),m(i));
    V_p(~logical(eye(size(V_p))))=r_p;
    
    V_g=zeros(m(i),m(i));
    for l=1:m(i)
        V_g(l,l)=h2_1*k^(l-1);
        for j=1:l-1
            V_g(j,l)=r_g*sqrt(V_g(l,l)*V_g(j,j));
            V_g(l,j)=V_g(j,l);
        end
    end
    
    A=V_p\V_g;
    [L,H] = eig(A);       
    [g2_max,I]=max(diag(H));
    [Lpca,Hpca] = eig(V_p);       
    [g2_max_pca,Ipca]=max(diag(Hpca));
    vpca=Lpca(:,Ipca);
    a_max = ((vpca'*V_g*vpca)./(vpca'*V_p*vpca))*pct_explained; 
    s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
    p_max_pca(i) = 1-ncx2cdf(15.136,1,s_max);
    if(g2_max>1)
        p_max(i) = 0;
        continue
    else
        v=L(:,I);
        a_max = ((v'*V_g*v)./(v'*V_p*v))*pct_explained; 
        s_max = n*a_max*sigma_x/(1-a_max*sigma_x);
        p_max(i) = 1-ncx2cdf(15.136,1,s_max);
    end
end


%% Plot %%

hold on
plot(m,p_max,':.','Color',ColOrd(5,:),'MarkerSize',30);
plot(m,p_max_pca,':.','Color',ColOrd(6,:),'MarkerSize',30);

figure(gcf);
xlabel('# of phenotypes');
ylabel('Power');
ylim([0 1]);
xlim([1 max(m)]);
set(gca,'Xtick',m);
name={};
for k = 1:2:6
   name{k} = ['k= ' num2str(prop(int8(k/2)))];
   name{k+1} = ['PCA k= ' num2str(prop(int8(k/2)))];
end
legend(gca, name,'Location','NorthWest');
grid on;
if printfig
    filename = 'power_positive_corr_pbg_unequal_g_vs_m';
    print('-depsc2',filename);
end