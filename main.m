%% Step 1 : Reading data : 
data1 = edfread('S002R03.edf');
info1 = edfinfo('S002R03.edf');
data2 = edfread('S001R04.edf');
info2 = edfinfo('S001R04.edf');
data3 = edfread('S001R03.edf');
info3 = edfinfo('S001R03.edf');
%% *********************************
%% Step : 2 Changing time table to table with all seconds 
T1=zeros([64,19680]);

for j=1:1:64
    A =[data1.(j)];
    Y=[A{1,:}];
    for i=2:1:123
        Y=[Y;A{i,:}];
    end
    T1(j,:)=Y;
end
T2=zeros([64,20000]);

for j=1:1:64
    A =[data2.(j)];
    K=[A{1,:}];
    for i=2:1:125
        K=[K;A{i,:}];
    end
    T2(j,:)=K;
end

T3=zeros([64,20000]);
for j=1:1:64
    A =[data3.(j)];
    O=[A{1,:}];
    for i=2:1:125
        O=[O;A{i,:}];
    end
    T3(j,:)=O;
end
%% ************************************
%%  Step 3 : Extracting annotations  
info1.Annotations
sample1=seconds(info1.Annotations.(2)).*160;
Annotation1 = string(info1.Annotations.(1));
time1 = [0;4.1;8.2;12.3;16.4;20.5 ;24.6;28.7;32.8;36.9;41 ;45.1;49.2;53.3;57.4 ;61.5;65.6 ;69.7;73.8;77.9;82 ;86.1;90.2;94.3 ;98.4;102.5;106.6 ;110.7 ;114.8 ;118.9].*160;
Label1 = [sample1 Annotation1  time1 ];

info2.Annotations
sample2=seconds(info2.Annotations.(2)).*160;
Annotation2 = string(info2.Annotations.(1));
time2 = [0 ;4.2;8.3;12.5;16.6;20.8;24.9;29.1;33.2;37.4;41.5;45.7;49.8;54;58.1;62.3;66.4;70.6;74.7;78.9;83 ;87.2;91.3;95.5;99.6;103.8;107.9;112.1 ;116.2 ;120.4 ]*160 ;
Label2 = [sample2 Annotation2  time2 ];

info3.Annotations
sample3=seconds(info3.Annotations.(2)).*160;
Annotation3 = string(info3.Annotations.(1));
time3 = [0 ;4.2;8.3;12.5;16.6;20.8;24.9;29.1;33.2;37.4;41.5;45.7;49.8;54;58.1;62.3;66.4;70.6;74.7;78.9;83 ;87.2;91.3;95.5;99.6;103.8;107.9;112.1 ;116.2 ;120.4 ]*160 ;
Label3 = [sample2 Annotation3  time2 ];
%% ***********************  Visualization *****************************
%% Spectrogram 
spectrogram(mean(T1),'yaxis')
pause 
spectrogram(mean(T2),'yaxis')
pause
spectrogram(mean(T3),'yaxis')
fs = 160;
ts = 0:1/fs:2;
M = 49;
L = 11;
g = bartlett(M);
Ndft = 1024;
%%
[s,f,t] = spectrogram(mean(T1),g,L,Ndft,fs);
waterplot(s,f,t)

%% %% ploting all the chanels 
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(T2(i,:))
    title(strcat("Canal ",info2.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% ***************************************
%%  computing spectre 1
Spec_Y1=zeros([64,19680]);
Pos1= zeros([64,9840]);

for i=1:1:64 
    Spec_Y1(i,:) = abs(fft(T1(i,:)));
    Pos1(i,:)= Spec_Y1(i,1:9840);
end 
%%  computing spectre 2
Spec_Y2=zeros([64,20000]);
Pos2= zeros([64,10000]);

for i=1:1:64 
    Spec_Y2(i,:) = abs(fft(T2(i,:)));
    Pos2(i,:)= Spec_Y2(i,1:10000);
end 
%%  computing spectre 2
Spec_Y3=zeros([64,20000]);
Pos3= zeros([64,10000]);

for i=1:1:64 
    Spec_Y3(i,:) = abs(fft(T3(i,:)));
    Pos3(i,:)= Spec_Y3(i,1:10000);
end 
%% ***************************************
%% ploting spectre 1
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(Pos1(i,:))
    title(strcat("Spectre of  ",info1.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% ploting spectre 2
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(Pos2(i,:))
    title(strcat("Spectre of  ",info2.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% ploting spectre 3
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(Pos3(i,:))
    title(strcat("Spectre of  ",info3.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% Mean of spectre 1
M = mean(Pos1,1);
plot(M)
hold on
%% Mean of spectre 2
M = mean(Pos2,1);
plot(M)
hold on
%% **************************************
%% Filtring 60 Hz 
fs = 160;             % sampling rate
f0 = 60;                % notch frequency
fn = fs/2;              % Nyquist frequency
freqRatio = f0/fn;      % ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       % width of the notch

% Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
notchPoles = (1-notchWidth) * notchZeros;



b = poly( notchZeros ); %  Get moving average filter coefficients
a = poly( notchPoles ); %  Get autoregressive filter coefficients



% Filtring de 1 Hz 
fs = 160;             % sampling rate
f0 = 0;                % notch frequency
fn = fs/2;              % Nyquist frequency
freqRatio = f0/fn;      % ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       % width of the notch

% Compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
notchPoles = (1-notchWidth) * notchZeros;


d = poly( notchZeros ); %  Get moving average filter coefficients
c = poly( notchPoles ); %  Get autoregressive filter coefficients



%% Loop for all canals 
Fil_T1 = zeros([64,19680]);
Fil_Spec_Y1=zeros([64,19680]);
Fil_Pos1= zeros([64,9840]);
for i=1:1:64 
   Fil_T1(i,:) = filter(b,a,T1(i,:));
   Fil_T1(i,:) = filter(d,c,Fil_T1(i,:));
   Fil_Spec_Y1(i,:) = abs(fft(Fil_T1(i,:)));
   Fil_Pos1(i,:)= Fil_Spec_Y1(i,1:9840);
end
Fil_T2 = zeros([64,20000]);
Fil_Spec_Y2=zeros([64,20000]);
Fil_Pos2= zeros([64,10000]);
for i=1:1:64 
   Fil_T2(i,:) = filter(b,a,T2(i,:));
   Fil_T2(i,:) = filter(d,c,Fil_T2(i,:));
   Fil_Spec_Y2(i,:) = abs(fft(Fil_T2(i,:)));
   Fil_Pos2(i,:)= Fil_Spec_Y2(i,1:10000);
end
%%
Fil_T3 = zeros([64,20000]);
Fil_Spec_Y3=zeros([64,20000]);
Fil_Pos3= zeros([64,10000]);
for i=1:1:64 
   Fil_T3(i,:) = filter(b,a,T3(i,:));
   Fil_T3(i,:) = filter(d,c,Fil_T3(i,:));
   Fil_Spec_Y3(i,:) = abs(fft(Fil_T3(i,:)));
   Fil_Pos3(i,:)= Fil_Spec_Y3(i,1:10000);
end
%% Ploting spectre 1 
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(Fil_Pos3(i,:))
    title(strcat("Spectre of  ",info1.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% Ploting spectre 2 
m = 8;
figure(1)
for i =1:64
    subplot(m,m,i)
    plot(Fil_Pos2(i,:))
    title(strcat("Spectre of  ",info2.SignalLabels(i,1)))
    ax= gca;
    ax.XTickLabel = {};
end
%% *********************************************
%% Spectrogram of the mean of preprocessed signal 
spectrogram(mean(Fil_T),'yaxis')
%% Water plot of the mean of preprocessed signal 
fs = 160;
ts = 0:1/fs:2;
M = 49;
L = 11;
g = bartlett(M);
Ndft = 1024;

[s,f,t] = spectrogram(mean(Fil_T),g,L,Ndft,fs);
waterplot(s,f,t)
%% Filtring alpha : 
alpha = bandpass(mean(Fil_T),[8 13],160);
[s,f,t] = spectrogram(alpha,g,L,Ndft,fs);
waterplot(s,f,t)
%%
spectrogram(alpha,'yaxis')
%% Filtring Betha 
Betha = bandpass(mean(Fil_T),[14 79],160);
[s,f,t] = spectrogram(Betha,g,L,Ndft,fs);
waterplot(s,f,t)
%% 
spectrogram(Betha,'yaxis')
%% Filtring Delta 
Delta = bandpass(mean(Fil_T),[0.5 4],160);
[s,f,t] = spectrogram(Delta,g,L,Ndft,fs);
waterplot(s,f,t)
%% 
spectrogram(Delta,'yaxis')
%% Filtring Theta
Theta = bandpass(mean(Fil_T),[4 7],160);
[s,f,t] = spectrogram(Theta,g,L,Ndft,fs);
waterplot(s,f,t)
%% 
spectrogram(Theta,'yaxis')
%% Filrting for alpha : 
%J = bandpass(Fil_T,[0.4 79],160);
% Filtering the signal
Fs = 160; % Sampling Frequency
fn = Fs/2; % Nyquist Frequency
start = 13; % Start Frequency
stop = 31; % Stop Frequency
W1 = start/fn;
W2 = stop/fn;
Wn=[W1 W2];
b=fir1(6,Wn,'bandpass'); %origine

%%
Fil_data_Spec_Y1=zeros([64,19680]);
Fil_data_Pos1= zeros([64,9840]);
J1=zeros([64,19680]);
for i=1:1:64 
    J1(i,:)=filter(b,1,Fil_T1(i,:))
    Fil_data_Spec_Y1(i,:) = abs(fft(J1(i,:)));
    Fil_data_Pos1(i,:)= Fil_data_Spec_Y1(i,1:9840);
end 
Fil_data_Spec_Y2=zeros([64,20000]);
Fil_data_Pos2= zeros([64,10000]);
J2=zeros([64,20000]);
for i=1:1:64 
    J2(i,:)=filter(b,1,Fil_T2(i,:))
    Fil_data_Spec_Y2(i,:) = abs(fft(J2(i,:)));
    Fil_data_Pos2(i,:)= Fil_data_Spec_Y2(i,1:10000);
end 
%%
Fil_data_Spec_Y3=zeros([64,20000]);
Fil_data_Pos3= zeros([64,10000]);
J3=zeros([64,20000]);
for i=1:1:64 
    J3(i,:)=filter(b,1,Fil_T3(i,:))
    Fil_data_Spec_Y3(i,:) = abs(fft(J3(i,:)));
    Fil_data_Pos3(i,:)= Fil_data_Spec_Y3(i,1:10000);
end 
%%
subplot(2,1,1)
plot(mean(Fil_data_Pos3))
subplot(2,1,2)
plot(mean(Fil_Pos3))
%%
subplot(2,1,1)
plot(mean(Fil_data_Pos2))
subplot(2,1,2)
plot(mean(Fil_Pos2))
%% *******************************************
%% feature vector applying PCA
q=21;
Data_PCA1=zeros([19680,21]);


for i =0:29
    l=(i)*656;
    k=(i+1)*656;
    [coeff1,Data_PCA1(l+1:k,:),latent1,tsquared1,explained1,mu1] = pca(J1(:,l+1:k)','NumComponents',q);
end
Data_PCA2=zeros([20000,21]);
Data_PCA3=zeros([20000,21]);
for i=0:28
    [coeff3,Data_PCA3(double(Label3(i+1,3))+1:double(Label3(i+2,3)),:),latent3,tsquared3,explained3,mu3] = pca(J3(:,double(Label3(i+1,3))+1:double(Label3(i+2,3)))','NumComponents',q);
    [coeff2,Data_PCA2(double(Label3(i+1,3))+1:double(Label3(i+2,3)),:),latent2,tsquared2,explained2,mu2] = pca(J2(:,double(Label3(i+1,3))+1:double(Label3(i+2,3)))','NumComponents',q) ; 
end
disp(strcat("Top ",string(q)," principle components explain ", string(sum(explained1(1:q))), " of variation 1"))
disp(strcat("Top ",string(q)," principle components explain ", string(sum(explained2(1:q))), " of variation 2"))
 disp(strcat("Top ",string(q)," principle components explain ", string(sum(explained3(1:q))), " of variation 3"))
%%  Training ICA model
Mdl1 = rica(Data_PCA1,q);
Mdl2 = rica(Data_PCA2,q);
Mdl3 = rica(Data_PCA3,q);
%% feature vector applying PCA and ICA
Data_ICA1 = transform(Mdl1, Data_PCA);
Data_ICA2 = transform(Mdl2, Data_PCA);
%%

%% plotting 
plotsPerCol = 7;
figure(2)
for i =1:q
    subplot(plotsPerCol,ceil(q/plotsPerCol),i)
    plot(Data_ICA1(:,i).^2)
    title(strcat("Component ",string(i)," Squared"))
    ax= gca;
    ax.XTickLabel = {};
end
%%
plotsPerCol = 7;
figure(2)
for i =1:q
    subplot(plotsPerCol,ceil(q/plotsPerCol),i)
    plot(Data_ICA2(:,i).^2)
    title(strcat("Component ",string(i)," Squared"))
    ax= gca;
    ax.XTickLabel = {};
end
%% **************************************
%% changing labels to doubles
Labels1 = zeros([1,30]);
Labels2 = zeros([1,30]);
Labels3 = zeros([1,30]);
for i=1:1:30
    if Label1(i,2)=='T0'
        Labels1(1,i)=0;
    end
    if Label1(i,2)=='T1'
        Labels1(1,i)=1;
    end
    if Label1(i,2)=='T2'
        Labels1(1,i)=2;
    end
    if Label2(i,2)=='T0'
        Labels2(1,i)=0;
    end
    if Label2(i,2)=='T1'
        Labels2(1,i)=1;
    end
    if Label2(i,2)=='T2'
        Labels2(1,i)=2;
    end
     if Label3(i,2)=='T0'
        Labels3(1,i)=0;
    end
    if Label3(i,2)=='T1'
        Labels3(1,i)=1;
    end
    if Label3(i,2)=='T2'
        Labels3(1,i)=2;
    end
   
end

 

%% ******************************************
%% Concatenating data 1 *
Conc_Data1= zeros([13776,15]);
    
for i=1:1:15
    P=[Data_PCA1((2*i-1)*656+1:2*i*656,1)];
    for j=2:1:21   
        P=[P;Data_PCA1((2*i-1)*656+1:2*i*656,j)];
    end
    Conc_Data1(:,i)=P;
end
%% Concatenating data 2 * 
Conc_Data2= zeros([13776,15]);
    
for i=1:1:15
    P=[Data_PCA2(i*672+(i-1)*656+1:i*(656+672),1)];
    for j=2:1:21   
        P=[P;Data_PCA2(i*672+(i-1)*656+1:i*(656+672),j)];
    end
    Conc_Data2(:,i)=P;
end

%% Concatenating data 3 * 
Conc_Data3= zeros([13776,15]);
    
for i=1:1:15
    P=[Data_PCA3(i*672+(i-1)*656+1:i*(656+672),1)];
    for j=2:1:21   
        P=[P;Data_PCA3(i*672+(i-1)*656+1:i*(656+672),j)];
    end
    Conc_Data3(:,i)=P;
end
%% Sub labels 
Sublabel1 = zeros([1,15]); 
Sublabel2 = zeros([1,15]); 
Sublabel3 = zeros([1,15]); 
for i=1:1:15
    Sublabel1(1,i)=Labels1(1,2*i);
    Sublabel2(1,i)=Labels2(1,2*i);
    Sublabel3(1,i)=Labels3(1,2*i);
end

%% ****************************
%% Concatenating data 2
Conc_Data2= zeros([13776,30]);
    
for i=0:1:29
     N=[Data_PCA2(double(Label2(i+1,3))+1:double(Label2(i+1,3))+656,1)];
    for j=2:1:21   
        N=[N;Data_PCA2(double(Label2(i+1,3))+1:double(Label2(i+1,3))+656,j)];
    end
    Conc_Data2(:,i+1)=N;
end
%% Concatenating data 3
Conc_Data3= zeros([13776,30]);
    
for i=0:1:29
     N=[Data_PCA3(double(Label3(i+1,3))+1:double(Label3(i+1,3))+656,1)];
    for j=2:1:21   
        N=[N;Data_PCA3(double(Label3(i+1,3))+1:double(Label3(i+1,3))+656,j)];
    end
    Conc_Data3(:,i+1)=N;
end
%% Making Input
Input1=[Conc_Data1 ; Sublabel1];
Input2=[Conc_Data2 ; Sublabel2];
Input3=[Conc_Data3 ; Sublabel3];
%%
Input =[ Input1 Input2(:,1:14) Input3(:,1:14)]'
%%

%% importing table 
writematrix(Input,'Input.csv');
%% Test
writematrix(Input2','Input2.csv');