readExactData = 1;
SoXdata = 1;

if readExactData
    T4days = readtable('220128_Raw_count_data_clonesize_update.xlsx','Sheet','1, 2, and 3 days (3 day ind)');
    T1months = readtable('220128_Raw_count_data_clonesize_update.xlsx','Sheet','1 Month clonal analysis');
    T3months = readtable('220128_Raw_count_data_clonesize_update.xlsx','Sheet','3 Month clonal analysis');
    T6months = readtable('220128_Raw_count_data_clonesize_update.xlsx','Sheet','6 Month clonal analysis');
    T12months = readtable('220128_Raw_count_data_clonesize_update.xlsx','Sheet','12 Month clonal analysis');
else
    T4days = readtable('211018_lineagetracing overview.xlsx','Sheet','1, 2, and 3 days (3 day ind)');
    T1months = readtable('211018_lineagetracing overview.xlsx','Sheet','1 Month clonal analysis');
    T3months = readtable('211018_lineagetracing overview.xlsx','Sheet','3 Month clonal analysis');
    T6months = readtable('211018_lineagetracing overview.xlsx','Sheet','6 Month clonal analysis');
    T12months = readtable('211018_lineagetracing overview.xlsx','Sheet','12 Month clonal analysis');
end
clear Td

Td{1}.data = [[1:5]' table2array(T4days(13:17,2:4))];
Td{1}.info = '1 day clonal count';
Td{2}.data = [[1:5]' table2array(T4days(13:17,5:7))];
Td{2}.info = '2 days clonal count';
Td{3}.data = [[1:5]' table2array(T4days(13:17,8:10))];
Td{3}.info = '3 days clonal count';
Td{4}.data = [[1:7]' table2array(T4days(2:8,2:4))];
Td{4}.info = '4 day clonal count';
Td{5}.data = [[1:7]' table2array(T4days(2:8,5:6))];
Td{5}.info = '5 days clonal count';
Td{6}.data = [[1:7]' table2array(T4days(2:8,8:10))];
Td{6}.info = '6 days clonal count';
Td{7}.data = table2array(T1months);
Td{7}.info = '1 month clonal analysis';
Td{8}.data = table2array(T3months(3:end,:));
Td{8}.info = '3 month clonal analysis';
Td{9}.data = table2array(T6months(3:end,:));
Td{9}.info = '6 month clonal analysis';
Td{10}.data = table2array(T12months(3:end,:));
Td{10}.info = '12 month clonal analysis';

for jj = 1:length(Td)
    for kk = 1:size(Td{jj}.data,1)
        for ll = 1:size(Td{jj}.data,2)
            if isnan(Td{jj}.data(kk,ll))
                if ll == 1
                    Td{jj}.data(kk,ll) = 2^(4+kk-32)+1;
                else
                    Td{jj}.data(kk,ll) = 0;
                end
            end
        end
    end
end

% Find average and sum clone sizes
for jj = 1:length(Td)
    Td{jj}.avg = mean(Td{jj}.data(:,2:end),2);
    Td{jj}.sum = sum(Td{jj}.data(:,2:end),2);
end

TdLong{1} = Td{7};
TdLong{2} = Td{8};
TdLong{3} = Td{9};
TdLong{4} = Td{10};


clear actualCloneSizeBins

for jj = 1:length(Td)
    % populate bins
    nbrBins = ceil(log2(Td{jj}.data(end,1)));
    bins = zeros(nbrBins,1);
    for kk = 1:size(Td{jj}.data,1)
        addToBin = floor(log2(Td{jj}.data(kk,1)))+1;
        bins(addToBin) = bins(addToBin)+Td{jj}.avg(kk);
    end
    if 0
        plotbins = [bins;bins(end)];
        figure;stairs(2.^[0:length(plotbins)-1],plotbins)
        set(gca,'XScale','log')
        set(gca,'XTick',2.^[1:length(plotbins)]);
        set(gca,'Xlim',[1 2^(length(plotbins)-1)]);
        set(gca,'Ylim',[0 max(plotbins)+2]);
    end
    actualCloneSizeBins{jj} = bins;
end





Td3table = readtable('211018_3m_lineage.xlsx','Sheet','220121_sorted ');
Td3.data = table2array(Td3table);
%Td3.data = Td3.data(3:end,:);
Td3.info = '3 month clonal analysis new';

Td3.data = [Td3.data(:,1:4) Td3.data(:,6:8) Td3.data(:,10:12)];
Td3.data = Td3.data(1:end-2,:);
Td3.data = [Td3.data(1:70,:);zeros(2,10);Td3.data(71,:);zeros(36,10);Td3.data(end,:)];
Td3.data(:,1) = [1:110]';

Td3.avg = mean(Td3.data(:,2:end)')';


Td15.data = [170;29;99;76;177;11;8;12;39;28;115;25;20;35;29;63;38;13;29;2;45;69;117;1;193;15;20];
Td15.info = '15 month clonal analysis';


TdSoxtable = readtable('230606_Sox2_3months_lineage.xlsx');
TdSox.data = table2array(TdSoxtable);
TdSox.info = '3 month clonal analysis Sox2';

TdSox.data = TdSox.data(3:end,:);
TdSox.avg = mean(TdSox.data(:,2:end)')';
TdSox.sum = sum(TdSox.data(:,2:end)')';

nbrBins = ceil(log2(TdSox.data(end,1)));
binsSox2 = zeros(nbrBins,1);
for kk = 1:size(TdSox.data,1)
    addToBin = floor(log2(TdSox.data(kk,1)))+1;
    binsSox2(addToBin) = binsSox2(addToBin)+TdSox.avg(kk);
end
