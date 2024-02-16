function [fC_train, fC_traindata, fC_trainlabel, fC_test, fC_testdata, fC_testlabel] = datapartition(feature, f_label, C, C_label)

if size(feature,1)==11
    PD = 0.8;
feature_label = [feature; f_label];
N = size(feature_label,2);
idx = randperm(N);
feature_train = feature_label(:,idx(1:round(N*PD)));
feature_traindata = feature_train(1:11,:)';
feature_trainlabel = feature_train(12,:)';
feature_test = feature_label(:,idx(round(N*PD)+1:end));
feature_testdata = feature_test(1:11,:)';
feature_testlabel = feature_test(12,:)';

Control_label = [C; C_label];
N = size(Control_label,2);
idx = randperm(N);
Control_train = Control_label(:,idx(1:round(N*PD)));
Control_traindata = Control_train(1:11,:)';
Control_trainlabel = Control_train(12,:)';
Control_test = Control_label(:,idx(round(N*PD)+1:end));
Control_testdata = Control_test(1:11,:)';
Control_testlabel = Control_test(12,:)';
else

PD = 0.8;
s1=size(feature,1);
feature_label = [feature; f_label];
N = size(feature_label,2);
idx = randperm(N);
feature_train = feature_label(:,idx(1:round(N*PD)));
feature_traindata = feature_train(1:s1,:)';
feature_trainlabel = feature_train(s1+1,:)';
feature_test = feature_label(:,idx(round(N*PD)+1:end));
feature_testdata = feature_test(1:s1,:)';
feature_testlabel = feature_test(s1+1,:)';

Control_label = [C; C_label];
N = size(Control_label,2);
idx = randperm(N);
Control_train = Control_label(:,idx(1:round(N*PD)));
Control_traindata = Control_train(1:s1,:)';
Control_trainlabel = Control_train(s1+1,:)';
Control_test = Control_label(:,idx(round(N*PD)+1:end));
Control_testdata = Control_test(1:s1,:)';
Control_testlabel = Control_test(s1+1,:)';
end

fC_train = [feature_train'; Control_train'];
fC_traindata = [feature_traindata; Control_traindata];
fC_trainlabel = [feature_trainlabel; Control_trainlabel];

fC_test = [feature_test'; Control_test'];
fC_testdata = [feature_testdata; Control_testdata]; 
fC_testlabel = [feature_testlabel; Control_testlabel];

end

