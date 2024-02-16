function [best_classifier,accuracy_svm, sens_svm, spec_svm, avg_acc, avg_sens, avg_spec,std_acc,std_sens,...
            std_spec,best_acc_indx] = kFoldCrossVal(train_x,train_y)
    k=10;
    indices = crossvalind('Kfold',56, k);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for s = 1:k
        test_indic = find(indices == s);
        train_indic = find(indices ~= s);
        % Var train-test
        train_datas =train_x(train_indic,:);% Find out training data and labels
        train_labels =train_y(train_indic,:);
        test_datas = train_x(test_indic,:);% Find out test data and labels
        test_labels{1,s}= train_y(test_indic,:);
        % Start svm multi-class training, fitcsvm is used for binary classification, fitcecoc is used for multiple classification,
         %pause(1)
        classifier{1,s} = fitcsvm(train_datas,train_labels,'KernelFunction','linear');
        %pause(1)
        save('file.mat', 'classifier');
        %pause(1)
        [predict_label{1,s},x{1,s}] = predict(classifier{1,s}, test_datas);% test
        %pause(1)
        accuracy_svm(s)= length(find(predict_label{1,s} == test_labels{1,s}))/length(test_labels{1,s});% accuracy
        
        %pause(1)
        tp = sum((predict_label{1,s} == -1) & (test_labels{1,s}==-1));
        %pause(1)
        fp = sum((predict_label{1,s} == -1) & (test_labels{1,s}==1));
        %pause(1)
        tn = sum((predict_label{1,s} == 1) & (test_labels{1,s}== 1));
        %pause(1)
        fn = sum((predict_label{1,s} == 1) & (test_labels{1,s}== -1));
        %pause(1)
        sens_svm(s)= tp/(tp+fp);% accuracy
        %pause(1)
        acc_svm(s) = (tp+tn)/(tp+fp+fn+tn);
        %pause(1)
        spec_svm(s) = tn/(tn+fn);% accuracy
        
    end
    %%%pause(1)
    [M,best_acc_indx] = max(accuracy_svm);
    %%pause(1)
    best_classifier = classifier{1,best_acc_indx};
    %%pause(1)
    avg_acc = mean(accuracy_svm,'omitnan');
    %%pause(1)
    std_acc = std(accuracy_svm,'omitnan');
    %%pause(1)
    avg_sens = mean(sens_svm,'omitnan');
    %%pause(1)
    std_sens = std(sens_svm,'omitnan');
    %%pause(1)
    avg_spec = mean(spec_svm,'omitnan');
    %%pause(1)
    std_spec  = std(spec_svm,'omitnan');
    
end



