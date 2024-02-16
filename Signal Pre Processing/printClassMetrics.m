function [F1,precision,recall,accuracy] = printClassMetrics (pred, y)
  verbose=1;
  a=find(pred=="burning");
  b=find(y=="burning");
  pred_val = zeros(size(pred));
  pred_val(a) =1;
  yval= zeros(size(y));
  yval(b) =1;
  accuracy = mean(double(pred_val == yval));
  acc_all0 = mean(double(0 == yval));

  actual_positives = sum(yval == 1);
  actual_negatives = sum(yval == 0);
  true_positives = sum((pred_val == 1) & (yval == 1));
  false_positives = sum((pred_val == 1) & (yval == 0));
  false_negatives = sum((pred_val == 0) & (yval == 1));
  precision = 0; 
  if ( (true_positives + false_positives) > 0)
    precision = true_positives / (true_positives + false_positives);
  end 

  recall = 0; 
  if ( (true_positives + false_negatives) > 0 )
    recall = true_positives / (true_positives + false_negatives);
  end

  F1 = 0; 
  if ( (precision + recall) > 0) 
    F1 = 2 * precision * recall / (precision + recall);
  end
 
  if (verbose) 
    sprintf("|-->  true_positives == %i  (actual positive =%i) \n",true_positives,actual_positives);
    sprintf("|-->  false_positives == %i \n",false_positives);
    sprintf("|-->  false_negatives == %i \n",false_negatives);
    sprintf("|-->  precision == %f \n",precision);
    sprintf("|-->  recall == %f \n",recall);
    sprintf("|-->  F1 == %f \n",F1);
  end
  
end
