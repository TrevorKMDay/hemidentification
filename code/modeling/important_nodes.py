import pickle

lda_file = "results/base/method-lda_outcome-hemi_test-A_hands-all_data-base.pickle"
svm_file = "results/base/method-svm_outcome-hemi_test-A_hands-all_data-base.pickle"
nn_file = "results/base/method-nn_outcome-hemi_test-A_hands-all_data-base.pickle"

with open(lda_file, 'rb') as f:
  lda = pickle.load(f)

with open(svm_file, "rb") as f:
  svm = pickle.load(f)
  
with open(nn_file, "rb") as f:
  nn = pickle.load(f)

  
  
