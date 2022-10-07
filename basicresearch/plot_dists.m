function plotdists
  close all

  dat = load('logs/qdist.mat');
  qdist_mat = dat.qdist_mat;

  qbin_mids = qdist_mat(:,1);
  qdist_e = qdist_mat(:,2);
  qdist_a = qdist_mat(:,3);
  qdist_b = qdist_mat(:,4);

  plot(qbin_mids,[qdist_e qdist_a qdist_b]);

end
